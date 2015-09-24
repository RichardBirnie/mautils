buildGraph = function(edgelist, all = FALSE, plotNetwork = FALSE) {
  #create a graph object edgelist - a two column matrix describing the
  #connections _from_ the first column _to_ the second column. The values are
  #the treatment numbers all - logical. Describes whether all studies of the
  #same treatments should be shown as separate connections in the graph or
  #collapsed into a single edge. If there are 3 studies comparing the same
  #treatments all=TRUE will show three lines, all=FALSE will only show one
  #(Default) plotNetwork - logical indicating whether to draw a plot of the
  #network or not

  #collapse the edges unless explicitly requested not to
  if (all != TRUE) {
    el = as.matrix(dplyr::distinct(edgelist))
  } else {
    el = as.matrix(edgelist)
  }
  g = igraph::graph_from_edgelist(el, directed = TRUE)

  if (plotNetwork == TRUE) {
    plot(g, layout = igraph::layout.fruchterman.reingold)
  }

  return(g)
}

getPathsByLength = function(g, len = 2) {
  #Function to return all combinations of 3 connected treatments identifies all
  #paths of length 2 (3 nodes, 2 edges) If there are multiple ways to connect
  #two treatments then all possible combinations are returned g - a graph object
  #as produced by igraph len - the length of path to be identified (Default=2).
  #Connecting 3 nodes has a path length of two
  sp = igraph::shortest.paths(g)
  sp[upper.tri(sp,TRUE)] = NA
  wp = which(sp == len, arr.ind = TRUE)
  mapply(function(a,b)
    igraph::get.all.shortest.paths(
      g, from = a, to = b, mode = 'all'
    ), wp[,1], wp[,2])
}

bucher = function(abTE, se.abTE, cbTE, se.cbTE, effect, model,
                  intervention, comparator, common, backtransf = FALSE,
                  ab.studies, cb.studies) {
  #abTE - treatment effect for a vs b. e.g. log OR, log HR, mean difference
  #se.abTE - standard error of the treatment effect for a vs b, e.g. se of log OR
  #cbTE - treatment effect for c vs b. e.g. log OR, log HR, mean difference
  #se.cbTE - standard error of the treatment effect for c vs b, e.g. se of log OR
  #effect - Character string describing the effect measure, e.g. 'Rate Ratio' or 'log Odds Ratio'
  #model - Character string indicating whether abTE and cbTE come from a fixed effect model or a
  #   random effect model
  #intervention - Character string. Name of the intervention treatment
  #comparator - Character string. Name of the comparator treatment
  #common - Character string. Name of the common comparator, e.g. placebo
  #backtransf - logical indicating whether the results should be exponentiated or not. If abTE and cbTE are
  #on the log scale set this to TRUE to return the exponentiated results. If TRUE this will return both the
  #log estimates and the exponentiated estimates
  #
  #effect, model, intervention, comparator and common are used for labels only.
  #The results depend only on abTE, cbTE and the respective standard errors.

  acTE = abTE - cbTE #treatment effect
  se.acTE = sqrt(se.abTE ^ 2 + se.cbTE ^ 2) #standard error - sqrt(sum of the two variances)
  log.lower = acTE - (1.96 * se.acTE) #calculate confidence intervals
  log.upper = acTE + (1.96 * se.acTE)
  df = data.frame(
    'log.TE.ind' = acTE, 'log.lower.ind' = log.lower, 'log.upper.ind' = log.upper, 'se.log.TE.ind' =
      se.acTE
  )

  #exponentiate if required
  if (backtransf == TRUE) {
    eTE = exp(df[,1:3])
    colnames(eTE) = gsub('log.', '', colnames(eTE), fixed = TRUE)
    df = cbind(df, eTE)
  }

  #identify and count the studies involved
  ab.studies = strsplit(as.character(ab.studies), split = ', ')[[1]]
  cb.studies = strsplit(as.character(cb.studies), split = ', ')[[1]]
  studies = c(ab.studies,cb.studies)
  studies = unique(studies)
  n.studies = length(studies)
  studies = paste0(studies, collapse = ', ')

  #tag on some labels
  df = data.frame(
    'Intervention' = intervention, 'Comparator' = comparator, 'Common' = common,
    'Effect' = effect, 'Model' = model, df, 'n.studies' = n.studies,
    'Studies' = studies
  )
}

doBucher = function(comparisons, direct, effectType = 'all', backtransf =
                      FALSE) {
  #comparisons - a data frame with four columns: StudyName, study, comparator, treatment. Describes the treatment comparisons
  # present in the dataset. study, comparator and treatment must be numbers. For example, study = 4, comparator=1, treatment=2
  # represents the compartison of treatment 2 vs treatment 1 in study 4
  #direct - a data frame containing the results of direct head-to-head meta-analysis for the treatment comparisons of interest
  # if only one study is available for a given comparison then the result of that study should be used. This data frame can be
  # created by using doDirectMeta and extractDirectRes in order
  #effectType - character string indicating what type of results are required. Default is all which will return both fixed effect
  # and random effect results. Alternatives are 'Fixed' or 'Random' (Case sensitive) if only one set of results is required

  #take information describing available comparisons from the direct MA dataset
  connections = dplyr::select(comparisons, from = treatment, to = comparator)
  g = buildGraph(connections)

  #get the set of possible BUcher comparisons
  bucherSet = getPathsByLength(g, 2)

  for (i in 1:ncol(bucherSet)) {
    #for each comparison get the list of possible alternatives
    res = bucherSet[,i]$res

    for (j in 1:length(res)) {
      #take the current triplet and get the edges connecting those treatments
      #the middle element of the triplet is the common comparator
      triplet = as.numeric(bucherSet[,i]$res[[j]])

      #identify the corresponding comparisons in the direct results
      #control which set of results we want. Default is all
      if (effectType == 'all') {
        mod = c('Fixed', 'Random')
      } else {
        mod = effectType
      }
      for (k in 1:length(mod)) {
        currComp = filter(
          direct, (direct$InterventionCode %in% triplet & direct$ComparatorCode %in% triplet),
          (Model == mod[k] | is.na(Model))
        )

        #if the common treatment is the intervention rather than the comparator
        #then we need to invert the treatment effect
        ix = currComp$InterventionCode == triplet[2]
        if (any(ix) == TRUE) {
          currComp$log.TE[ix] = -currComp$log.TE[ix]
          currComp[ix, c('InterventionCode', 'ComparatorCode')] = currComp[ix, c('ComparatorCode', 'InterventionCode')]
          currComp[ix, c('Intervention', 'Comparator')] = currComp[ix, c('Comparator', 'Intervention')]
        }

        #run the comparison
        e = currComp$Effect[1]
        m = currComp$Model[!is.na(currComp$Model)]
        abTE = currComp$log.TE[currComp$InterventionCode == triplet[1]]
        se.abTE = currComp$seTE.log[currComp$InterventionCode == triplet[1]]
        cbTE = currComp$log.TE[currComp$InterventionCode == triplet[3]]
        se.cbTE = currComp$seTE.log[currComp$InterventionCode == triplet[3]]
        int = currComp$Intervention[currComp$InterventionCode == triplet[1]]
        com = currComp$Intervention[currComp$InterventionCode == triplet[3]]
        indMA = bucher(
          abTE = abTE, se.abTE = se.abTE,
          cbTE = cbTE, se.cbTE = se.cbTE,
          backtransf = backtransf, effect = e,
          model = ifelse(length(m) == 0, NA, m),
          intervention = int, comparator = com,
          common = currComp$Comparator[1],
          ab.studies = currComp$studies[1],
          cb.studies = currComp$studies[2]
        )

        #collect FE and RE results
        if (k == 1) {
          fere = indMA
        } else {
          fere = rbind(fere, indMA)
        }

      }
      #assemble the results
      if (j == 1) {
        ind.df = fere
      } else {
        ind.df = rbind(ind.df, fere)
      }
    }
    if (i == 1) {
      df = ind.df
    } else {
      df = rbind(df, ind.df)
    }
  }
  #calculate the inverse comparisons
  df_inv = df
  df_inv[, c('Intervention', 'Comparator')] = df[, c('Comparator', 'Intervention')]
  df_inv$log.TE.ind = -df$log.TE.ind
  df_inv$log.lower.ind = df_inv$log.TE.ind - (1.96 * df_inv$se.log.TE.ind)
  df_inv$log.upper.ind = df_inv$log.TE.ind + (1.96 * df_inv$se.log.TE.ind)
  df_inv[,10:12] = exp(df_inv[,6:8])

  #combine results
  df = rbind(df, df_inv)
  df = dplyr::arrange(df, dplyr::desc(Intervention), dplyr::desc(Comparator))

  #fix variable types
  df = rbutils::factorToCharacter(df)

}
