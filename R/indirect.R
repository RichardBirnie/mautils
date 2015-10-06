#' Run simple indirect meta-analysis for all possible pairwise comparisons
#' in a dataset
#'
#' @param df A \code{data.frame} This should be in one of two formats. Arm level
#'   data must contain the columns 'study' and 'treatment' where study is a
#'   study id number (1, 2, 3 ...) and treatment is a treatment id number. If
#'   the data are binary then the data frame should also contain columns
#'   'responders' and 'sampleSize' reporting the total number of events and
#'   total number analysed respectively for each arm. Relative effect data (e.g.
#'   log odds ratio, log rate ratio) must contain the same study and treatment
#'   columns plus the columns 'diff' and 'std.err'. Set diff=NA for the baseline
#'   arm. The column std.err should be the standard error of the relative effect
#'   estimate. For trials with more than two arms set std.err as the standard
#'   error of the baseline arm. This determines the covariance which is used to
#'   adjust for the correlation in multiarm studies. This is used to identify
#'   the set of comparisons available in the data set
#' @param data_type A character string specifying which type of data has been
#'   provided. Currently only 'treatment difference' or 'binary' are supported
#' @param direct_results A data frame containing the results of direct
#'   meta-analysis as returned by \code{runDirect}. These results are required
#'   to provide the inputs for the indirect comparisons
#' @param effect_type A character string indicating what kind of analysis is
#'   required. Set to 'Fixed' for fixed effect, 'Random' for random effects or
#'   'all' to get both (Default).
#' @param back_calc A logical indicating whether results should be back transformed.
#'   This is used to set the corresponding \code{backtransf} argument of the
#'   underlying functions from the \code{meta} package. If
#'   \code{backtransf=TRUE} then log odds ratios (or hazard ratios etc) will
#'   be converted to odds ratios on plots and print outs. Default is FALSE
#' @param order_treatments An optional argument to specify the order in which
#'   treatment comparisons are sorted in the output. The default is NA in
#'   which case comparisons will be sorted alphabetically by intervention. If
#'   a specific order is required then this should be provided as a data frame
#'   with two columns named 'description' and 'Order'. Note that column
#'   headers are specific and case sensitive. The description column should
#'   contain the names of the treatments \emph{exactly} as they are specified
#'   in the data set. The id column should contain the numbered order of
#'   treatments required.
#'
#' @details This function performs indirect meta-analysis using the Bucher
#'   method for all possible comparisons in a given data set. This function
#'   takes a set of treatment comparisons from one or more studies and
#'   identifies all possible indirect comparisons where two treatments can be
#'   connected via a common comparator. If there is more than one way to
#'   connect two treatments then all possible variations are calculated. This
#'   function calls \code{doBucher} internally to calculate the
#'   treatment effects
#'
#'   The inputs for this function will usually be the results from direct
#'   meta-analysis for a given set of treatments. The recommended workflow is
#'   to use \code{\link{runDirect}} to perform head to head meta-analysis for
#'   a given set of treatments then use the resulting data frame to provide
#'   the inputs for this function.
#'
#' @return A data frame containing the results of all possible indirect
#'   comparisons in the data set. The help page for \code{\link{doBucher}}
#'   provides a detailed description of the columns in the output
#'
#' @seealso \code{\link{doBucher}}, \code{\link{runDirect}}
runIndirect = function(df, data_type, direct_results, effect_type = 'all',
                       back_calc = FALSE, order_treatments = NA) {
  message('Run indirect (Bucher) meta-analysis')

  #convert data to a suitable format
  reDataDirect = formatDataToDirectMA(df, dataType = data_type)

  indMA = doBucher(
    comparisons = reDataDirect[,1:4], direct = direct_results,
    effectType = effect_type, backtransf = back_calc
  )
  if(is.data.frame(order_treatments)){
    #set the order of results
    indMA = dplyr::left_join(indMA, order_treatments, by = c('Intervention' = 'description'))
    indMA = dplyr::arrange(indMA, Order, Comparator, Common, Model)
  } else {
    indMA = dplyr::arrange(indMA, Intervention, Model)
  }
}

#'Create an igraph object
#' @param edgelist A two column matrix describing the connections _from_ the first
#'   column _to_ the second column. The values are the treatment numbers
#' @param Logical. Describes whether all studies of the same treatments should be
#'   shown as separate connections in the graph or collapsed into a single edge.
#'   If there are 3 studies comparing the same treatments all=TRUE will show three
#'   lines, all=FALSE will only show one(Default)
#' @param plotNetwork Logical indicating whether to draw a plot of the
#'   network or not
#'
#' @return An igraph object
.buildGraph = function(edgelist, all = FALSE, plotNetwork = FALSE) {
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

#' Return all combinations of nodes with a set path length
#'
#' @param g An igraph object
#' @param len Integer. The desired path length. Default = 2.
#'
#' @details A function to search an igraph object and return
#'  all combinations of nodes connected by a certain path
#'  length. The default is to return all combinations of three
#'  nodes that can be connected in exactly two steps.
#'
#' @return A data frame
#'
.getPathsByLength = function(g, len = 2) {
  sp = igraph::shortest.paths(g)
  sp[upper.tri(sp,TRUE)] = NA
  wp = which(sp == len, arr.ind = TRUE)
  mapply(function(a,b)
    igraph::get.all.shortest.paths(
      g, from = a, to = b, mode = 'all'
    ), wp[,1], wp[,2])
}

#' Simple indirect (Bucher) meta-analysis
#'
#' @param abTE Numeric. The treatment effect for a vs b. e.g. log OR, log HR,
#'   mean difference
#' @param se.abTE Numeric. The standard error of the treatment effect for
#'   a vs b, e.g. se of log OR
#' @param cbTE Numeric. The treatment effect for c vs b. e.g. log OR, log HR,
#'   mean difference
#' @param se.cbTE Numeric. The standard error of the treatment effect for c vs
#'   b, e.g. se of log OR
#' @param effect Character string describing the effect measure, e.g. 'Rate
#'   Ratio' or 'log Odds Ratio'
#' @param model Character string indicating whether abTE and cbTE come from a
#'   fixed effect model or a random effect model
#' @param intervention Character string. Name of the intervention treatment
#' @param comparator Character string. Name of the comparator treatment
#' @param common Character string. Name of the common comparator that links the
#'   indirect comparison, e.g. placebo
#' @param backtransf Logical indicating whether the results should be
#'   exponentiated or not. If abTE and cbTE are on the log scale (e.g. log
#'   hazard ratio) set this to TRUE to return the exponentiated results (e.g.
#'   hazard ratio). If TRUE this will return both the log estimates and the
#'   exponentiated estimates
#'
#' @details Calculate an indirect estimate of the relative effect of two
#'   treatments using the Bucher method
#'
#'   The inputs are the relative treatment effects for two pairs of treatments
#'   linked by a common comparator. For example, if you have the relative
#'   effects (e.g. log hazard ratio) for treatment A vs placebo and treatment C
#'   vs placebo then this function will return the relative effect of treatment
#'   A compared to treatment C.
#'
#'   This function is mainly intended to be called from \code{\link{doBucher}} but can
#'   be used directly if required.
#'
#'   \code{effect}, \code{model}, \code{intervention}, \code{comparator} and
#'   \code{common} are used for labels only. The results depend only on
#'   \code{abTE}, \code{cbTE} and the respective standard errors.
#'
#'   @return A data frame with the following columns:
#' \itemize{
#'  \item \code{Intervention} The name of the intervention
#'  \item \code{Comparator} The name of the comparator
#'  \item \code{Common} The name of the common treatment linking intervention
#'  and comparator
#'  \item \code{Effect} The type of effect measure. Takes the
#'    value of the \code{effect} argument
#'  \item \code{Model} The type of model. Fixed effect or Random Effects.
#'  \item \code{log.TE.ind} The treatment effect on log scale, e.g. log OR
#'  \item \code{log.lower.ind}, \code{log.upper.ind} The upper and lower 95\% confidence
#'    intervals for the log treatment effect
#'  \item \code{se.log.TE.ind} The standard error for the log treatment effect

#'  \item \code{TE.ind}, \code{lower.ind}, \code{upper.ind} The treatment effect with lower
#'    and upper confidence intervals backtransformed to a linear scale
#'  \item \code{n.studies} The number of studies included in the analysis
#'  \item \code{Studies} The names of the studies included in the analysis
#' }
#'
#'   @seealso \code{\link{doBucher}}

bucher = function(abTE, se.abTE, cbTE, se.cbTE, effect, model,
                  intervention, comparator, common, backtransf = FALSE,
                  ab.studies, cb.studies) {

  acTE = abTE - cbTE #indirect treatment effect
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

#' Perform all possible indirect comparisons for a given data set
#'
#' @param comparisons A data frame with four columns: StudyName, study,
#'   comparator, treatment. Describes the treatment comparisons present in the
#'   dataset. study, comparator and treatment must be numbers. For example,
#'   study = 4, comparator=1, treatment=2 represents the comparison of treatment
#'   2 vs treatment 1 in study 4
#'@param direct A data frame containing the results of direct head-to-head
#'  meta-analysis for the treatment comparisons of interest if only one study is
#'  available for a given comparison then the result of that study should be
#'  used. This data frame can be created by using \code{doDirectMeta} and
#'  \code{extractDirectRes} in that order
#' @param effectType Character string indicating what type of results are
#'   required. Default is 'all' which will return both fixed effect and random
#'   effect results. Alternatives are 'Fixed' or 'Random' (Case sensitive) if
#'   only one set of results is required
#'
#' @details This function performs indirect meta-analysis for all possible
#'   comparisons in a given data set. This function takes a set of treatment
#'   comparisons from one or more studies and identifies all possible indirect
#'   comparisons where two treatments can be connected via a common comparator.
#'   If there is more than one way to connect two treatments then all possible
#'   variations are calculated. This function calls \code{\link{bucher}} internally
#'   to calculate the treatment effects
#'
#'   The inputs for this function will usually be the results from direct
#'   meta-analysis for a given set of treatments. The recommended workflow is
#'   to use \code{\link{doDirectMeta}} to perform head to head meta-analysis for
#'   a given set of treatments, extract the results as a data frame using
#'   \code{\link{extractDirectRes}} then use that data frame to provide the inputs
#'   for this function.
#'
#'   @return A data frame with the following columns:
#' \itemize{
#'  \item \code{Intervention} The name of the intervention
#'  \item \code{Comparator} The name of the comparator
#'  \item \code{Common} The name of the common treatment linking intervention
#'  and comparator
#'  \item \code{Effect} The type of effect measure. Takes the
#'    value of the \code{effect} argument
#'  \item \code{Model} The type of model. Fixed effect or Random Effects.
#'  \item \code{log.TE.ind} The treatment effect on log scale, e.g. log OR
#'  \item \code{log.lower.ind}, \code{log.upper.ind} The upper and lower 95\% confidence
#'    intervals for the log treatment effect
#'  \item \code{se.log.TE.ind} The standard error for the log treatment effect
#'  \item \code{TE.ind}, \code{lower.ind}, \code{upper.ind} The treatment effect with lower
#'    and upper confidence intervals backtransformed to a linear scale
#'  \item \code{n.studies} The number of studies included in the analysis
#'  \item \code{Studies} The names of the studies included in the analysis
#' }
#' @seealso \code{\link{bucher}}, \code{\link{doDirectMeta}},
#'   \code{\link{extractDirectRes}}
doBucher = function(comparisons, direct, effectType = 'all',
                    backtransf = FALSE) {

  #take information describing available comparisons from the direct MA dataset
  connections = dplyr::select(comparisons, from = treatment, to = comparator)
  g = .buildGraph(connections)

  #get the set of possible Bucher comparisons
  bucherSet = .getPathsByLength(g, 2)

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

  #rearrange the columns into a more report friendly format
  df = dplyr::select(
    df, Effect, Intervention, Comparator, Common, TE.ind, lower.ind, upper.ind,
    n.studies, Studies, Model, log.TE.ind, log.lower.ind, log.upper.ind,
    se.log.TE.ind
  )

}
