reDataToDirectMA = function(input.df, dataType) {
  #Rearrange data from gemtc input format into a format suitable for direct meta-analysis
  #dataType - Character string indicating what type of data has been provided.
  # Must be one of: 'treatment difference', 'binary'

  #get a list of unique study IDs and count how many
  studyID = unique(input.df$study)

  for (s in studyID) {
    #pull out the current study
    study = dplyr::filter(input.df, study == s)

    #identify the set of all possible pairwise comparisons in this study
    comparisons = combn(study$treatment, 2)

    #rearrange treatment difference data
    if (dataType == 'treatment difference') {
      #set up a temporary data frame
      df = data_frame(
        StudyName = NA, study = NA, comparator = NA, treatment = NA, diff = NA,
        std.err = NA, NumberAnalysedComparator = NA, NumberAnalysedTreatment = NA,
        ComparatorName = NA, TreatmentName = NA
      )

      #loop through the set of comparisons and rearrange the data
      for (i in 1:ncol(comparisons)) {
        comp = dplyr::filter(study, study$treatment %in% comparisons[,i])

        #only report the direct comparisons as reported in the data
        if (any(is.na(comp$diff))) {
          df[i, 'StudyName'] = as.character(comp$StudyName[1])
          df[i,'study'] = as.integer(s)
          df[i,'treatment'] = as.integer(comp$treatment[2])
          df[i,'comparator'] = as.integer(comp$treatment[1])
          df[i,'diff'] = comp$diff[2]
          df[i,'std.err'] = comp$std.err[2]
          df[i,'NumberAnalysedComparator'] = as.integer(comp$NumberAnalysed[1])
          df[i,'NumberAnalysedTreatment'] = as.integer(comp$NumberAnalysed[2])
          df[i,'ComparatorName'] = as.character(comp$TreatmentName[1])
          df[i,'TreatmentName'] = as.character(comp$TreatmentName[2])
        }
      }

      if (s == 1) {
        directData = df
      } else {
        directData = dplyr::bind_rows(directData, df)
      }
    }

    #rearrange binary data
    if (dataType == 'binary') {
      #set up a temporary data frame
      df = dplyr::data_frame(
        StudyName = NA, study = NA, comparator = NA, treatment = NA,
        NumberEventsComparator = NA, NumberAnalysedComparator = NA,
        NumberEventsTreatment = NA, NumberAnalysedTreatment = NA,
        ComparatorName = NA, TreatmentName = NA
      )

      #loop through the set of comparisons and rearrange the data
      for (i in 1:ncol(comparisons)) {
        comp = dplyr::filter(study, study$treatment %in% comparisons[,i])
        df[i, 'StudyName'] = as.character(comp$StudyName[1])
        df[i,'study'] = as.integer(s)
        df[i,'treatment'] = as.integer(comp$treatment[2])
        df[i,'comparator'] = as.integer(comp$treatment[1])
        df[i,'NumberEventsComparator'] = as.integer(comp$responders[1])
        df[i,'NumberAnalysedComparator'] = as.integer(comp$sampleSize[1])
        df[i,'NumberEventsTreatment'] = as.integer(comp$responders[2])
        df[i,'NumberAnalysedTreatment'] = as.integer(comp$sampleSize[2])
        df[i,'ComparatorName'] = as.character(comp$TreatmentName[1])
        df[i,'TreatmentName'] = as.character(comp$TreatmentName[2])
      }

      if (s == 1) {
        directData = df
      } else {
        directData = dplyr::bind_rows(directData, df)
      }
    }
  }
  return(directData)
}

doDirectMeta = function(df, effectCode, dataType, backtransf = FALSE) {
  #dataType - Character string indicating what type of data has been provided.
  # Must be one of: 'treatment difference', 'binary'

  #create a list object to store the results
  resList = list()

  #identify the set of treatment comparisons present in the data
  comparisons = dplyr::distinct(df[,3:4])

  for (i in 1:nrow(comparisons)) {
    #get data for the first comparison
    comp = dplyr::filter(df, comparator == comparisons$comparator[i],
                  treatment == comparisons$treatment[i])

    #run the analysis for different data types
    if (dataType == 'treatment difference') {
      #Generic inverse variance method for treatment differences
      #backtransf=TRUE converts log effect estimates (e.g log OR) back to linear scale
      directRes = meta::metagen(
        TE = comp$diff, seTE = comp$std.err,sm = effectCode, backtransf = backtransf,
        studlab = comp$StudyName, n.e = comp$NumberAnalysedTreatment,
        n.c = comp$NumberAnalysedComparator, label.e = comp$TreatmentName,
        label.c = comp$ComparatorName
      )
    }
    if (dataType == 'binary') {
      #Analysis of binary data provided as n/N
      directRes = meta::metabin(
        event.e = comp$NumberEventsTreatment, n.e = comp$NumberAnalysedTreatment,
        event.c = comp$NumberEventsComparator, n.c = comp$NumberAnalysedComparator,
        sm = effectCode, backtransf = backtransf, studlab = comp$StudyName,
        label.e = comp$TreatmentName, label.c = comp$ComparatorName,
        method = ifelse(nrow(comp) > 1, "MH", "Inverse")
      )
    }
    #add the treatment codes to the results object. These will be needed later
    directRes$e.code = comparisons$treatment[i]
    directRes$c.code = comparisons$comparator[i]

    #compile the results into a list
    if (length(resList) == 0) {
      resList[[1]] = directRes
    } else {
      resList[[length(resList) + 1]] = directRes
    }
  }
  return(resList)
}

backtransform = function(df) {
  #simple function to convert log OR (or HR or RR) back to a linear scale
  #df - a data frame derived from a metagen summary object

  #relabel column names to preserve the log results
  colnames(df)[c(1,3:4)] = paste0('log.', colnames(df)[c(1,3:4)])
  colnames(df)[2] = paste0(colnames(df)[2], '.log')
  #exponentiate the effect estimate and CI
  df$TE = exp(df$log.TE)
  df$lower = exp(df$log.lower)
  df$upper = exp(df$log.upper)
  df
}

extractDirectRes = function(metaRes, effect, intervention = 'Int', comparator = 'Con',
                            interventionCode = NA, comparatorCode = NA,
                            backtransf = FALSE) {
  #metaRes - an object of class c("metagen", "meta") as returned by the function
  #metagen in the package meta effect - a character string describing the effect
  #estimate, e.g. 'Rate Ratio', 'Odds Ratio', 'Hazard Ratio' intervention - name
  #of the intervention treatment comparator - name of the comparator treatment
  #backtransf - logical indicating whether the results should be exponentiated
  #or not. If the results of the meta-analysis are log odds ratios set this to
  #TRUE to return the odds ratios. If TRUE this will return both the log
  #estimates and the exponentiated estimates

  #create a summary of the results, extract fixed and random
  res = summary(metaRes)
  fixed = as.data.frame(res$fixed)
  random = as.data.frame(res$random)[1:7]
  if (backtransf == TRUE) {
    #exponentiate the effect estimates if required
    fixed = backtransform(fixed)
    random = backtransform(random)
  }

  #if more than one study then this must be a meta-analysis and both fixed and
  #random are expected if there is only one study then use the fixed results as
  #all results are just the result of the original study
  if (res$k > 1) {
    df = rbind(fixed, random)
    model = c('Fixed', 'Random')
    df = data.frame('Model' = model, df, stringsAsFactors = FALSE)
  } else {
    df = fixed
    df = data.frame('Model' = NA, df, stringsAsFactors = FALSE)
  }
  studies = paste0(metaRes$studlab, collapse = ', ')
  df = data.frame(
    'Intervention' = intervention, 'InterventionCode' = interventionCode,
    'Comparator' = comparator, 'ComparatorCode' = comparatorCode,
    'Effect' = effect, df,'Tau.sq' = res$tau,
    'method.tau' = res$method.tau, 'I.sq' = res$I2, 'n.studies' = res$k,
    'studies' = studies, stringsAsFactors = FALSE
  )

}

drawForest = function(meta, showFixed = TRUE, showRandom = TRUE, ...) {
  #meta - an object of class c("metagen", "meta") as returned by the function
  #metagen in the package meta showFixed, showRandom - logical indicating
  #whether pooled estimates from fixed or random effects models should be be
  #shown on the forest plot

  #work out sensible values for the x-axis limits
  limits = c(
    meta$lower, meta$upper, meta$lower.fixed, meta$upper.fixed,
    meta$lower.random, meta$upper.random
  )
  limits = range(exp(limits))
  xlower = ifelse(limits[1] < 0.2, round(limits[1], 1), 0.2)
  xupper = ifelse(limits[2] > 5, round(limits[2]), 5)
  xlimits = c(xlower, xupper)

  #don't show the pooled estimate if there is only one study
  if (meta$k == 1) {
    showFixed = FALSE
    showRandom = FALSE
  }

  #forest plot
  meta::forest(
    meta, hetlab = NULL, text.I2 = 'I-sq', text.tau2 = 'tau-sq', xlim = xlimits,
    comb.fixed = showFixed, comb.random = showRandom, lty.fixed = 0,
    lty.random = 0, ...
  )
}
