#' Run direct head to head meta-analysis for all possible pairwise comparisons
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
#'   adjust for the correlation in multiarm studies.
#' @param treatments An optional (but recommended) data frame describing the
#'   treatments. If present then this data frame should contain the columns
#'   'id', 'description', 'shortname'. id gives the ID number of each treatment,
#'   description gives the full name of each treatment, shortname gives a
#'   shortened version of the name for each treatment. shortname is used for
#'   labels and filenames in forest plots which is the main purpose of this
#'   argument. If this is not provided then forest plots will be named
#'   'Comparison1', Comparison2 etc.
#' @param file A character string specifying the directory where the results
#'   will be saved. This function will create two new subdirectories
#'   'Results/Direct' and 'Results/Direct/Figures' to store the output. These
#'   are only created if they do not already exist
#' @param data_type A character string specifying which type of data has been
#'   provided. Must be one of 'treatment difference', 'binary' or 'continuous'
#' @param effect_code A character string indicating the underlying effect
#'   estimate. This is used to set the \code{sm} argument of the underlying
#'   analysis functions from the \code{meta} package. Acceptable values are
#'   'RD', 'RR', 'OR', 'HR', 'ASD', 'MD', 'SMD'
#' @param outcome A character string indicating what outcome we are currently
#'   looking at. This is mainly used to automatically generate meaningful file
#'   names for output so keep it short.
#' @param effect_measure A character string indicating what type of effect
#'   measure is used, e.g. 'Rate Ratio', 'Odds Ratio' etc. This is used as a
#'   label in forest plots so keep it short.
#' @param back_calc A logical indicating whether results should be back transformed.
#'   This is used to set the corresponding \code{backtransf} argument of the
#'   underlying functions from the \code{meta} package. If
#'   \code{backtransf=TRUE} then log odds ratios (or hazard ratios etc) will
#'   be converted to odds ratios on plots and print outs. Default is FALSE
#' @param forest_plot A logical indicating whether a forest plot should be
#'   produced for each pairwise comparison.
#' @param show_fixed,show_random Logical indicating whether fixed effect
#'   and/or random effects results should be shown on the forest plot. By
#'   default both are TRUE and both results are shown. Set the appropriate
#'   argument to FALSE if you want to exclude that result from the plot.
#' @param method A character string indicating which method is used for pooling
#'   studies. This argument is only applicable if data_type='binary'. Must be one
#'   of 'MH' (Default, Mantel-Haenszel method), 'inverse' or 'Peto'
#'
#' @details This function runs direct head to head meta-analysis for all
#'   possible pairwise comparisons in a dataset. The objective is to provide a
#'   minimal set of options and abstract as much of the technical detail as
#'   reasonably possible. This function calls a series of internal functions
#'   to do the work (linked below). These functions may accessed directly if
#'   additional flexibility is required but the intention is that this wrapper
#'   should cover the majority of common applications. The actual
#'   meta-analysis depends on functions from the \code{meta} package as does
#'   drawing forest plots. Currently only treatment differences (OR, HR, RR
#'   etc) and binary data are supported. Additional data types will be added
#'   in the future.
#'
#'   The basic order of events is:
#'   \itemize{
#'    \item Run all possible pairwise comparisons in the dataset.
#'    \item Draw forest plots if requested (see \code{forest_plot} above).
#'    These are saved as jpg files in a subfolder named figures.
#'    \item Extract the key information from the meta-analysis output and save
#'    this as an excel file
#'   }
#'
#' @return A data frame containing the results of the meta-analysis
#'
#' @seealso \code{\link{formatDataToDirectMA}}, \code{\link{doDirectMeta}},
#'   \code{\link{drawForest}}, \code{\link{extractDirectRes}}
#' @export
runDirect = function(df, treatments=NULL, file, data_type, effect_code, outcome, effect_measure,
                     back_calc = FALSE, forest_plot = TRUE,
                     show_fixed = TRUE, show_random = TRUE,
                     method='MH') {

  message('Run direct meta-analysis')
  #set up folders for the results and figures. No need to edit this
  directResultDir = file.path(baseFile, 'Results','Direct')
  if (!file.exists(directResultDir)) {
    dir.create(directResultDir, recursive = TRUE)
  }
  directFigDir = file.path(directResultDir, 'Figures')
  if (!file.exists(directFigDir)) {
    dir.create(directFigDir, recursive = TRUE)
  }

  #convert data to format suitable for direct meta-analysis using the meta packages
  reDataDirect = formatDataToDirectMA(df, dataType = data_type)

  #run the meta-analysis
  directRes = doDirectMeta(
    df = reDataDirect, dataType = data_type, effectCode = effect_code,
    backtransf = back_calc, method = method
  )

  #draw the forest plots if requested
  if(forest_plot) {
    message('Drawing forest plots')
    for (i in 1:length(directRes)) {

      #set up meaningful file names and treatment labels if possible
      if('shortname' %in% colnames(treatments)) {
        int = treatments$shortname[treatments$id == directRes[[i]]$e.code]
        com = treatments$shortname[treatments$id == directRes[[i]]$c.code]
        f = paste0(outcome, '_', int, '_', com, '.jpg')
        directRes[[i]]$e.name = int
        directRes[[i]]$c.name = com
      } else {
        #fall back if nothing more helpful is provided
        f = paste0(outcome, '_Comparison', i, '.jpg')
        directRes[[i]]$e.name = directRes[[i]]$label.e
        directRes[[i]]$c.name = directRes[[i]]$label.c
      }

      figFile = file.path(directFigDir, f)
      jpeg(
        file = figFile, width = 25, height = 15, units = 'cm', res = 300,
        quality = 100
      )
      drawForest(
        directRes[[i]], col.square = 'red', col.diamond = 'black',
        smlab = effect_measure, showFixed = show_fixed,
        showRandom = show_random
      )
      graphics.off()
    }
  }

  #extract the results
  #get results in a useful format
  message('Extract results')
  for (i in 1:length(directRes)) {
    res = directRes[[i]]
    res = extractDirectRes(
      metaRes = res, effect = effect_measure, backtransf = back_calc,
      intervention = res$label.e[1], comparator = res$label.c[1],
      interventionCode = res$e.code, comparatorCode = res$c.code
    )
    if (i == 1) {
      dirMA = res
    } else {
      dirMA = dplyr::bind_rows(dirMA, res)
    }
  }
  #return the table of results
  return(dirMA)
}

#' Rearrange data from gemtc input format to a format suitable for direct
#' meta-analysis
#'
#' @param input.df \code{data.frame} This should be in one of two formats. Arm level
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
#'   adjust for the correlation in multiarm studies.
#' @param dataType A character string specifying which type of data has been
#'   provided. Must be one of 'treatment difference', 'binary' or 'continuous'
#'
#' @return A data frame
#'
#' @seealso \code{\link[gemtc]{mtc.network}}
#' @export
formatDataToDirectMA = function(input.df, dataType) {

  #if data is treatment difference then change column name for convenience
  if(dataType == 'treatment difference') {colnames(input.df)[2] = 'sampleSize'}
  #get a list of unique study IDs
  studyID = unique(input.df$study)

  for (s in studyID) {
    #pull out the current study
    study = dplyr::filter(input.df, study == s)

    #identify the set of all possible pairwise comparisons in this study
    comparisons = combn(study$treatment, 2)

    df = dplyr::data_frame(
      StudyName = NA, study = NA, comparator = NA, treatment = NA,
      NumberAnalysedComparator = NA, NumberAnalysedTreatment = NA,
      ComparatorName = NA, TreatmentName = NA
    )

    #loop through the set of comparisons and rearrange the data
    for (i in 1:ncol(comparisons)) {
      #handle generic columns first. Same for all data types
      comp = dplyr::filter(study, study$treatment %in% comparisons[,i])
      if (dataType == 'treatment difference' && !anyNA(comp$diff)) {
        next #exception for three arm studies. Usually only two effects reported
      } else {
        df[i,'StudyName'] = as.character(comp$StudyName[1])
        df[i,'study'] = as.integer(s)
        df[i,'treatment'] = as.integer(comp$treatment[2])
        df[i,'comparator'] = as.integer(comp$treatment[1])
        df[i,'NumberAnalysedComparator'] = as.integer(comp$sampleSize[1])
        df[i,'NumberAnalysedTreatment'] = as.integer(comp$sampleSize[2])
        df[i,'ComparatorName'] = as.character(comp$TreatmentName[1])
        df[i,'TreatmentName'] = as.character(comp$TreatmentName[2])
        if (dataType == 'treatment difference') {
          #only report the direct comparisons as reported in the data
          # if (any(is.na(comp$diff))) {
          df[i,'diff'] = comp$diff[2]
          df[i,'std.err'] = comp$std.err[2]
          df = dplyr::select(
            df, 1:4, diff, std.err,
            NumberAnalysedComparator, NumberAnalysedTreatment,
            ComparatorName, TreatmentName
          )
        }
        if (dataType == 'binary') {
          #binary data only
          df[i,'NumberEventsComparator'] = as.integer(comp$responders[1])
          df[i,'NumberEventsTreatment'] = as.integer(comp$responders[2])
          #rearrange the column order
          df = dplyr::select(
            df, 1:4, NumberEventsComparator, NumberAnalysedComparator,
            NumberEventsTreatment, NumberAnalysedTreatment,
            ComparatorName, TreatmentName
          )
        }
        if (dataType == 'continuous') {
          #continuous data only
          df[i,'MeanComparator'] = as.numeric(comp$mean[1])
          df[i,'SDComparator'] = as.numeric(comp$std.dev[1])
          df[i,'MeanTreatment'] = as.numeric(comp$mean[2])
          df[i,'SDTreatment'] = as.numeric(comp$std.dev[2])
          #rearrange the column order
          df = dplyr::select(
            df, 1:4, MeanComparator, SDComparator, NumberAnalysedComparator,
            MeanTreatment, SDTreatment, NumberAnalysedTreatment,
            ComparatorName, TreatmentName
          )
        }
      }
    }
    if (s == 1) {
      directData = df
    } else {
      directData = dplyr::bind_rows(directData, df)
    }
  }
  return(directData)
}

#' Perform multiple direct head to head meta-analyses from a single data frame
#'
#' @param df A data frame as returned by \code{formatDataToDirectMA}
#' @param effectCode A character string indicating the underlying effect
#'   estimate. This is used to set the \code{sm} argument of the underlying
#'   analysis functions from the \code{meta} package. Acceptable values are
#'   'RD', 'RR', 'OR', 'HR', 'ASD', 'MD', 'SMD'
#' @param dataType A character string specifying which type of data has been
#'   provided. Must be one of 'treatment difference', 'binary' or 'continuous'
#' @param backtransf A logical indicating whether results should be back
#'   transformed. This is used to set the corresponding \code{backtransf}
#'   argument of the underlying functions from the \code{meta} package. If
#'   \code{backtransf=TRUE} then log odds ratios (or hazard ratios etc) will be
#'   converted to odds ratios on plots and print outs
#' @param method A character string indicating which method is used for pooling
#'   studies. This argument is only applicable if data_type='binary'. Must be one
#'   of 'MH' (Default, Mantel-Haenszel method), 'inverse' or 'Peto'
#' @details This function provides a wrapper around the \code{metagen} or
#'   \code{metabin} functions from the \code{meta} package to one or more
#'   analyses to be carried out from a single data frame. This is most useful
#'   when direct meta-analysis is required to support all pairwise comparisons
#'   in a network meta-analyis or effect estimates are required for multiple
#'   pairs of treatments to perform simple indirect meta-analysis using the
#'   Bucher method
#'
#' @return A data frame
#'
#' @seealso \code{\link{formatDataToDirectMA}}, \code{\link[meta]{metagen}},
#'   \code{\link[meta]{metabin}}
#' @export
doDirectMeta = function(df, effectCode, dataType, backtransf = FALSE, method =
                          'MH', ...) {
  #identify the set of treatment comparisons present in the data
  com = dplyr::distinct(df[,3:4]) %>%
    tidyr::unite(contrast, comparator, treatment, remove = FALSE)
  comparisons = data.frame(contrast = NA, comparator = 0, treatment = 0)
  for (i in 1:nrow(com)) {
    if (!com$contrast[i] %in% comparisons$contrast &&
        !str_reverse(com$contrast[i]) %in% comparisons$contrast) {
      comparisons[i,] = com[i,]
    }
  }
  comparisons = comparisons %>% dplyr::filter(complete.cases(.))
  for (i in 1:nrow(comparisons)) {
    d = dplyr::filter(df, comparator == comparisons$comparator[i],
                      treatment == comparisons$treatment[i])
    inv_comp = dplyr::filter(df, comparator == comparisons$treatment[i],
                             treatment == comparisons$comparator[i])
    if (nrow(inv_comp) > 0) {
      inv_comp$diff = -inv_comp$diff
      inv_comp[, c('comparator', 'treatment')] = inv_comp[, c('treatment', 'comparator')]
      inv_comp[, c('NumberAnalysedComparator', 'NumberAnalysedTreatment')] = inv_comp[, c('NumberAnalysedTreatment', 'NumberAnalysedComparator')]
      inv_comp[, c('ComparatorName', 'TreatmentName')] = inv_comp[, c('TreatmentName', 'ComparatorName')]
      d = dplyr::bind_rows(d, inv_comp)
    }
    if (i == 1) {
      comp = d
    } else {
      comp = dplyr::bind_rows(comp, d)
    }

  }
  comp = comp %>%
    #add a column that defines the comparisons present in the dataset
    #group dataset by treatment comparison
    #nest the data. Produces a data frame with one row per treament comparison
    #"data" column is a list of data.frames. One per treatment comparison
    #map combined with map_meta applies the appropriate meta-analysis function from
    #the meta package according to the data_type argument
    tidyr::unite(comparison, comparator, treatment, remove = FALSE) %>%
    dplyr::group_by(comparison) %>%
    tidyr::nest() %>%
    dplyr::mutate(
      metaRes = purrr::map(
        data, map_meta, data_type = dataType, method =
          method, effectCode = effectCode, backtransf = backtransf
      )
    )
  return(comp$metaRes)
}


map_meta = function(df, data_type, method, effectCode, backtransf, ...) {
  if (data_type == 'continuous') {
    directRes = meta::metacont(
      n.e = df$NumberAnalysedTreatment, mean.e = df$MeanTreatment, sd.e = df$SDTreatment,
      n.c = df$NumberAnalysedComparator, mean.c = df$MeanComparator, sd.c = df$SDComparator,
      sm = effectCode, studlab = df$StudyName, label.e = df$TreatmentName,
      label.c = df$ComparatorName
    )
  }
  if (data_type == 'binary') {
    directRes = meta::metabin(
      event.e = df$NumberEventsTreatment, n.e = df$NumberAnalysedTreatment,
      event.c = df$NumberEventsComparator, n.c = df$NumberAnalysedComparator,
      sm = effectCode, backtransf = backtransf, studlab = df$StudyName,
      label.e = df$TreatmentName, label.c = df$ComparatorName,
      method = ifelse(nrow(df) > 1, method, "Inverse")
    )
  }
  if (data_type == 'treatment difference') {
    directRes = meta::metagen(
      TE = df$diff, seTE = df$std.err,sm = effectCode, backtransf = backtransf,
      studlab = df$StudyName, n.e = df$NumberAnalysedTreatment,
      n.c = df$NumberAnalysedComparator, label.e = df$TreatmentName,
      label.c = df$ComparatorName
    )
  }

  directRes$e.code = df$treatment[1]
  directRes$c.code = df$comparator[1]
  directRes
}

.backtransform = function(df) {
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

#' Extract summary results for a direct meta-analysis
#'
#' @param metaRes An object of class \code{c("metagen", "meta")} or c("metabin",
#'   "meta"). These objects are lists containing the results of direct
#'   meta-analysis. See \code{\link[meta]{metagen}}, \code{\link[meta]{metabin}}
#'   for a description of exactly what is contained in the list
#' @param effect a character string describing the effect estimate, e.g. 'Rate
#'   Ratio', 'Odds Ratio', 'Hazard Ratio'
#' @param intervention A character string describing the name of the
#'   intervention. Defaults to 'Int' if not provided
#' @param comparator A character string describing the name of the comparator.
#'   Defaults to 'Con' if not provided
#' @param backtransf A logical indicating whether results should be back
#'   transformed. If \code{backtransf=TRUE} then log odds ratios (or hazard
#'   ratios etc) will be converted to odds ratios on plots and print outs.
#'
#' @details This function extracts the results of a meta-analysis from a
#'   list-type object produced by the \code{meta} package and returns them as a
#'   data frame. If there is only one study comparing two treatments then no
#'   meta-anlaysis is performed but the results of the primary study are
#'   included in the output. In this case the fields \code{Model},
#'   \code{Tau.sq}, \code{method.tau} and \code{I.sq} in the output will be
#'   blank as these fields have no meaning for a single study. This is useful if
#'   you want to use these results to perform simple indirect (Bucher) analyses.
#'
#'   If more than one study is available then both fixed effect and random
#'   effects results will be returned
#'
#' @return A data frame with the following columns:
#' \itemize{
#'  \item \code{Intervention} The name of the intervention
#'  \item \code{InterventionCode} The ID number of the intervention in the
#'  current set of analyses. NA if not provided
#'  \item \code{Comparator} The name of the comparator
#'  \item \code{ComparatorCode} The ID number of the comparator in the current
#'    set of analyses
#'  \item \code{Effect} The type of effect measure. Takes the
#'    value of the \code{effect} argument
#'  \item \code{Model} The type of model. Fixed effect or Random Effects.
#'    Blank if there is only one study
#'  \item \code{log.TE} The treatment effect on log scale, e.g. log OR
#'  \item \code{seTE.log} The standard error for the log treatment effect
#'  \item \code{log.lower}, \code{log.upper} The upper and lower 95\% confidence
#'    intervals for the log treatment effect
#'  \item \code{z}, \code{p} The z-value and corresponding p-value for the test
#'    of effect
#'  \item \code{level} The level for the confidence intervals. Defaults to 0.95
#'    for a 95\% confidence interval
#'  \item \code{TE}, \code{lower}, \code{upper} The treatment effect with lower
#'    and upper confidence intervals backtransformed to a linear scale
#'  \item \code{Tau.sq} The heterogeneity variance
#'  \item \code{I.sq}, \code{I.sq.lower}, \code{I.sq.upper} The heterogeneity statistic
#'    I-squared with upper and lower confidence intervals.
#' }
#'
#' @seealso \code{\link[meta]{metagen}}, \code{\link[meta]{metabin}}
#' @export
extractDirectRes = function(metaRes, effect, intervention = 'Int',
                            comparator = 'Con', interventionCode = NA,
                            comparatorCode = NA, backtransf = FALSE) {

  #create a summary of the results, extract fixed and random
  res = summary(metaRes)
  fixed = as.data.frame(res$fixed)
  random = as.data.frame(res$random)[1:7]
  if (backtransf == TRUE) {
    #exponentiate the effect estimates if required
    fixed = .backtransform(fixed)
    random = .backtransform(random)
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

  #re-arrange the output columns into a more report friendly order
  #for future versions consider adapting this to use dplyr and select by
  #name instead of position.
  if (effect != 'Mean Difference') {
    df = df[,c(5, 1, 3, 14:16, 22:23, 6, 2, 4, 7:13, 17:21)]
  } else {
    df = dplyr::select(
      df, Effect, Intervention, Comparator, TE, lower, upper,
      n.studies, studies, Model, InterventionCode, ComparatorCode,
      seTE, z, p, level, Tau.sq, method.tau,
      I.sq.TE, I.sq.lower, I.sq.upper
    )
  }

}

#' Draw a forest plot
#'
#' @param meta an object of class c("metagen", "meta"), c("metacont", "meta") or
#'   c("metabin", "meta") as returned by the functions \code{metagen},
#'   \code{metabin} or \code{metacont} in the package meta
#' @param showFixed,showRandom Logical indicating whether fixed effect and/or
#'   random effects results should be shown on the forest plot. By default both
#'   are shown set the appropriate argument to FALSE if you want to exclude that
#'   result from the plot.
#' @param ... additional arguments to be passed to \code{forest}
#'   e.g. \code{col.square='red'}, \code{col.diamond='black'}, \code{smlab='Odds
#'   Ratio'}
#'
#' @details This function provides a very simple wrapper around \code{forest}
#'   from the \code{meta} package. The main purpose is to allow some defaults
#'   to be set for the arguments to \code{forest} and to work out a reasonable
#'   scale for the x-axis automatically.
#'
#'   By default pooled estimates for both fixed and random effects models will
#'   be shown. If there is only one study comparing a given pair of treatments
#'   then the results of that study are shown but no pooled estimates are
#'   displayed.
#'
#'   @return NULL
#'
#' @seealso \code{\link[meta]{forest}}
#' @export
drawForest = function(meta, showFixed = TRUE, showRandom = TRUE, ...) {

  #work out sensible values for the x-axis limits
  limits = c(
    meta$lower, meta$upper, meta$lower.fixed, meta$upper.fixed,
    meta$lower.random, meta$upper.random
  )
  if (!'MD' %in% meta$sm) {
    #if data are not continuous set limits appropriate for
    #ratio measures
    limits = range(exp(limits))
    xlower = ifelse(limits[1] < 0.2, round(limits[1], 1), 0.2)
    #check that the lower limit does not end up as zero
    xlower = ifelse(xlower == 0, 0.01, xlower)
    xupper = ifelse(limits[2] > 5, round(limits[2]), 5)
  } else {
    #if data are continuous then set limits accordingly
    limits = range(limits)
    xlower = ifelse(limits[1] < -2, round(limits[1]), -2)
    xupper = ifelse(limits[2] > 2, round(limits[2]), 2)
  }

  xlimits = c(xlower, xupper)

  #check and set sensible treatment names if available
  if('e.name' %in% names(meta)){
    meta$label.e = meta$e.name
    meta$label.c = meta$c.name
  }

  #don't show the pooled estimate if there is only one study
  if (meta$k == 1) {
    showFixed = FALSE
    showRandom = FALSE
  }

  #forest plot
  meta::forest(
    meta, hetlab = NULL, text.I2 = 'I-sq', text.tau2 = 'tau-sq', xlim = xlimits,
    comb.fixed = showFixed, comb.random = showRandom, lty.fixed = 0,
    lty.random = 0, just.studlab = 'right', fontsize = 10, ...
  )
}
