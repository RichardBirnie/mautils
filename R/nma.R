#' Run a mixed treatment comparison (MTC, a.k.a network meta-analysis)
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
#' @param file A character string specifying the directory where the results
#'   will be saved. This function will create a subdirectory 'Results/MTC'
#'   within this directory will be one or more further directories
#'   'Results/MTC/FixedEffects' and 'Results/MTC/RamdomEffects' depending on
#'   the values of \code{doFixed} and \code{doRandom}. These are only created
#'   if they are required and do not exist already. Each results directory
#'   will also contain a 'Figures' subdirectory.
#' @param data_type A character string specifying which type of data has been
#'   provided. Currently only 'treatment difference' or 'binary' are supported
#' @param treatmentID A data frame with columns 'description' defining the
#'   treatment names and 'id' defining the treatment ID numbers. This is used
#'   to set the \code{treatments} argument of the \code{mtc.network} function
#'   from \code{gemtc}. An optional third column 'Order' may also be provided.
#'   If present this should be a column of integers indicating the order in
#'   which the interventions should be presented in output tables and figures
#' @param effect_measure A character string indicating what type of effect
#'   measure is used, e.g. 'Rate Ratio', 'Odds Ratio' etc. This is used as a
#'   label in forest plots so keep it short.
#' @param toi A named vector specifying the treatments of interest. The names
#'   should be the treatment names. The values should be the ID number of the
#'   treatments in the network.
#' @param analysis_set A character string indicating which analysis set this
#'   is. This is useful if there are multiple different sets of comparators in
#'   a single project. If nothing is provided this will be set to the string
#'   'Default'. This is used to set the \code{description} argument of the
#'   underlying \code{mtc.network} function from the \code{gemtc} package
#' @param outcome A character string indicating which outcome is being analysed.
#' @param analysis_case A character string indicating what the current analysis
#'   case is e.g. 'Base Case'
#' @param doFixed,doRandom Logical indicating which model is required. Fixed
#'   effect, random effect or both. The default is to run both. Set the
#'   appropriate argument to FALSE if you do not want to run that analysis.
#' @param model_type A character string describing the type of model. This is
#'   passed to the \code{type} argument of the underlying \code{mtc.model}
#'   function from \code{gemtc}. The underlying function accepts the values
#'   "consistency", "nodesplit", "ume", or "use". Only consistency has been
#'   tested in the context of this function.
#' @param max_val Used to set \code{om.scale} and create vague priors for the
#'   between trials heterogeneity. For the log odds-ratio, values between 2
#'   and 5 are considered reasonable. Default = 5.
#' @param prior A function call to \code{mtc.hy.prior} which specifies the type and
#'   parameters of the prior distribution for the between trials
#'   heterogeneity. This is passed directly to the \code{hy.prior} argument of
#'   the underlying \code{mtc.model} function from the \code{gemtc} package.
#'   The default is a uniform distribution on the between trials standard
#'   deviation with a range from 0-5. This is adequate for log odds ratios,
#'   log hazard ratios etc but should be changed if you are working with
#'   continuous outcomes such as mean difference. For more details see
#'   \code{\link[gemtc]{mtc.model}} and \code{\link[gemtc]{mtc.hy.prior}}
#' @param burn_in The number of iterations for the MCMC burn in period
#'   (Default=10000)
#' @param iterations The number of iterations from the MCMC simulation to be
#'   retained to estimate the results (Default=20000)
#' @param save_convergence Logical indicating whether run and keep plots to
#'   assess the convergence of the MCMC chains. If TRUE (Default) then this
#'   plots will be produced for the trace of each chain, autocorrelation
#'   within each chain and the posterior distribution of each model parameter.
#' @param back_calc A logical indicating whether results should be back transformed.
#'   If set to TRUE then log odds ratios (or hazard ratios etc) will
#'   be converted to odds ratios on plots and print outs. Default is FALSE
#' @param includes_placebo Logical indicating whether the network includes
#'   placebo. This requires slightly different handling internally as placebo
#'   is always the comparator and never the intervention. Default is FALSE. If
#'   this is set to TRUE then it is essential to also set the
#'   \code{placebo_code} argument appropriately.
#' @param placebo_code An integer specifying the ID number of the placebo
#'   treatment in the network if a placebo is included. Often placebo is given
#'   the number 1 but not always.
#' @param report_order A character string indicating whether the treatments
#'   should be reported in the order the order they were provided (default) or
#'   if a custom order is required. Acceptable values are 'default' or
#'   'custom'. If this is set to 'custom' then \code{treatmentID} must contain a
#'   column named 'Order', see above.
#' @param check_inconsistency A logical indicating whether checks for
#'   inconsistency should be performed. If set to \code{TRUE} (default) then the
#'   \code{mtc.nodesplit} function from gemtc is used to run all possible
#'   nodesplitting models for the current network.
#'
#' @details This function provides an interface to run MTC analyses using the
#'   gemtc package. Although the function has a large number of arguments the
#'   default values should produce reasonable results in the majority of cases.
#'   The aim is to abstract away as much of the technical detail as possible so
#'   that reasonable values are chosen based on the type of data provided.
#'
#'   The output of the analysis is automatically rearranged and saved in a
#'   'report ready' form to save time in downstream processing of the results
#'
#'
#' @seealso \code{\link[gemtc]{mtc.network}}, \code{\link[gemtc]{mtc.model}},
#'   \code{\link[gemtc]{mtc.run}}, \code{\link[gemtc]{mtc.run}},
#'   \code{\link{extractMTCResults}}, \code{\link{plotEstimates}}

runMTC = function(df, file, data_type, treatmentID, effect_measure, toi,
                  analysis_set = 'Default', outcome, analysis_case,
                  doFixed = TRUE, doRandom = TRUE, model_type = 'consistency',
                  max_val = 5, prior = mtc.hy.prior("std.dev", "dunif", 0, max_val),
                  burn_in = 10000, iterations = 20000, save_convergence = TRUE,
                  back_calc = FALSE, includes_placebo = FALSE, placebo_code,
                  report_order = 'default', check_inconsistency = TRUE) {
  message('Start MTC')
  #create an mtc.network object
  #note different data argument for treatment differences versus per arm data
  #types
  if (data_type == 'treatment difference') {
    network = gemtc::mtc.network(
      data.re = df, treatments = treatmentID[, 1:2], description = analysis_set
    )
    #specify likelihood and link appropriate to the data type
    likelihood = 'normal'
    link = 'identity'
  }
  if (data_type == 'binary') {
    network = gemtc::mtc.network(
      data.ab = df, treatments = treatmentID[, 1:2], description = analysis_set
    )
    #specify likelihood and link appropriate to the data type
    likelihood = 'binom'
    link = 'logit'
  }

  #keep the input data to avoid overwriting
  inputData = df

  #set which models are required
  EffectsModel = c('fixed', 'random')[c(doFixed, doRandom)]

  for (i in 1:length(EffectsModel)) {
    message('Run MTC: ', EffectsModel[i], ' effects')
    #set up folders for the results and figures. No need to edit this
    f = paste0(outcome, '_', analysis_case, '.xlsx')
    MTCresultsFile = file.path(baseFile, 'Results', 'MTC',
                               paste0(capwords(EffectsModel[i]), 'Effects'), f)
    MTCresDir = dirname(MTCresultsFile)
    if (!dir.exists(MTCresDir)) {
      dir.create(MTCresDir, recursive = TRUE)
    }
    MTCfigDir = file.path(MTCresDir, 'Figures')
    if (!file.exists(MTCfigDir)) {
      dir.create(MTCfigDir, recursive = TRUE)
    }
    #save the treatment codes and the input data with the results file
    rbutils::saveXLSX(
      as.data.frame(treatmentID), file = MTCresultsFile, sheetName = 'Code',
      showNA = FALSE, row.names = FALSE, append = TRUE
    )
    rbutils::saveXLSX(
      as.data.frame(inputData), file = MTCresultsFile, sheetName = 'Data',
      showNA = FALSE, row.names = FALSE, append = TRUE
    )

    #create the model object. See ?mtc.model for details
    message('Run MCMC')
    if (EffectsModel[i] == 'random' && doRandom == TRUE) {
      #Random effects model
      #hy.prior = the prior for the between trials standard deviation
      #set as a uniform distribution. Upper limit is set from om.scale
      model = suppressWarnings(
        gemtc::mtc.model(
          network, type = model_type, n.chain = 3, likelihood = likelihood,
          link = link, linearModel = 'random',
          om.scale = max_val, hy.prior = prior
        )
      )
    } else {
      #Fixed effects model
      model = suppressWarnings(
        gemtc::mtc.model(
          network, type = model_type, n.chain = 3, likelihood = likelihood,
          link = link, linearModel = 'fixed',
          om.scale = max_val
        )
      )
    }
    #plot(model, layout=igraph::layout.fruchterman.reingold)

    #run the MCMC simulation
    mtcResults = gemtc::mtc.run(model, n.adapt = burn_in, n.iter = iterations,
                                thin = 1)

    #Save the JAGS (or BUGS) code that was generated for the model
    #set up a file to save the model. If necessary create the directory
    message('Save model file')
    m = paste0(
      mtcResults$model$linearModel, '_effects_', mtcResults$model$type,'_',
      '_model', '.txt'
    )
    modelFile = file.path(MTCresDir, m)
    saveModelCode(mtcResults, modelFile)

    if (save_convergence) {
      message('Save convergence plots')
      diagDir = file.path(MTCresDir, 'ConvergenceDiagnostics')
      saveDiagnostics(mtcResults, directory = diagDir)
      rm(diagDir)
    }

    #extract the results for the treatment comparisons
    #this function also saves the results to the excel file specified by resultsFile
    message('Extract and save results')
    pairwiseResults = extractMTCResults(
      mtcResults, resultsFile = MTCresultsFile, expon = back_calc,
      includesPlacebo = includes_placebo, coding = treatmentID,
      reportOrder = report_order
    )

    #slice out results for treatments of interest and save these as separate
    #sheets for convenience everything vs placebo ALWAYS APPEND=TRUE OR YOU WILL
    #OVERWRITE THE EXISTING RESULTS
    message('Drawing plots')
    if (includes_placebo) {
      placebo = extractTOI(
        pairwiseResults, toi = placebo_code, treatments = treatmentID,
        intervention = FALSE, reportOrder = report_order
      )
      saveXLSX(
        as.data.frame(placebo), file = MTCresultsFile, sheetName = 'Placebo',
        showNA = FALSE, row.names = FALSE, append = TRUE
      )

      #make a plot of all treatments vs placebo
      p = plotEstimates(
        df = placebo, yvar = 'nameB', xlabel = effect_measure,
        reportOrder = report_order
      )
      #save the plot
      figFile = file.path(MTCfigDir, 'AllVsPlacebo.jpg')
      jpeg(
        file = figFile, width = 22, height = 15, units = 'cm', res = 300,
        quality = 100
      )
      suppressWarnings(print(p))
      graphics.off()
      rm(placebo)
    }

    #slice out results for treatments of interest and save these as separate
    #sheets for convenience
    for (j in 1:length(toi)) {
      n = names(toi)[j]
      tr = extractTOI(
        pairwiseResults, toi = toi[j], treatments = treatmentID,
        intervention = TRUE, reportOrder = report_order
      )
      #ALWAYS APPEND=TRUE OR YOU WILL OVERWRITE THE EXISTING RESULTS
      saveXLSX(
        as.data.frame(tr), file = MTCresultsFile, sheetName = n, showNA = FALSE,
        row.names = FALSE, append = TRUE
      )

      #make a plot for each treatment of interest vs all other treatments
      p = plotEstimates(
        df = tr, yvar = 'nameA', xlabel = effect_measure,
        reportOrder = report_order
      )
      #save the plot
      f = paste0(outcome, ' ', n, 'VsAll.jpg')
      figFile = file.path(MTCfigDir, f)
      jpeg(
        file = figFile, width = 22, height = 15, units = 'cm', res = 300,
        quality = 100
      )
      suppressWarnings(print(p))
      graphics.off()
    }
    rm(tr)

    #check for presence of inconsistency in the network using node splitting
    #node splitting function provided by gemtc
    if(check_inconsistency){
      splitComp = gemtc::mtc.nodesplit.comparisons(network)
      if(nrow(splitComp) > 0) {
        splitRes = gemtc::mtc.nodesplit(
          network = network, comparisons = splitComp,
          n.chain = 3, likelihood = likelihood, link = link,
          linearModel = EffectsModel[i], hy.prior = prior, om.scale = max_val,
          n.adapt = burn_in, n.iter = iterations
        )

        #format the node splitting results as a more useful data frame
        #and save the output
        nsRes = extractNodesplit(splitRes, treatments = treatmentID[,1:2], backtransf = back_calc)
        incDir = file.path(MTCresDir, 'Inconsistency')
        if (!dir.exists(incDir)) {
          dir.create(incDir, showWarnings = FALSE, recursive = TRUE)
        }
        incFile = file.path(incDir, paste0(outcome, '_', analysis_case, '_inconsistency.xlsx'))
        rbutils::saveXLSX(
          as.data.frame(nsRes), file = incFile, sheetName = 'Inconsistency', showNA = FALSE,
          row.names = FALSE, append = TRUE
        )

        #plot the ratio of effect estimates to visualise inconsistency
        df = as.data.frame(nsRes)
        #create a column to use for the y-axis
        df$comparison = paste0(df$Intervention, '\n', df$Comparator)
        p = plotEstimates(
          df, yvar = 'comparison', xvar = 'RoR', lowLimit = 'RoR.lower.ci',
          hiLimit = 'RoR.upper.ci', xlabel = 'RoR'
        )
        f = paste0(outcome, '_', analysis_case, '_inconsistency.jpg')
        incfig = file.path(incDir, f)
        jpeg(
          file = incfig, width = 22, height = 15, units = 'cm', res = 300,
          quality = 100
        )
        suppressWarnings(print(p))
        graphics.off()
      }
    }
  }
}

#' Extract treatment comparison information from mtc summary output
#'
#' @param df A data frame as returned by \code{calcAllPairs}
#'
#' @details In the summary output from a network meta-analysis treatment
#'   comparisons are commonly labelled as 'd.1.2' for the relative effect of
#'   treatment 2 compared to treatment 1. This function removes the 'd' and
#'   splits '1' and '2' into separate columns which is more useful for
#'   downstream reporting. This function is intended to work on the output from
#'   \code{calcAllPairs}. This functionality will usually be accessed via
#'   extractMTCResults and should only be used directly if you understand what
#'   you are doing.
#'
#' @return A data frame with additional columns
#'
#' @seealso \code{\link{calcAllPairs}}, \code{\link{extractMTCResults}}
extractComparison = function(df) {
  #extract the information of which comparison is on each row from the
  #first column of the winBUGS output and store it in individual columns
  treatments = as.data.frame(stringr::str_extract_all(df$node, pattern = '\\d+',
                                             simplify = TRUE))
  treatments[,1] = as.integer(as.character(treatments[,1]))
  treatments[,2] = as.integer(as.character(treatments[,2]))
  colnames(treatments) = c('tA', 'tB')
  df = dplyr::tbl_df(dplyr::bind_cols(treatments, df[,2:ncol(df)]))
}

#' Match treatment names to ID numbers
#'
#' @param results A data frame as returned by \code{extractComparison}
#' @param coding A data frame with two columns 'id' and 'description'. 'id' must
#'   be the treatment id numbers corresponding to the way the treatments were
#'   coded in the network. 'description' should be the name of the treatment.
#'
#' @details This function matches the coded treatment id numbers in the network
#'   to the corresponding human readable names. The mapping from id number to
#'   name should be provided as a two column data frame via the \code{coding}
#'   argument. This function is intended to work on the data frame generated as
#'   the output from \code{extractComparison}. This function will mainly be
#'   called via \code{extractMTCResults} and should only be used directly if you
#'   understand what you are doing. The general work flow is \code{mtc.run} >
#'   \code{calcAllPairs} > \code{extractComparison} > \code{nameTreatments} >
#'   \code{makeTab}. \code{extractMTCResults} will handle the last four steps
#'   for you.
#'
#' @return The same data frame with the treatment names appended
#'
#' @seealso \code{\link{extractComparison}}, \code{\link{calcAllPairs}}, \code{\link{extractMTCResults}}
nameTreatments = function(results, coding, ...) {

  #drop order column from coding. Makes no sense here
  coding = coding[,1:2]

  #match intervention codes to their names
  results = dplyr::right_join(coding, results, by = c('id' = 'tB'))
  colnames(results)[1:2] = c('TreatmentB', 'nameB')

  #match comparator codes to their names
  results = dplyr::right_join(coding, results, by = c('id' = 'tA'))
  colnames(results)[1:2] = c('TreatmentA', 'nameA')
  return(results)
}

#' Generate a square grid reporting all pairwise comparisons of treatments in
#' the network
#'
#' @param results A data frame as returned by \code{nameTreatments}
#' @param coding  A data frame with two columns 'id' and 'description'. 'id'
#'   must be the treatment id numbers corresponding to the way the treatments
#'   were coded in the network. 'description' should be the name of the
#'   treatment. An optional third column named 'Order' may also be provided. If
#'   present this controls the order in which treatments are presented in the
#'   output. The values should be a sequence of numbers indicating the order in
#'   which the treatments should be sorted. If the 'Order' column is present
#'   then \code{reportOrder}  should be set to 'custom'
#' @param rounding The number of decimal places to be included in the output.
#'   The default is two.
#' @param reportOrder A character string indicating whether the treatments
#'   should be reported in the order the order they were provided (default) or
#'   if a custom order is required. Acceptable values are 'default' or 'custom'.
#'   If this is set to 'custom' then \code{coding} must contain a column named
#'   'Order'. See above.
#' @param ... Arguments passed from other functions, primarily
#'   \code{extractMTCResults}
#'
#' @details This function takes a data frame as returned by
#'   \code{nameTreatments} and returns a square grid containing all pairwise
#'   treatment comparisons in the network. The results are presented for column
#'   versus row. This is a common approach to reporting MTC results. This
#'   function should mainly be used via \code{extractMTCResults}
#'
#' @return A data frame
#'
#' @seealso \code{\link{extractMTCResults}}, \code{\link{nameTreatments}}
makeTab = function(results, coding, rounding = 2, reportOrder = 'default', ...) {
  #create summary string
  results[,'effect'] = paste0(
    round(results$median, rounding), ' (',
    round(results$CrI_lower,rounding),' to ',
    round(results$CrI_upper,rounding), ')'
  )

  #The order of treatments in the table can be controlled manually
  #using the argument reportOrder.
  #The default is to order the table in the same order the treatments are coded
  #in the network
  if (reportOrder == 'custom') {
    coding = dplyr::arrange(coding, Order)
  } else {
    #order by treatment coding in the network
    coding$Order = coding$id
  }

  #Create a simple data frame to add output
  reportTab = data.frame('Treatment' = coding$id)

  for (i in 1:nrow(coding)) {
    comparator = coding$id[i]
    comp = dplyr::filter(results, TreatmentA == comparator)
    comp = dplyr::select(comp, TreatmentB, effect)
    colnames(comp)[ncol(comp)] = coding$id[i]
    reportTab = dplyr::full_join(reportTab, comp,
                          by = c('Treatment' = 'TreatmentB')
                          )
  }
  rownames(reportTab) = coding$description
  reportTab = dplyr::select(reportTab,-Treatment)
  colnames(reportTab) = coding$description

  reportTab = as.data.frame(t(reportTab))
}

#internal function used to access the S3 summary() method for mtc.result objects
.summariseMTC = getS3method('summary', 'mtc.result')

#' Calculate all pairwise treatment effects in a network meta-analysis
#'
#' @param mtcRes An object of class \code{mtc.results} as returned by mtc.run
#' @param expon Logical indicating whether the results should be exponentiated
#'   or not. If your analysis is based on log hazard ratio or log odds ratio etc
#'   then set this to TRUE to get the results as the corresponding hazard ratio
#'   or odds ratio.
#' @param ... Additional arguments passed down from extractMTCResults
#'
#' @details This function takes a \code{mtc.results} object and extracts summary
#'   statistics for all pairwise treatment comparisons in the network using
#'   \code{relative.effect} from gemtc. This function should mainly be used via
#'   \code{extractMTCresults} but can be called directly on any
#'   \code{mtc.results} object. Only results for the treatment comparisons are
#'   returned, additional nodes such as hetereogeneity parameters are not
#'   returned. Model comparisons statistics and heterogeneity parameters can be
#'   obtained using \code{extractModelFit}
#'
#' @return A data frame reporting summary statistics (mean, median etc) for the
#'   posterior distribution of each treatment comparison. Each row in the data
#'   frame represents one treatment comparison.
#'
#' @seealso \code{\link[gemtc]{mtc.run}}, \code{\link[gemtc]{relative.effect}}
#'   \code{\link{extractMTCResults}}, \code{\link{extractModelFit}}
calcAllPairs = function(mtcRes, expon = FALSE, ...) {
  tid = as.integer(mtcRes$model$network$treatments$id)
  for (t in 1:length(tid)) {
    re = suppressWarnings(.summariseMTC(
      gemtc::relative.effect(mtcRes, t1 = tid[t], preserve.extra = FALSE)
    ))

    stats = re$summaries$statistics
    stats = data.frame(node = rownames(stats), stats, row.names = NULL)
    quan = re$summaries$quantiles
    quan = data.frame(node = rownames(quan), quan[,c(1,3,5)],
                      row.names = NULL)
    colnames(quan) = c('node','CrI_lower', 'median', 'CrI_upper')
    out = dplyr::full_join(stats, quan, by = 'node')

    if (t == 1) {
      output = out
    } else {
      output = suppressWarnings(bind_rows(output, out))
    }
  }

  if (expon == TRUE) {
    output[,c(2,6:8)] = exp(output[,c(2,6:8)])
    output$SD = (output$CrI_upper - output$CrI_lower) / 3.92
  }
  output$Sample = nrow(mtcRes$samples[[1]])
  output = output[,c(1,7,6,8,2:5,9)]
}

#' Extract model comparison statistics from a \code{mtc.result} object
#'
#' @param mtcRes an object of class \code{mtc.result} as returned by
#'   \code{mtc.run} from gemtc
#'
#' @details This function returns model comparison statistics (DIC, pD, Dbar)
#'   from a \code{mtc.result} object. If the results object is derived from a
#'   random effects model then the between trials standard deviation will also
#'   be provided as sd.d. This function should mainly be used via
#'   \code{extractMTCResults} but can be used directly on a \code{mtc.result}
#'   object.
#'
#' @return A data frame containing the model comparison statistics.
#'
#' @seealso \code{\link{extractMTCResults}}, \code{\link[gemtc]{mtc.run}}
extractModelFit = function(mtcRes) {
  modelSummary = .summariseMTC(mtcRes)
  dic = data.frame(
    'Mean' = modelSummary$DIC,
    'SD' = NA,
    row.names = c('Dbar', 'pD', 'DIC')
  )
  if (mtcRes$model$linearModel == 'random') {
    modelSD = modelSummary$summaries$statistics
    ix = rownames(modelSD) == 'sd.d'
    msd = data.frame(
      'Mean' = modelSD[ix,1],
      'SD' = modelSD[ix,2],
      row.names = rownames(modelSD)[ix]
    )
    dic = rbind(msd, dic)
  }
  return(dic)
}

#' Extract the results of a network meta-analysis from a \code{mtc.results} object
#'
#' @param res A \code{mtc.results} object as returned by \code{mtc.run}
#' @param resultsFile A character string indicating the path to the excel file
#'   where the results should be saved. This is passed on directly to
#'   \code{\link[rbutils]{saveXLSX}}
#' @param includesPlacebo Logical indicating whether the network includes
#'   placebo. This requires slightly different handling internally as placebo is
#'   always the comparator and never the intervention
#' @param ... Additional arguments passed to underlying extractor functions,
#'   particularly calcAllPairs, nameTreatments and makeTab
#'
#' @details This function takes an object of class \code{mtc.result} and
#'   extracts the results for all pairwise treatment comparisons in the network.
#'   The results are returned as a data frame. This function is a wrapper that
#'   calls a series of other functions in the right order to do the actual work.
#'   \itemize{
#'    \item calcAllPairs calculates all pairwise treament comparisons
#'    in the network using \code{relative.effect} from the gemtc package
#'    \item extractComparison separates out the which pair of treatments
#'    are being compared
#'    \item nameTreatments matches the numbered treatment ID used by mtc.run to
#'    human readable treatment names
#'    \item makeTab produces a square grid of all pairwise treatment comparisons
#'    in the network. The results are presented for column versus row. This is a
#'    common approach to reporting MTC results.
#'    \item extractModelFit Extracts and returns model comparison statistics such
#'    as DIC, pD and Dbar
#'   }
#'
#'   The underlying functions can be used directly provided this order is
#'   maintained however this should not usually be necessary. If you take this
#'   approach you will also need to save all the respective outputs
#'   appropriately.
#'
#'   The \code{...} can be used to pass arguments to the underlying functions.
#'   The following arguments in particular will be required in almost all cases:
#'   \itemize{
#'   \item \code{expon} Logical indicating whether the results should be exponentiated
#'   or not. If your analysis is based on log hazard ratio or log odds ratio etc
#'   then set this to TRUE to get the results as the corresponding hazard ratio
#'   or odds ratio. Passed to \code{calcAllPairs}
#'   \item \code{coding} A data frame with two columns 'id' and 'description'. 'id'
#'   must be the treatment id numbers corresponding to the way the treatments
#'   were coded in the network. 'description' should be the name of the
#'   treatment. An optional third column named 'Order' may also be provided. If
#'   present this controls the order in which treatments are presented in the
#'   output. The values should be a sequence of numbers indicating the order in
#'   which the treatments should be sorted. If the 'Order' column is present
#'   then \code{reportOrder}  should be set to 'custom'. Passed to \code{makeTab}
#'   \item \code{reportOrder} A character string indicating whether the treatments
#'   should be reported in the order the order they were provided (default) or
#'   if a custom order is required. Acceptable values are 'default' or 'custom'.
#'   If this is set to 'custom' then \code{coding} must contain a column named
#'   'Order'. See above. Passed to \code{makeTab}
#'   }
#'
#' @seealso \code{\link[gemtc]{mtc.run}}, \code{\link[rbutils]{saveXLSX}},
#'   \code{\link{calcAllPairs}}, \code{\link{extractComparison}},
#'   \code{\link{nameTreatments}}, \code{\link{makeTab}},
#'   \code{\link{extractModelFit}}
#'
extractMTCResults = function(res, resultsFile, includesPlacebo = FALSE, ...) {
  #calculate all pairwise effects
  pairwiseResults = calcAllPairs(res, ...)
  rbutils::saveXLSX(
    as.data.frame(pairwiseResults), file = resultsFile, sheetName = 'Raw',
    showNA = FALSE, row.names = FALSE, append = TRUE
  )

  #extract the comparison info
  #map treatment names to numbers
  pairwiseResults = extractComparison(pairwiseResults)
  pairwiseResults = nameTreatments(pairwiseResults, ...)
  rbutils::saveXLSX(
    as.data.frame(pairwiseResults), file = resultsFile,
    sheetName = 'ProcessedAll', showNA = FALSE,
    row.names = FALSE, append = TRUE
  )

  #make a table of all pairwise comparisons and save it as an excel file
  #ALWAYS APPEND=TRUE OR YOU WILL OVERWRITE THE EXISTING RESULTS
  reportTab = makeTab(results = pairwiseResults, ...)
  if (includesPlacebo) {
    reportTab = dplyr::select(reportTab,-matches('Placebo')) #drop the placebo column
  }
  rt = as.data.frame(reportTab, check.names = FALSE)
  rbutils::saveXLSX(
    rt, file = resultsFile, sheetName = 'Report', row.names = TRUE,
    showNA = FALSE, append = TRUE
  )

  #extract and save the model fit information
  modelFit = extractModelFit(res)
  rbutils::saveXLSX(
    modelFit, file = resultsFile, sheetName = 'DIC', showNA = FALSE,
    row.names = TRUE, append = TRUE
  )

  return(pairwiseResults)
}

#' Extract all treatment comparisons involving a specific treatment if interest
#'
#' @param df A data frame as produced by \code{extractMTCResults}
#' @param toi An integer specifying the id number of the treatment of interest
#'   in the network
#' @param treatments A data frame with two columns 'id' and 'description'. 'id'
#'   must be the treatment id numbers corresponding to the way the treatments
#'   were coded in the network. 'description' should be the name of the
#'   treatment. An optional third column named 'Order' may also be provided. If
#'   present this controls the order in which treatments are presented in the
#'   output. The values should be a sequence of numbers indicating the order in
#'   which the treatments should be sorted. If the 'Order' column is present
#'   then set \code{orderResults=TRUE}
#' @param intervention Logical indicating whether the treatment of interest is
#'   an intervention or a comparator (e.g. placebo). This controls which column
#'   of \code{df} is searched to identify relevant treatment comparisons
#' @param reportOrder A character string indicating whether the treatments
#'   should be reported in the order the order they were provided (default) or
#'   if a custom order is required. Acceptable values are 'default' or
#'   'custom'. If this is set to 'custom' then \code{treatments} must contain a
#'   column named 'Order', see above.
#'
#' @details This is a simple filtering function to return all treatment
#'   comparisons involving a specific treatment of interest. This function is
#'   designed to be used with the output of \code{extractMTCResults}. There are
#'   two intended use cases:
#'   \enumerate{
#'   \item To return results for a given intervention of interest compared to
#'   all other treatments in the network
#'   \item To return results for all interventions compared to a given control
#'   treatment, e.g. placebo.
#'   }
#'
#' @return A data frame containing all treatment comparisons involving the
#'   specified treatment of interest.
#'
#' @seealso \code{\link{extractMTCResults}}
extractTOI = function(df, treatments, toi, intervention = TRUE,
                      reportOrder = 'default') {

  #extract the treatment of interest
  if (intervention == TRUE) {
    df = dplyr::filter(df, TreatmentB == toi)
  } else {
    df = dplyr::filter(df, TreatmentA == toi)
  }

  #sort if required
  if (intervention == TRUE && reportOrder == 'custom') {
    df = dplyr::left_join(df, treatments[,2:3], by = c('nameA' = 'description'))
    df = dplyr::arrange(df, Order)
    df$Order = 1:nrow(df)
  }
  if (intervention == FALSE && reportOrder == 'custom') {
    df = dplyr::left_join(df, treatments[,2:3], by = c('nameB' = 'description'))
    df = dplyr::arrange(df, Order)
    df$Order = 1:nrow(df)
  }

  return(df)
}

#' Extract node splitting results
#'
#' @param ns.res An object of class \code{mtc.nodesplit}
#' @param treatments A data frame with columns 'description' defining the
#'   treatment names and 'id' defining the treatment ID numbers.
#' @param back_calc A logical indicating whether results should be back transformed.
#'   If set to TRUE then log odds ratios (or hazard ratios etc) will
#'   be converted to odds ratios (or hazard ratios etc).
#'
#' @details This function takes the output from the function mtc.nodesplit and
#'   returns a data frame of the results for each treatment comparison that can
#'   be split. The results are presented as one row per treatment comparison.
#'
#'
#' @return A data frame. Column headers starting 'dir.' are the treatment effect
#'   and corresponding 95\% confidence intervals for the direct comparison. Column
#'   headers starting 'ind.' report the same information for the indirect
#'   comparison. Column headers starting 'cons.' report the pooled estimate
#'   using both direct and indirect evidence under a standard consistency model.
#'   Columns labelled RoR report the ratio of direct:indirect effect estimates
#'   and the corresponding 95\% confidence interval
#'
#' @seealso \code{\link[gemtc]{mtc.nodesplit}}
extractNodesplit = function(ns.res, treatments, backtransf) {
  #get summary and retrieve direct, indirect and consistency effect estimates
  ns.res = summary(ns.res)
  splitDir = ns.res$dir.effect
  colnames(splitDir)[3:5] = paste0('dir.', colnames(splitDir)[3:5])
  splitInd = ns.res$ind.effect
  colnames(splitInd)[3:5] = paste0('ind.', colnames(splitInd)[3:5])
  splitCons = ns.res$cons.effect
  colnames(splitCons)[3:5] = paste0('cons.', colnames(splitCons)[3:5])

  #put all the results together in a sensible table
  nsRes = dplyr::left_join(splitDir, splitInd, by = c("t1", "t2")) %>%
    dplyr::left_join(splitCons, by = c("t1", "t2")) %>%
    dplyr::left_join(ns.res$p.value, by = c("t1", "t2"))
  #calculate the inconsistency factor
  #difference of log odds ratios (or rate ratio, hazard ratio etc)
  nsRes$RoR = nsRes$dir.pe - nsRes$ind.pe
  #calculate variance of inconsistency factor
  #sum of Var(direct) + Var(indirect)
  varIF = ((nsRes$dir.ci.u - nsRes$dir.ci.l) / 3.92) ^ 2 +
    ((nsRes$ind.ci.u - nsRes$ind.ci.l) / 3.92) ^ 2
  #Std Error of inconsistency factor
  seIF = sqrt(varIF)
  #CI for inconsistency factor
  nsRes$RoR.lower.ci = nsRes$RoR - (1.96 * seIF)
  nsRes$RoR.upper.ci = nsRes$RoR + (1.96 * seIF)

  #rearrange column order
  nsRes = dplyr::select(nsRes, 1:11, 13:ncol(nsRes), 12)
  nsRes$t1 = as.numeric(nsRes$t1)
  nsRes$t2 = as.numeric(nsRes$t2)

  #convert log values to linear scale if required
  if(backtransf) {
    nsRes[,3:14] = exp(nsRes[,3:14])
  }
  nsRes[,3:15] = round(nsRes[,3:15], 3)
  #add treatment names
  nsRes = dplyr::left_join(nsRes, treatments, by = c('t1' = 'id')) %>%
    dplyr::left_join(treatments[,1:2], by = c('t2' = 'id'))

  colnames(nsRes)[16:17] = c('Comparator', 'Intervention')
  #Rearrange the column order and return the result
  dplyr::select(
    nsRes, 2:1, 17:16, starts_with('dir'),
    starts_with('ind'), starts_with('cons'),
    contains('RoR'), p
  )
}

#' Extract JAGS model code from a \code{mtc.results} object
#'
#' @param mtcRes An object of class \code{mtc.result} as returned by mtc.run
#'   from the gemtc package
#' @param modelFile A character string specify the file path where the code will
#'   be saved.
#'
#' @details The \code{mtc.run} function automatically generates JAGS code to
#'  specify a network meta-analysis model. This function extracts the JAGS code
#'  from the \code{mtc.results} object and saves it to a file to have a record of what
#'  was done.
#'
#' @seealso \code{\link[gemtc]{mtc.run}}
saveModelCode = function(mtcRes, modelFile) {
  #create some basic descriptive text from the model object
  des = paste0('#Description: ', mtcRes$model$network$description, '\n')
  cat(des, file = modelFile)
  modelType = paste0('#Model Type: ', mtcRes$model$linearModel, ' effects\n')
  cat(modelType, file = modelFile, append = TRUE)
  consistency = paste0('#Consistency Assumption: ', mtcRes$model$type, '\n')
  cat(consistency, file = modelFile, append = TRUE)
  like = paste0('#Likelihood: ', mtcRes$model$likelihood, '\n')
  cat(like, file = modelFile, append = TRUE)
  link = paste0('#Link: ', mtcRes$model$link, '\n')
  cat(link, file = modelFile, append = TRUE)
  nchains = paste0('#Number of chains: ', mtcRes$model$n.chain,'\n')
  cat(nchains, file = modelFile, append = TRUE)

  #lastly save the actual model code
  cat(mtcRes$model$code, file = modelFile, append = TRUE)
}

#' Save convergence diagnostics
#' @param mtc An object of class \code{mtc.result} as returned by \code{mtc.run}
#'   in the \code{gemtc} package
#' @param directory A character string indicating the directory to save the
#'   results. By default a directory named ConvergenceDiagnostics will be
#'   created as a subfolder of the current directory
#' @details Create pdf files of some basic diagnostic plots to check the
#'   convergence of MCMC models. Plots are produced by the ggmcmc package for
#'   the trace of the MCMC chain(s), the posterior density of the model
#'   parameters and the autocorrelation within chains
#'
#' @seealso \code{\link[ggmcmc]{ggmcmc}}
saveDiagnostics = function(mtc, directory='./ConvergenceDiagnostics') {

  #if there is no directory to save the plot files then create one
  if (!dir.exists(directory)) {
    dir.create(directory, recursive = TRUE)
  }
  diagData = ggmcmc::ggs(as.mcmc.list(mtc)) #convert data to ggplot friendly format
  f = file.path(directory, 'Traceplot.pdf') #set file name and make traceplots
  ggmcmc::ggmcmc(
    diagData, file = f, plot = 'traceplot', param_page = 3,
    simplify_traceplot = 0.25
  )
  #set file name and make density plots
  f = file.path(directory, 'Density.pdf')
  ggmcmc::ggmcmc(diagData, file = f, plot = 'density', param_page = 3)
  #set file name and make autocorrelation plots
  f = file.path(directory, 'Autocorrelation.pdf')
  ggmcmc::ggmcmc(diagData, file = f, plot = 'autocorrelation', param_page = 3)
  #set file name and make gelman plots
  #this uses the plot from the coda package. The version in ggmcmc is not
  #very useful
  message('Saving BGR')
  pdf(file=file.path(directory, 'BGR.pdf'), paper='a4')
  coda::gelman.plot(coda::as.mcmc.list(mtc), ask=FALSE, lty='solid')
  graphics.off()
}

#' Plot the results of a network meta-analysis
#'
#' @param df A data frame containing summary statistics for the posterior
#'   distribution of the parameters in a network meta-analysis.
#' @param yvar A character string specifying the name of the column in \code{df}
#'   that contains the variable to be used on the y-axis. Usually this will be
#'   the names of the treatments in the network.
#' @param xvar A character string specifying the name of the column in \code{df}
#'   that contains the variable to be used on the x-axis. Usually this will be
#'   the median of the posterior distribution of the treatment effect
#' @param lowLimit,hiLimit Character strings specifying the names of the
#'   columns in \code{df} that specify the lower and upper limits of the 95\%
#'   credible intervals. If you have followed the work flow in this package
#'   these will be called 'CrI_lower' and 'CrI_upper' respectively.
#' @param xlabel A character string to be used as a label on the x-axis. This is
#'   usually the type of effect measure, e.g. 'Odds ratio', 'Rate Ratio' etc.
#' @param noEffectLine An integer indicating where to draw the line of no
#'   effect. For ratio measures this is usually 1. For mean differences this is
#'   usually 0.
#' @param report_order A character string indicating whether the treatments
#'   should be reported in the order the order they were provided (default) or
#'   if a custom order is required. Acceptable values are 'default' or
#'   'custom'.
#'
#' @details This function should most commonly be used with the output from
#'   \code{extractTOI} to plot the results for a particular treatment of interest. This
#'   function uses ggplot to do the actual plotting and returns a ggplot object
#'   which may be further modified if required or just saved to a file.
#'
#' @return A ggplot object
#'
#' @seealso \code{\link{extractTOI}}, \code{\link[ggplot2]{geom_point}},
#'   \code{\link[ggplot2]{geom_errorbarh}}
plotEstimates = function(df, yvar, xvar = 'median', lowLimit = 'CrI_lower',
                         hiLimit = 'CrI_upper', xlabel = 'Effect Estimate',
                         noEffectLine = 1, reportOrder = 'default') {
  df = as.data.frame(df)
  if (reportOrder == 'default') {
    df[,yvar] = factor(df[,yvar], levels = sort(unique(df[,yvar])))
    p = ggplot2::ggplot(df) + ggplot2::geom_point(ggplot2::aes_string(x = xvar, y = yvar), size = 4)
    p = p + ggplot2::geom_errorbarh(
      ggplot2::aes_string(x = xvar, y = yvar, xmax = hiLimit, xmin = lowLimit),
      height = 0.15
      )
    p = p + ggplot2::scale_y_discrete(limits = rev(levels(df[,yvar])))
  }
  if(reportOrder == 'custom') {
    p = ggplot2::ggplot(df) + ggplot2::geom_point(ggplot2::aes_string(x = xvar, y = 'Order'), size = 4)
    p = p + ggplot2::geom_errorbarh(
      ggplot2::aes_string(x = xvar, y = 'Order', xmax = hiLimit, xmin = lowLimit),
      height = 0.15
      )
    p = p + ggplot2::scale_y_reverse(breaks = df[,'Order'], labels = df[,yvar])
  }

  #set sensible scale for x-axis
  #default is for ratio measures (OR, HR etc).
  #Mean difference needs special handling
  xrange = range(df[,lowLimit], df[,hiLimit])
  if (xlabel != 'Mean Difference' | xlabel != 'Effect Estimate') {
    xrange[1] = ifelse(xrange[1] > 0, 0, NA)
    xrange[2] = ifelse(xrange[2] < 2, 2, NA)
  } else {
    xrange = c(NA,NA)
  }

  p = p + ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n = 10), limits = xrange)
  p = p + ggplot2::labs(x = xlabel)
  p = p + ggplot2::geom_vline(xintercept = noEffectLine)
  p = p + ggplot2::theme_bw()
  p = p + ggplot2::theme(
    axis.title.y = ggplot2::element_blank(), panel.grid = ggplot2::element_blank(),
    axis.text = ggplot2::element_text(size = 12),
    axis.title = ggplot2::element_text(size = 12, vjust = 0),
    panel.border = ggplot2::element_rect(colour = 'black')
  )
}
#close but not quite
#
# yvar = 'comparison'
# xvar = 'effect'
# lowLimit = 'lower'
# hiLimit = 'upper'
#
# df = as.data.frame(nsRes)
# #create a column to use for the y-axis
# df$comparison = paste0(df$Intervention, '\n', df$Comparator)
#
# #re-arrange the data
# dir = dplyr::select(df, 3:7, comparison)
# colnames(dir)[3:5] = c('effect', 'lower', 'upper')
# dir$method = 'direct'
# ind = dplyr::select(df, 3:4, 8:10, comparison)
# colnames(ind)[3:5] = c('effect', 'lower', 'upper')
# ind$method = 'indirect'
# cons = dplyr::select(df, 3:4, 11:13, comparison)
# colnames(cons)[3:5] = c('effect', 'lower', 'upper')
# cons$method = 'pooled'
# pd = bind_rows(dir, ind, cons)
#
# pd = as.data.frame(pd)
# pd[,'comparison'] = factor(pd[,'comparison'], levels = sort(unique(pd[,'comparison'])))
#
# dodge = ggplot2::position_dodge(width=1)
# p = ggplot2::ggplot(pd) +
#   ggplot2::geom_point(ggplot2::aes_string(x = yvar, y = xvar, colour = 'method'), size = 4,
#                       position = dodge)
#
# p = p + ggplot2::geom_errorbar(
#   ggplot2::aes_string(x = yvar, y = xvar, ymax = hiLimit, ymin = lowLimit, colour = 'method'),
#   height = 0.25, position = dodge
# )
# p = p + ggplot2::scale_x_discrete(limits = rev(levels(pd[,yvar])), expand = c(0, 1))
#
# p = p + ggplot2::coord_flip()
#
#
