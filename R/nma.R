#' Extract treatment comparison information from mtc summary output
#'
#' @param A data frame as returned by \code{calcAllPairs}
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
  }
  #order by treatment coding in the network
  coding$Order = coding$id

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
    rt, file = resultsFile, sheetName = 'Report', showNA = FALSE, append = TRUE
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
#' @param orderResults Logical indicating whether the results are returned in a
#'   specific order or in the order they were given. The default is
#'   \code{FALSE}. If this is set to \code{TRUE} then \code{treatments} must
#'   contain a column named 'Order', see above.
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
                      orderResults = FALSE) {

  #extract the treatment of interest
  if (intervention == TRUE) {
    df = dplyr::filter(df, TreatmentB == toi)
  } else {
    df = dplyr::filter(df, TreatmentA == toi)
  }

  #sort if required
  if (intervention == TRUE && orderResults == TRUE) {
    df = dplyr::left_join(df, treatments[,2:3], by = c('nameA' = 'description'))
    df = dplyr::arrange(df, Order)
    df$Order = 1:nrow(df)
  }
  if (intervention == FALSE && orderResults == TRUE) {
    df = dplyr::left_join(df, treatments[,2:3], by = c('nameB' = 'description'))
    df = dplyr::arrange(df, Order)
    df$Order = 1:nrow(df)
  }

  return(df)
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
#' @param lowLimit, hiLimit Character strings specifying the names of the
#'   columns in \code{df} that specify the lower and upper limits of the 95%
#'   credible intervals. If you have followed the work flow in this package
#'   these will be called 'CrI_lower' and 'CrI_upper' respectively.
#' @param xlabel A character string to be used as a label on the x-axis. This is
#'   usually the type of effect measure, e.g. 'Odds ratio', 'Rate Ratio' etc.
#' @param noEffectLine An integer indicating where to draw the line of no
#'   effect. For ratio measures this is usually 1. For mean differences this is
#'   usually 0.
#' @param yOrder A character string giving the name of the column in \code{df}
#'   that specifies the order of treatments on the y-axis. The default is NA in
#'   which case the treatments will be displayed on the y-axis alphabetically
#'   from top to bottom. If a specific order is required this is best specified
#'   in a column named 'Order' containing a series of numbers indicating the
#'   plotting order from top to bottom.
#'
#' @details This function should most commonly be used with the output from
#'   extractTOI to plot the results for a particular treatment of interest. This
#'   function uses ggplot to do the actual plotting and returns a ggplot object
#'   which may be further modified if required or just saved to a file.
#'
#' @return A ggplot object
#'
#' @seealso \code{\link{extractTOI}}, \code{\link[ggplot2]{geom_point}},
#'   \code{\link[ggplot2]{geom_errorbarh}}
plotEstimates = function(df, yvar, xvar = 'median', lowLimit = 'CrI_lower',
                         hiLimit = 'CrI_upper', xlabel = 'Effect Estimate',
                         noEffectLine = 1, yOrder = NA) {
  df = as.data.frame(df)
  if (is.na(yOrder)) {
    df[,yvar] = factor(df[,yvar], levels = sort(unique(df[,yvar])))
    p = ggplot2::ggplot(df) + ggplot2::geom_point(ggplot2::aes_string(x = xvar, y = yvar), size = 4)
    p = p + ggplot2::geom_errorbarh(
      ggplot2::aes_string(x = xvar, y = yvar, xmax = hiLimit, xmin = lowLimit),
      height = 0.15
      )
    p = p + ggplot2::scale_y_discrete(limits = rev(levels(df[,yvar])))
  } else {
    p = ggplot2::ggplot(df) + ggplot2::geom_point(ggplot2::aes_string(x = xvar, y = yOrder), size = 4)
    p = p + ggplot2::geom_errorbarh(
      ggplot2::aes_string(x = xvar, y = yOrder, xmax = hiLimit, xmin = lowLimit),
      height = 0.15
      )
    p = p + ggplot2::scale_y_reverse(breaks = df[,yOrder], labels = df[,yvar])
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