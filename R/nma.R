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
#' @export
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
#' @export
nameTreatments = function(results, coding, ...) {

  #drop order column from coding. Makes no sense here
  coding = coding[,1:2]
  coding$id = as.integer(coding$id)

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
#' @export
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
#' @export
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
      output = suppressWarnings(dplyr::bind_rows(output, out))
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
#' @export
extractModelFit = function(mtcRes) {
  modelSummary = .summariseMTC(mtcRes)
  dic = data.frame(
    'Mean' = modelSummary$DIC[1:3],
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
#'   where the results should be saved.
#' @param reference Name of the reference treatment for this analysis. This
#'   should match the name used in the input data file exactly
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
#' @seealso \code{\link[gemtc]{mtc.run}},
#'   \code{\link{calcAllPairs}}, \code{\link{extractComparison}},
#'   \code{\link{nameTreatments}}, \code{\link{makeTab}},
#'   \code{\link{extractModelFit}}
#'
#' @export
extractMTCResults = function(res, resultsFile, reference, ...) {
  #calculate all pairwise effects
  pairwiseResults = calcAllPairs(res, ...)

  #extract the comparison info
  #map treatment names to numbers
  pairwiseResults = extractComparison(pairwiseResults)
  pairwiseResults = nameTreatments(pairwiseResults, ...)
  XLConnect::writeWorksheetToFile(file = resultsFile, data = as.data.frame(pairwiseResults),
                                  sheet = 'AllComparisons', clearSheets = TRUE)

  #make a table of all pairwise comparisons and save it as an excel file
  reportTab = makeTab(results = pairwiseResults, ...)
  reportTab = dplyr::select(reportTab,-matches(reference)) #drop the reference column

  rt = as.data.frame(reportTab, check.names = FALSE) %>%
    tibble::rownames_to_column(var = 'Comparator')

  XLConnect::writeWorksheetToFile(file = resultsFile, data = as.data.frame(rt),
                                  sheet = 'Report', clearSheets = TRUE)

  #extract and save the model fit information
  modelFit = extractModelFit(res) %>%
    tibble::rownames_to_column(var = 'Parameter')
  XLConnect::writeWorksheetToFile(file = resultsFile, data = as.data.frame(modelFit),
                                  sheet = 'DIC', clearSheets = TRUE)

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
#' @export
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
#' @export
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
  #can also be difference of mean differences for continuous
  nsRes$IF = nsRes$dir.pe - nsRes$ind.pe
  #calculate variance of inconsistency factor
  #sum of Var(direct) + Var(indirect)
  varIF = ((nsRes$dir.ci.u - nsRes$dir.ci.l) / 3.92) ^ 2 +
    ((nsRes$ind.ci.u - nsRes$ind.ci.l) / 3.92) ^ 2
  #Std Error of inconsistency factor
  seIF = sqrt(varIF)
  #CI for inconsistency factor
  nsRes$IF.lower.ci = nsRes$IF - (1.96 * seIF)
  nsRes$IF.upper.ci = nsRes$IF + (1.96 * seIF)

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
    nsRes, 2:1, 17:16, dplyr::starts_with('dir'),
    dplyr::starts_with('ind'), dplyr::starts_with('cons'),
    dplyr::contains('IF'), p
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
#' @export
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

#'Extract ranking probabilities
#'
#' @param ranks An object of class \code{mtc.rank.probability} as returned by
#'   the \code{probability} function in the \code{gemtc} package
#' @param treatments A data frame with columns 'description' defining the
#'   treatment names and 'id' defining the treatment ID numbers.
#' @details This function takes the ranking probablilities returned by \code{rank.probability}, matches the treatment names to the id numbers and returns the results as a data frame
#'
#' @return A data frame
#'
#' @seealso \code{\link[gemtc]{rank.probability}}
#' @export
extractRanks = function(ranks, treatments) {
  treatments$id = as.integer(treatments$id)

  class(ranks) = 'matrix'
  ranks = as.data.frame(ranks)
  colnames(ranks) = paste0('Rank', 1:ncol(ranks))
  ranks$id = as.integer(rownames(ranks))

  ranks = dplyr::left_join(treatments[,1:2], ranks, by = 'id')
  ranks = dplyr::arrange(ranks, dplyr::desc(Rank1))
}

#' Extract the coda from an mtc.result object
#'
#' @param mtcResults An \code{mtc.result} object as returned by \code{mtc.run}
#'   from the gemtc package
#' @param summarise logical. If TRUE (default) and the mtc.result object
#'   includes more than one chain then the mean across all chains is calculated
#'   for each variable at each iteration of the chains
#'
#' @details This function takes the output from running an MTC model using
#'   mtc.run, extracts the coda from the \code{mtc.result} object and returns a
#'   data frame. If the summarise argument is TRUE and multiple chains were used
#'   for the MCMC then this function will return the mean across all chains for
#'   each variable in the model at each iteration of the chain. This is the
#'   default behaviour. If the summarise argument is FALSE then results for each
#'   individual chain are preserved.
#'
#' @return If summarise = TRUE:
#' \itemize{
#' \item A data frame with one column per variable plus one column indicating
#' the iteration of the chain. The values for each variable will be the mean
#' across all chains for each iteration of the chain
#' }
#' If summarise = FALSE:
#' \itemize{
#' \item A data frame with one column per variable plus one column indicating
#' the iteration of each chain and one column indicating the chain numbered 1,
#' 2, 3 etc. The values for each variable will be the value of the variable at
#' the indicated iteration of the indicated chain.
#' }
#'
#' @seealso \code{\link[gemtc]{mtc.run}}
#' @export
extractCoda = function(mtcResults, summarise=TRUE) {
  #starts with the output from mtc.run
  #extract coda from results object
  #coerce to a matrix and then to a data frame
  mtcCoda = coda::as.mcmc.list(mtcResults)
  mtcCoda = as.data.frame(as.matrix(mtcCoda, iters = TRUE, chains = TRUE))

  if(summarise & max(unique(mtcCoda$CHAIN)) >1) {
    #Get the mean across all three chains for each variable separately at each
    #iteration of the chains.
    #results are on log scale
    mtcCoda = dplyr::group_by(mtcCoda, ITER) %>%
      dplyr::summarise_each(funs(mean),-CHAIN)
  }
  return(mtcCoda)
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
#' @export
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
#' @export
plotEstimates = function(df, yvar, xvar = 'median', lowLimit = 'CrI_lower',
                         hiLimit = 'CrI_upper', xlabel = 'Effect Estimate',
                         noEffectLine = 1) {
  df = as.data.frame(df)
  if (!'Order' %in% colnames(df)) {
    df = df[order(df[,yvar]),]
    df$Order = 1:nrow(df)
  }
  df$values = paste0(sprintf("%.2f", df[,xvar]), ' (', sprintf("%.2f",df[,lowLimit]), ' to ', sprintf("%.2f", df[,hiLimit]), ')')

  #build the basic plot
  p = ggplot2::ggplot(df) + ggplot2::geom_point(ggplot2::aes_string(x = xvar, y = 'Order'), size = 2)
  p = p + ggplot2::geom_errorbarh(ggplot2::aes_string(
    x = xvar, y = 'Order', xmax = hiLimit, xmin = lowLimit
  ),
  height = 0.15)
  p = p + ggplot2::scale_y_reverse(breaks = df[,'Order'], labels = df[,yvar], expand = c(0, 0.2))


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

  #add annotations - values on right side
  p = p + ggplot2::theme(plot.margin = grid::unit(c(1,7,1,1), 'lines'))
  p = p + ggplot2::geom_text(ggplot2::aes_string(label = 'values',x = Inf, y = 'Order'),
                             hjust = -0.1, size = 3)

  #tidy up the final appearance
  p = p + ggplot2::theme(
    axis.title.y = ggplot2::element_blank(), panel.grid = ggplot2::element_blank(),
    axis.text = ggplot2::element_text(size = 8),
    axis.title = ggplot2::element_text(size = 8, vjust = 0),
    panel.border = ggplot2::element_rect(colour = 'black')
  )

  #Code to override clipping so that annotation is visible
  gt <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(p))
  gt$layout$clip[gt$layout$name == "panel"] <- "off"
  return(gt)
  # grid::grid.draw(gt)
}

#' Plot rank probabilities
#'
#' @param ranks A data frame of ranking probabilities for each treatment as
#'   returned by \code{extractRanks}
#' @details This function takes a data frame of ranking probabilities for each
#'   treatment in a network and constructs a 'rankogram'; i.e. a bar chart
#'   showing the probability of being ranked first, second, third etc. for each
#'   treatment
#' @return A ggplot object which can then be saved using an appropriate graphics
#'   device, e.g. jpeg, png, pdf etc.
#' @seealso \code{\link{extractRanks}}, \code{\link[gemtc]{rank.probability}}
#' @export
plotRanks = function(ranks) {
  ranks = tidyr::gather(ranks[,2:ncol(ranks)], description, Probability)
  colnames(ranks)[2] = 'Rank'
  p = ggplot2::ggplot(ranks) +
    ggplot2::geom_bar(
      ggplot2::aes(x = description, y = Probability, fill = Rank),
      position = 'dodge', stat = 'identity'
    )
  p = p + ggplot2::scale_fill_brewer(palette = 'Dark2')
  p = p + ggplot2::theme_bw()
  p = p + ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    axis.text = ggplot2::element_text(size = 12),
    axis.text.x = ggplot2::element_text(
      angle = 45, hjust = 1, vjust = 1
    ),
    axis.title = ggplot2::element_text(size = 12, vjust = 0),
    axis.title.x = ggplot2::element_blank(),
    panel.border = ggplot2::element_rect(colour = 'black')
  )
}
