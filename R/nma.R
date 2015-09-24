#load packages required by functions in this file
#library(tools, quietly = TRUE)
#require(gridExtra)

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

nameTreatments = function(results, coding, ...) {
  #map the treatment names to the numerical codes
  #results = data frame of winbugs output that has already been cleaned up with
  #extract comparisons
  #coding = data frame that maps treatment names to their numerical codes

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

summariseMTC = getS3method('summary', 'mtc.result')

calcAllPairs = function(mtcRes, expon = FALSE, ...) {
  #This function takes two arguments
  #mtcRes is an object of class mtc.result, i.e. the output from mtc.run in the
  #gemtc package
  #expon controls whether the output is exponentiated or not. If the input is
  #log OR (or HR or RR)
  # set this to TRUE to get output on the linear scale

  tid = as.integer(mtcRes$model$network$treatments$id)
  for (t in 1:length(tid)) {
    re = suppressWarnings(summariseMTC(gemtc::relative.effect(
      mtcRes, t1 = tid[t], preserve.extra = FALSE
    )))

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

extractModelFit = function(mtcRes) {
  #mtcRes - an mtc.result object as returned by mtc.run from the gemtc package
  modelSummary = summariseMTC(mtcRes)
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

extractResults = function(res, resultsFile, includesPlacebo = FALSE, ...) {
  #calculate all pairwise effects
  pairwiseResults = calcAllPairs(res, ...)
  rbutils::saveXLSX(
    as.data.frame(pairwiseResults), file = resultsFile, sheetName = 'Raw',
    showNA = FALSE, row.names = FALSE, append = TRUE
  )

  #extract the comparison info from the first column
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

extractTOI = function(df, treatments, toi, intervention = TRUE,
                      orderResults = FALSE) {
  #df   - a data frame of results as produced after running extractComparison
  #and nameTreatments
  #toi  - the code number of the treatment of interest (toi) in the network
  #treatments - a table of treatments with columns 'id', 'description', 'Order'
  #The 'Order' column is optional. If present it should be a column of integers
  #describing the order the results should be returned in. Note that the name is
  #case sensitive. If you wish to use thise feature set orderResults=TRUE
  #intervention - specifies whether the treatment of interest is an intervention
  #or a comparator (i.e. placebo)
  #orderResults - controls whether the results are returned in a specific order
  #or in the order they were given.
  #If present this should be a vector if integers

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

saveDiagnostics = function(mtc, directory) {
  #Save some convergence diagnostics as pdf files. Plots are based on ggmcmc package
  #mtc - an object of class mtc.result as returned but mtc.run in the gemtc package
  #directory - a file path indicating the folder to save the results

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

plotEstimates = function(df, yvar, xvar = 'median', lowLimit = 'CrI_lower',
                         hiLimit = 'CrI_upper', xlabel = 'Effect Estimate',
                         noEffectLine = 1, yOrder = NA) {
  df = as.data.frame(df)
  if (is.na(yOrder)) {
    df[,yvar] = factor(df[,yvar], levels = sort(unique(df[,yvar])))
    p = ggplot2::ggplot(df) + ggplot2::geom_point(aes_string(x = xvar, y = yvar), size = 4)
    p = p + ggplot2::geom_errorbarh(
      aes_string(x = xvar, y = yvar, xmax = hiLimit, xmin = lowLimit),
      height = 0.15
      )
    p = p + ggplot2::scale_y_discrete(limits = rev(levels(df[,yvar])))
  } else {
    p = ggplot2::ggplot(df) + ggplot2::geom_point(aes_string(x = xvar, y = yOrder), size = 4)
    p = p + ggplot2::geom_errorbarh(
      aes_string(x = xvar, y = yOrder, xmax = hiLimit, xmin = lowLimit),
      height = 0.15
      )
    p = p + ggplot::scale_y_reverse(breaks = df[,yOrder], labels = df[,yvar])
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

  p = p + ggplot2::scale_x_continuous(breaks = pretty_breaks(n = 10), limits = xrange)
  p = p + ggplot::labs(x = xlabel)
  p = p + ggplot2::geom_vline(xintercept = noEffectLine)
  p = p + ggplot2::theme_bw()
  p = p + ggplot2::theme(
    axis.title.y = element_blank(), panel.grid = element_blank(),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12, vjust = 0),
    panel.border = element_rect(colour = 'black')
  )
}