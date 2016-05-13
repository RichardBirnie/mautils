#' Perform meta-analysis of diagnostic accuracy data using the Reitsma bivariate
#' model
#'
#' @param df A data frame with a column of study IDs and integer variables for
#'   observed frequencies of true positives, false negatives, false positives
#'   and true negatives. The names of the columns should be StudyID, TP, FN, FP
#'   and TN.
#' @param file A character string specifying the base directory for the
#'   analysis. This function will create a subdirectory 'Results' within this
#'   directory will be a 'Figures' subdirectory.
#' @param analysis_case A character string indicating a short name that
#'   identifies this analysis. This is used in filenames so keep it short.
#' @param test_id A character string indicating the identity of the diagnostic
#'   test

runDAcc = function(df, file, test_id, analysis_case = 'Default', ...){
  #set up folders for the results and figures.
  resDir = file.path(file, 'Results')
  if (!dir.exists(resDir)) {
    dir.create(resDir, recursive = TRUE)
  }
  figDir = file.path(resDir, 'Figures')
  if (!file.exists(figDir)) {
    dir.create(figDir, recursive = TRUE)
  }
  resFile = file.path(resDir, paste0(test_id, '_', analysis_case, '.xlsx'))
  #save the input data with the results file
  rbutils::saveXLSX(
    as.data.frame(df), file = resFile, sheetName = 'Data',
    showNA = FALSE, row.names = FALSE, append = TRUE
  )

  #use madad to calculate descriptive statistics
  #includes sensitivity, specificity, false positive rate, DOR, positive and
  #negative likelihood ratios, tests for heterogeneity
  des = mada::madad(df, ...)
  des$studlab = df$StudyID
  dacc = extractDACC(des = des, input = df)
  daccSummary = reportDACC(dacc)

  #save the descriptive results
  rbutils::saveXLSX(
    as.data.frame(dacc), file = resFile, sheetName = 'DAccPerStudy',
    showNA = FALSE, row.names = FALSE, append = TRUE
  )
  rbutils::saveXLSX(
    as.data.frame(daccSummary), file = resFile, sheetName = 'DAccPerStudyReport',
    showNA = FALSE, row.names = FALSE, append = TRUE
  )

  #I-squared
  i2 = i2dacc(df)
  rbutils::saveXLSX(
    as.data.frame(i2), file = resFile, sheetName = 'Univariate I2',
    showNA = FALSE, row.names = FALSE, append = TRUE
  )

  #make some forest plots of sensitivity and specificity for the individual
  #studies for descriptive purposes only
  jpeg(
    filename = file.path(figDir, 'Sens_Spec.jpg'), width = 35, height = 12,
    units = 'cm', quality = 100, res = 300
  )
  par(mfcol = c(1, 2))
  forestmada2(x = des$sens$sens, ci = des$sens$sens.ci,
              main = NULL, xlab='Sensitivity',
              snames = paste0('Study ', des$studlab), xlim = c(0,1))
  forestmada2(x = des$spec$spec, ci = des$spec$spec.ci,
              main = NULL, xlab='Specificity',
              snames = paste0('Study ', des$studlab), xlim = c(0,1))
  graphics.off()

  #Bivariate model using the reitsma method
  fit.r = mada::reitsma(df, method='reml', ...)

  #extract summary for the key parameters and save
  summaryRes = extractSummary(fit.reitsma = fit.r)
  rbutils::saveXLSX(
    as.data.frame(summaryRes), file = resFile, sheetName = 'Bivariate',
    showNA = FALSE, row.names = FALSE, append = TRUE
  )

  #plot the SROC curve
  jpeg(
    filename = file.path(figDir, 'SROC.jpg'),
    width = 15, height = 15, units = 'cm',
    quality = 100, res = 300
  )
  plotSROC(fit.r)
  graphics.off()
}

#' Extract descriptive accuracy results from the output of \code{madad}
#'
#' @param des An object of class \code{madad} as returned by the function
#'   \code{madad} in the \code{mada} package
#' @param input A data frame of 2x2 diagnostic accuracy data. The following
#'   columns must be present 'StudyID, 'TP', 'TN', 'FP', 'FN' corresponding to
#'   the number of True Positives, True Negatives, False Postives and False
#'   Negatives respectively.
#'
#' @return A data frame containing the input data plus summary estimates from
#'   each study for sensitvity, specificity, positive likelihood ratio, negative
#'   likelihood ratio and diagnostic odds ratio
#' @seealso \code{\link[mada]{madad}}
extractDACC = function(des, input){
  dacc = dplyr::data_frame(
    'Sens' = des$sens$sens,
    'Sens_Lower_CI' = des$sens$sens.ci[,1],
    'Sens_Upper_CI' = des$sens$sens.ci[,2],
    'Spec' = des$spec$spec,
    'Spec_Lower_CI' = des$spec$spec.ci[,1],
    'Spec_Upper_CI' = des$spec$spec.ci[,2],
    'PosLR' = des$posLR$posLR,
    'PosLR_Lower_CI' = des$posLR$posLR.ci[,1],
    'PosLR_Upper_CI' = des$posLR$posLR.ci[,2],
    'NegLR' = des$negLR$negLR,
    'NegLR_Lower_CI' = des$negLR$negLR.ci[,1],
    'NegLR_Upper_CI' = des$negLR$negLR.ci[,2],
    'DOR' = des$DOR$DOR,
    'DOR_Lower_CI' = des$DOR$DOR.ci[,1],
    'DOR_Upper_CI' = des$DOR$DOR.ci[,2]
  )
  #add the total N to the data
  input = dplyr::group_by(input, StudyID) %>%
    dplyr::mutate('TotalN' = TP + FN + FP + TN) %>%
    dplyr::ungroup()
  #combine data and accuracy
  dacc = dplyr::bind_cols(input, dacc)
  dacc[,7:12] = round(dacc[,7:12]*100, 2)
  dacc[,13:ncol(dacc)] = round(dacc[,13:ncol(dacc)], 3)
  return(dacc)
}

#' Summarise diagnostic accuracy data in a convenient format
#'
#' @param df A data frame of diagnostic accuracy data as returned by
#'   \code{extractDACC}
#'
#' @return A data frame summarising the diagnostic accuracy information for each
#'   study
#' @seealso \code{extractDACC}, \code{\link[mada]{madad}}
reportDACC = function(df){
  data_frame(
    'StudyID' = df$StudyID,
    'TP' = df$TP,
    'FN' = df$FN,
    'FP' = df$FP,
    'TN' = df$TN,
    'TotalN' = df$TotalN,
    'Sensitivity' = paste0(sprintf('%.2f', df$Sens),
                           ' (', sprintf('%.2f', df$Sens_Lower_CI), ' to ',
                           sprintf('%.2f', df$Sens_Upper_CI), ')'),
    'Specificity' = paste0(sprintf('%.2f', df$Spec),
                           ' (', sprintf('%.2f', df$Spec_Lower_CI), ' to ',
                           sprintf('%.2f', df$Spec_Upper_CI), ')'),
    'Positive LR' = paste0(sprintf('%.2f', df$PosLR),
                           ' (', sprintf('%.2f', df$PosLR_Lower_CI), ' to ',
                           sprintf('%.2f', df$PosLR_Upper_CI), ')'),
    'Negative LR' = paste0(sprintf('%.2f', df$NegLR),
                           ' (', sprintf('%.2f', df$NegLR_Lower_CI), ' to ',
                           sprintf('%.2f', df$NegLR_Upper_CI), ')'),
    'DOR' = paste0(sprintf('%.2f', df$DOR),
                   ' (', sprintf('%.2f', df$DOR_Lower_CI), ' to ',
                   sprintf('%.2f', df$DOR_Upper_CI), ')')
  )
}

#' Extract essential summary estimates after fitting the Reitsma bivariate model
#'
#' @param fit.reitsma An object of class \code{reitsma} returned by fitting the
#'   Reitsma bivariate model using the function named \code{reitsma} in the mada
#'   package
#'
#' @details This function takes the output from a fit of the Reitsma bivariate
#'   model, extracts the essential summary estimates and returns them in a
#'   convenient data frame. The Positive/Negative likelihood ratios and the
#'   diagnostic odds ratio are calculated using the method of Zwindermann and
#'   Bossuyt. This method uses an MCMC sampling based approach to calculate the
#'   results based on a fit of the bivariate model.
#'
#' @return A data frame containing the mean and 95% confidence intervals for the
#'   summary estimates of sensitivity, specificity, positive LR, negative LR and
#'   DOR.
#'
#' @seealso \code{\link[mada]{reitsma}}, \code{\link[mada]{SummaryPts}}
extractSummary = function(fit.reitsma) {
  #get summary from reitsma object
  #convert the essential information into a data frame
  summ.fit.reitsma = summary(fit.reitsma)
  summEst = as.data.frame(summ.fit.reitsma$coefficients)[3:4, c(1,5:6)]
  colnames(summEst) = c('Mean', 'Lower CI', 'Upper CI')
  #convert FPR to Spec
  rownames(summEst)[2] = 'specificity'
  summEst[2,] = 1 - summEst[2,]
  summEst[2, c('Upper CI', 'Lower CI')] = summEst[2, c('Lower CI', 'Upper CI')]
  #report Sens and Spec as percentages
  summEst = summEst * 100
  summEst = data.frame(
    'Parameter' = rownames(summEst), summEst,
    row.names = NULL, stringsAsFactors = FALSE
  )

  #get Likelihood ratios and DOR based on Zwindermann and Bossuyt method
  mcmc = mada::SummaryPts(fit.reitsma)
  summMCMC = as.data.frame(summary(mcmc))
  colnames(summMCMC)[3:4] = c('Lower CI', 'Upper CI')
  summMCMC = summMCMC[c(1:2,4), c(1,3:4)]
  summMCMC = data.frame(
    'Parameter' = rownames(summMCMC), summMCMC,
    row.names = NULL, stringsAsFactors = FALSE
  )

  #combine the two sets of results, use better names and return
  summEst = rbind(summEst, summMCMC)
  summEst$Parameter = c('Sensitivity', 'Specificity',
                        'Positive LR', 'Negative LR', 'DOR')
  return(summEst)
}

#' Calculate I-squared for diagnostic accuracy parameters
#'
#' @param df A data frame of 2x2 diagnostic accuracy data. The following
#'   columns must be present 'StudyID, 'TP', 'TN', 'FP', 'FN' corresponding to
#'   the number of True Positives, True Negatives, False Postives and False
#'   Negatives respectively.
#'
#' @details This function will calculate the I-squared heterogeneity statistic
#'   with 95% confidence intervals for each of the following diagnostic accuracy
#'   parameters: Sensitivity, Specificity, Diagnostic Odds Ratio and both
#'   positive and negative likelihood ratios.
#'
#'   The I-squared values are calculated using functions from the \code{meta}
#'   package but no pooling of the study estimates are performed. These
#'   functions are used only for convenience to provide the I-squared values.
#'   Sensitivity and specificity use \code{metaprop} for proportions. DOR,
#'   Positive LR and negative LR use \code{metagen} for generic inverse
#'   variance.
#'
#' @return A data frame of I-squared values with confidence intervals for each
#'   of the parameters.
#'
#' @seealso \code{\link[meta]{metaprop}}, \code{\link[meta]{metagen}}
i2dacc = function(df) {
  #calculate sens, spec etc for each study
  des = mada::madad(df, correction.control = 'single')

  #calculate I-squared based on results for each study
  #use metaprop for univariate sens and spec using logit transform.
  #There is no meta-analysis here, this is just a convenient way to get
  #an I-squared value
  #sensitivity
  sensN = df$TP + df$FN
  sens = meta::metaprop(
    df$TP, sensN, comb.fixed = FALSE, comb.random = FALSE,
    sm = 'PLOGIT'
  )
  i2Sens = data.frame(
    'I2' = sens$I2, 'I2.lower' = sens$lower.I2,
    'I2.upper' = sens$upper.I2
  ) * 100
  #specificity
  specN = df$FP + df$TN
  spec = meta::metaprop(
    df$TN, specN, comb.fixed = FALSE, comb.random = FALSE,
    sm = 'PLOGIT'
  )
  i2Spec = data.frame(
    'I2' = spec$I2, 'I2.lower' = spec$lower.I2,
    'I2.upper' = spec$upper.I2
  ) * 100

  #Use metagen for DOR
  #There is no meta-analysis here, this is just a convenient way to get
  #an I-squared value
  dor = meta::metagen(
    TE = log(des$DOR$DOR), seTE = des$DOR$se.lnDOR,
    comb.fixed = FALSE, comb.random = FALSE, sm = 'OR'
  )
  i2DOR = data.frame(
    'I2' = dor$I2, 'I2.lower' = dor$lower.I2,
    'I2.upper' = dor$upper.I2
  ) * 100

  #Use metagen for Positive LR
  #There is no meta-analysis here, this is just a convenient way to get
  #an I-squared value
  posLR = meta::metagen(
    TE = log(des$posLR$posLR), seTE = des$posLR$se.lnposLR,
    comb.fixed = FALSE, comb.random = FALSE, sm = 'OR'
  )
  i2posLR = data.frame(
    'I2' = posLR$I2, 'I2.lower' = posLR$lower.I2,
    'I2.upper' = posLR$upper.I2
  ) * 100
  #Use metagen for Negative LR
  #There is no meta-analysis here, this is just a convenient way to get
  #an I-squared value
  negLR = meta::metagen(
    TE = log(des$negLR$negLR), seTE = des$negLR$se.lnnegLR,
    comb.fixed = FALSE, comb.random = FALSE, sm = 'OR'
  )
  i2negLR = data.frame(
    'I2' = negLR$I2, 'I2.lower' = negLR$lower.I2,
    'I2.upper' = negLR$upper.I2
  ) * 100

  #combine the results and return
  i2 = dplyr::bind_rows(
        i2Sens,
        i2Spec,
        i2posLR,
        i2negLR,
        i2DOR
      )
  i2$Parameter = c('Sensitivity', 'Specificity', 'Positive LR', 'Negative LR', 'DOR')
  dplyr::select(i2, Parameter, contains('I2'))
}

#' Extend the forestmada function with more flexible axes
#'
#' @details This function extends the forestmada function from the mada package
#'   to allow the user to control the range of the x-axes. The primary purpose
#'   is to allow the user to set \code{xlim = c(0,1)} in forest plots of
#'   sensitivity or specificity
#'
#' @seealso \code{\link[mada]{forest}} for details of the other arguments
forestmada2 = function (x, ci, plotci = TRUE, main = "Forest plot", xlab = NULL,
          digits = 2L, snames = NULL, subset = NULL, pch = 15, cex = 1,
          cipoly = NULL, polycol = NA, xlim, ...)
{
  stopifnot(length(x) == dim(ci)[1], all(!is.na(c(ci, x))),
            is.logical(plotci))
  if (!is.null(snames)) {
    stopifnot(length(snames) == length(x))
  }
  if (is.null(snames)) {
    snames <- paste("Study", 1:length(x))
  }
  if (!is.null(subset)) {
    stopifnot(length(subset) > 0, is.integer(subset))
  }
  if (is.null(subset)) {
    subset <- 1:length(x)
  }
  if (!is.null(cipoly)) {
    stopifnot(length(cipoly) == length(x))
    cireg <- which(!cipoly)
    cipoly <- which(cipoly)
  }
  if (is.null(cipoly)) {
    cireg <- 1:length(x)
  }
  x <- x[subset]
  lb <- ci[subset, 1]
  ub <- ci[subset, 2]
  snames <- snames[subset]
  if (is.null(xlab)) {
    xlab <- ""
  }
  N <- length(x)
  plotrange <- max(ub) - min(lb)
  if(missing(xlim)) {
    if (plotci) {
      xlim <- c(min(lb) - plotrange * 1.2, max(ub) + plotrange *
                  1.2)
    }
    else {
      xlim <- c(min(lb) - plotrange * 1.2, max(ub) + plotrange *
                  0.4)
    }
    default.xlim = TRUE
  } else {
    xlim = xlim
    default.xlim = FALSE
  }

  ylim <- c(0.5, N + 1)
  if(!default.xlim) { par(mar = c(5,8,1,8))}
  plot(
    NA, NA, xlim = xlim, ylim = ylim, xlab = xlab, ylab = "",
    yaxt = "n", xaxt = "n", xaxs = "i", bty = "n", main = main,
    ...
  )
  abline(h = ylim[2], ...)
  for (i in cireg) {
    points(x[i], (N:1)[i], pch = pch)
    arrows(
      lb[i], (N:1)[i], ub[i], angle = 90, code = 3,
      length = 0.05, ...
    )
  }
  for (i in cipoly) {
    polygon(
      x = c(lb[i], x[i], ub[i], x[i]), y = c((N:1)[i],
                                             (N:1)[i] + 0.25, (N:1)[i], (N:1)[i] - 0.25), col = polycol
    )
  }
  if(default.xlim) {
    text(
      x = xlim[1], N:1, labels = snames[1:N], pos = 4, cex = cex
    )
  } else {
    mtext(snames[1:N], side = 2, at = N:1, line = 1, las = 1)
  }
  if (plotci) {
    citext <- format(round(cbind(x, ci), digits = digits),
                     nsmall = digits)
    citext <- paste(citext[, 1], " [", citext[, 2], ", ",
                    citext[, 3], "]", sep = "")
    if(default.xlim) {
      text(
        x = xlim[2], N:1, labels = citext, pos = 2, cex = cex
      )
    } else {
      mtext(citext, side = 4, at = N:1, line = 1, las = 1)
    }
  }

  if(default.xlim) {
    axis(1, at = round(seq(
      from = min(lb), to = max(ub), length.out = 5
    ),
    digits), cex = cex)
  } else {
    axis(1, at = seq(from = min(xlim), to = max(xlim), by = 0.1),
    cex = cex)
  }

  return(invisible(NULL))
}


#' Extend the plot method for reitsma objects
#'
#' @details This function extends the plot method for reitsma objects from the mada package
#'   to allow the user to set the label for the x-axis.
#'
#' @seealso \code{\link[mada]{reitsma-class}} for details of the other arguments
plot.reitsma2 <- function(x, extrapolate = FALSE, plotsumm = TRUE, level = 0.95,
                          ylim = c(0,1), xlim = c(0,1), pch = 1,
                          sroclty = 1, sroclwd = 1,
                          predict = FALSE, predlty = 3, predlwd = 1,
                          type = "ruttergatsonis", xlab,
                          ...)
{
  plot(c(2,2), ylim = ylim, xlim = xlim,
       xlab = ifelse(missing(xlab), "False Positive Rate", xlab),
       ylab = "Sensitivity", ...)
  if(length(coef(x)) == 2){
    FP <- x$freqdata$FP
    negatives <- FP + x$freqdata$TN
    FPR <- FP/negatives

    if(extrapolate){bound = c(0,1)}
    if(!extrapolate){bound = c(min(FPR), max(FPR))}
    srocmat <- mada::sroc(x, type = type)
    lines(srocmat[cut(srocmat[,1],bound, "withinbound") == "withinbound",],
          lty = sroclty, lwd = sroclwd)
  }else{
    warning("Not plotting any SROC for meta-regression")
  }
  if(plotsumm){mada::ROCellipse(x, level = level, add = TRUE, pch = pch, ...)}

  if(predict){
    alpha.sens <- x$alphasens
    alpha.fpr <- x$alphafpr
    mu <- x$coefficients["(Intercept)",]
    Sigma <- x$Psi + vcov(x)
    talphaellipse <- ellipse::ellipse(Sigma, centre = mu, level = level)
    predellipse <- matrix(0, ncol = 2, nrow = nrow(talphaellipse))
    predellipse[,1] <- mada:::inv.trafo(alpha.fpr, talphaellipse[,2])
    predellipse[,2] <- mada:::inv.trafo(alpha.sens, talphaellipse[,1])
    lines(predellipse, lty = predlty, lwd = predlwd)
  }

  return(invisible(NULL))
}

#' Draw a Summary ROC curve based on the results of the bivariate model
#'
#' @param fit.reitsma An object of class \code{reitsma} returned by fitting the
#'   Reitsma bivariate model using the function named \code{reitsma} in the mada
#'   package
#' @param show.data A logical indicating whether or not to show the mean
#'   sensitivity and specificity of the individual studies
#' @param show.prediction A logical indicating whether or not draw a prediction region
#' @param summary.legend.pos A character string indicating the position of the
#'   legend for the summary point and data points. This is passed to the \code{legend} function
#'   and must be one of 'bottomleft', 'bottomright' (default), 'topleft',
#'   'topright'.
#' @param sroc.legend.pos A character string indicating the position of the
#'   legend for the SROC curve. This is passed to the \code{legend} function
#'   and must be one of 'bottomleft' (default), 'bottomright', 'topleft',
#'   'topright'.
#' @param ... Additional arguments passed to the underlying plot method
#'
#' @details This function will draw a summary ROC curve based on a fit of the
#'   Reitsma bivariate model. The plot will always include a summary point with
#'   a confidence ellipse and a SROC curve. By default the plot will include the
#'   data points from the individual studies and a 95% prediction region around
#'   the summary point. These can be excluded by setting the arguments show.data
#'   or show.prediction to \code{FALSE}
#'
#'   The argument extrapolate may also be set to TRUE or FALSE to indicate
#'   whether to extrapolate the SROC curve beyond the observed data points. The
#'   default is FALSE indicating no extrapolation
#'
#' @seealso \code{\link[mada]{reitsma-class}}
plotSROC = function(fit.reitsma, show.data = TRUE, show.prediction = TRUE,
                    summary.legend.pos = 'bottomright',
                    sroc.legend.pos = 'bottomleft', ...) {
  plot.reitsma2(fit.reitsma, sroclwd = 2,
                predict = show.prediction, predlty = 3,
                pch = 15, las = 1, cex = 1.5,
                xlab = '1 - Specificity', ...
                )
  #control whether or not to show data from individual studies
  if(show.data) {
    points(
      mada::fpr(fit.reitsma$data), mada::sens(fit.reitsma$data),
      pch = 19, col = 'grey'
    )
    legend(summary.legend.pos, c("data", "summary estimate"),
      pch = c(19,15), col = c('grey', 'black')
    )
  } else {
    legend(summary.legend.pos, c("summary estimate"),
           pch = 15, col = c('black')
    )
  }
  #modify the legend depending on whether or not the prediction region
  #is shown
  if(show.prediction) {
    legend(sroc.legend.pos,
      c("SROC", "95% confidence region", "95% prediction region"),
      lwd = c(2,1,1), lty = c(1,1,3)
    )
  } else {
    legend(sroc.legend.pos,
           c("SROC", "95% confidence region"),
           lwd = c(2,1), lty = c(1,1)
    )
  }
}
