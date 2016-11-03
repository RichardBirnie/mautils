library(mautils)
context('Indirect meta-analysis')

test_that('Indirect meta-analysis of ratio measures is correct for single cases', {
  #Input data are log odds ratios for binary data
  #This test should cover generic inverse variance pooling for any log ratio measure
  #Split RE and FE results and test both
  direct_input = read.table('../data/smoking_anno_OR_bucher_test.txt', sep='\t', header=TRUE)
  fe_input = direct_input[direct_input$Model=='Fixed',]
  re_input = direct_input[direct_input$Model=='Random',]

  #Fixed Effect
  fe_ind = bucher(
    abTE = fe_input$log.TE[1], se.abTE = fe_input$seTE.log[1],
    cbTE = fe_input$log.TE[2], se.cbTE = fe_input$seTE.log[2],
    backtransf = TRUE, effect = fe_input$Effect[1],
    model = fe_input$Model[1], intervention = fe_input$Intervention[1],
    comparator = fe_input$Intervention[2], common = fe_input$Comparator[1],
    ab.studies = fe_input$studies[1],
    cb.studies = fe_input$studies[2],
    continuous = FALSE
  )
  #test that the object has the right class and the expected content
  expect_is(fe_ind, 'data.frame')
  expect_that(fe_ind, equals(dget('../data/smoking_anno_OR_bucher_test_resFE.txt')))

  #Random Effects
  re_ind = bucher(
    abTE = re_input$log.TE[1], se.abTE = re_input$seTE.log[1],
    cbTE = re_input$log.TE[2], se.cbTE = re_input$seTE.log[2],
    backtransf = TRUE, effect = re_input$Effect[1],
    model = re_input$Model[1], intervention = re_input$Intervention[1],
    comparator = re_input$Intervention[2], common = re_input$Comparator[1],
    ab.studies = re_input$studies[1],
    cb.studies = re_input$studies[2],
    continuous = FALSE
  )
  #test that the object has the right class and the expected content
  expect_is(re_ind, 'data.frame')
  expect_that(re_ind, equals(dget('../data/smoking_anno_OR_bucher_test_resRE.txt')))

})

test_that('Indirect meta-analysis of continuous measures is correct for single cases', {
  #Input data are mean differences for continuous data
  #Split RE and FE results and test both
  direct_input = read.table('../data/parkinsons_anno_MD_bucher_test.txt', sep='\t', header=TRUE)
  fe_input = direct_input[direct_input$Model=='Fixed',]
  re_input = direct_input[direct_input$Model=='Random',]

  #Fixed Effect
  fe_ind = bucher(
    abTE = fe_input$TE[1], se.abTE = fe_input$seTE[1],
    cbTE = fe_input$TE[2], se.cbTE = fe_input$seTE[2],
    backtransf = FALSE, effect = fe_input$Effect[1],
    model = fe_input$Model[1], intervention = fe_input$Intervention[1],
    comparator = fe_input$Intervention[2], common = fe_input$Comparator[1],
    ab.studies = fe_input$studies[1],
    cb.studies = fe_input$studies[2],
    continuous = TRUE
  )
  #test that the object has the right class and the expected content
  expect_is(fe_ind, 'data.frame')
  expect_that(fe_ind, equals(dget('../data/parkinsons_anno_MD_bucher_test_resFE.txt')))

  #Random Effects
  re_ind = bucher(
    abTE = re_input$TE[1], se.abTE = re_input$seTE[1],
    cbTE = re_input$TE[2], se.cbTE = re_input$seTE[2],
    backtransf = FALSE, effect = re_input$Effect[1],
    model = re_input$Model[1], intervention = re_input$Intervention[1],
    comparator = re_input$Intervention[2], common = re_input$Comparator[1],
    ab.studies = re_input$studies[1],
    cb.studies = re_input$studies[2],
    continuous = TRUE
  )
  #test that the object has the right class and the expected content
  expect_is(re_ind, 'data.frame')
  expect_that(re_ind, equals(dget('../data/parkinsons_anno_MD_bucher_test_resRE.txt')))
})

test_that('Indirect meta-analysis of binary arm-level data is correct', {
  #load test data
  direct_input = dget('../data/smoking_anno_arm_direct.txt')
  #run direct meta-analysis and extract results. This should be valid because this has
  #it's own set of tests
  direct_res = doDirectMeta(df=direct_input, effectCode='OR',
                            dataType = 'binary', backtransf = TRUE)
  for (i in 1:length(direct_res)) {
    res = direct_res[[i]]
    df = extractDirectRes(
      metaRes = res, effect = 'Odds Ratio', backtransf = TRUE,
      intervention = res$label.e[1], comparator = res$label.c[1],
      interventionCode = res$e.code, comparatorCode = res$c.code
    )
    if (i == 1) {
      dirMA = df
    } else {
      dirMA = dplyr::bind_rows(dirMA, df)
    }
  }
  #run bucher analysis for the standard case
  indMA = doBucher(
    comparisons = direct_input[,1:4], direct = dirMA,
    effectType = 'all', continuous = FALSE,
    backtransf = TRUE
  )
  #test that the object has the right class and the expected content
  expect_is(indMA, 'data.frame')
  expect_that(indMA, equals(dget('../data/smoking_anno_arm_indirect_res_all.txt')))

  #test effect_type argument
  indMA = doBucher(
    comparisons = direct_input[,1:4], direct = dirMA,
    effectType = 'Fixed', continuous = FALSE,
    backtransf = TRUE
  )
  expect_false('Random' %in% indMA$Model)
  indMA = doBucher(
    comparisons = direct_input[,1:4], direct = dirMA,
    effectType = 'Random', continuous = FALSE,
    backtransf = TRUE
  )
  expect_false('Fixed' %in% indMA$Model)
})

test_that('Indirect meta-analysis of binary contrast data is correct', {
  #load test data
  direct_input = dget('../data/smoking_anno_diff_direct.txt')
  #run direct meta-analysis and extract results. This should be valid because this has
  #it's own set of tests
  direct_res = doDirectMeta(df=direct_input, effectCode='OR',
                            dataType = 'treatment difference', backtransf = TRUE)
  for (i in 1:length(direct_res)) {
    res = direct_res[[i]]
    df = extractDirectRes(
      metaRes = res, effect = 'Odds Ratio', backtransf = TRUE,
      intervention = res$label.e[1], comparator = res$label.c[1],
      interventionCode = res$e.code, comparatorCode = res$c.code
    )
    if (i == 1) {
      dirMA = df
    } else {
      dirMA = dplyr::bind_rows(dirMA, df)
    }
  }
  #run bucher analysis for the standard case
  indMA = doBucher(
    comparisons = direct_input[,1:4], direct = dirMA,
    effectType = 'all', continuous = FALSE,
    backtransf = TRUE
  )
  #test that the object has the right class and the expected content
  expect_is(indMA, 'data.frame')
  expect_that(indMA, equals(dget('../data/smoking_anno_diff_indirect_res_all.txt')))

  #test effect_type argument
  indMA = doBucher(
    comparisons = direct_input[,1:4], direct = dirMA,
    effectType = 'Fixed', continuous = FALSE,
    backtransf = TRUE
  )
  expect_false('Random' %in% indMA$Model)
  indMA = doBucher(
    comparisons = direct_input[,1:4], direct = dirMA,
    effectType = 'Random', continuous = FALSE,
    backtransf = TRUE
  )
  expect_false('Fixed' %in% indMA$Model)
})

test_that('Indirect meta-analysis of continuous arm-level data is correct', {
  #load test data
  direct_input = dget('../data/parkinsons_anno_continuous_arm_direct.txt')
  #run direct meta-analysis and extract results. This should be valid because this has
  #it's own set of tests
  direct_res = doDirectMeta(df=direct_input, effectCode='MD',
                            dataType = 'continuous', backtransf = FALSE)
  for (i in 1:length(direct_res)) {
    res = direct_res[[i]]
    df = extractDirectRes(
      metaRes = res, effect = 'Mean Difference', backtransf = FALSE,
      intervention = res$label.e[1], comparator = res$label.c[1],
      interventionCode = res$e.code, comparatorCode = res$c.code
    )
    if (i == 1) {
      dirMA = df
    } else {
      dirMA = dplyr::bind_rows(dirMA, df)
    }
  }
  #run bucher analysis for the standard case
  indMA = doBucher(
    comparisons = direct_input[,1:4], direct = dirMA,
    effectType = 'all', continuous = TRUE,
    backtransf = FALSE
  )
  #test that the object has the right class and the expected content
  expect_is(indMA, 'data.frame')
  expect_that(indMA, equals(dget('../data/parkinsons_anno_continuous_arm_indirect_res_all.txt')))

  #test effect_type argument
  indMA = doBucher(
    comparisons = direct_input[,1:4], direct = dirMA,
    effectType = 'Fixed', continuous = TRUE,
    backtransf = FALSE
  )
  expect_false('Random' %in% indMA$Model)
  indMA = doBucher(
    comparisons = direct_input[,1:4], direct = dirMA,
    effectType = 'Random', continuous = TRUE,
    backtransf = FALSE
  )
  expect_false('Fixed' %in% indMA$Model)
})

test_that('Indirect meta-analysis of continuous contrast data is correct', {
  #load test data
  direct_input = dget('../data/parkinsons_anno_continuous_diff_direct.txt')
  #run direct meta-analysis and extract results. This should be valid because this has
  #it's own set of tests
  direct_res = doDirectMeta(df=direct_input, effectCode='MD',
                            dataType = 'treatment difference', backtransf = FALSE)
  for (i in 1:length(direct_res)) {
    res = direct_res[[i]]
    df = extractDirectRes(
      metaRes = res, effect = 'Mean Difference', backtransf = FALSE,
      intervention = res$label.e[1], comparator = res$label.c[1],
      interventionCode = res$e.code, comparatorCode = res$c.code
    )
    if (i == 1) {
      dirMA = df
    } else {
      dirMA = dplyr::bind_rows(dirMA, df)
    }
  }
  #run bucher analysis for the standard case
  indMA = doBucher(
    comparisons = direct_input[,1:4], direct = dirMA,
    effectType = 'all', continuous = TRUE,
    backtransf = FALSE
  )
  #test that the object has the right class and the expected content
  expect_is(indMA, 'data.frame')
  expect_that(indMA, equals(dget('../data/parkinsons_anno_continuous_diff_indirect_res_all.txt')))

  #test effect_type argument
  indMA = doBucher(
    comparisons = direct_input[,1:4], direct = dirMA,
    effectType = 'Fixed', continuous = TRUE,
    backtransf = FALSE
  )
  expect_false('Random' %in% indMA$Model)
  indMA = doBucher(
    comparisons = direct_input[,1:4], direct = dirMA,
    effectType = 'Random', continuous = TRUE,
    backtransf = FALSE
  )
  expect_false('Fixed' %in% indMA$Model)
})