library(mautils)
context('Indirect meta-analysis')

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