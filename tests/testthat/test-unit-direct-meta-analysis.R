library(mautils)
context('Direct meta-analysis')

test_that('Direct meta-analysis of binary arm level data is correct',{
  #do direct meta-analysis of test data.
  res = doDirectMeta(df=dget('../data/smoking_anno_arm_direct.txt'), effectCode='OR',
             dataType = 'binary', backtransf = TRUE)
  #check the result is a list and it is not empty
  expect_is(res, 'list')
  expect_true(length(res) > 0)
  #each element of res should be a list
  expect_true(all(unlist(lapply(res, is.list))))
  #All results should be of class metabin
  for (i in 1:length(res)) {
    expect_is(res[[i]], 'metabin')
  }

  #test that the extraction and reformatting of results
  res = res[[1]]
  df = extractDirectRes(
    metaRes = res, effect = 'Odds Ratio', backtransf = TRUE,
    intervention = res$label.e[1], comparator = res$label.c[1],
    interventionCode = res$e.code, comparatorCode = res$c.code
  )
  #test that the object has the right class and the expected content
  expect_is(df, 'data.frame')
  expect_that(df, equals(dget('../data/smoking_anno_arm_direct_res.txt')))
})

test_that('Direct meta-analysis of binary contrast data is correct',{
  #do direct meta-analysis of test data.
  res = doDirectMeta(df=dget('../data/smoking_anno_diff_direct.txt'), effectCode='OR',
                     dataType = 'treatment difference', backtransf = TRUE)
  #check the result is a list and it is not empty
  expect_is(res, 'list')
  expect_true(length(res) > 0)
  #each element of res should be a list
  expect_true(all(unlist(lapply(res, is.list))))
  #All results should be of class metagen
  for (i in 1:length(res)) {
    expect_is(res[[i]], 'metagen')
  }

  #test that the extraction and reformatting of results
  res = res[[1]]
  df = extractDirectRes(
    metaRes = res, effect = 'Odds Ratio', backtransf = TRUE,
    intervention = res$label.e[1], comparator = res$label.c[1],
    interventionCode = res$e.code, comparatorCode = res$c.code
  )
  #test that the object has the right class and the expected content
  expect_is(df, 'data.frame')
  expect_that(df, equals(dget('../data/smoking_anno_diff_direct_res.txt')))
})

test_that('Direct meta-analysis of continuous arm level data is correct',{
  #do direct meta-analysis of test data.
  res = doDirectMeta(df=dget('../data/parkinsons_anno_continuous_arm_direct.txt'), effectCode='MD',
                     dataType = 'continuous', backtransf = FALSE)
  #check the result is a list and it is not empty
  expect_is(res, 'list')
  expect_true(length(res) > 0)
  #each element of res should be a list
  expect_true(all(unlist(lapply(res, is.list))))
  #All results should be of class metacont
  for (i in 1:length(res)) {
    expect_is(res[[i]], 'metacont')
  }

  #test that the extraction and reformatting of results
  res = res[[1]]
  df = extractDirectRes(
    metaRes = res, effect = 'Mean Difference', backtransf = FALSE,
    intervention = res$label.e[1], comparator = res$label.c[1],
    interventionCode = res$e.code, comparatorCode = res$c.code
  )
  #test that the object has the right class and the expected content
  expect_is(df, 'data.frame')
  expect_that(df, equals(dget('../data/parkinsons_anno_continuous_arm_direct_res.txt')))
})

test_that('Direct meta-analysis of continuous contrast data is correct',{
  #do direct meta-analysis of test data.
  res = doDirectMeta(df=dget('../data/parkinsons_anno_continuous_diff_direct.txt'), effectCode='MD',
                     dataType = 'treatment difference', backtransf = FALSE)
  #check the result is a list and it is not empty
  expect_is(res, 'list')
  expect_true(length(res) > 0)
  #each element of res should be a list
  expect_true(all(unlist(lapply(res, is.list))))
  #All results should be of class metagen
  for (i in 1:length(res)) {
    expect_is(res[[i]], 'metagen')
  }

  #test that the extraction and reformatting of results
  res = res[[1]]
  df = extractDirectRes(
    metaRes = res, effect = 'Mean Difference', backtransf = FALSE,
    intervention = res$label.e[1], comparator = res$label.c[1],
    interventionCode = res$e.code, comparatorCode = res$c.code
  )
  #test that the object has the right class and the expected content
  expect_is(df, 'data.frame')
  expect_that(df, equals(dget('../data/parkinsons_anno_continuous_diff_direct_res.txt')))
})