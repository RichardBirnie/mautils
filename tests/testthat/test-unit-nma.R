library(mautils)
context('Network meta-analysis')

test_that('Ranking of NMA results is correct', {
  #load test data. One example should be sufficient.
  #This should behave the same for any mtc.result object
  fe_res = readRDS('../data/smoking_anno_arm_nma_resFE.RData')

  #calculate ranks and test class of the result
  r = gemtc::rank.probability(fe_res, preferredDirection = 1)
  expect_is(r, 'mtc.rank.probability')

  #convert to a better format
  #check class, dimensions and content of the result
  ranks = extractRanks(ranks = r, treatments = fe_res$model$network$treatments)
  expect_is(ranks, 'data.frame')
  expect_true(nrow(ranks) == nrow(fe_res$model$network$treatments))
  expect_true(ncol(ranks) == (ncol(r) + ncol(fe_res$model$network$treatments)))
  expect_equal(ranks, dget('../data/smoking_anno_arm_nma_resFE_ranks.txt'))
})

test_that('Code is extracted correctly',{
  #load test data. One example should be sufficient.
  #This should behave the same for any mtc.result object
  fe_res = readRDS('../data/smoking_anno_arm_nma_resFE.RData')
  coda = extractCoda(mtcResults = fe_res, summarise = FALSE)
  expect_is(coda, 'data.frame')
  expect_equal(coda[1:2000,], dget('../data/smoking_anno_arm_nma_resFE_coda.txt'))
})

test_that('JAGS model code is annotated and exported correctly',{
  #load test data. One example should be sufficient.
  #This should behave the same for any mtc.result object
  fe_res = readRDS('../data/smoking_anno_arm_nma_resFE.RData')
  saveModelCode(fe_res, modelFile = '../data/smoking_anno_arm_nma_resFE_new_model.txt')
  ref = readLines('../data/smoking_anno_arm_nma_resFE_model.txt')
  new = readLines('../data/smoking_anno_arm_nma_resFE_new_model.txt')
  expect_identical(ref, new)
  file.remove('../data/smoking_anno_arm_nma_resFE_new_model.txt')
})

test_that('Extraction and formatting of NMA results is correct',{
  fe_res = readRDS('../data/smoking_anno_arm_nma_resFE.RData')
  pairwise_res = calcAllPairs(fe_res, expon = TRUE)
  expect_is(pairwise_res, 'data.frame')
  expect_equal(pairwise_res, dget('../data/smoking_anno_arm_nma_fe_pairwise.txt'))
})

