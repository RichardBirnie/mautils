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

  pairwise_res = extractComparison(pairwise_res)
  expect_is(pairwise_res, 'data.frame')
  expect_is(pairwise_res$tA, 'integer')
  expect_is(pairwise_res$tB, 'integer')

  pairwise_res = nameTreatments(pairwise_res, coding = fe_res$model$network$treatments)
  expect_is(pairwise_res, 'data.frame')
  expect_is(pairwise_res$TreatmentA, 'integer')
  expect_is(pairwise_res$TreatmentB, 'integer')
  expect_is(pairwise_res$nameA, 'character')
  expect_is(pairwise_res$nameB, 'character')
  sapply(pairwise_res[,5:11], expect_is, 'numeric')
  expect_is(pairwise_res$Sample, 'integer')

  #get treatment coding from results object
  #simulate a custom order to test sorting of results
  treatments = fe_res$model$network$treatments
  treatments$id = as.integer(treatments$id)
  treatments$Order = rev(treatments$id)
  fnames = treatments$description
  rnames = rev(fnames)
  #run function and test output
  tab = makeTab(results = pairwise_res, coding = treatments, reportOrder = 'default')
  #check class and dimensions of output
  expect_is(tab, 'data.frame')
  expect_true(nrow(tab) == nrow(treatments))
  expect_true(ncol(tab) == nrow(treatments))

  #check ordering of output
  expect_true(all(rownames(tab) == fnames))
  expect_true(all(colnames(tab) == fnames))

  #test custom ordering
  tab = makeTab(results = pairwise_res, coding = treatments, reportOrder = 'custom')
  #check class and dimensions of output
  expect_is(tab, 'data.frame')
  expect_true(nrow(tab) == nrow(treatments))
  expect_true(ncol(tab) == nrow(treatments))

  #check ordering of output
  expect_true(all(rownames(tab) == rnames))
  expect_true(all(colnames(tab) == rnames))
})

test_that('Model comparison stats are correctly extracted',{
  #load test data. One FE and one RE as these have diff model fit stats
  fe_res = readRDS('../data/smoking_anno_arm_nma_resFE.RData')
  re_res = readRDS('../data/smoking_anno_arm_nma_resRE.RData')
  #test for FE model
  fe_fit = extractModelFit(fe_res)
  expect_is(fe_fit, 'data.frame')
  expect_true(nrow(fe_fit) == 3)
  expect_true(ncol(fe_fit) == 2)
  expect_true(all(is.na(fe_fit[,2])))
  #test for RE model
  re_fit = extractModelFit(re_res)
  expect_is(re_fit, 'data.frame')
  expect_true(nrow(re_fit) == 4)
  expect_true(ncol(re_fit) == 2)
  expect_false(all(is.na(re_fit[,2])))
})

