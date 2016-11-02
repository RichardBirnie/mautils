library(mautils)
context('formatDataToDirectMA')

test_that('Annotated binary arm level data is correctly reformatted',{
  #call format data function on test data
  input = formatDataToDirectMA(input.df = read.table(
    file = '../data/smoking_anno_arm.txt', header = TRUE, sep = '\t'
  ), dataType = 'binary')
  #compared to pre-stored output
  expect_that(input, equals(dget('../data/smoking_anno_arm_direct.txt')))
})

test_that('Annotated binary contrast data is correctly reformatted',{
  #call format data function on test data
  input = formatDataToDirectMA(input.df = read.table(
    file = '../data/smoking_anno_diff.txt', header = TRUE, sep = '\t'
  ), dataType = 'treatment difference')
  #compared to pre-stored output
  expect_that(input, equals(dget('../data/smoking_anno_diff_direct.txt')))
})

test_that('Annotated continuous arm level data is correctly reformatted',{
  #call format data function on test data
  input = formatDataToDirectMA(input.df = read.table(
    file = '../data/parkinsons_anno_continuous_arm.txt', header = TRUE, sep = '\t'
  ), dataType = 'continuous')
  #compared to pre-stored output
  expect_that(input, equals(dget('../data/parkinsons_anno_continuous_arm_direct.txt')))
})

test_that('Annotated continuous contrast data is correctly reformatted',{
  #call format data function on test data
  input = formatDataToDirectMA(input.df = read.table(
    file = '../data/parkinsons_anno_continuous_diff.txt', header = TRUE, sep = '\t'
  ), dataType = 'treatment difference')
  #compared to pre-stored output
  expect_that(input, equals(dget('../data/parkinsons_anno_continuous_diff_direct.txt')))
})