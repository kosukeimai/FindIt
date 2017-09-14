rm(list=ls())
library(FindIt)
library(testthat)
context("tests FindIt")

test_that("tests FindIt on CausalANOVA", {
  # set random seed
  set.seed(12345)

  # load the data
  data(Carlson)
  # Specify the order of each factor
  Carlson$newRecordF<- factor(Carlson$newRecordF,ordered=TRUE,
                              levels=c("YesLC", "YesDis","YesMP",
                                       "noLC","noDis","noMP","noBusi"))
  Carlson$promise <- factor(Carlson$promise,ordered=TRUE,levels=c("jobs","clinic","education"))
  Carlson$coeth_voting <- factor(Carlson$coeth_voting,ordered=FALSE,levels=c("0","1"))
  Carlson$relevantdegree <- factor(Carlson$relevantdegree,ordered=FALSE,levels=c("0","1"))
  
  ######################################### 
  # Without Screening and Collapsing
  ######################################### 
  #################### only AMEs ####################
  fit1 <- CausalANOVA(formula=won ~ newRecordF + promise + coeth_voting + relevantdegree,
                      data=Carlson, pair.id=Carlson$contestresp, diff=TRUE,
                      cluster=Carlson$respcodeS, nway=1)
  x <- summary(fit1)
  expect_that(length(x), is_equivalent_to(5))
  expect_true("range.name" %in% names(x))
  expect_equal(x$range["coeth_voting", "range"], 0.065, tolerance = 0.001)
  expect_equal(x$AME["newRecordF3", "AME"], 0.093, tolerance = 0.001)
  expect_true(is.null(x$AMIE2))
  expect_true(is.null(x$AMIE3))
  
  plot(fit1)
})  
  


test_that("tests FindIt on ConditionalEffect", {
  # set random seed
  set.seed(12345)
  
  data(Carlson)
  # Specify the order of each factor
  Carlson$newRecordF<- factor(Carlson$newRecordF,ordered=TRUE,
                              levels=c("YesLC", "YesDis","YesMP",
                                       "noLC","noDis","noMP","noBusi"))
  Carlson$promise <- factor(Carlson$promise,ordered=TRUE,levels=c("jobs","clinic","education"))
  Carlson$coeth_voting <- factor(Carlson$coeth_voting,ordered=FALSE,levels=c("0","1"))
  Carlson$relevantdegree <- factor(Carlson$relevantdegree,ordered=FALSE,levels=c("0","1"))
  
  ######################################### 
  # Collapsing Without Screening
  ######################################### 
  #################### AMEs and two-way AMIEs ####################
  # We show a very small example for illustration.
  # Recommended to use cv.collapse.cost=c(0.1,0.3,0.5) and nfolds=10 in practice.
  x <- cv.CausalANOVA(formula=won ~ newRecordF + promise + coeth_voting + relevantdegree,
                           int2.formula = ~ newRecordF:coeth_voting,
                           data=Carlson, pair.id=Carlson$contestresp,diff=TRUE,
                           cv.collapse.cost=c(0.1,0.3), nfolds=2,
                           cluster=Carlson$respcodeS, nway=2)
  expect_that(length(x), is_equivalent_to(5))
  expect_true("cv.each.mat" %in% names(x))
  expect_equal(x$cv.min, 0.3, tolerance = 0.05)
  expect_equal(x$cv.1Std, 0.1, tolerance = 0.05)
  expect_equal(as.numeric(x$cv.error[1]), 0.433, tolerance = 0.001)
  expect_equal(as.numeric(x$cv.each.mat[1, "0.1"]), 0.422, tolerance = 0.001)
  expect_equal(as.numeric(x$cv.cost[1]), 0.1, tolerance = 0.05)
})  

test_that("tests FindIt on SVMHet", {
  # set random seed
  set.seed(12345)
  
  ################################################### 
  # Example 1: Treatment-Covariate Interaction
  ################################################### 
  data(LaLonde)
  
  # The model includes a treatment variable, 
  # nine covariates to be interacted with the treatment variable,
  # and the same nine covariates to be adjusted.
  
  # Run to find the LASSO parameters
  x  <-FindIt(model.treat= outcome ~ treat,
               model.main= ~ age+educ+black+hisp+white+
               marr+nodegr+log.re75+u75,
               model.int= ~ age+educ+black+hisp+white+
               marr+nodegr+log.re75+u75,
               data = LaLonde,  
               type="binary",
               treat.type="single") 
  expect_that(length(x), is_equivalent_to(40))
  expect_true("Treatment.version" %in% names(x))
  expect_equal(x$coefs[8], -0.147, tolerance = 0.001)
  expect_equal(x$ATE, 0.079, tolerance = 0.001)
  expect_equal(as.numeric(x$frame.meanC["marr1"]), 0.162, tolerance = 0.001)

  # set random seed
  set.seed(12345)

  # Fit with uncovered lambda parameters.
  F1  <-FindIt(model.treat= outcome ~ treat,
               model.main= ~ age+educ+black+hisp+white+
               marr+nodegr+log.re75+u75,
               model.int= ~ age+educ+black+hisp+white+
               marr+nodegr+log.re75+u75,
               data = LaLonde, 
               type="binary",
               treat.type="single",
               search.lambdas=FALSE,
               lambdas = c(-3.8760,-4.0025) )
  
  x <- summary(F1)
  expect_that(length(x), is_equivalent_to(3))
  expect_true("GCV" %in% names(x))
  expect_equal(as.numeric(x$coefficients["white1"]), 0.00577, tolerance = 0.00001)
  expect_equal(x$GCV[1], 0.9857783, tolerance = 0.00001)
  expect_equal(x$misclass[2], 0.3185596, tolerance = 0.00001)
  
  # Returns all the estimated treatment effects. 
  x <- predict(F1)
  expect_that(length(x), is_equivalent_to(5))
  expect_true("coefs" %in% names(x))
  expect_equal(x$ATE, 0.07903784, tolerance = 0.00001)
  expect_equal(x$data["10", "Treatment.effect"], 0.1787417, tolerance = 0.00001)
  expect_equal(x$data["567", "log.re75"], 10.52945, tolerance = 0.01)
  expect_equal(as.numeric(x$coefs["log.re75.2"]), -0.3272249, tolerance = 0.001)
  expect_equal(as.numeric(x$orig.coef["u751:age"]), 0.007013081, tolerance = 0.00001)
})  



test_that("tests FindIt on test.CausalANOVA", {
  # set random seed
  set.seed(12345)
  
  ######################################## 
  # With Screening and Collapsing
  ########################################
  data(Carlson)
  # Specify the order of each factor
  Carlson$newRecordF<- factor(Carlson$newRecordF,ordered=TRUE,
                              levels=c("YesLC", "YesDis","YesMP",
                                       "noLC","noDis","noMP","noBusi"))
  Carlson$promise <- factor(Carlson$promise,ordered=TRUE,levels=c("jobs","clinic","education"))
  Carlson$coeth_voting <- factor(Carlson$coeth_voting,ordered=FALSE,levels=c("0","1"))
  Carlson$relevantdegree <- factor(Carlson$relevantdegree,ordered=FALSE,levels=c("0","1"))

  # Sample Splitting
  train.ind <- sample(unique(Carlson$respcodeS), 272, replace=FALSE)
  test.ind <- setdiff(unique(Carlson$respcodeS), train.ind)
  Carlson.train <- Carlson[is.element(Carlson$respcodeS,train.ind), ]
  Carlson.test <- Carlson[is.element(Carlson$respcodeS,test.ind), ]
   
  #################### AMEs and two-way AMIEs ####################
  fit.r2 <- CausalANOVA(formula=won ~ newRecordF + promise + coeth_voting + relevantdegree,
                        data=Carlson.train, pair.id=Carlson.train$contestresp,diff=TRUE,
                        screen=TRUE, collapse=TRUE,
                        cluster=Carlson.train$respcodeS, nway=2)
  x <- summary(fit.r2)
  expect_that(length(x), is_equivalent_to(5))
  expect_true("range.name" %in% names(x))
  expect_equal(x$range["coeth_voting", "range"], 0.065, tolerance = 0.001)
  expect_equal(x$AME["newRecordF3", "AME"], 0.1351, tolerance = 0.001)
  expect_equal(x$AMIE2[10, "AMIE"], 0.0624, tolerance = 0.001)
  expect_true(is.null(x$AMIE3))
  
  # refit with test.CausalANOVA
  fit.r2.new <- test.CausalANOVA(fit.r2, newdata=Carlson.test, diff=TRUE,
                                 pair.id=Carlson.test$contestresp, cluster=Carlson.test$respcodeS)
  
  x <- summary(fit.r2.new)
  expect_that(length(x), is_equivalent_to(5))
  expect_true("range.name" %in% names(x))
  expect_equal(x$range["coeth_voting", "range"], 0.024, tolerance = 0.001)
  expect_equal(x$AME["newRecordF3", "AME"], -0.0247, tolerance = 0.001)
  expect_equal(x$AMIE2[10, "AMIE"], -0.0651, tolerance = 0.001)
  expect_true(is.null(x$AMIE3))
  
  plot(fit.r2.new)
  plot(fit.r2.new, type="ConditionalEffect", fac.name=c("newRecordF","coeth_voting"))
})  
