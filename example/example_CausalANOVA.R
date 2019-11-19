data(Carlson)
## Specify the order of each factor
Carlson$newRecordF<- factor(Carlson$newRecordF,ordered=TRUE,
                            levels=c("YesLC", "YesDis","YesMP",
                                     "noLC","noDis","noMP","noBusi"))
Carlson$promise <- factor(Carlson$promise,ordered=TRUE,levels=c("jobs","clinic","education"))
Carlson$coeth_voting <- factor(Carlson$coeth_voting,ordered=FALSE,levels=c("0","1"))
Carlson$relevantdegree <- factor(Carlson$relevantdegree,ordered=FALSE,levels=c("0","1"))

## #######################################
## Without Screening and Collapsing
## #######################################
#################### only AMEs ####################
fit1 <- CausalANOVA(formula=won ~ newRecordF + promise + coeth_voting + relevantdegree,
                    data=Carlson, pair.id=Carlson$contestresp, diff=TRUE,
                    cluster=Carlson$respcodeS, nway=1)
summary(fit1) ## OK

#################### AMEs and two-way AMIEs ####################
fit2 <- CausalANOVA(formula=won ~ newRecordF + promise + coeth_voting + relevantdegree,
                    int2.formula = ~ newRecordF:coeth_voting,
                    data=Carlson, pair.id=Carlson$contestresp,diff=TRUE,
                    cluster=Carlson$respcodeS, nway=2)
summary(fit2) ## OK 
plot(fit2, type="ConditionalEffect", fac.name=c("newRecordF", "coeth_voting"))
ConditionalEffect(fit2, treat.fac="newRecordF", cond.fac="coeth_voting")

\dontrun{
#################### AMEs and two-way and three-way AMIEs ####################
## Note: All pairs within thee-way interactions should show up in int2.formula (Strong Hierarchy).
fit3 <- CausalANOVA(formula=won ~ newRecordF + promise + coeth_voting + relevantdegree,
                    int2.formula = ~ newRecordF:promise + newRecordF:coeth_voting
                                       + promise:coeth_voting,
                    int3.formula = ~ newRecordF:promise:coeth_voting,
                    data=Carlson, pair.id=Carlson$contestresp,diff=TRUE,
                    cluster=Carlson$respcodeS, nway=3)
summary(fit3) ## OK
plot(fit3, type="AMIE", fac.name=c("newRecordF","promise", "coeth_voting"),space=25,adj.p=2.2)
}

## #######################################
## With Screening and Collapsing
## #######################################
## Sample Splitting
train.ind <- unique(Carlson$respcodeS)[1:272]
test.ind <- setdiff(unique(Carlson$respcodeS), train.ind)
Carlson.train <- Carlson[is.element(Carlson$respcodeS,train.ind), ]
Carlson.test <- Carlson[is.element(Carlson$respcodeS,test.ind), ]

#################### AMEs and two-way AMIEs ####################
fit.r2 <- CausalANOVA(formula=won ~ newRecordF + promise + coeth_voting + relevantdegree,
                      data=Carlson.train, pair.id=Carlson.train$contestresp,diff=TRUE,
                      screen = TRUE, collapse = TRUE,
                      cluster=Carlson.train$respcodeS, nway=2)
summary(fit.r2) ## OK 

## refit with test.CausalANOVA
fit.r2.new <- test.CausalANOVA(fit.r2, newdata=Carlson.test, diff=TRUE,
                               pair.id=Carlson.test$contestresp, cluster=Carlson.test$respcodeS)

summary(fit.r2.new) ## OK 
plot(fit.r2.new)
plot(fit.r2.new, type="ConditionalEffect", fac.name=c("newRecordF","coeth_voting"))
ConditionalEffect(fit.r2.new, treat.fac="newRecordF", cond.fac="coeth_voting")