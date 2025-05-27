#!/usr/bin/env Rscript
knitr::opts_chunk$set(echo = TRUE)
library(tidyverse)
library(readr)
library(rpart)
library(partykit)
library(rDNAse)
library(ROCR)
df <- read_csv("Anabaena7120.csv")

df$Annotation <- df$`Locus Name`
fox1 <- read_csv("FoxGenePositives.csv")
fox1 <- fox1 %>% mutate(Annotation = GENES)

df_tree2 <- df #as.data.frame(cbind(df, kmer, kmer2, kmer4))

negative <- read_csv("FoxGeneNegatives.csv")
negative$Annotation <- negative$GENES

unlabeled <- read_csv("UnlabeledFoxGenes.csv")
unlabeled$Annotation <- unlabeled$GENES

negative <- anti_join(negative, unlabeled, by = "Annotation")



f <- semi_join(df_tree2, fox1, by = "Annotation")
f2 <- semi_join(df_tree2, negative, by = "Annotation")
f$FOX <- 1L
f2$FOX <- 0L



df_tree2 <- bind_rows(f, f2)


is.na(df_tree2) <- sapply(df_tree2, is.infinite)
is.na(df_tree2) <- sapply(df_tree2, is.nan)




df_tree2 <- df_tree2 %>%
  group_by(Annotation) %>%
  dplyr::summarise(across(where(is.numeric), sum, na.rm = T))







df_tree2 <- df_tree2 %>% ungroup  %>% select(Annotation, where(is.numeric))

df_tree3 <- df_tree2 %>% select(-Annotation)


library(rpart)
library(partykit)

library(caret)
set.seed(123)
 is.na(df_tree3) <- sapply(df_tree3, is.infinite)
 is.na(df_tree3) <- sapply(df_tree3, is.nan)

 
 df_tree3 <- as.data.frame(df_tree3)
 for(i in 1:ncol(df_tree3)) {
  df_tree3[ , i][is.na(df_tree3[ , i])] <- median(df_tree3[ , i], na.rm = TRUE)
}



# create a list of 80% of the rows in the original dataset we can use for training
validation_index <- createDataPartition(df_tree3$FOX, p= 0.8, list=FALSE)
# select 20% of the data for validation
validation <- df_tree3[-validation_index,]

# use the remaining 80% of data to training and testing the models
dataset <- df_tree3[validation_index,]
library(mlbench)
# Run algorithms using 10-fold cross validation


control <- trainControl(method = "repeatedcv",
                          number = 10,
                          repeats = 3,
                          classProbs = TRUE,
                          sampling = "smote",
                          summaryFunction = twoClassSummary,
                          search = "grid")
library(e1071)

dataset$Annotation <- NULL
dataset$FOX <- ifelse(dataset$FOX == 0, "FOX0", "FOX1")
validation$FOX <- ifelse(validation$FOX == 0, "FOX0", "FOX1")

set.seed(7)
fit.glm <- train(as.factor(FOX)~ `Chromosome region end`    + `Fold Change in RPKM 0 to 21 hours`, data=dataset, method='glm')
summary(fit.glm)

#dataset$FOX <- ifelse(dataset$FOX == 1, "FOX", "NOT")
set.seed(7)
fit.lda <- train(as.factor(FOX)~., data=dataset, method="lda",  trControl=control)
# CART
set.seed(7)
fit.cart <- train(as.factor(FOX)~ ., data=dataset, method="ctree",  trControl=control)


set.seed(7)
fit.mars <- train(as.factor(FOX)~., data=dataset, method="earth",  trControl=control)


 

# Random Forest
set.seed(7)
fit.rf <- train(as.factor(FOX)~., data=dataset, method="rf",  trControl=control)
#

set.seed(7)
fit.xgb <- caret::train(as.factor(FOX)~ ., data=dataset, method='xgbTree',  trControl=control)






library(party)
plot(fit.cart$finalModel, type="simple",           # no terminal plots
  inner_panel=node_inner(fit.cart$finalModel,
       abbreviate = TRUE,            # short variable names
       pval = FALSE,                 # no p-values
       id = FALSE),                  # no id of node
  terminal_panel=node_terminal(fit.cart$finalModel, 
       abbreviate = TRUE,
       digits = 1,                   # few digits on numbers
       fill = c("white"),            # make box white not grey
       id = FALSE)
   )

summary(fit.glm$finalModel)
summary(fit.mars$finalModel)

library(ROCR)
mars_tr <- prediction(predict(fit.mars, newdata = dataset, type = "prob")[,2], dataset$FOX)
mars_tr_auc = ROCR::performance(mars_tr, "auc")
mars_ts <- prediction(predict(fit.mars, newdata = validation, type = "prob")[,2], validation$FOX)
mars_ts_auc = ROCR::performance(mars_ts, "auc")

lda_tr <- prediction(predict(fit.lda, newdata = dataset, type = "prob")[,2], dataset$FOX)
lda_tr_auc = ROCR::performance(lda_tr, "auc")
lda_ts <- prediction(predict(fit.lda, newdata = validation, type = "prob")[,2], validation$FOX)
lda_ts_auc = ROCR::performance(lda_ts, "auc")

cart_tr <- prediction(predict(fit.cart, newdata = dataset, type = "prob")[,2], dataset$FOX)
cart_tr_auc = ROCR::performance(cart_tr, "auc")
cart_ts <- prediction(predict(fit.cart, newdata = validation, type = "prob")[,2], validation$FOX)
cart_ts_auc = ROCR::performance(cart_ts, "auc")

rf_tr <- prediction(predict(fit.rf, newdata = dataset, type = "prob")[,2], dataset$FOX)
rf_tr_auc = ROCR::performance(rf_tr, "auc")
rf_ts <- prediction(predict(fit.rf, newdata = validation, type = "prob")[,2], validation$FOX)
rf_ts_auc = ROCR::performance(rf_ts, "auc")

xgb_tr <- prediction(predict(fit.xgb, newdata = dataset, type = "prob")[,2], dataset$FOX)
xgb_tr_auc = ROCR::performance(xgb_tr, "auc")
xgb_ts <- prediction(predict(fit.xgb, newdata = validation, type = "prob")[,2], validation$FOX)
xgb_ts_auc = ROCR::performance(xgb_ts, "auc")




glm_tr <- prediction(predict(fit.glm, newdata = dataset, type = "prob")[,2], dataset$FOX)
glm_tr_auc = ROCR::performance(glm_tr, "auc")
glm_ts <- prediction(predict(fit.glm, newdata = validation, type = "prob")[,2], validation$FOX)
glm_ts_auc = ROCR::performance(glm_ts, "auc")

ensemble_tr <- (predict(fit.xgb, newdata = dataset, type = "prob")[,2] + predict(fit.rf, newdata = dataset, type = "prob")[,2] + predict(fit.glm, newdata = dataset, type = "prob")[,2])/3

ensemble_ts <- ( predict(fit.xgb, newdata = validation, type = "prob")[,2] + predict(fit.rf, newdata = validation, type = "prob")[,2] + predict(fit.glm, newdata = validation, type = "prob")[,2]  )/3

ensemble_auc <- prediction(ensemble_tr, dataset$FOX)
ens_tr_auc = ROCR::performance(ensemble_auc, "auc")

ensemble_auc <- prediction(ensemble_ts, validation$FOX)
ens_ts_auc = ROCR::performance(ensemble_auc, "auc")

mars <- ifelse(predict(fit.mars, newdata = validation, type = "prob")[,2] < 0.5, "FOX0", "FOX1")
lda <- ifelse(predict(fit.lda, newdata = validation, type = "prob")[,2] < 0.5, "FOX0", "FOX1")
cart <- ifelse(predict(fit.cart, newdata =validation, type = "prob")[,2] < 0.5,"FOX0", "FOX1")
rf <- ifelse(predict(fit.rf, newdata = validation, type = "prob")[,2] < 0.5, "FOX0", "FOX1")
xgb <- ifelse(predict(fit.xgb, newdata = validation, type = "prob")[,2] < 0.5, "FOX0", "FOX1")
glm <- ifelse(predict(fit.glm, newdata = validation, type = "prob")[,2] < 0.5, "FOX0", "FOX1")
ens <- ifelse(ensemble_ts < 0.5, "FOX0", "FOX1")

library(caret)

mars_precision <- posPredValue(as.factor(mars), as.factor(validation$FOX), positive= "FOX1")
mars_recall <- sensitivity(as.factor(mars), as.factor(validation$FOX), positive="FOX1")
mars_F1 <- (2 * mars_precision * mars_recall) / (mars_precision + mars_recall)

lda_precision <- posPredValue(as.factor(lda), as.factor(validation$FOX), positive="FOX1")
lda_recall <- sensitivity(as.factor(lda), as.factor(validation$FOX), positive="FOX1")
lda_F1 <- (2 * lda_precision * lda_recall) / (lda_precision + lda_recall)

rf_precision <- posPredValue(as.factor(rf), as.factor(validation$FOX), positive="FOX1")
rf_recall <- sensitivity(as.factor(rf), as.factor(validation$FOX), positive="FOX1")
rf_F1 <- (2 * rf_precision * rf_recall) / (rf_precision + rf_recall)

cart_precision <- posPredValue(as.factor(cart), as.factor(validation$FOX), positive="FOX1")
cart_recall <- sensitivity(as.factor(cart), as.factor(validation$FOX), positive="FOX1")
cart_F1 <- (2 * cart_precision * cart_recall) / (cart_precision + cart_recall)

xgb_precision <- posPredValue(as.factor(xgb), as.factor(validation$FOX), positive="FOX1")
xgb_recall <- sensitivity(as.factor(xgb), as.factor(validation$FOX), positive="FOX1")
xgb_F1 <- (2 * xgb_precision * xgb_recall) / (xgb_precision + xgb_recall)




glm_precision <- posPredValue(as.factor(glm), as.factor(validation$FOX), positive="FOX1")
glm_recall <- sensitivity(as.factor(glm), as.factor(validation$FOX), positive="FOX1")
glm_F1 <- (2 * glm_precision * glm_recall) / (glm_precision + glm_recall)


ens_precision <- posPredValue(as.factor(ens), as.factor(validation$FOX), positive="FOX1")
ens_recall <- sensitivity(as.factor(ens), as.factor(validation$FOX), positive="FOX1")
ens_F1 <- (2 * ens_precision * ens_recall) / (ens_precision + ens_recall)

Train_AUC = rbind(mars_tr_auc@y.values[[1]], cart_tr_auc@y.values[[1]], lda_tr_auc@y.values[[1]], rf_tr_auc@y.values[[1]], xgb_tr_auc@y.values[[1]],  glm_tr_auc@y.values[[1]], ens_tr_auc@y.values[[1]])

Test_AUC = rbind(mars_ts_auc@y.values[[1]], cart_ts_auc@y.values[[1]], lda_ts_auc@y.values[[1]], rf_ts_auc@y.values[[1]], xgb_ts_auc@y.values[[1]], glm_ts_auc@y.values[[1]], ens_ts_auc@y.values[[1]])

Test_Recall <- rbind(mars_recall, cart_recall, lda_recall, rf_recall, xgb_recall, glm_recall, ens_recall)

Test_precision <- rbind(mars_precision, cart_precision, lda_precision, rf_precision, xgb_precision, glm_precision, ens_precision)

Test_F1 <- rbind(mars_F1, cart_F1, lda_F1, rf_F1, xgb_F1, glm_F1, ens_F1)

Model <- c("MARS", "CART", "LDA", "RF", "XGBoost",  "GLM", "ENS")

performance <- as.data.frame(cbind(Train_AUC, Test_AUC, Test_Recall, Test_precision, Test_F1))
performance[,1:5] <- round(performance[,1:5], 2)

perf <- cbind(Model, performance)
colnames(perf)[2:6] <- c("Train_AUC", "Test_AUC", "Test_Recall", "Test_precision", "Test_F1")
perf
pred2 <- predict(fit.lda, newdata = validation, type = "prob")[,2]
pred3 <- predict(fit.rf, newdata = validation, type = "prob")[,2]
pred4 <- predict(fit.xgb, newdata = validation, type = "prob")[,2]
pred5 <- predict(fit.mars, newdata = validation, type = "prob")[,2]
pred6 <- predict(fit.cart, newdata = validation, type = "prob")[,2]
pred7 <- predict(fit.glm, newdata = validation, type = "prob")[,2]

#predict(fit.knn, dataset, probability = T)




#z <-cbind(pred[,2], df_tree2$FOX)
z <-cbind(pred4, as.numeric(validation$FOX))
z <- as.data.frame(z)

library(gains)
par(mfrow=c(1,1.5))
par(mar=c(5, 4, 4, 6) + 0.1)
dt2 <- gains( z$V2, z$pred, groups=10, optimal=F)

plot(dt2$depth, dt2$cume.lift, type="l", ylab="Cumulative lift", xlab="Rank Buckets", main = "Cumulative Lift and Response: XGBoost Model")
par(new = TRUE)
plot(dt2$depth, dt2$cume.pct.of.total, type = "l",col="red", axes = FALSE, bty = "n", xlab = "", ylab = "")
axis(side=4, at = pretty(range(dt2$cume.pct.of.total)))
mtext("Cumulative Response",side=4,col="red",line=2) 
legend("right",legend=c("Cumulative \n Lift","Cumulative \n Response"),
       text.col=c("black","red"),pch=c(16,16), col=c("black","red"))
lines(dt2$depth, dt2$cume.pct.of.total, type="l",col="red")
dt2
library(ROCR)
summary(fit.knn)

p2 <- predict(fit.xgb, dataset, type = 'prob')[,2]

p2 <- as.numeric(p2)
p3 <- predict(fit.xgb, validation, type = 'prob')[,1]

p3 <- as.numeric(p3)
pred1 <- prediction(p2, dataset$FOX)
roc.perf = ROCR::performance(pred1, measure = "tpr", x.measure = "fpr")
pred2 <- prediction(p3, validation$FOX)
roc.perf2 = ROCR::performance(pred2, measure = "tpr", x.measure = "fpr")
plot(roc.perf,col='red', lty=1, lwd=3, main = "Random Forest Model")
abline(a=0, b= 1)
plot(roc.perf2, add=TRUE, lty=1, lwd=3)
#roc.perf
auctest <- ROCR::performance(pred2,"auc")
auctrain <- ROCR::performance(pred1,"auc")
# now converting S4 class to vector
auctest <- unlist(slot(auctest, "y.values"))
auctrain <- unlist(slot(auctrain, "y.values"))
# adding min and max ROC AUC to the center of the plot
auctest<-mean(round(auctest, digits = 3))
auctrain<-mean(round(auctrain, digits = 3))
minauct <- paste(c("Train (AUC)  = "), auctrain,sep="")
maxauct <- paste(c("Test (AUC) = "),auctest,sep="")
legend(0.62,0.6,c(maxauct),cex=1.2,box.col = "white", text.col = "black")
legend(0.61,0.4,c(minauct),cex=1.2,box.col = "white", text.col = "red")

library(xgboost)
#xgb.importance(colnames(dataset, do.NULL = TRUE, prefix = "col"), model = fit.xgb)

 importance <- varImp(fit.xgb, scale=FALSE)
# # summarize importance
 print(importance)
# # plot importance
 plot(importance)

summary(fit.mars)
plot(fit.cart$finalModel)
library(readr)
library(tidyverse)

df <- read_csv("Anabaena7120.csv")
df$Annotation <- df$`Locus Name`
negative <- read_csv("FoxGeneNegatives.csv")
negative$Annotation <- negative$GENES

unlabeled <- read_csv("UnlabeledFoxGenes.csv")
unlabeled$Annotation <- unlabeled$GENES

df_tree2 <- semi_join( negative, unlabeled, by = "Annotation")

df_tree2 <- semi_join( df, df_tree2, by = "Annotation")

is.na(df_tree2) <- sapply(df_tree2, is.infinite)
is.na(df_tree2) <- sapply(df_tree2, is.nan)





df_tree2 <- df_tree2 %>%
  group_by(Annotation) %>%
  dplyr::summarise(across(where(is.numeric), sum, na.rm = T))







df_tree2 <- df_tree2 %>% ungroup  %>% select(Annotation, where(is.numeric))

df_tree3 <- df_tree2

library(caret)
set.seed(123)
 is.na(df_tree3) <- sapply(df_tree3, is.infinite)
 is.na(df_tree3) <- sapply(df_tree3, is.nan)
#  df_tree3[is.na(df_tree3)] <- 0
 
 df_tree3 <- as.data.frame(df_tree3)
 for(i in 2:ncol(df_tree3)) {
  df_tree3[ , i][is.na(df_tree3[ , i])] <- median(df_tree3[ , i], na.rm = TRUE)
}

 
 
 library(e1071)


# 
df_tree3$LDA_PRED <- predict(fit.lda, newdata = df_tree3, type = "prob")[,2]
df_tree3$RF_PRED <- predict(fit.rf, newdata = df_tree3, type = "prob")[,2]
df_tree3$XGB_PRED <- predict(fit.xgb, newdata = df_tree3, type = "prob")[,2]
df_tree3$MARS_PRED <- predict(fit.mars, newdata = df_tree3, type = "prob")[,2]
df_tree3$CART_PRED <- predict(fit.cart, newdata = df_tree3, type = "prob")[,2]
df_tree3$GLM_PRED <- predict(fit.glm, newdata = df_tree3, type = "prob")[,2]

df_tree3$ENS_PRED <- (predict(fit.glm, newdata = df_tree3, type = "prob")[,2] + predict(fit.xgb, newdata = df_tree3, type = "prob")[,2] + predict(fit.rf, newdata = df_tree3, type = "prob")[,2])/3

write_csv(df_tree3, "Predictions.csv")

