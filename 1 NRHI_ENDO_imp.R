#Project: Endotype Project
#
# Purpose: Imputation of lung function with XGBoost
# Version: 1.1 
# Date: 26/05/22
# Author:HPF
# Intro
# This is a kernel aims for minimum data wrangling. 
# The data preparation workflow is summarized as following:
# 1. Fix typos, problematic names and values
# 2. Create new features
# 3. Perform knn imputation on missing values (both numerical and categorical variables)
# 4. One-hot encoding for categorical variables
# 5. Normalization by BoxCox transformation to reduce skewness
# 6. Standarize (center and scale)
# 7. Remove zero variance variables
# 8. Remove highly correlated variables
# 9. Log transform the outcome variables

# Input: MWDataBIO_impLungF_406x49_260522.csv
#       
# Output: All ppFVC Impute 
#
# Add:
#
# Dependencies: F/:NHRI 2 Imp FVC 250522
#
# Notes:

# =    1 working space =========================================================

setwd("F:/")
#if (! dir.exists("Project_IPF")) dir.create("Project_IPF")
list.files("F:/NHRI 3 imp Bio 260522")
setwd("F:/NHRI 3 imp Bio 260522")
list.files("F:/NHRI 3 imp Bio 260522")

# = 2 packages need ============================================================

#analytical pack
require(data.table)
require(skimr)
require(plyr) #revalue
require(caret)
require(caretEnsemble)
require(xgboost)
require(kernlab)
require(Matrix)
require(recipes)
library(tidyverse)
library(Metrics)
#Plotting and color options packages
library(gplots)
library(ggplot2)
library(scales)
library(RColorBrewer)
library(grid)
library(gridExtra)
library(ggpubr)

# 3 open dataset ===============================================================

MWD<- read.csv("MWDataBIO_impLungF_463x49_070622.csv",header=TRUE, sep=",", stringsAsFactors=FALSE, na.strings=c("","NA"))
str(MWD)
dim(MWD)
MWD$X=NULL
rownames(MWD)<-MWD$patientid


# 4 select data for CYFRA211=====================================================

MWD_1<-MWD%>%select(PROC3,PROC6,REC1M,C3M,C6M,PROC28,XFIB,PROC4,CYFRA211)

MWD_1_C<-MWD_1%>%drop_na(CYFRA211)
MWD_1_I<-MWD_1%>%filter(is.na(CYFRA211))

# = 5 detail of training and testing ===========================================

set.seed(2000)
ind <- sample(2, nrow(MWD_1_C), replace = TRUE, prob = c(0.7, 0.3))
CYFRA211TrainA_C <-MWD_1_C[ind==1,]
CYFRA211TestA_CC <-MWD_1_C[ind==2,]
CYFRA211TestA_CI <-MWD_1_I


# = 6 Data preprocessing using recipes =========================================

set.seed(2019)
rec_obj <- recipe(CYFRA211 ~ ., data =CYFRA211TrainA_C) %>%
  step_impute_knn(all_predictors()) %>%
  step_dummy(all_predictors(), -all_numeric()) %>%
  step_BoxCox(all_predictors()) %>%
  step_center(all_predictors())  %>%
  step_scale(all_predictors()) %>%
  step_zv(all_predictors()) %>%
  step_corr(all_predictors(), threshold = .90) %>%
  step_log(all_outcomes()) %>%
  check_missing(all_predictors())

rec_obj

trained_rec <- prep(rec_obj, training = CYFRA211TrainA_C)
train <- bake(trained_rec, new_data = CYFRA211TrainA_C)
rownames(CYFRA211TrainA_C)->rownames(train)
test <- bake(trained_rec, new_data = CYFRA211TestA_CC)
rownames(CYFRA211TestA_CC)->rownames(test)
test2<- bake(trained_rec, new_data = CYFRA211TestA_CI)
rownames(CYFRA211TestA_CI)->rownames(test2)

# = 7 Check the data after preprocessing =======================================

anyNA(subset(train, select=CYFRA211))
qqnorm(train$CYFRA211)

# = 8 Modelling xgboost turnning ===============================================

set.seed(2019)
trControl <- trainControl(method='cv',
                          savePredictions="final",
                          index = createFolds(train$CYFRA211,
                                              k = 10, 
                                              returnTrain = TRUE),
                          allowParallel =TRUE, 
                          verboseIter = TRUE)



grid_tune <- expand.grid(nrounds = seq(100,300, length.out = 100),             # number of trees
                         max_depth = c(2,3,4),                                 # tree complexity max 6 max overfitt
                         eta = c(0.02),                                        # Learning rate increase speed but no converese innacuare results
                         gamma = c(0,5),                                       # pruning --> Should be tuned. i.e inf overfitt 
                         colsample_bytree = c(0.285,0.65),                     # subsample ratio of columns for tree % low overfit 
                         min_child_weight = c(1,2),                            # the larger, the more conservative the model hig regualiser use mean 
                         subsample =c(0.6,0.75))                               # prevent overfitting by sampling X% training - 1 overfit



xgb_tune <- train(x = subset(train, select=-c(CYFRA211)),
                  y = train$CYFRA211,
                  trControl = trControl,
                  tuneGrid = grid_tune,
                  method= "xgbTree",
                  objective="reg:squarederror",
                  verbose = F,
                  verbosity = 0)

xgb_tune

# = 9 Modelling xgboost nround and overfitting   ===============================


dtrain <- xgb.DMatrix(data = sparse.model.matrix(CYFRA211~ ., train), label= train$CYFRA211)
dtest <- xgb.DMatrix(data = sparse.model.matrix(CYFRA211~ ., test), label= test$CYFRA211)

watchlist <- list(train = dtrain, test = dtest)


Opt_list<-list(max_depth =3,
               eta = 0.02,
               gamma = 5,
               colsample_bytree = 0.65,
               subsample = 0.6,
               min_child_weight = 3)

set.seed(2019)
xgbcv <- xgb.cv( params = Opt_list, 
                 data = dtrain, 
                 label = train$CYFRA211,
                 watchlist = watchlist,
                 nrounds = 2000, 
                 nfold = 10, 
                 showsd = F, 
                 stratified = T, 
                 print_every_n = 25, 
                 early_stopping_rounds = 50, 
                 maximize = F)
xgbcv

# =10 Training & test error plot ===============================================

e <- data.frame(xgbcv$evaluation_log)

e1<- e[1:2]
e1<-e1%>%mutate(Data="Train")
names(e1)[2]<-paste("RMSE")
e2<-data.frame(e[1],e[4])
names(e2)[2]<-paste("RMSE")
e2<-e2%>%mutate(Data="Test")
e<-rbind(e1,e2)

ntree <- xgbcv$best_iteration
ntree

F1<-ggplot(e,aes(x =iter, y = RMSE, colour = as.factor(Data))) +
  geom_line(size = 1, alpha=0.2)+ geom_point(size = 0.5, alpha=0.2) + theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ 
  theme(legend.position = "none")+
  labs(y="RMSE ",x="Iter", title="",caption = " ")+
  geom_hline(yintercept= 0.6661667, linetype="dashed", color = "darkred", size=0.5)+
  geom_hline(yintercept= 0.6131968, linetype="dashed", color = "royalblue", size=0.5)+
  geom_vline(xintercept= ntree, linetype="dashed", color = "red", size=0.5)+
  theme(plot.title = element_text(size=20, hjust = 0.5, face="bold"),
        axis.title.x = element_text(face="bold", size=15, hjust = 0.5),
        axis.title.y = element_text(face="bold", size=15, hjust = 0.5)) +
  theme(axis.text.x = element_text(face="bold", color="#993333",size=12, angle=0),
        axis.text.y = element_text(face="bold", color="#993333",size=12, angle=0))+
  scale_color_manual(name= " " ,values=c("darkred","royalblue"))
F1

# = 11 most important predictors ===============================================

xgb_mod <- xgb.train(data = dtrain, params=Opt_list, nrounds = ntree)
importance <- xgb.importance(feature_names = xgb_mod$feature_names, model = xgb_mod)

IMPPredictors<-ggplot(importance[1:40,])+
  geom_bar(aes(x=reorder(Feature, Gain), y=Gain), stat='identity', fill='blue')+
  xlab(label = "Top Predictors")+ coord_flip()

IMPPredictors

# = 12  ensemble model training by caret =======================================

set.seed(2019)
trControl <- trainControl(method='cv',
                          savePredictions="final",
                          index = createFolds(train$CYFRA211,
                                              k = 10, 
                                              returnTrain = TRUE),
                          allowParallel =TRUE, 
                          verboseIter = TRUE)



xgbTreeGrid <- expand.grid(nrounds = c(ntree),
                           max_depth = c(2,3,4), 
                           eta =  c(0.02), 
                           gamma = c(5),
                           colsample_bytree =c(0.65),  
                           subsample = c(0.6), 
                           min_child_weight = c(2,3,4))


modelList <<- caretList(x = subset(train, select=-c(CYFRA211)),
                        y = train$CYFRA211,
                        trControl=trControl,
                        metric="RMSE",
                        tuneList=list(xgbTree = caretModelSpec(method="xgbTree",tuneGrid = xgbTreeGrid)))
                                     

modelList$xgbTree

# = 13 visualize the tuning result ==============================================

dev.off()

plot(modelList$xgbTree)
plot(modelList$glmnet)
plot(modelList$svm)


dev.off()
dpi=100
tiff("PredictionCYFRA211pt310222.tiff",  res=dpi, height=8*dpi, width=15*dpi)
bwplot(resamples(modelList),metric="RMSE")
dev.off()

modelCor(resamples(modelList))


greedyEnsemble <- caretEnsemble(   ### equals caretStack (method='glm')
  modelList, 
  metric="RMSE",
  trControl=trainControl(number=10, method = "repeatedcv", repeats=3))
summary(greedyEnsemble)


# = 14 Visualize the reiduals  =================================================

pred.train <- predict(greedyEnsemble, newdata=subset(train, select=-c(CYFRA211)))
dt.plot <- CYFRA211TrainA_C
dt.plot$pred <- pred.train
dt.plot<-dt.plot%>%mutate(residual= log(CYFRA211)-pred)%>%mutate(patientid=rownames(dt.plot))


outlier.res <-dt.plot%>%filter((residual>0.4))

Residual<-ggplot(dt.plot, aes(x = pred, y = residual)) + 
  geom_pointrange(aes(ymin = 0, ymax = residual)) + 
  geom_hline(yintercept = 0, linetype = 3)+
  ggtitle("Residuals vs. model prediction") +
  xlab("prediction") +
  ylab("residual") +
  theme(text = element_text(size=9)) 

Residual


# = 15 Make final prediction  ==================================================

pred1 <- predict(greedyEnsemble, newdata=subset(test, select=-c(CYFRA211)))
test$pred<-pred1
test<-test%>%mutate(SourceCYFRA211="Test_Gboox")
rownames(CYFRA211TestA_CC)->rownames(test)

pred2 <- predict(greedyEnsemble, newdata=subset(test2, select=-c(CYFRA211)))
test2$pred<-pred2
test2<-test2%>%mutate(SourceCYFRA211="Pred_Gboox")
rownames(CYFRA211TestA_CI)->rownames(test2)

predA <- predict(greedyEnsemble, newdata=subset(train, select=-c(CYFRA211)))
train$pred<-predA
train<-train%>%mutate(SourceCYFRA211="Train_Gboox")
rownames(CYFRA211TrainA_C)->rownames(train)

MDataF1<-rbind(train,test,test2)
MDataF1$patientid<-rownames(MDataF1)
MDataFT<-MDataF1%>%select(patientid,pred,SourceCYFRA211)
MDataFT<-MDataFT%>%mutate(pred=(exp(1)^pred))


names(MDataFT)[2]<-paste("CYFRA211pred")


MData<-join(MWD,MDataFT, by="patientid")

trainCYFRA211<-MData%>%filter(SourceCYFRA211=="Train_Gboox")
testCYFRA211<-MData%>%filter(SourceCYFRA211=="Test_Gboox")


# = 15 error calculation  ======================================================

TrainCYFRA211<-trainCYFRA211 %>% dplyr:: summarise(n = n(),max=max(CYFRA211),min=min(CYFRA211),sd= sd(CYFRA211))
TestCYFRA211<-testCYFRA211 %>% dplyr:: summarise(n = n(),max=max(CYFRA211, na.rm=TRUE),min=min(CYFRA211,na.rm=TRUE),sd= sd(CYFRA211, na.rm=TRUE))


FVCTrain_RMSECYFRA211<-(round(rmse(trainCYFRA211$CYFRA211pred,trainCYFRA211$CYFRA211)/(TrainCYFRA211$sd),2))
FVCTest_RMSECYFRA211<-(round(rmse(testCYFRA211$CYFRA211pred,testCYFRA211$CYFRA211)/(TestCYFRA211$sd),2))




CYFRA211g1<-ggscatter(MData, x = "CYFRA211", y = "CYFRA211pred", 
                     add = "reg.line", conf.int = TRUE, 
                     cor.coef = TRUE, cor.method = "pearson")



CYFRA211g1

MData<-MData%>%mutate(CYFRA2111= ifelse(is.na(CYFRA211), CYFRA211pred, CYFRA211)) 
MData$CYFRA211=NULL
dev.off()
dpi=100
tiff("EndoPredGbooxCYFRA211_310522.tiff",  res=dpi, height=10*dpi, width=15*dpi)
grid.arrange(IMPPredictors,Residual,CYFRA211g1,
             F1,ncol=2)
dev.off()

write.csv(MData, file="MWDataBIO_Imp1_463x51_070622.csv", na="")
