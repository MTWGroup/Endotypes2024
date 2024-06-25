# =    1 working space =========================================================
setwd("G:/")
#if (! dir.exists("Project_IPF")) dir.create("Project_IPF")
list.files("G:/")
setwd("G:/NHRI 11 validation predition 160223")
list.files("G:/NHRI 11 validation predition 160223")

# = 2 packages need ============================================================

#analytical pack
library(dplyr)
library(tidyr)
library(plyr)
library(tidymodels)

library(recipes)
library(caret)
require(caretEnsemble)
require(xgboost)
require(kernlab)
require(Matrix)
library(SHAPforxgboost)

library(dendextend)
library(flextable)
webshot::install_phantomjs()
library(grid)
library(gridExtra)
library(lattice)
library(ggpubr)
library(pROC)

# step one only decline and ppFVC00D ===========================================
# 3 open dataset ===============================================================

MWD<- read.csv("MWDataBIO_ImpFinal5_462x59_070622.csv",header=TRUE, sep=",", stringsAsFactors=FALSE, na.strings=c("","NA"))
str(MWD)
dim(MWD)
rownames(MWD)<-MWD$patientid
MWD_B<-MWD%>%select(patientid,MMP7,SPD,CA125,PROC3,PROC4,PROC6,REC1M,XFIB)

MWD_C<-read.csv("PROFILEEndoClustOutPut210223_457x24.csv",header=TRUE, sep=",", stringsAsFactors=FALSE, na.strings=c("","NA"))
str(MWD_C)
dim(MWD)
rownames(MWD_C)<-MWD_C$patientid
MWD_C<-MWD_C%>%select(patientid,ClusterC2)

MWD_B<-join(MWD_B,MWD_C, by="patientid")

MWD_A<-read.csv("MWD_validation_all_bioYuben_260123.csv",header=TRUE, sep=",", stringsAsFactors=FALSE, na.strings=c("","NA"))
str(MWD_A)
dim(MWD_A)
rownames(MWD_A)<-MWD_A$patientid
MWD_A<-MWD_A%>%select(MMP7,SPD,CA125,PROC3,PROC4,PROC6,REC1M,XFIB)
MWD_A<-MWD_A%>%mutate(Class=" ")

summary(MWD_B)
summary(MWD_A)

MWD_B<-MWD_B%>%mutate(Class=case_when(ClusterC2==1 ~ 1, ClusterC2==2 ~ 0, ClusterC2==3 ~ 0))
MWD_B$ClusterC2=NULL
MWD_B<-na.omit(MWD_B)

# = 4 split dataset as a trainning and testing =================================

set.seed(2000)
ind <- sample(2, nrow(MWD_B), replace = TRUE, prob = c(0.7, 0.3))
Train <-MWD_B[ind==1,]
Test <-MWD_B[ind==2,]

rownames(Train)<-Train$patientid
Train$patientid=NULL
rownames(Test)<-Test$patientid
Test$patientid=NULL

# = 5 Data preprocessing using recipes =========================================

set.seed(2019)
rec_obj <- recipe(Class ~ ., data =Train) %>%
  step_impute_knn(all_predictors()) %>%
#  step_center(all_predictors()) %>%
#  step_scale(all_predictors()) %>%
  step_dummy(all_predictors(), -all_numeric()) %>%
  check_missing(all_predictors())

rec_obj

trained_rec <- prep(rec_obj, training = Train)
train <- bake(trained_rec, new_data = Train)
test <- bake(trained_rec, new_data = Test)
test_A <- bake(trained_rec, new_data = MWD_A)

train$Class <-  as.factor(train$Class)
rownames(Train)->rownames(train)
test$Class <- as.factor(test$Class)
rownames(Test)->rownames(test)
test_A$Class<-as.factor(test_A$Class)
rownames(MWD_A)->rownames(test_A)

str(train)
str(test)
str(test_A)

table(train$Class)
table(test$Class)


barplot(prop.table(table(train$Class)),
        col = rainbow(2),
        ylim = c(0, 0.7),
        main = "Class Distribution")

barplot(prop.table(table(test$Class)),
        col = rainbow(2),
        ylim = c(0, 0.7),
        main = "Class Distribution")

set.seed(2022)
train<- upSample(x = train[, -ncol(train)],y = train$Class)

set.seed(2022)
test<- upSample(x = test[, -ncol(train)],y = test$Class)

barplot(prop.table(table(train$Class)),
        col = rainbow(2),
        ylim = c(0, 0.7),
        main = "Class Distribution")


str(train)

train$Class <-  (as.integer(train$Class)-1)


# = 8 Modelling 1)xgboost ======================================================

final_param<-list( max_depth = 4,
                   eta = 0.02,
                   gamma = 0,
                   colsample_bytree = 0.65,
                   subsample = 1,
                   min_child_weight = 3,
                   eval_metric = list("rmse","auc"),
                   objective = "binary:logistic")

dtrain <- xgb.DMatrix(data = sparse.model.matrix(Class~ ., train), label= train$Class)

set.seed(2019)
xgbcv <- xgb.cv( params = final_param, 
                 data = dtrain, 
                 label = train$Class, 
                 nrounds = 2000, 
                 nfold = 10,
                 showsd = F, 
                 stratified = T, 
                 print_every_n = 25, 
                 early_stopping_rounds = 50, 
                 maximize = F)

ntree <- xgbcv$best_iteration

e<- data.frame(xgbcv$evaluation_log)
e1<- e[1:2]
e1<-e1%>%mutate(Data="Train")
names(e1)[2]<-paste("RMSE")
e2<-data.frame(e[1],e[4])
names(e2)[2]<-paste("RMSE")
e2<-e2%>%mutate(Data="Test")

e<-rbind(e1,e2)
MT<-e%>%group_by(Data)%>%summarise(min=min(RMSE))

F1<-ggplot(e,aes(x =iter, y = RMSE, colour = as.factor(Data))) +
  geom_line(size = 1, alpha=0.2)+ geom_point(size = 0.5, alpha=0.2) + theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ 
  theme(legend.position = "none")+
  labs(y="RMSE",x="Iter", title="",caption = " ")+
  geom_hline(yintercept= MT$min[2], linetype="dashed", color = "royalblue", size=0.5)+
  geom_hline(yintercept= MT$min[1], linetype="dashed", color = "red", size=0.5)+
  geom_vline(xintercept= ntree, linetype="dashed", color = "red", size=0.5)+
  theme(plot.title = element_text(size=20, hjust = 0.5, face="bold"),
        axis.title.x = element_text(face="bold", size=15, hjust = 0.5),
        axis.title.y = element_text(face="bold", size=15, hjust = 0.5)) +
  theme(axis.text.x = element_text(face="bold", color="#993333",size=12, angle=0),
        axis.text.y = element_text(face="bold", color="#993333",size=12, angle=0))+
  scale_color_manual(name= " " ,values=c("darkred","royalblue"))
F1

# = 9 most important predictors ================================================

xgb_mod <- xgb.train(data = dtrain, params=final_param, nrounds = ntree)
importance <- xgb.importance(feature_names = xgb_mod$feature_names, model = xgb_mod)


IMPPredictors<-ggplot(importance[1:8,])+
  geom_bar(aes(x=reorder(Feature, Gain), y=Gain), stat='identity', fill='royalblue')+
  xlab(label = "Top Predictors")+
  theme_bw() +
  labs(x="Top Biomarkers",y="Gini coefficient", title="Distribution of Biomarkers Cluster 1 ",caption = " ")+
  theme(plot.title = element_text(size=15, hjust = 0.5, face="bold"),
        axis.title.x = element_text(face="bold", size=17, hjust = 0.5),
        axis.title.y = element_text(face="bold", size=17, hjust = 0.5)) +
  theme(axis.text.x = element_text(face="bold", color="#993333",size=15, angle=0),
        axis.text.y = element_text(face="bold", color="#993333",size=15, angle=0))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ 
  coord_flip()

IMPPredictors

str(train)
str(test)
str(test_A)

# = 10  ensemble model training by caret =======================================

train_2<-train
train_2$Class <- as.factor(train_2$Class)
levels(train_2$Class) <- c("X0", "X1")
str(train_2)

test_2<-test
test_2$Class <- as.factor(test_2$Class)
levels(test_2$Class) <- c("X0", "X1")
str(test_2)

test_A2<-test_A
test_A2$Class <- as.factor(test_A2$Class)
levels(test_A2$Class) <- c("X0", "X1")
str(test_A2)

set.seed(2019)
trControl <- trainControl(method='cv',
                          savePredictions="final",
                          index = createFolds(train_2$Class,
                                              k = 3, 
                                              returnTrain = TRUE),
                          allowParallel =TRUE, 
                          verboseIter = TRUE,
                          classProbs=TRUE)


xgbTreeGrid <- expand.grid(nrounds = ntree, 
                           max_depth = c(3,4,5), 
                           eta =c(0, 0.05, 0.1), 
                           gamma = c(0,0.05,1), 
                           colsample_bytree = c(0.4,0.6,1),  
                           subsample = c(0.5,0.6,1), 
                           min_child_weight = c(3,4,5))


glmnetGrid <- expand.grid(alpha = 1, lambda = seq(0.00001,0.01,by = 0.0001))

svmGrid <- expand.grid(sigma= 2^seq(-11, -16, -0.5), C= 2^seq(4,9,1))


modelList <<- caretList(x = subset(train_2, select=-c(Class)),
                        y = train_2$Class,
                        trControl=trControl,
                        metric="accuracy",
                        tuneList=list(xgbTree = caretModelSpec(method="xgbTree",tuneGrid = xgbTreeGrid),
                                      glmnet=caretModelSpec(method="glmnet", tuneGrid = glmnetGrid),
                                      svmRadial = caretModelSpec(method="svmRadial", tuneGrid = svmGrid, preProcess=c("nzv", "pca"))))

modelCor(resamples(modelList))

greedyEnsemble <- caretEnsemble(   ### equals caretStack (method='glm')
  modelList, 
  metric  = "Accuracy",
  trControl=trainControl(number=10, method = "repeatedcv", repeats=3))

results<-data.frame(summary(greedyEnsemble))


param_list <- list(max_depth = 5, 
                   eta = 0.05, 
                   gamma = 0.05, 
                   colsample_bytree= 0.6, 
                   min_child_weight = 3, 
                   subsample = 1,
                   eval_metric = list("auc"))


set.seed(2019)
mod <- xgboost::xgboost(data = as.matrix(subset(train, select=-c(Class))),
                        label = as.matrix(train$Class),
                        params = param_list, 
                        verbose = T, 
                        nrounds = ntree,
                        nthread = parallel::detectCores() - 2,
                        early_stopping_rounds = 50,
                        objective = "binary:logistic")


shap_values <- shap.values(xgb_model = mod, X_train = as.matrix(subset(train, select=-c(Class))))
shap_values$mean_shap_score

shap_long <- shap.prep(xgb_model = mod, X_train =as.matrix(subset(train, select=-c(Class))))

Shap1<-shap.plot.summary(shap_long)
Shap1<-Shap1 +  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

Shap1

# = 12 Clustrer Prediction  ====================================================

Prob<-predict(greedyEnsemble, newdata=subset(test_2, select=-c(Class)),type="prob")
test_2$Prob<-Prob

Pred<-predict(greedyEnsemble, newdata=subset(test, select=-c(Class)))
test_2$Pred<-Pred
test_2<-test_2%>%mutate(ProbCl1=(1-Prob))

test_2<- test_2%>% mutate(TrendRF1 = case_when((Class =="X1" &  Pred=="X1") ~ "True_Poss",
                                           (Class =="X0" &  Pred=="X0") ~ "True_Neg",
                                           (Class =="X1" &  Pred=="X0") ~ "False_Poss",
                                           (Class =="X0" &  Pred=="X1") ~ "False_Neg"))

ConfM<-data.frame(test_2%>% group_by(TrendRF1) %>% dplyr:: summarise(n = n()))
ConfM1<-data.frame(c(" "," "))
ConfM2<-data.frame(c(" "," "))
ConfM3<-cbind(ConfM1,ConfM2)
names(ConfM3)[1:2]<-paste(c("TrendRF1","n"))
ConfM<-rbind(ConfM,ConfM3)
names(ConfM)[1:2]<-paste(c("Confusion","Patients"))

Conf<-confusionMatrix(table(test_2$Pred,test_2$Class))
Conf$overall[1]
Stat2<-data.frame(c(Conf$overall[1]))

Blue_roc <- roc(response = test_2$Class,
               predictor = test_2$Prob,
               levels = c('X0', 'X1'))

ggroc(Blue_roc,legacy.axes = TRUE)
auc <- round(auc(test_2$Class, test_2$Prob),4)

F1<-ggroc(Blue_roc, colour = 'steelblue', size = 2 ,legacy.axes = TRUE) +
  labs(x = 'False-positive rate', y = 'True-positive rate', title = ' ROC curve Cluster 1 - Blue')+
  ggtitle(paste0('ROC Curve Cluster 1 (Blue)', ' (AUC = ', auc, ')'))+
  geom_abline(intercept = 0, slope = 1, colour="grey", size = 1.5)+
  theme_bw()+ theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(plot.title = element_text(size=10, hjust = 0, face="bold"),
        axis.title.x = element_text(face="bold", size=15, hjust = 0.5),
        axis.title.y = element_text(face="bold", size=15, hjust = 0.5)) +
  theme(axis.text.x = element_text(face="bold", color="#993333",size=12, angle=0),
        axis.text.y = element_text(face="bold", color="#993333",size=12, angle=0))
F1

ProbA<-predict(greedyEnsemble, newdata=subset(test_A2, select=-c(Class)),type="prob")
test_A$Prob<-ProbA

Pred<-predict(greedyEnsemble, newdata=subset(test_A2, select=-c(Class)))
test_A$Pred<-Pred
test_A<-test_A%>%mutate(ProbCl1=(1-Prob))

Conf<-data.frame(as.matrix(Conf, what = "classes"))
Conf_RF1<-data.frame(t(Conf))
Conf_RF1<-Conf_RF1%>%mutate(Balanced.Accuracy=auc)
Conf_RF1<-Conf_RF1%>%select(Balanced.Accuracy,Sensitivity,Specificity,Pos.Pred.Value,Neg.Pred.Value)
names(Conf_RF1)[1:5]<-paste (c("AUC","Sensitivity","Specificity","NPV","PPV"))
names(Stat2)[1]<-paste("Accuracy")
Conf_RF1<-cbind(Stat2,Conf_RF1)
Conf_RF1<-data.frame(t(Conf_RF1))
names(Conf_RF1)[1]<-paste("Results")
Conf_RF1<-cbind(Conf_RF1,ConfM)

Conf_RF1<-Conf_RF1%>%mutate(Parameters=rownames(Conf_RF1))
Conf_RF1<-Conf_RF1%>%select(Parameters,Results,Confusion,Patients)%>%mutate(Results=as.numeric(Results))
Conf_RF1$Results<-round(Conf_RF1$Results,3)
Conf_RF1$ConfusionR=NULL

ft <- flextable((Conf_RF1))
ft <- color(ft, part = "footer", color = "#666666")
ft<-theme_vanilla(ft)
#ft <- set_caption(ft, "Table 1: Comparison of baseline clinical characteristics between PROFILE centres",  style = "Table Caption")
ft <- align(ft, align = "center")
ft

ftR<-as_raster(ft)
gflRF <- ggplot() + 
  theme_void() + 
  annotation_custom(rasterGrob(ftR), xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)
gflRF

dev.off()
dpi=300
tiff("Prediction Cluster1.tiff", res=dpi, height=15*dpi, width=18*dpi)
grid.newpage()
# Create layout : nrow = 2, ncol = 2
pushViewport(viewport(layout = grid.layout(6,8)))

define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 
# Arrange the plots

print(IMPPredictors, vp=define_region(1:3,1:4))
print(gflRF, vp=define_region(1:3,5:8))
print(F1, vp=define_region(4:6,1:4))
print(Shap1, vp=define_region(4:6,5:8))
dev.off()

rownames(MWD_A)->rownames(test_A)
rownames(test_A)->test_A$patientid
test_A<-test_A%>%mutate(Pred=case_when(Pred=="X1"~"Cl_1", Pred=="X0"~"Other_Cl"))

write.table(test_A, "Australia Cluster 1 prediction.csv", row.names = F, sep = ",")
