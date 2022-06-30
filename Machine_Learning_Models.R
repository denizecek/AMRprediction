
#Machine Learning

library(randomForest)
library(mlbench)
library(caret)
library(tidyverse)
library(preprocessCore)
library(foreach)
library(reshape2)
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(plotly)

combine_final_data <- "#combine final data"


# get indices for 80% train 10% test 10% validation of the data set

intrain <- createDataPartition(y = combine_final_data$antibiotic, p= 0.80)[[1]]

# seperate test,validation and training sets

training <- combine_final_data[intrain,]

testing_and_validation <- combine_final_prs_bin[-intrain,]

intrain_2 <- createDataPartition(y = testing_and_validation$antibiotic, p= 0.50)[[1]]

testing <- testing_and_validation[intrain_2,]

validation <- testing_and_validation[-intrain_2,]

######################################################################################################################################################

# 1.Random Forest

# 1.1 Random Forest Default 


#Create a train control object with 10-fold CV
trainControl <- trainControl(method="cv", number = 10,
                             savePredictions = "final",
                             search = "random")


#mtry: Number of variables randomly sampled as candidates at each split.
#ntree: Number of trees to grow.


#Using the recommend defaults for each parameter and mtry=floor(sqrt(ncol(x))) 

x <- training["all columns"]

y <- training["antibiotic column"]

# Create model with default parameters

control <- trainControl(method="repeatedcv", number=10, repeats=3,classProbs = TRUE )

mtry <- sqrt(ncol(x))

tunegrid <- expand.grid(.mtry=mtry)

rf_default <- train(antibiotic~., data=training, method="rf", metric=metric, tuneGrid=tunegrid, trControl=control)

#Results

rf_default_result <- rf_default$results


#density plot for CV folds
lapply(rf_default$control$index, function(index) training[index,]) %>%
  bind_rows(.id = "Fold") %>%
  ggplot(aes(antibiotic, col = Fold)) +
  geom_density()+
  ggtitle("Density Plot for CV Folds")

# Importance and top ten feature plots

plot(varImp(rf_default))

plot(varImp(rf_default), top=10)

# training accuracy

class.res.rfdefault_training=predict(rf_default,training)

class.res.rfdefault_training

# test accuracy 

class.res.rfdefault_test=predict(rf_default,testing)

class.res.rfdefault_test

# Confusion Matrix


confusionMatrix((as.factor(testing["antibiotic column"])),(as.factor(class.res.rfdefault_test)))

# F1 Score


library(MLmetrics)

y_pred = class.res.rfdefault_test

y_actual = testing["antibiotic column"]

res = F1_Score(y_pred,y_actual)


# 1.2 Random Search

#Random search is a method in which random combinations of hyperparameters are selected and used to train a model. The best random hyperparameter combinations are used. 

control <- trainControl(method="repeatedcv", number=10, repeats=3, search="random", classProbs = TRUE )

set.seed(seed)

mtry <- sqrt(ncol(x))

rf_random <- train(antibiotic~., data=training, method="rf", metric=metric, tuneLength=15, trControl=control)


plot(rf_random)

#Results

rf_random_result <- rf_random$results


#density plot for CV folds

lapply(rf_random$control$index, function(index) training[index,]) %>%
  bind_rows(.id = "Fold") %>%
  ggplot(aes(antibiotic, col = Fold)) +
  geom_density()+
  ggtitle("Density Plot for CV Folds")

# Importance and top ten feature plots

plot(varImp(rf_random),top=10)

plot(varImp(rf_random))


# training accuracy 

class.res.rfrandom_training=predict(rf_random,training)

class.res.rfrandom_training

# training accuracy 

class.res.rfrandom_validation=predict(rf_random,validation)

class.res.rfrandom_validation


# test accuracy 

class.res.rfrandom_test=predict(rf_random,testing)

class.res.rfrandom_test

cm_random = table(testing["antibiotic column"], class.res.rfrandom_test)


#Confusion Matrix and Statistics

confusionMatrix((as.factor(testing["antibiotic column"])),(as.factor(class.res.rfrandom_test)))


# F1 Score

y_pred = class.res.rfdefault_test

y_actual = testing["antibiotic column"]

res = F1_Score(y_pred,y_actual)


# 1.3 Random Forest Grid Search


#Random search allowed us to narrow down the range for each hyperparameter, now with grid search specify every combination of settings to try.

control <- trainControl(method="repeatedcv", number=10, repeats=3, search="grid",classProbs = TRUE)

set.seed(seed)

tunegrid <- expand.grid(.mtry=c("best parameters"))

rf_gridsearch <- train(antibiotic~., data=training, method="rf", metric=metric, tuneGrid=tunegrid, trControl=control)

#Results

rf_gridsearch_result <- rf_gridsearch$results


#density plot for CV folds

lapply(rf_gridsearch$control$index, function(index) training[index,]) %>%
  bind_rows(.id = "Fold") %>%
  ggplot(aes(antibiotic, col = Fold)) +
  geom_density()+
  ggtitle("Density Plot for CV Folds")

# Importance and top ten feature plots

plot(varImp(rf_gridsearch),top=10)

plot(varImp(rf_gridsearch))


# training accuracy 

class.res.rfgridsearch_training=predict(rf_gridsearch,training)

# validation accuracy 

class.res.rfgridsearch_val=predict(rf_gridsearch,validation)

cm_gridsearch = table(validation["antibiotic column"], class.res.rfgridsearch_val)

#Confusion Matrix and Statistics

confusionMatrix((as.factor(validation["antibiotic column"])),(as.factor(class.res.rfgridsearch_val)))


# test accuracy 

class.res.rfgridsearch_test=predict(rf_gridsearch,testing)

class.res.rfgridsearch_test

cm_gridsearch = table(testing["antibiotic column"], class.res.rfgridsearch_test)


#Confusion Matrix and Statistics

confusionMatrix((as.factor(testing["antibiotic column"])),(as.factor(class.res.rfgridsearch_test)))

# Plots of model

trellis.par.set(caretTheme())

plot(rf_gridsearch)  

trellis.par.set(caretTheme())

plot(rf_gridsearch, metric = "Kappa")

ggplot(rf_gridsearch) 

# F1 Score

y_pred = class.res.rfgridsearch_test

y_actual = testing["antibiotic column"]

res = F1_Score(y_pred,y_actual)

#####################################################################################################################################################

# 2.Support Vector Machine

# 2.1 Support Vector Machines with Linear

# Fit the model 

svm1 <- train(antibiotic ~., data = training, method = "svmLinear", trControl = trainControl)

#Results

svm1_model_result <- svm1$results

# Importance and top ten feature plots

plot(varImp(svm1),top=10)

plot(varImp(svm1))

# training accuracy 
class.res.svm1_training=predict(svm1,training)

class.res.svm1_training

# test accuracy 
class.res.svm1_test=predict(svm1,testing)

class.res.svm1_test

cm_svm1_model = table(testing["antibiotic column"], class.res.svm1_test)

#Confusion Matrix and Statistics

confusionMatrix((as.factor(testing["antibiotic column"])),(as.factor(class.res.svm1_test)))

# F1 Score

y_pred = class.res.svm1_test

y_actual = testing["antibiotic column"]

res = F1_Score(y_pred,y_actual)

# Plots of the model

trellis.par.set(caretTheme())

plot(svm1)  

trellis.par.set(caretTheme())

plot(svm1, metric = "Kappa")

ggplot(svm1) 

# 2.2 Support Vector Machines with Linear w/ choice of cost  

svm2 <- train(antibiotic ~., data = training, method = "svmLinear", 
              trControl = trainControl, 
              tuneGrid = expand.grid(C = seq(0, 2, length = 20)))

# Print the best tuning parameter

svm2$bestTune

#Results

svm2_model_result <- svm2$results

# Importance and top ten feature plots

plot(varImp(svm2),top=10)

plot(varImp(svm2))

# training accuracy 

class.res.svm2_training=predict(svm2,training)

class.res.svm2_training

# test accuracy 
class.res.svm2_test=predict(svm2,testing)

class.res.svm2_test

cm_svm2_model = table(testing["antibiotic column"], class.res.svm2_test)

#Confusion Matrix and Statistics

confusionMatrix((as.factor(testing["antibiotic column"])),(as.factor(class.res.svm2_test)))

# F1 Score

y_pred = class.res.svm2_test

y_actual = testing["antibiotic column"]

res = F1_Score(y_pred,y_actual)

# Plots of the model

trellis.par.set(caretTheme())

plot(svm2)  

trellis.par.set(caretTheme())

plot(svm2, metric = "Kappa")

ggplot(svm2)

# 2.3 Support Vector Machines with Radial Basis Function Kernel 


# Fit the model 
svm3 <- train(antibiotic ~., data = training, method = "svmRadial", trControl = trainControl, tuneLength = 10)

# Print the best tuning parameter

svm3$bestTune

#Results

svm3_model_result <- svm3$results

# Importance and top ten feature plots

plot(varImp(svm3),top=10)

plot(varImp(svm3))

# training accuracy 

class.res.svm3_training=predict(svm3,training)

class.res.svm3_training

# test accuracy 

class.res.svm3_test=predict(svm3,testing)

class.res.svm3_test

cm_svm3_model = table(testing["antibiotic column"], class.res.svm3_test)

#Confusion Matrix and Statistics

confusionMatrix((as.factor(testing["antibiotic column"])),(as.factor(class.res.svm3_test)))

# F1 Score

y_pred = class.res.svm3_test

y_actual = testing["antibiotic column"]

res = F1_Score(y_pred,y_actual)

# Plots of the model

trellis.par.set(caretTheme())

plot(svm3)  

trellis.par.set(caretTheme())

plot(svm3, metric = "Kappa")

ggplot(svm3)

# 2.4 Support Vector Machines with Polynomial Kernel

# Fit the model 

svm4 <- train(antibiotic~., data = training, method = "svmPoly", trControl = trainControl, tuneLength = 4)

# Print the best tuning parameter sigma and C that maximizes model accuracy

svm4$bestTune

#Results

svm4_model_result <- svm4$results

# Importance and top ten feature plots

plot(varImp(svm4),top=10)

plot(varImp(svm4))

# training accuracy 

class.res.svm4_training=predict(svm4,training)


# test accuracy 

class.res.svm4_test=predict(svm4,testing)

cm_svm4_model = table(testing["antibiotic column"], class.res.svm4_test)

#Confusion Matrix and Statistics

confusionMatrix((as.factor(testing["antibiotic column"])),(as.factor(class.res.svm4_test)))

# F1 Score

y_pred = class.res.svm4_test

y_actual = testing["antibiotic column"]

res = F1_Score(y_pred,y_actual)

# Plots of the model

trellis.par.set(caretTheme())

plot(svm4)  

trellis.par.set(caretTheme())

plot(svm4, metric = "Kappa")

ggplot(svm4)

#####################################################################################################################################################

# 3. Stochastic Gradient Boosting 

# 3.1 Gradient Boosting with Default Parameters 

library(gbm)

gbm_model <- train(antibiotic~., data = training, 
                   method = "gbm",
                   trControl= trainControl)

# Print the best tuning parameter

gbm_model$bestTune

# Importance and top ten feature plots

plot(varImp(gbm_model),top=10)

plot(varImp(gbm_model))

#Results

gbm_model_result <- gbm_model$results

#density plot for CV folds

lapply(gbm_model$control$index, function(index) training[index,]) %>%
  bind_rows(.id = "Fold") %>%
  ggplot(aes(antibiotic, col = Fold)) +
  geom_density()+
  ggtitle("Density Plot for CV Folds")

# training accuracy 

class.res.gbmmodel_training=predict(gbm_model,training)


# test accuracy 

class.res.gbmmodel_test=predict(gbm_model,testing)

cm_gbm_model = table(testing["antibiotic column"], class.res.gbmmodel_test)

#Confusion Matrix and Statistics

confusionMatrix((as.factor(testing["antibiotic column"])),(as.factor(class.res.gbmmodel_test)))

# Plots of the model

trellis.par.set(caretTheme())

plot(gbm_model)  

trellis.par.set(caretTheme())

plot(gbm_model, metric = "Kappa")

trellis.par.set(caretTheme())

plot(gbm_model, metric = "Kappa", plotType = "level",
     scales = list(x = list(rot = 90)))

ggplot(gbm_model) 

# F1 Score

y_pred = class.res.gbmmodel_test

y_actual = testing["antibiotic column"]

res = F1_Score(y_pred,y_actual)

# 3.2 Gradient Boosting with 10-fold CV

#gridversion

fitControl <- trainControl(## 10-fold CV
  method = "repeatedcv",
  number = 10,
  repeats = 2)

gbmFit1 <- train(antibiotic ~ ., data = training, 
                 method = "gbm", 
                 trControl = fitControl,
                 verbose = FALSE)

# Print the best tuning parameter

gbmFit1$bestTune

# Importance and top ten feature plots

plot(varImp(gbmFit1),top=10)

plot(varImp(gbmFit1))

#Results

gbmFit1_result <- gbmFit1$results

# training accuracy

class.res.gbmFit1_training=predict(gbmFit1,training)

# test accuracy 

class.res.gbmFit1_test=predict(gbmFit1,testing)

cm_gbmFit1 = table(testing["antibiotic column"], class.res.gbmFit1_test)

#Confusion Matrix and Statistics

confusionMatrix((as.factor(testing["antibiotic column"])),(as.factor(class.res.gbmFit1_test)))

# F1 Score

y_pred = class.res.gbmFit1_test

y_actual = testing["antibiotic column"]

res = F1_Score(y_pred,y_actual)

# Plots of the model

trellis.par.set(caretTheme())

plot(gbmFit1)  

trellis.par.set(caretTheme())

plot(gbmFit1, metric = "Kappa")

trellis.par.set(caretTheme())

plot(gbmFit1, plotType = "level",
     scales = list(x = list(rot = 90)))

ggplot(gbmFit1) 


## 3.3 Gradient Boosting with Grid Version

gbmGrid <-  expand.grid(interaction.depth = c("new parameters"), 
                        n.trees = ("new parameters"), 
                        shrinkage = "new parameters",
                        n.minobsinnode = "new parameters")


gbmFit2 <- train(antibiotic ~ ., data = training, 
                 method = "gbm", 
                 trControl = fitControl, 
                 verbose = FALSE, 
                 tuneGrid = gbmGrid)


# Importance and top ten feature plots

plot(varImp(gbmFit2),top=10)

plot(varImp(gbmFit2))

# training accuracy 

class.res.gbmfit2_training=predict(gbmFit2,training)

# validation accuracy 

class.res.gbmfit2_val=predict(gbmFit2,validation)

cm_gridsearch = table(validation["antibiotic column"], class.res.gbmFit2_val)

# test accuracy 

class.res.gbmfit2_test=predict(gbmFit2,testing)

cm_gbmfit2_model = table(testing["antibiotic column"], class.res.gbmfit2_test)

#Confusion Matrix and Statistics

confusionMatrix((as.factor(testing["antibiotic column"])),(as.factor(class.res.gbmfit2_test)))

# F1 Score

y_pred = class.res.gbmfit2_test

y_actual = testing["antibiotic column"]

res = F1_Score(y_pred,y_actual)

# Plots of the model

trellis.par.set(caretTheme())

plot(gbmFit2)  

trellis.par.set(caretTheme())

plot(gbmFit2, metric = "Kappa")

trellis.par.set(caretTheme())

plot(gbmFit2, plotType = "level",
     scales = list(x = list(rot = 90)))

ggplot(gbmFit2)


####################################################################################################################################################


# 4. eXtreme Gradient Boosting 

# 4.1 eXtreme Gradient Boosting w Default Parameters

library(xgboost)

xgb_model <- train(antibiotic~., data = training, 
                   trControl= trainControl,
                   method = "xgbTree")

# Print the best tuning parameter

xgb_model$bestTune

# Importance and top ten feature plots

plot(varImp(xgb_model),top=10)

plot(varImp(xgb_model))

#Results

xgb_model_result <- xgb_model$results

#density plot for CV folds

lapply(xgb_model$control$index, function(index) training[index,]) %>%
  bind_rows(.id = "Fold") %>%
  ggplot(aes(antibiotic, col = Fold)) +
  geom_density()+
  ggtitle("Density Plot for CV Folds")

# training accuracy 

class.res.xgb_training=predict(xgb_model,training)

# test accuracy 

class.res.xgb_test=predict(xgb_model,testing)

cm_xgb_model = table(testing["antibiotic column"], class.res.xgb_test)

#Confusion Matrix and Statistics

confusionMatrix((as.factor(testing["antibiotic column"])),(as.factor(class.res.xgb_test)))

# Plots of the model

trellis.par.set(caretTheme())

plot(xgb_model)  

trellis.par.set(caretTheme())

plot(xgb_model, metric = "Kappa")

trellis.par.set(caretTheme())

plot(xgb_model, plotType = "level",
     scales = list(x = list(rot = 90)))

ggplot(xgb_model) 

# F1 Score

y_pred = class.res.xgb_test

y_actual = testing["antibiotic column"]

res = F1_Score(y_pred,y_actual)

# 4.2 eXtreme Gradient Boosting w 10-fold CV

#gridversion
fitControl <- trainControl(## 10-fold CV
  method = "repeatedcv",
  number = 10,
  ## repeated ten times
  repeats = 2)

xgbFit1 <- train(antibiotic ~ ., data = training, 
                 method = "xgbTree", 
                 trControl = fitControl,
                 ## This last option is actually one
                 ## for gbm() that passes through
                 verbose = FALSE)

# Print the best tuning parameter

xgbFit1$bestTune

# Importance and top ten feature plots

plot(varImp(xgbFit1),top=10)

plot(varImp(xgbFit1))

#Results

xgbFit1_result <- xgbFit1$results

#density plot for CV folds

lapply(xgbFit1$control$index, function(index) training[index,]) %>%
  bind_rows(.id = "Fold") %>%
  ggplot(aes(antibiotic, col = Fold)) +
  geom_density()+
  ggtitle("Density Plot for CV Folds")

# training accuracy 

class.res.xgbFit1_training=predict(xgbFit1,training)

class.res.xgbFit1_training

# test accuracy 

class.res.xgbFit1_test=predict(xgbFit1,testing)

class.res.xgbFit1_test

cm_xgbFit1 = table(testing["antibiotic column"], class.res.xgbFit1_test)

#Confusion Matrix and Statistics

confusionMatrix((as.factor(testing["antibiotic column"])),(as.factor(class.res.xgbFit1_test)))
 
# Plots of the model

trellis.par.set(caretTheme())

plot(xgbFit1)  

trellis.par.set(caretTheme())

plot(xgbFit1, metric = "Kappa")

trellis.par.set(caretTheme())

plot(xgbFit1, plotType = "level",
     scales = list(x = list(rot = 90)))

# F1 Score

y_pred = class.res.xgbFit1_test

y_actual = testing["antibiotic column"]

res = F1_Score(y_pred,y_actual)

## 4.3 eXtreme Gradient Boosting w Grid version


xgbGrid <- expand.grid(nrounds = c("new parameters"),  
                       max_depth = c("new parameters"),
                       colsample_bytree = seq("new parameters", length.out = "new parameters"),
                       eta = "new parameters",
                       gamma="new parameters",
                       min_child_weight =c("new parameters") ,
                       subsample =c("new parameters"))


xgbFit2 <- train(antibiotic ~ ., data = training, 
                 method = "xgbTree", 
                 trControl = fitControl, 
                 verbose = FALSE, 
                 ## Now specify the exact models 
                 ## to evaluate:
                 tuneGrid = xgbGrid)

# Print the best tuning parameter

xgbFit2$bestTune

# Importance and top ten feature plots

plot(varImp(xgbFit2),top=10)

plot(varImp(xgbFit2))

#Results

xgbFit2_result <- xgbFit2$results

#density plot for CV folds

lapply(xgbFit2$control$index, function(index) training[index,]) %>%
  bind_rows(.id = "Fold") %>%
  ggplot(aes(antibiotic, col = Fold)) +
  geom_density()+
  ggtitle("Density Plot for CV Folds")

# training accuracy 

class.res.xgbFit2_training=predict(xgbFit2,training)


# validation accuracy 

class.res.xgbfit2_val=predict(xgbFit2,validation)

cm_gridsearch = table(validation["antibiotic column"], class.res.xgbFit2_val)

# test accuracy 

class.res.xgbFit2_test=predict(xgbFit2,testing)

cm_xgbFit2 = table(testing["antibiotic column"], class.res.xgbFit2_test)

#Confusion Matrix and Statistics

confusionMatrix((as.factor(testing["antibiotic column"])),(as.factor(class.res.xgbFit2_test)))

# F1 Score

y_pred = class.res.xgbFit2_test

y_actual = testing["antibiotic column"]

res = F1_Score(y_pred,y_actual)

# Plots of the model

trellis.par.set(caretTheme())

plot(xgbFit2)  

trellis.par.set(caretTheme())

plot(xgbFit2, metric = "Kappa")

trellis.par.set(caretTheme())

plot(xgbFit2, plotType = "level",
     scales = list(x = list(rot = 90)))

ggplot(xgbFit2) 

####################################################################################################################################################



