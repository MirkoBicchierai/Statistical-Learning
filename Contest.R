options(java.parameters = "-Xmx16g")
library(caret)
library(rpart)
library(partykit)
library(randomForest)
library(ggplot2)
library(bartMachine)


data <- read.csv("C:\\Users\\mirko\\Desktop\\games.csv")

summary(data)
x <- data[c(3,6:11,27:31,52:56)]

# Select the target variable
y <- data[, 5]  # Use column index to extract y_train as a vector

n <- nrow(x)
set.seed(111)

id.test <- sample(n, floor(n/3))

x_train <- x[-id.test,]
y_train <- y[-id.test]

x_test <- x[id.test,]
y_test <- y[id.test]


#plot histogram basic variable  T1

data_new <- data.frame(data["winner"], data["t1_towerKills"])
means <- aggregate(. ~ winner, data = data_new, FUN = mean)
ggplot(means, aes(x = factor(winner), y = t1_towerKills)) +
  geom_bar(stat = "identity", fill = "steelblue", alpha = 0.85, width = NULL) +  # Remove width
  labs(title = "Tower Kills T1 - Winners",
       x = "Winner",
       y = "Average Tower Kills") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~., ncol = 2) +  # Change nrow to ncol for vertical layout
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(data_new$winner)  # Setting x-axis ticks


data_new <- data.frame(data["winner"], data["t1_inhibitorKills"])
means <- aggregate(. ~ winner, data = data_new, FUN = mean)
ggplot(means, aes(x = factor(winner), y = t1_inhibitorKills)) +
  geom_bar(stat = "identity", fill = "steelblue", alpha = 0.85, width = NULL) +  # Remove width
  labs(title = "t1_inhibitorKills - Winners",
       x = "Winner",
       y = "Average Inhibitor Kills") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~., ncol = 2) +  # Change nrow to ncol for vertical layout
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(data_new$winner)  # Setting x-axis ticks


data_new <- data.frame(data["winner"], data["t1_baronKills"])
means <- aggregate(. ~ winner, data = data_new, FUN = mean)
ggplot(means, aes(x = factor(winner), y = t1_baronKills)) +
  geom_bar(stat = "identity", fill = "steelblue", alpha = 0.85, width = NULL) +  # Remove width
  labs(title = "t1_baronKills - Winners",
       x = "Winner",
       y = "Average Baron Kills") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~., ncol = 2) +  # Change nrow to ncol for vertical layout
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(data_new$winner)  # Setting x-axis ticks


data_new <- data.frame(data["winner"], data["t1_dragonKills"])
means <- aggregate(. ~ winner, data = data_new, FUN = mean)
ggplot(means, aes(x = factor(winner), y = t1_dragonKills)) +
  geom_bar(stat = "identity", fill = "steelblue", alpha = 0.85, width = NULL) +  # Remove width
  labs(title = "t1_dragonKills - Winners",
       x = "Winner",
       y = "Average Dragon Kills") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~., ncol = 2) +  # Change nrow to ncol for vertical layout
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(data_new$winner)  # Setting x-axis ticks


#plot histogram basic variable  T2

data_new <- data.frame(data["winner"], data["t2_towerKills"])
means <- aggregate(. ~ winner, data = data_new, FUN = mean)
ggplot(means, aes(x = factor(winner), y = t2_towerKills)) +
  geom_bar(stat = "identity", fill = "steelblue", alpha = 0.85, width = NULL) +  # Remove width
  labs(title = "Tower Kills T2 - Winners",
       x = "Winner",
       y = "Average Tower Kills") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~., ncol = 2) +  # Change nrow to ncol for vertical layout
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(data_new$winner)  # Setting x-axis ticks


data_new <- data.frame(data["winner"], data["t2_inhibitorKills"])
means <- aggregate(. ~ winner, data = data_new, FUN = mean)
ggplot(means, aes(x = factor(winner), y = t2_inhibitorKills)) +
  geom_bar(stat = "identity", fill = "steelblue", alpha = 0.85, width = NULL) +  # Remove width
  labs(title = "t2_inhibitorKills - Winners",
       x = "Winner",
       y = "Average Inhibitor Kills") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~., ncol = 2) +  # Change nrow to ncol for vertical layout
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(data_new$winner)  # Setting x-axis ticks


data_new <- data.frame(data["winner"], data["t2_baronKills"])
means <- aggregate(. ~ winner, data = data_new, FUN = mean)
ggplot(means, aes(x = factor(winner), y = t2_baronKills)) +
  geom_bar(stat = "identity", fill = "steelblue", alpha = 0.85, width = NULL) +  # Remove width
  labs(title = "t2_baronKills - Winners",
       x = "Winner",
       y = "Average Baron Kills") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~., ncol = 2) +  # Change nrow to ncol for vertical layout
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(data_new$winner)  # Setting x-axis ticks


data_new <- data.frame(data["winner"], data["t2_dragonKills"])
means <- aggregate(. ~ winner, data = data_new, FUN = mean)
ggplot(means, aes(x = factor(winner), y = t2_dragonKills)) +
  geom_bar(stat = "identity", fill = "steelblue", alpha = 0.85, width = NULL) +  # Remove width
  labs(title = "t2_dragonKills - Winners",
       x = "Winner",
       y = "Average Dragon Kills") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~., ncol = 2) +  # Change nrow to ncol for vertical layout
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(data_new$winner)  # Setting x-axis ticks








#Check if the training set is balanced
mean(y_train)
# Combine x_train and y_train into a single data frame for modeling
train_data <- cbind(x_train, y_train)


# Fit a decision tree model using the Gini index as a loss function
fit <- rpart(y_train ~ ., method = "class", data = train_data, control = rpart.control(cp = 0, minsplit=20))
printcp(fit)

plotcp(fit)
minCP <- fit$cptable[which.min(fit$cptable[,"xerror"]),"CP"]
minCP

# Prune the tree
fit.pr <- prune(fit, cp = minCP)

rpart.plot::rpart.plot(fit.pr, main = "Base")



# Fit a decision tree model using the Information gain index as a loss function
fit_entropy <- rpart(y_train ~ ., data = train_data, method = "class", 
                     parms = list(split = "information"),
                     control = rpart.control(cp = 0, minsplit = 20))
printcp(fit_entropy)

plotcp(fit_entropy)
minCP_e <- fit_entropy$cptable[which.min(fit_entropy$cptable[,"xerror"]),"CP"]
minCP_e
# Plot the pruned decision tree
fit_entropy.pr <- prune(fit_entropy, cp = minCP_e)

rpart.plot::rpart.plot(fit_entropy.pr, main = "Entropy")



#Confusion matrix for the predictions on the training set
predetto.g <- predict(fit.pr,  type = "class")
predetto.e <- predict(fit_entropy.pr,  type = "class")

caret::confusionMatrix(table(predicted = predetto.g, actual = y_train), positive = "1")
caret::confusionMatrix(table(predicted = predetto.e, actual = y_train), positive = "1")


#Confusion matrix for the predictions on the testing set
predetto.gtest <- predict(fit.pr,  type = "class", newdata = x_test)
predetto.etest <- predict(fit_entropy.pr,  type = "class", newdata = x_test)

caret::confusionMatrix(table(predicted = predetto.gtest, actual = y_test), positive = "1")
caret::confusionMatrix(table(predicted = predetto.etest, actual = y_test), positive = "1")


# Variable Importance, find best ntree
y_train <- as.factor(y_train)
RF <- randomForest(y_train ~ . , data = x_train, ntree=500)
plot(RF, col="#A20045", main="Random forest")
RF

# Variable Importance, with best ntree
RF <- randomForest(y_train ~ . , data = x_train, ntree=200, mtry=5)
varImpPlot(RF, main="Variable importance", pch = 19, color='blue')



#BART CV
y_train = factor(y_train)  
num_trees_values <- c(50, 100, 150, 200)
cv_results <- bartMachineCV(
  X = x_train,        # Predictors
  y = y_train,      # Response variable (classification example)
  num_tree_cvs = num_trees_values,  # Values of num_trees to tune
  k_folds = 4,
  k_cvs = c(2,3)
)
print(cv_results)

cv_results$num_trees # = 100
cv_results$k # = 2
cv_results$confusion_matrix


#We have to interpret those results
cv_results$cv_stats 



#Re-train the selected best model with ntree = 100 and k = 2, using 8 cores(it can crash tho)
set_bart_machine_num_cores(8)
bestBART <- bartMachine(x_train, y_train,num_trees = 100, k=2, serialize = TRUE)

#Summary of the trained model with confusion matrix on train data
summary(bestBART)



#Prediction on test data
testy <- factor(y_test)
pred_test <- bart_predict_for_test_data(bestBART, x_test, testy)


#Print the confusion matrix on test data and compute test accuracy

pred_test$confusion_matrix
cm_test <- pred_test$confusion_matrix[-3,-3]
test_accuracy <- (cm_test[1,1] + cm_test[2,2])/(cm_test[1,1] + cm_test[1,2] + cm_test[2,1] + cm_test[2,2])
print(paste("Test Accuracy:", round(test_accuracy * 100,2), "%"))


investigate_var_importance(bestBART, plot= TRUE, bottom_margin=8, type="splits")
investigate_var_importance(bestBART, plot= TRUE, bottom_margin=8, type="trees") 




