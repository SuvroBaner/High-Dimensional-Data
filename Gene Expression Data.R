################ Gene expression data ###################

# We'll now examine the Khan data set, which consists of a number of tissue samples 
# corresponding to four distinct types of small round blue cell tumors.

# For each tissue sample, gene expression measurements are available.

# The data set consists of training data, xtrain and ytrain
# and testing data xtest and ytest.


install.packages("ISLR")
install.packages("e1071")
install.packages("LiblineaR") # useful for very large linear problems.

library(ISLR)
library(e1071)
library(LiblineaR)


names(Khan)  

dim(Khan$xtrain)
# 63 tissue samples and 2308 predictors(regarding gene expression measurements)

dim(Khan$xtest)
# 20 tissue samples and 2308 predictors(regarding gene expression measurements)

length(Khan$ytrain)
# 63 observations in a vector

length(Khan$ytest)
# 20 observations ina vector

# So, the data set consists of expression measurements of 2,308 genes.

table(Khan$ytrain)
#1  2  3  4 
#8 23 12 20

# Note here, 1:4 are the four distinct types of cell tumors and their distribution is shown in the train set.

table(Khan$ytest)
#1 2 3 4 
#3 6 6 5

# Here the distribution of the four distinct cell tumors are shown in the test set.

# Now, we'll use a Support Vector approach to predict cancer subtype using the gene expression measurements.

# Note: In this data set, there are a very large number of features relative to the number of observations.
# This suggests that we should use a linear kernel, because the additional flexibility that will result from
# using a polynomial or radial kernel is unnecessary.

dat = data.frame(x = Khan$xtrain, y = as.factor(Khan$ytrain))

out = svm(y ~ . , data = dat, kernel = "linear", cost = 10)

summary(out)

#Call:
#  svm(formula = y ~ ., data = dat, kernel = "linear", cost = 10)
#
#
#Parameters:
#  SVM-Type:  C-classification 
#SVM-Kernel:  linear 
#cost:  10 
#gamma:  0.0004332756 
#
#Number of Support Vectors:  58
#
#( 20 20 11 7 )
#
#
#Number of Classes:  4 
#
#Levels: 
#  1 2 3 4

names(out)

length(out$fitted)  # 63 observations (of four different classes of tumors)

table(pred = out$fitted, truth = dat$y)

#      truth
#pred  1  2  3  4
#   1  8  0  0  0
#   2  0 23  0  0
#   3  0  0 12  0
#   4  0  0  0 20

# We see that there is no training errors. In fact, this is not surprising, because the large number of variables
# relative to the number of observations implies that it is easy to find hyperplanes that fully separate the classes.

# But as usual we are more interested in knowing the SVM's performance on the test data set.

data.test = data.frame(x = Khan$xtest, y = as.factor(Khan$ytest))

pred.test = predict(out, newdata = data.test)

table(pred = pred.test, truth = data.test$y)


#     truth
#pred 1 2 3 4
#   1 3 0 0 0
#   2 0 6 2 0
#   3 0 0 4 0
#   4 0 0 0 5

mean(pred.test == data.test$y)

# So, the correct prediction rate on the test data set is 90 % using a linear kernel and cost = 10