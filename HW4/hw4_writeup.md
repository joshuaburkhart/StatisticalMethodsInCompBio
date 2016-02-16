>Joshua Burkhart  
>2/18/2016  
>BMI 651
# HW4: Classification Trees

(100 pts) 

For this assignment, you will utilize the spam data set  (58 columns total, with 57 being features and last a factor : email or spam, and 4601 rows - representing emails ) 

You can download data at: http://statweb.stanford.edu/~tibs/ElemStatLearn/datasets/spam.data) 

See details at: http://statweb.stanford.edu/~tibs/ElemStatLearn/datasets/spam.info.txt

Develop an R script that does the following:

##1.

- Allows you to assess classifier you develop to a constant classifier (naïve case- always predicts the same class, no matter what the input features).

- For the constant classifier, it should  output what fraction are actually spam, what the constant classifier predicts and the error rate of the constant classifier

- Must provide assessment of which approaches below out-performs the constant classifier

##2.

- Splits the data set at random into a training set and testing set (1:1 split)

- checks to ensure they are non-overlapping

- split requested is correct

- compares percentage of spam emails in each set. 

##3.

- Classification tree for training data that you prune by cross-validation.

-  Script must generate plots for CV error versus tree size

- plot best tree and provide its error rate on the testing data. 

- It should also summarize the variables on the tree. 

##4.

- Bagging for ensemble of 100 trees for training data.

- Script must generate a plot of the importance of the variables, (for ensemble)

-  provide error rate (ensemble) on the testing data.

##5.

- Boosting for ensemble of 100 trees for training data. 

- Script must generate a plot of  importance of the variables

- provide error rate (ensemble) on the testing data.