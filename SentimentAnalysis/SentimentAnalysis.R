#=============================================================================#
# Twitter Data Sentiment Analysis                                             #
#=============================================================================#

#=============================================================================#
# Import Libraries:

# if not already installed, uncomment and run these commands
install.packages("e1071")
install.packages("text2vec")
install.packages("glmnet")
install.packages("randomForest")
install.packages("tm")

# change the following to set directory to preferred destination
setwd("~/R") # twitter data required to be in a directory called "R"
library(e1071)
library(text2vec)
library(glmnet)
library(randomForest)
library(tm)

# loading the twitter train, development and test datasets
YTrain = read.table("train-labels.txt", sep = "\t", header = F)
YTrain = YTrain[,2]
XTrain = readLines("train-tweets.txt")
YDev = read.table("dev-labels.txt", sep = "\t", header = F)
YDev = YDev[,2]
XDev = readLines("dev-tweets.txt")
XTest = readLines("test-tweets.txt")

#=============================================================================#
# Preprocessing of the data: Vectorisation and term document matrix build

prepFunc = tolower
tokFunc = word_tokenizer 

TokenTrain = XTrain %>%
  prepFunc %>%
  tokFunc

TokTrain = itoken(TokenTrain, ids = YTrain, progressbar = F)
vocabTrain = create_vocabulary(TokTrain)

vectoriser  = vocab_vectorizer(vocabTrain)
termMatTrain = create_dtm(TokTrain, vectoriser)

TokenDev = XDev %>%
  prepFunc %>%
  tokFunc
TokDev = itoken(TokenDev, ids = YDev, progressbar = F)
termMatDev = create_dtm(TokDev, vectoriser)

#=============================================================================#
# Training and testing a LASSO-penalised GLM model:

M1.1 = cv.glmnet(termMatTrain, YTrain, family="multinomial",
                    type.multinomial = "grouped", parallel = TRUE,
                    nfolds = 5, type.measure = "deviance")
plot(M1.1)

M1.1Pred = predict(M1.1, termMatDev,  s = "lambda.min",type = "class")
glm1Error = sum(M1Pred != YDev) / length(YDev)
glm1Error

# removing terms with a low frequency. 
prunedVocab = prune_vocabulary(vocabTrain, doc_proportion_min = 0.01)
vectoriser2 = vocab_vectorizer(prunedVocab)
# create dtm_train with new pruned vocabulary vectoriser
termMatTrain2  = create_dtm(TokTrain, vectoriser2)
termMatDev2 = create_dtm(TokDev, vectoriser2)
dim(termMatTrain)
dim(termMatTrain2)

M1.2 = cv.glmnet(termMatTrain2, YTrain, family="multinomial",
               type.multinomial = "grouped", parallel = TRUE,
               nfolds = 5, type.measure = "deviance")

M1.2Pred = predict(M1.2, termMatDev2,  s = "lambda.min",type = "class")
glm2Error = sum(M1.2Pred != YDev) / length(YDev)
glm2Error
plot(M1.2)

#=============================================================================#
# Train and test a naive bayes model

M2 = naiveBayes(as.factor(YTrain), as.matrix(TermMatTrain) ,laplace = 3)

M2Pred = predict(M2, as.matrix(TermMatDev), type = "class")
nbError = sum(M2Pred != YDev) / length(YDev)
nbError

#=============================================================================#
# Train and test a random forest model

M3 = randomForest(as.factor(YTrain), TermMatTrain)
M3Pred = predict(M3, as.matrix(TermMatDev), type = "class")
rfError = sum(M3Pred != YDev) / length(YDev)
rfError
