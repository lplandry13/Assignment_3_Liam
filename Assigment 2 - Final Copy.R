####Assignment 2 - V2.0-Final####
####Supervised machine learning####
####Part 1: required R packages####

#Testedit

#install.packages("tidyverse")
library(tidyverse)

#install.packages("stringi")
library(stringi)

#install.packages("ape")
library(ape)

#install.packages("randomForest")
library(randomForest)

#install.packages("rentrez")
library(rentrez)

#install.packages("seqinr")
library(seqinr)

#install.packages("BiocManager")
library(BiocManager)
library(Biostrings)

#install.packages("caret")
library(caret)

#install.packages("tree")
library(tree)

#install.packages("e1071")
library(e1071)

####Part 2: NCBI data extraction####

#This first section of code is for searching the NCBI database for all sequences of both the COX1 and COI gene in the subfamily Leuciscinae. Since we do not know the sequence length of the genes, we will first run without the SLEN search term in order to generate a histogram of all sequence lengths. This will allows us to determine what to set the range of our SLEN search to.

#First we will run a search using entrez for our subfamily and first gene. We are using the use_history set to true to pull all possible hits for this search term and save it to the web in order to save space, since it will be a large file.

cytb_search_test <- entrez_search(db = "nuccore", term = "Leuciscinae[ORGN] AND CYTB[Gene]", use_history = T)


#Now repeat for our second gene.

coi_search_test <- entrez_search(db = "nuccore", term = "Leuciscinae[ORGN] AND COI[Gene]", use_history = T)

#For this next section, will be using the Entrez_Functions.R special functions made by Jacqueline May for extracting web history. The FetchFastaFiles function should be run into a separate working directory for each gene as the mergeFastaFiles function will merge all fasta files it detects in a folder.

source("Entrez_Functions.R")

FetchFastaFiles(searchTerm = "Leuciscinae[ORGN] AND CYTB[Gene]", seqsPerFile = 100, fastaFileName = "Leuciscinae_CYTB_test")

dfCYTB_test <- MergeFastaFiles(filePattern = "Leuciscinae_CYTB_test*")

FetchFastaFiles(searchTerm = "Leuciscinae[ORGN] AND COI[Gene]", seqsPerFile = 100, fastaFileName = "Leuciscinae_COI_test")

dfCOI_test <- MergeFastaFiles(filePattern = "Leuciscinae_COI_test*")

#Next is altering the data frame to make a new column for species names

dfCYTB_test$Species_Name <- word(dfCYTB_test$Title, 2L, 3L)

dfCYTB_test <- dfCYTB_test[, c("Title", "Species_Name", "Sequence")]
view(dfCYTB_test)

dfCOI_test$Species_Name <- word(dfCOI_test$Title, 2L, 3L)

dfCOI_test <- dfCOI_test[, c("Title", "Species_Name", "Sequence")]
view(dfCOI_test)

#Here is the code to make histogram to look at spread of sequence length

hist(nchar(dfCYTB_test$Sequence), xlab = "Sequence Length", main =  "Histogram of CYTB Sequence Length")

#This histogram shows a high frequency at the 0:2500 range, and a very low frequency above 15000. The high frequency is likely the CYTB gene, with the low frequency being the full mitochondrial genome.

#Now check for the COI gene.

hist(nchar(dfCOI_test$Sequence), xlab = "Sequence Length", main = "Histogram of COI Sequence Length")

#This shows a very similar result to the COX1 histogram, with the high frequency range a little smaller, at 0:1000. Since there are no sequences between the high and low frequency, we will narrow our search by setting a sequence length of 0:2500 for both genes, to ensure that we grab all instances of the genes without the full mitochondrial genome.
#To do this, we will repeat all of the code from before, with the addition of the sequence length search term.

cytb_search_main <- entrez_search(db = "nuccore", term = "Leuciscinae[ORGN] AND CYTB[Gene] AND 0:2500 [SLEN]", use_history = T)

FetchFastaFiles(searchTerm = "Leuciscinae[ORGN] AND CYTB[Gene] AND 0:2500 [SLEN]", seqsPerFile = 100, fastaFileName = "Leuciscinae_CYTB_main")

dfCYTB_main <- MergeFastaFiles(filePattern = "Leuciscinae_CYTB_main*")

coi_search_main <- entrez_search(db = "nuccore", term = "Leuciscinae[ORGN] AND COI[Gene] AND 0:2500 [SLEN]", use_history = T)

FetchFastaFiles(searchTerm = "Leuciscinae[ORGN] AND COI[Gene] AND 0:2500 [SLEN]", seqsPerFile = 100, fastaFileName = "Leuciscinae_COI_main")

dfCOI_main <- MergeFastaFiles(filePattern = "Leuciscinae_COI_main*")

#Now to clean up the new data frames like we did before.

dfCYTB_main$Species_Name <- word(dfCYTB_main$Title, 2L, 3L)

#This time, let's also make a new column for the marker code.

dfCYTB_main$Marker_Code <- c("CYTB")

#And order the columns as before.

dfCYTB_main <- dfCYTB_main[, c("Title", "Species_Name", "Marker_Code", "Sequence")]
view(dfCYTB_main)

#Repeat for the COI gene.

dfCOI_main$Species_Name <- word(dfCOI_main$Title, 2L, 3L)

dfCOI_main$Marker_Code <- c("COI")

dfCOI_main <- dfCOI_main[, c("Title", "Species_Name", "Marker_Code", "Sequence")]

view(dfCOI_main)

#Now we generate the histograms again to ensure that the full genomes were filtered out. 

hist(nchar(dfCYTB_main$Sequence), xlab = "Sequence Length", main = "Histogram of CYTB Sequence Length")

hist(nchar(dfCOI_main$Sequence), xlab = "Sequence Length", main = "Histogram of COI Sequence Length")

#Zach - Could Consider including some summary statistics of the sequence lengths here. Since K-mer could be indirectly influenced by sequence length, this could influence classifier bias.

#Calculate in variables
mean_COI <- mean(nchar(dfCOI_main$Sequence))
sd_COI <- sd(nchar(dfCOI_main$Sequence))

mean_CYTB <- mean(nchar(dfCYTB_main$Sequence))
sd_CYTB <- sd(nchar(dfCYTB_main$Sequence))

#Print as a summary data frame
summ_table <- data.frame(
  Gene = c("COI", "CYTB"),
  Mean_Length = c(mean_COI, mean_CYTB),
  SD_Length = c(sd_COI, sd_CYTB)
)

print(summ_table)

#Zach - There is a lot more variation in CYTB, which is consistent with the literature. It would be interesting to use a subset of the data with matching sequence lengths to see if sequence length is an important contributor to the classifier through it's indirect effect on k-mer frequency.


#The filtering by sequence length has worked. We can now move forward knowing that the data we are looking at includes only the genes we are analyzing. Before moving on, remove the test data as it won't be necessary for further analysis.

rm(coi_search_test, cytb_search_test, dfCOI_test, dfCYTB_test)

#Now we will merge the two data frames for use in our classification tests.

dfAllSeqs <- rbind(dfCYTB_main, dfCOI_main)

view(dfAllSeqs)

#Use the sum function to check if any NAs are present. Sum should be zero.

sum(is.na(dfAllSeqs$Sequence))

#No NAs are present so we can move forward with the analysis.

#Check the class of the Sequence column.

class(dfAllSeqs$Sequence)

#This data is currently in character form, so we need to convert to DNA string set for use in our analysis.

dfAllSeqs <- as.data.frame(dfAllSeqs)

dfAllSeqs$Sequence <- DNAStringSet(dfAllSeqs$Sequence)

#Check the class again to verify that the conversion worked.

class(dfAllSeqs$Sequence)

####Part 3: random Forest machine learning####

#The first step for this part will be to calculate the frequency of each nucleotide in the sequences and attach those values into the data frame using the functon cbind.

dfAllSeqs <- cbind(dfAllSeqs, as.data.frame(letterFrequency(dfAllSeqs$Sequence, letters = c("A", "C","G", "T"))))

#Check the data frame to ensure that the new columns have been added.

view(dfAllSeqs)

#Now we will add the proportions of A, G, and T in relation to the total nucleotides (excluding C since it will be the remainder).

dfAllSeqs$Aprop <- (dfAllSeqs$A) / (dfAllSeqs$A + dfAllSeqs$T + dfAllSeqs$C + dfAllSeqs$G)

dfAllSeqs$Tprop <- (dfAllSeqs$T) / (dfAllSeqs$A + dfAllSeqs$T + dfAllSeqs$C + dfAllSeqs$G)

dfAllSeqs$Gprop <- (dfAllSeqs$G) / (dfAllSeqs$A + dfAllSeqs$T + dfAllSeqs$C + dfAllSeqs$G)

#Now check the data frame again to ensure the columns have been added.

view(dfAllSeqs)

#For further analysis, we will add nucleotide frequencies of k-mer lengths of 2, 3, and 4. Make sure to check data frame after each addition.

dfAllSeqs <- cbind(dfAllSeqs, as.data.frame(dinucleotideFrequency(dfAllSeqs$Sequence, as.prob = TRUE)))

dfAllSeqs <- cbind(dfAllSeqs, as.data.frame(trinucleotideFrequency(dfAllSeqs$Sequence, as.prob = TRUE)))

dfAllSeqs <- cbind(dfAllSeqs, as.data.frame(oligonucleotideFrequency(x = dfAllSeqs$Sequence, width = 4, as.prob = TRUE)))

#Now that we have all of this nucleotide data, we can covert the sequences back to character data for use with the tidyverse functions.

dfAllSeqs$Sequence <- as.character(dfAllSeqs$Sequence)

#Check count by marker code.

table(dfAllSeqs$Marker_Code)

#COI has the smaller sequence size, so we will use that as the basis for setting sample sizes for our analysis. First, we will write code to set the COI sequence amount size as a new variable smaller_sample for use in later code.

smaller_sample <- min(table(dfAllSeqs$Marker_Code))
smaller_sample

#For the random Forest analysis, we typically want to set aside about 20% of the data to use as a validation set. To do this, we will use set.seed and sample_n to randomly sample sequences equal to 20% of the value of smaller_sample and assign them to a new data frame that will be kept separate.

set.seed(001)

dfSequence_validation <- dfAllSeqs %>%
  group_by(Marker_Code) %>%
  sample_n(floor(0.2 * smaller_sample))

#Check to ensure it worked.

table(dfSequence_validation$Marker_Code)

#Next we will create the training data set, which will not overlap with the validation data set and will contain sequences equal to 80% of smaller_sample for each gene. We are using the same sample size for both genes to ensure class balance.

set.seed(002)

dfSequence_training <- dfAllSeqs %>%
  filter(!Title %in% dfSequence_validation$Title) %>%
  group_by(Marker_Code) %>%
  sample_n(ceiling(0.8 * smaller_sample))

#Check that it worked.

table(dfSequence_training$Marker_Code)

#Now we can create our gene classifier using random forest. We will first test it using only the A, T and G proportion columns (should be columns 9:11). We will first test with a ntree value of 500

gene_classifier <- randomForest::randomForest(x = dfSequence_training[, 9:11], y = as.factor(dfSequence_training$Marker_Code), ntree = 500, importance = TRUE)

#Time to check the results.

gene_classifier

#1.38% is alright, but lets keep trying things until it reaches close to 0%.

gene_classifier <- randomForest::randomForest(x = dfSequence_training[, 9:27], y = as.factor(dfSequence_training$Marker_Code), ntree = 500, importance = TRUE)

gene_classifier

#The addition of 2-mer data greatly improved the error rate. Let's add the 3-mer data to ensure the classifier works perfectly.

gene_classifier <- randomForest::randomForest(x = dfSequence_training[, 9:91], y = as.factor(dfSequence_training$Marker_Code), ntree = 500, importance = TRUE)

gene_classifier

#Close, but not quite there. Finally, let's add the 4-mer data to the classifier and see what we get.

gene_classifier <- randomForest::randomForest(x = dfSequence_training[, 9:347], y = as.factor(dfSequence_training$Marker_Code), ntree = 500, importance = TRUE)

gene_classifier

#The addition of the 4-mers did not improve the error rate any further, so instead let's increase the number of trees to 1000.

gene_classifier <- randomForest::randomForest(x = dfSequence_training[, 9:347], y = as.factor(dfSequence_training$Marker_Code), ntree = 1000, importance = TRUE)

gene_classifier

#still no improvement! Let's check at 1500 trees.

gene_classifier <- randomForest::randomForest(x = dfSequence_training[, 9:347], y = as.factor(dfSequence_training$Marker_Code), ntree = 1500, importance = TRUE)

gene_classifier

#The error rate was worse in that case, so let's use all proportion data and an ntree of 1000 as our best measure of classification. Run that line of code again before continuing. Now let's use the classifier on our validation data to ensure that it works as intended.

predict_validation <- predict(gene_classifier, dfSequence_validation[, c(3, 9:347)])

#Let's look at the results.

predict_validation
class(predict_validation)
length(predict_validation)

#To verify that the classifier worked, create a confusion matrix like the one output by gene_classifier.

table(observed = dfSequence_validation$Marker_Code, predicted = predict_validation)

#Perfect results. The random forest method of classification of unknown genes works as intended. 

####Part 4: Gradient Boosting####

#To verify if random forest is a good method of classification of unknown genes, we need to compare it to another method. In this case, we will be using a classification tree. This method is very similar to random Forest, but as you may guess, it generates one tree instead of many. This means that for lots of data, Forest is likely the better choice. But how about our data set? Let's test it.

#To ensure ease of running the data, let's simplify our data frame by removing all columns that we don't need.

dfAllSeqsTree <- dfAllSeqs %>%
  select(-Title, -Species_Name, -Sequence, -A, -C, -G, -T,)

#For prediction and the confusion matrix at the end to work, we need our Marker_Code data set as factor data.

dfAllSeqsTree$Marker_Code <- as.factor(dfAllSeqsTree$Marker_Code)

#We will now separate our training and validation data again, but this time we will use a different method that uses the createDataPartition function. This simplifies the code a bit from what we did for the random Forest. 

set.seed(1001)

intrain <- createDataPartition(y = dfAllSeqsTree$Marker_Code, p = 0.8, list = FALSE)

train <- dfAllSeqsTree[intrain, ]

test <- dfAllSeqsTree[-intrain, ]

#The first function we will use is the tree function, which will generate our initial tree with Marker_Code as the formula and our training data set as the data.

treemod <- tree(Marker_Code ~ ., data = train, )

#Next, we will use the function cv.tree to determine the number of misclassifications and then trim them from the tree using the function prune.tree.

cv.trees <- cv.tree(treemod, FUN = prune.tree)

#Let's plot the result to determine how to improve our classification.

plot(cv.trees)

#This plot shows us that we have the least deviances starting at size 3 of trees. So let's use the prune.tree function to narrow our classification down to improve the results.

prune.trees <- prune.tree(treemod, best = 3)

plot(prune.trees)

text(prune.trees, pretty = 0)

#This plot shows us the number of branches it took to classify all of the genes. As we can see, due to the amount of k-mer data we had available, the tree only need GGAA and AAA probabilities to classify the genes. We should now be able to run our prediction data through the tree.

treepred <- predict(prune.trees, newdata = test, type = "class")

#Now we can generate the confusion matrix, using a different line of code using the caret package.

confusionMatrix(treepred, test$Marker_Code)

#This method had 99.8% accuracy with all of the k-mer data and without any additions to run size in the tree function, so for a data set of this size, classification trees seem to be a good alternative to random Forest!