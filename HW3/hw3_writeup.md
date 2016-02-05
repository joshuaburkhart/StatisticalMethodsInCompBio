write [#] references and fill in XXX numbers

> Input datasets are as below:
> 
> - a: raw intensity values from the original CEL files arranged in a matrix layout, where each column represents one hybridization, and rows stand for individual array features
> - x: RMA normalized dataset in the assayData and annotation in the phenoData
> - xq: single-cell gene expression levels measured by qPCR
> - xql: single-cell gene expression measured by qPCR, with cells facing the blastocyst cavity labelled fluorescently
> 
> Let us assume our final goals are to build 1) a classifier that can be generalized for the microarray platform & probeset used in x and 2) a classifier that can be generalized for the qPCR assay used in xq that can predict the developmental stage better than each dataset's class distribution[1] and to compare their performance. We'd like to not only correctly label test data as "E3.25" or not in each experiment but build classifiers whose results may be interpreted by experts in the field to gain knowledge about any biological interactions that occur between selected features and developmental stage. We'll use the RMA normalized (x) dataset for 1 and the xq dataset for 2. The xql dataset would be of interest if it contained our target developmental stage, which it does not[2].
> 
> Though we're ignoring EDA and potential batch effects for the moment, we'll still reject metadata such as microarray date, number of cells, or index in the original dataset as feature candidates due to our noted constraint of creating a model that may increase biological insight into developmental stage. Additionally, we'll restrict our feature selection strategies to only those who avoid transforms, like PCA, and we'll restrict our models to avoid those whose results are notoriously difficult  to interpret, like ANN.
> 
> Our null hypothesis shall be: All probe intensities in dataset x and gene expression values in dataset xq are evenly and randomly distributed between developmental stage E3.25 and the other developmental stages.
> 
> Analysis Plan:
> 
> 1. Extract feature candidates, probe intensities, gene expression levels, etc., from their original data sources and combine them into a single data table. This will allow for simpler handling for the remainder of the analysis. 
>
> 2. Scan & address missing data. Because we have so few samples and so many total features, we'll prefer to drop features than drop samples. We must scan probe intensities, gene expression levels, etc..
>
> 3. Convert feature datatypes as nessecary. For example, genotype is stored as one of two values: "FGF4-KO" and "WT". We'll change "FGF4-KO" to 1 and "WT" to 0.
>
> 4. Extract development stage and convert it to a numeric value. Currently, Embryonic.day is stored as one of three values: "E3.25", "E3.5", and "E4.5". We'll copy this column to a separate vector and store our positive class, "E3.25", as 1 and our negative class, "E3.5" and "E4.5", as 0.
>
> 5. Randomly split data into training and test sets. We will assure similar proportians of genotype and class are split between training and test sets. Also, seeing as we have so few samples, we'll use an 80%/20% training/test split.
>
> 6. Randomly split training data to prepare for 10-fold cross validation. This will require a 90/10 validation/training split.
>
> 7. Z score transform all the training feature values, storing the means and standard deviations for all columns so validation & test data can be transformed similarly later on. Putting features on the same scale allows some optimization algorithms, such as gradient descent, to converge more quickly and generally allows for more direct comparisons between feature distributions.
>
> 8. We'd like to explain a lot of variance without using too many features, as that can make things confusing. Normally we'd favor PCA as a logical next step to a learning system but, keeping in mind we'd like to create a model that's interpretable to biologists later on, it will be better to avoid it. Instead, we'll use another nonparametric multivariate filter that avoids feature transforms, a Markov Blanket[3]. The difficulty with this technique is in constructing the nessecary bayesian net from which to extract the blanket. For the x dataset we'll need to filter the features apriori.
> 
> 9. While the Markov Blanket feature selection approach will work well for the xq dataset containing only XXX features, the construction of a bayesian net for the x dataset (containing XXX features) requires unavailable computational resources. Thus, we'll use a nonparametric univariate filter, the Wilcoxon Rank Sum, to reduce our feature set before hand.
> 
> 10. As there has been skepticism regarding the use of the Wilcoxon Rank Sum test to select features[4], we'll remember that an accepted parametric univariate filter method, linear model fitting, relies on the OLS optimization which assumes normally distributed feature values, while the Wilcoxon Rank Sum test makes no such assumption about the distribution of feature values as it simply relies on the ranking of feature values to discover features whose values for one class are greater or less than those of another. Also, we'll show that linear model fitting finds a similar set of features as the Wilcoxon Rank Sum approach, providing at least some sense of validation. 
>
> 11. After extracting our Markov Blankets for each dataset, we'll train a support vector machine (SVM) with a linear kernel. This model selection was made because SVM's with linear kernels are thought to be rather straight forward to interpret, which is a main concern.
>
> 12. The trained SVM will be run against the previously separated validation data and its error rate, a simple misclassification rate, will be recorded for each of the 10 cross-validation folds, along with the generated feature set and SVM.
> 
> 13. Following the cross-validation, the feature sets and SVM's from the fold with the lowest error rate will be selected as our 'best model'.
> 
> 14. To test our model, we'll apply Z score transformations to test data using stored means and standard deviations from training data, select the features used in our most accurate fold, run the SVM from the best fold, and report our misclassification rate.

## (2) Write R script for performing EDA on the data set and then perform EDA

```{r global_options, echo=FALSE, include=FALSE, error=FALSE}
knitr::opts_chunk$set(fig.path = "Figs/",
                      message = FALSE,
                      warning = FALSE,
                      include = TRUE,
                      echo = TRUE,
                      error = TRUE,
                      fig.width = 11,
                      comment = NA)
```

```{r, echo=FALSE, include=FALSE}
library(MASS)
library(plyr) #this must be loaded before dplyr 
library(dplyr)
library(ggplot2)
library(broom)
library(knitr)
library(magrittr)
library(reshape2)
library(infotheo)
library(stats)
library(ggbiplot)
library(car)
library("Hiiragi2013")
set.seed(2013) #does this need to be set to 2013? why?
```

### Data Source Description

> a: raw intensity values from the original CEL files arranged in a matrix layout, where each column represents one hybridization, and rows stand for individual array features  
  x: RMA normalized dataset in the assayData and annotation in the phenoData  
  xq: single-cell gene expression levels measured by qPCR  
  xql: single-cell gene expression measured by qPCR, with cells facing the blastocyst cavity labelled fluorescently
  
### Load Data
  
```{r}
data("x")
data("xq")
data("xql")

#feature candidates
hw3A.x.genotypes <- x@phenoData@data$genotype
hw3A.x.probe_intensities <- data.frame(assayDataElement(x@assayData,'exprs'))
hw3A.xq.cell_type <- xq@phenoData@data$Cell.type
hw3A.xq.gene_expressions <- data.frame(assayDataElement(xq@assayData,'exprs'))
hw3A.xql.label <- xql@phenoData@data$Label
hw3A.xql.gene_expressions <- data.frame(assayDataElement(xql@assayData,'exprs'))

#classes
hw3A.x.classes <- x@phenoData@data$Embryonic.day
hw3A.xq.classes <- xq@phenoData@data$Embryonic.day
hw3A.xql.classes <- xql@phenoData@data$Embryonic.day

#additional data we'll consider when testing for batch effects 
hw3A.x.pheno_dates <- as.Date(x@phenoData@data$ScanDate)
hw3A.x.proto_dates <- as.Date(x@protocolData@data$ScanDate)
hw3A.x.num_cells <- as.numeric(x@phenoData@data$Total.number.of.cells)
```

### Scan for missing Data

```{r, eval=FALSE}
sum(is.na(hw3A.x.genotypes))
sum(is.na(hw3A.x.probe_intensities))
sum(is.na(hw3A.xq.cell_type))
sum(is.na(hw3A.xq.gene_expressions))
sum(is.na(hw3A.xql.label))
sum(is.na(hw3A.xql.gene_expressions))

sum(is.na(hw3A.x.classes))
sum(is.na(hw3A.xq.classes))
sum(is.na(hw3A.xql.classes))
```

> (Output omitted) 302 na's reported in the xq gene expression data but nowhere else. The na's are from Spp1 (40), Prdm14 (40), Sdc4 (40), Morc1 (67), Tbpl1 (2), and Zp3 (113). For now we'll remove them. We may want to replace them with dummy values later as we have so few features for the xq dataset (38 gene expression levels and cell type).

```{r, eval=FALSE}
t(hw3A.xq.gene_expressions) %>% summary()
hw3A.xq.gene_expressions <- na.omit(hw3A.xq.gene_expressions)
```

### Checking Distributions

> Check outcome ratios. None show a uniform distribution and xql doesn't contain E3.25 developmental stage results. We can safely remove it from our analysis going forward.

```{r, echo=FALSE}
plot(x=hw3A.x.classes,
     xlab="Developmental Stage",
     ylab="Abundance",
     main="Distribution of Developmental Stages in Dataset x")
plot(x=as.factor(hw3A.xq.classes),     
     xlab="Developmental Stage",
     ylab="Abundance",
     main="Distribution of Developmental Stages in Dataset xq")
plot(x=as.factor(hw3A.xql.classes),
     xlab="Developmental Stage",
     ylab="Abundance",
     main="Distribution of Developmental Stages in Dataset xql")
```

> Check feature vectors. Genotypes are skewed. The probe intensities overall are slightly skewed while the probe intensity distributions all appear uniformly distributed with only two having noticable differences. It's odd to have noticably different distributions as they were all supposed to have been normalized in the x dataset. The gene expressions show wide variation among samples.

```{r, echo=FALSE}
#x
plot(x=hw3A.x.genotypes,
     xlab="Genotype",
     ylab="Abundance",
     main="Distribution of Genotypes in Dataset x")

#probe intensities contain too much data for plot()
range(hw3A.x.probe_intensities) %>% invisible()
rdiff <- round(range(hw3A.x.probe_intensities)[2] - range(hw3A.x.probe_intensities)[1])
bins <- seq(1:rdiff)
Boxplot(.bincode(hw3A.x.probe_intensities[,1],bins),
     ylab="Intensity",
     main="Distribution of Probe Intensities in Dataset x") %>% invisible()
mm <- melt(hw3A.x.probe_intensities)
ggplot(mm) +
  geom_boxplot(
    aes(y=round(value),x=variable,colour=floor(value))) +
  labs(x="Sample",
       y="Probe Intensity Distribution",
       title="Probe Intensity Distributions Across Samples")

#xq
plot(as.factor(hw3A.xq.cell_type),
     xlab="Cell Type",
     ylab="Abundance",
     main="Distribution of Cell Type in Dataset xq")

#gene expressions too much data for plot()
range(hw3A.xq.gene_expressions) %>% invisible()
rdiff <- round(range(hw3A.xq.gene_expressions)[2] - range(hw3A.xq.gene_expressions)[1])
bins <- seq(1:rdiff)
Boxplot(.bincode(hw3A.xq.gene_expressions[,1],bins),
     ylab="Expression Level",
     main="Distribution of Gene Expression Levels in Dataset xq") %>% invisible()
mm2 <- melt(hw3A.xq.gene_expressions)
ggplot(mm2) +
  geom_boxplot(
    aes(y=round(value),x=variable,colour=floor(value))) +
  labs(x="Sample",
       y="Gene Expression Distribution",
       title="Gene Expression Distributions Across Samples")
```

> Test distributions vs outcomes. These may not be visualized well but genotypes and cell types can be tested easily using Chi-squared tests. Genotypes are fairly well balanced among developmental stages but cell types are skewed.

```{r}
#x
table(x=hw3A.x.genotypes,y=hw3A.x.classes) %>% chisq.test()
```

```{r, echo=FALSE}
smush <- t(data.frame(hw3A.x.probe_intensities))
mm3 <- melt(smush)
mm3$class <- hw3A.x.classes
ggplot(mm3) +
  geom_boxplot(aes(y=round(value), x=Var1, colour=floor(value))) +
  facet_wrap(~ class) +
  labs(x="Developmental Stage",
       y="Probe Intensity Distributions Across Samples",
       title="Probe Intensity Distributions Across Samples Across Developmental Stages in Dataset x")
```

```{r}
#xq
table(hw3A.xq.cell_type,hw3A.xq.classes) %>% chisq.test()
```

```{r, echo=FALSE}
smush <- t(data.frame(hw3A.xq.gene_expressions))
mm4 <- melt(smush)
mm4$class <- hw3A.xq.classes
ggplot(mm4) +
  geom_boxplot(aes(y=round(value), x=Var1, colour=floor(value))) +
  facet_wrap(~ class) +
  labs(x="Developmental Stage",
       y="Gene Expression Distributions Across Samples",
       title="Gene Expression Distributions Across Samples Across Developmental Stages in Dataset xq")
```

### Test for Batch Effects

> Plotting array index vs developmental stage. Note xq doesn't store Embryonic.day as a factor like x does. We should consider the regularity of the positions in these datasets when separating our data into training and test sets. Correlation coefficient agrees in the latter case, though visual inspection of the first plot also shows obvious bias.

```{r, echo=FALSE}
# x
x.idx <- seq(1,length(hw3A.x.classes))
plot(x=x.idx,
     y=hw3A.x.classes,
     xlab="Index",
     ylab="Developmental Stage",
     main="Developmental Stages Across Index in Dataset x")
```

```{r}
cor(x.idx,as.numeric(hw3A.x.classes))
```
```{r, echo=FALSE}
# xq
xq.idx <- seq(1,length(hw3A.xq.classes))
plot(x=xq.idx,
     y=as.factor(hw3A.xq.classes),
     xlab="Index",
     ylab="Developmental Stage",
     main="Developmental Stages Across Index in Dataset xq")
```

```{r}
cor(xq.idx,as.numeric(as.factor(hw3A.xq.classes)))
```

> Plotting dates vs developmental stage for dataset x. Unfortunately, developmental stage appears highly correlated with both date variables over several years. This indicates unknown and unrecorded variables could be influencing the results of our data. Chi-squared tests agree in both cases, shown by extreemly low p-values.

```{r, echo=FALSE}
# x phenotypic scan dates
plot(x=hw3A.x.pheno_dates,
     y=hw3A.x.classes,
     xlab="Phenotypic Scan Date",
     ylab="Developmental Stage",
     main="Phenotypic Scan Dates Across Developmental Stages in Dataset x")
```

```{r}
table(hw3A.x.pheno_dates,hw3A.x.classes) %>% chisq.test()
```

```{r, echo=FALSE}
# x protocol scan dates
plot(x=hw3A.x.proto_dates,
     y=hw3A.x.classes,
     xlab="Protocol Scan Date",
     ylab="Developmental Stage",
     main="Protocol Scan Dates Across Developmental Stages in Dataset x")
```

```{r}
table(hw3A.x.proto_dates,hw3A.x.classes) %>% chisq.test()
```

> Plotting number of cells vs developmental stage. Like the dates above, the number of cells appears to correlate well with developmental stage. And again, a Chi-squared test agrees.

```{r, echo=FALSE}
# x number of cells
plot(x=hw3A.x.num_cells,
     y=hw3A.x.classes,
     xlab="Number of Cells",
     ylab="Developmental Stage",
     main="Number of Cells Across Developmental Stages in Dataset x")
```

```{r}
table(hw3A.x.num_cells,hw3A.x.classes) %>% chisq.test()
```

## (3) Create version 2 of analysis plan that includes EDA and any revisions (if needed) based on EDA 

> Considering the serious batch effects discovered during EDA, we should include a step to address this. Also, it was previously thought that the datasets x, xq, and xql contained the same samples and that their probe intensities, gene expression, etc. could be combined into a single table. After examining the available data, there appears no way to map data from one dataset to another, and xql doesn't have any "E3.25" developmental stage results, thus two separate learning experiments will have to be performed. We'll still use decision tree, logistic regression, and support vector machine for each dataset and still use a training test data split of 80/20.

> 0. Address batch effects of both scan dates and numbers of cells.  
> 1-10 (as above)

## (4) Write R script for analyzing the data 

### Transform positive (E3.25) and negative (E3.5 and E4.5) values 1 and 0, respectively

```{r}
hw3A.x.classes <- as.character(hw3A.x.classes)
hw3A.x.classes[hw3A.x.classes == "E3.25"] = 1
hw3A.x.classes[hw3A.x.classes == "E3.5"] = 0
hw3A.x.classes[hw3A.x.classes == "E4.5"] = 0
hw3A.x.classes <- as.factor(hw3A.x.classes)

hw3A.xq.classes <- as.character(hw3A.xq.classes)
hw3A.xq.classes[hw3A.xq.classes == "E3.25"] = 1
hw3A.xq.classes[hw3A.xq.classes == "E3.5"] = 0
hw3A.xq.classes[hw3A.xq.classes == "E4.5"] = 0
hw3A.xq.classes <- as.factor(hw3A.xq.classes)
```

### Combine Features for x and xq datasets

```{r}
hw3A.x.features <- rbind(hw3A.x.genotypes,hw3A.x.probe_intensities)
hw3A.xq.features <- rbind(as.factor(hw3A.xq.cell_type),hw3A.xq.gene_expressions)
```

### Split into Training & Test sets

```{r}
# x
index <- sample(1:ncol(hw3A.x.features),round(0.8*ncol(hw3A.x.features)))
hw3A.x.features_train <- hw3A.x.features[,index]
hw3A.x.classes_train <- hw3A.x.classes[index]
hw3A.x.features_test <- hw3A.x.features[,-index]
hw3A.x.classes_test <- hw3A.x.classes[-index]

# xq
index <- sample(1:ncol(hw3A.xq.features),round(0.8*ncol(hw3A.xq.features)))
hw3A.xq.features_train <- hw3A.xq.features[,index]
hw3A.xq.classes_train <- hw3A.xq.classes[index]
hw3A.xq.features_test <- hw3A.xq.features[,-index]
hw3A.xq.classes_test <- hw3A.xq.classes[-index]
```

### Z-score training set

```{r}
# x
x.train_sd = matrix()
x.train_mean = matrix()
hw3A.x.features_train_normalized <- hw3A.x.features_train
for(i in 2:nrow(hw3A.x.features_train)) # each row is a feature, first is a factor (genotype)
  x.train_sd[i] = sd(hw3A.x.features_train[i,])
  x.train_mean[i] = apply(hw3A.x.features_train[i,],1,mean)
  for(j in 1:ncol(hw3A.x.features_train)) # each column is a sample
    hw3A.x.features_train_normalized[i,j] =
    (hw3A.x.features_train[i,j] - x.train_mean[i]) / x.train_sd[i]

# xq  
xq.train_sd = matrix()
xq.train_mean = matrix()
hw3A.xq.features_train_normalized <- hw3A.xq.features_train
for(i in 2:nrow(hw3A.xq.features_train)) # each row is a feature, first is a factor (cell_type)
  xq.train_sd[i] = sd(hw3A.xq.features_train[i,])
  xq.train_mean[i] = apply(hw3A.xq.features_train[i,],1,mean)
  for(j in 1:ncol(hw3A.xq.features_train)) # each column is a sample
    hw3A.xq.features_train_normalized[i,j] =
    (hw3A.xq.features_train[i,j] - xq.train_mean[i]) / xq.train_sd[i]
```

> With our datasets now split and Z scored, we are now ready to move on to feature selection (in the next assignment).

## References

> [4] Dr. Shannon McWeeney, Personal Communication.
