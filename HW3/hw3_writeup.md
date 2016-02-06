write [R#] references, fill in XXX numbers, add [A#] appendix information

> Input datasets are as below:
> 
> - a: raw intensity values from the original CEL files arranged in a matrix layout, where each column represents one hybridization, and rows stand for individual array features
> - x: RMA normalized dataset in the assayData and annotation in the phenoData
> - xq: single-cell gene expression levels measured by qPCR
> - xql: single-cell gene expression measured by qPCR, with cells facing the blastocyst cavity labelled fluorescently
> 
> Let us assume our final goals are to build 1) a classifier that can be generalized for the microarray platform & probeset used in x and 2) a classifier that can be generalized for the qPCR assay used in xq that can predict the developmental stage better than each dataset's class distribution[A1] and to compare their performance. We'd like to not only correctly label test data as "E3.25" or not in each experiment but build classifiers whose results may be interpreted by experts in the field to gain knowledge about any biological interactions that occur between selected features and developmental stage. We'll use the RMA normalized (x) dataset for 1 and the xq dataset for 2. The xql dataset would be of interest if it contained our target developmental stage, which it does not[A2].
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
> 8. We'd like to explain a lot of variance without using too many features, as that can make things confusing. Normally we'd favor PCA as a logical next step to a learning system but, keeping in mind we'd like to create a model that's interpretable to biologists later on, it will be better to avoid it. Instead, we'll use another nonparametric multivariate filter that avoids feature transforms, a Markov Blanket[R1]. The difficulty with this technique is in constructing the nessecary bayesian net from which to extract the blanket. For the x dataset we'll need to filter the features apriori.
> 
> 9. While the Markov Blanket feature selection approach will work well for the xq dataset containing only XXX features, the construction of a bayesian net for the x dataset (containing XXX features) requires unavailable computational resources. Thus, we'll use a nonparametric univariate filter, the Wilcoxon Rank Sum, to reduce our feature set before hand.
> 
> 10. As there has been skepticism regarding the use of the Wilcoxon Rank Sum test to select features[R2], we'll remember that an accepted parametric univariate filter method, linear model fitting, relies on the OLS optimization which assumes normally distributed feature values, while the Wilcoxon Rank Sum test makes no such assumption about the distribution of feature values as it simply relies on the ranking of feature values to discover features whose values for one class are greater or less than those of another. Also, we'll show that linear model fitting finds a similar set of features as the Wilcoxon Rank Sum approach, providing at least some sense of validation. 
>
> 11. After extracting our Markov Blankets for each dataset, we'll train a support vector machine (SVM) with a linear kernel. This model selection was made because SVM's with linear kernels are thought to be rather straight forward to interpret, which is a main concern.
>
> 12. The trained SVM will be run against the previously separated validation data and its error rate, a simple misclassification rate, will be recorded for each of the 10 cross-validation folds, along with the generated feature set and SVM.
> 
> 13. Following the cross-validation, the feature sets and SVM's from the fold with the lowest error rate will be selected as our 'best model'.
> 
> 14. To test our model, we'll apply Z score transformations to test data using stored means and standard deviations from training data, select the features used in our most accurate fold, run the SVM from the best fold, and report our misclassification rate.

### Missing Data

> 302 na's reported[A3] in the xq gene expression data but nowhere else. The na's are from Spp1 (40), Prdm14 (40), Sdc4 (40), Morc1 (67), Tbpl1 (2), and Zp3 (113). We'll remove them[A4].

### Class & Feature Distributions

> Developmental stage (class) distributions are not uniform[A5] but with so few samples, we won't remove any.
>
> Regarding feature vectors, genotypes are skewed[A6]. The probe intensities overall are slightly skewed[A7] while the probe intensity distributions all appear uniformly distributed with only two having noticable differences[A8]. It's odd to have noticably different distributions as they were all supposed to have been normalized in the x dataset. The gene expressions show wide variation among samples[A9].
>
> Class distributions among features may not be visualized well[A8] but genotypes and cell types can be tested easily using Chi-squared tests[A9]. Genotypes are fairly well balanced among developmental stages but cell types are skewed.

### Batch Effects

> Plotting array index vs developmental stage[A10]. We should consider the regularity of the positions in these datasets when separating our data into training and test sets, which can be addressed by randomly sampling when separating our training, test, and validation sets. Correlation coefficient agrees in the latter case[A11], though visual inspection of the first plot also shows obvious bias[A12].
>
> Plotting dates vs developmental stage for dataset x[A13]. Unfortunately, developmental stage appears highly correlated with both date variables over several years. This indicates unknown and unrecorded variables could be influencing the results of our data. Chi-squared tests agree in both cases, shown by extreemly low p-values[A14].
>
> Plotting number of cells vs developmental stage[A15]. Like the dates above, the number of cells appears to correlate well with developmental stage. And again, a Chi-squared test agrees[A16].
>
> Sadly, aside from randomly sampling when selecting our training, test, and validation sets, we're currently unable to adjust for the batch effects we discovered.

### Cross Validation

### Model Selection & Validation

### Model Assessment & Results

## References

> [] Dr. Shannon McWeeney, Personal Communication, Feb. 4, 2016.

## APPENDIX

> []
