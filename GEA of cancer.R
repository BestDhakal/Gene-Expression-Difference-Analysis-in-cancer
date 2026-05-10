# gene expression differences analysis in cancer using R
                #Prabhat_Dhakal_2026 May
  #Learned from Hamidreza Bolhasani Youtube Channel

# Step 1 :: installing and loading packages  -----
#if the below packages are not install in R, you need to install before loading the library
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!require("BiocGenerics", quietly = TRUE))
    install.packages("BiocGenerics")
 #Installing BiocManager package
BiocManager::install("GEOquery")
  
#loading the installed packages,
library("BiocManager")   #to insatall Bioconductor packages
library("reshape2")  # for data reshaping
library("ggplot2")   #for data visualization
library("GEOquery")   #for downloading GEO data
library("limma")   # for differential expression analysis

# Step 2 :: Data Acquisition ----
  #here GSE is the local database where the gene dataset can be downloaded,
  # creating objects to store datasets
gseid <- "GSE72094"
gsedata <- getGEO(gseid, GSEMatrix = TRUE)
if (length(gsedata) == 0) {
  stop("Dataset not found.")
}
exprdata <- exprs(gsedata[[1]])
sampleinfo <- pData(gsedata[[1]])

          ## Step 3 ##

# Step 3 : Data preprocessing ----
    #reorderd sampleinfo using match() function to mtach the column names of exprdata 
sampleinfo <- sampleinfo[match(colnames(exprdata), rownames(sampleinfo)), ]
  # checking for the presence of the TP53 status ; if not found ,show error message 
if (!"tp53_status:ch1" %in% colnames(sampleinfo)) {
  stop("TP53 status column not found. Check column names.")
  }
  # defining condition based on TP54_status Labeled samples as Tumor , if 
condition <- ifelse(sampleinfo$'tp53_status:ch1' == "Mut", "Tumor", "Normal")
  #creating a sample table data frame with condition as a fcator
sample_table <- data.frame(condition = factor(condition))
  #setting row names to sample table , to match the column names of exprdata
row.names(sample_table) <- colnames(exprdata)

      ## Step 4 ##

# Step 4 : Differential Expression Analysis ----
  #creating design matrix using model.matrix(), to represent the relationship condition and the expression data
design <- model.matrix(~ condition, data = sample_table)
  #applying linear modeling with lmfit(), to fit the expression data to the design matrix
fit <- lmFit(exprdata, design)
  #empirical bayes moderation on the fitted model using eBayes()
fit <- eBayes(fit)
  #generating results table , for top 10 features using topTable()
results <- topTable(fit, coef = 2, number = 10)
  #printing or views few top rows using head()
print(head(results))

# Step 5 : Data visulization ----
 
 #A. using volcano plot to display significance and magnitude of changes in gene expression
volcanoplot(fit, coef = 2, main = "Volcano Plot of TP53 Expression", highlight = 10)
topgenes <- topTable(fit, coef = 2, number = 10)
topgeneids <- rownames(topgenes)
topexprdata <- exprdata[topgeneids, ]
exprmelted <- melt(topexprdata)
exprmelted$condition <- rep(sample_table$condition, each = nrow(topexprdata))

  #B. creating box plot of to 10 differentially expressed genes
ggplot(exprmelted, aes(x = Var2, y = value, fill = condition)) +
  geom_boxplot() +
  labs(x = "Sample", y = "Expression Level", title = "Boxplot of Top 10 Differentially Expressed Genes") +
  facet_wrap(~ Var1, scales = "free_y") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), 
        axis.text = element_text(size = 10))
colorvector <- ifelse(sample_table$condition == "Tumor", "red", "blue")

  #C. use of boxplot to compare TP53 expression levels betwn normal and cancer cell
boxplot(exprdata[topgeneids, ],
        main = "Gene Expression Boxplot", 
        xlab = "Samples", 
        ylab = "Expression Levels", 
        col = colorvector,
        las = 2)
legend("topright", legend = c("Tumor", "Normal"), fill = c("red", "blue"))


# Step 6 : Data Interpretation
 # A. Volcano plot analysis
  ## log2 fold change ranges from -1.5 to 1.0
  ## Statistically significant changes occur above 1.3 on the y axis
  ## Transcripts 37404 and 34113 showed the greatest significance, with a balanced expression of up and downregulated genes

  # B. Box plot analysis
  ## Expression levels ranges btween 4 abd 12 log scale for both lung cancer and normal tissue
  ## Tumor and normal samples show overlapping expression ranges, with outlier suggesting potential biomarkers
  ## subtle differences in expression patterns hint at underlying regulatory mechanism in tumor tissues

# Conclusion : Lung cancer shows significant dyresulation of tp53 related transcripts , with both up and down regulated genes
