---
title: "Gene-level differential expression analysis using DESeq2"
author: "Meeta Mistry"
date: "Friday February 19, 2016"
---

Approximate time: 2.5 hours

## Learning Objectives 

* Understanding the different components of differential expression analysis in the context of DESeq2
* Exploring different objects in DESeq2 
* Summarizing and filtering results to find significant DEGs
* Using data visualization to look at results


## DESeq2: Differential expression analysis

### Getting setup

Let's get started by opening RStudio and opening up the project that we created last lesson. 

1. Go to the File menu and select 'Open project ...'
2. Navigate to `~/Desktop/DEanalysis/` and double click on the `DEanalysis.Rproj` file

You should see your environment become populated with all of the variables created last lesson. The only thing that we will need to do is reload the required libraries:

```
library(ggplot2)
library(RColorBrewer)
library(DESeq2)
library(pheatmap)
```


### Running DESeq2

To run the differential expression pipeline on the raw counts in DESeq2, we use a **single call to the function `DESeq()`**. The required input is the `DESeqDataSet` object that we created in the last lesson. By re-assigning the results of the function back to the same variable name, we can continue to fill in the `slots` of our `DESeqDataSet` object.

	##Run analysis
	dds <- DESeq(dds)
 
This function will print out a message for the various steps it performs: 

```
estimating size factors
estimating dispersions
gene-wise dispersion estimates
mean-dispersion relationship
final dispersion estimates
fitting model and testing
``` 

<img src="../img/slide16+33_DGE.png" width="400">


**Everything from normalization to linear modeling was carried out by the use of a single function!** The results of each step were inserted into the object that you initialized.

![deseq1](../img/deseq_obj2.png)


> *NOTE:* There are individual functions available in DESeq2 that would allow us to carry out each step in the workflow in a step-wise manner, rather than a single call. We demonstrated one example when generating size factors to create a normalized matrix. By calling `DESeq()`, the individual functions for each step are run for you.


## Normalization

To normalize the count data DESeq2 calculates size factors for each sample, using the *median of ratios method*. Let's take a quick look at size factor values we have for each sample:

```
> sizeFactors(dds)
Mov10_kd_2 Mov10_kd_3 Mov10_oe_1 Mov10_oe_2 Mov10_oe_3 Irrel_kd_1 Irrel_kd_2 Irrel_kd_3 
 1.5646728  0.9351760  1.2016082  1.1205912  0.6534987  1.1224020  0.9625632  0.7477715 
 
```
 
These numbers should be identical to those we generated initially when we had run the function `estimateSizeFactors(dds)`. Take a look at the total number of reads for each sample using `colSums(counts(dds))`. *How do the numbers correlate with the size factor?*

> *NOTE:* it can be advantageous to calculate gene-specific normalization factors (size factors) to account for further sources of technical biases such as differing dependence on GC content, gene length or the like, and these can be supplied to DESeq2 instead of using the median of ratios method.


## Dispersion estimates

In our model, the **within group variability** is accounted for using the dispersion parameter. Dispersion estimates are computed **per gene**, because different genes naturally have a different scale of biological variability. DESeq2 does a first pass estimate on dispersion for each gene (using maximum-likelihood estimate), but with such small sample sizes we will make very bad estimates of gene-wise dispersion unless we **share information across genes**. The next step is therefore taking information from all gene dispersion estimates to shrink them to more reasonable values.

Let's take a look at the dispersion estimates for our data:

	# Plot dispersion estimates
	plotDispEsts(dds)
	

<img src="../img/plotDispersion.png"">
 

The black dots are the original estimates for each gene. The red smooth curve provides an accurate estimate for the expected dispersion value for genes of a given expression strength. The blue dots represent shrunken estimates. The circles indicate outliers, where we don't perform shrinkage. 

We use an empirical Bayes approach which lets the strength of shrinkage depend (i) on an estimate of how close true dispersion values tend to be to the fit and (ii) on the degrees of freedom. **Since we have a small sample size, for many genes we see quite a bit of shrinkage.**


## Identifying gene expression changes

We have three sample classes so we can make three possible pairwise comparisons:

1. Control vs. Mov10 overexpression
2. Control vs. Mov10 knockdown
3. Mov10 knockdown vs. Mov10 overexpression

We are really only interested in #1 and #2 from above. Using the design formula we provided `~sampletype`, DESeq 2 internally created the following design matrix:

```
   	      Intercept	sampletypecontrol sampletypeMOV10_knockdown	sampletypeMOV10_overexpression
Mov10_kd_2	 1		0		1		0
Mov10_kd_3	 1		0		1		0
Mov10_oe_1   1		0		0		1
Mov10_oe_2   1		0		0		1
Mov10_oe_3   1		0		0		1
Irrel_kd_1	 1		1		0		0
Irrel_kd_2	 1		1		0		0
Irrel_kd_3	 1		1		0		0	

```
This design matrix is now used to setup the contrasts to request the comparisons we want to make.


### Hypothesis testing: Wald test

To build a results table, we use the `results()` function on the `dds` object. Additionally we need to specify **which comparisons we are interested in** looking at. 

The comparisons are provided to DESeq2 in the form of **contrasts**, in one of three different ways. In this lesson we will demonstrate the method that is most intuitive. By providing contrasts we are telling DESeq2 **which coefficients to use for the hypothesis testing** procedure; this also corresponds to the headers in your design matrix. To find out how the coefficients are named we can use the `resultsNames()` function:

	# Find names of coefficients
	resultsNames(dds)

To specify the specific contrasts we are interested in, we need to provide the column names from the coefficients table as a list of 2 character vectors:

	## Define contrasts
	contrast_oe <- list( "sampletypeMOV10_overexpression", "sampletypecontrol")

**The order of the names, determines the direction of fold change that is reported.** The name provided in the second element is the level that is used to baseline. So for example, if we observe a log2 fold change of -2 this would mean the gene expression is lower in Mov10_oe relative to the control. Pass the contrast vector as an argument to the `results()` function:

	# Extract results table
	res_tableOE <- results(dds, contrast=contrast_oe)


Let's take a look at what information is stored in the results:

	head(res_tableOE)

```
log2 fold change (MAP): sampletype MOV10_overexpression vs control 
Wald test p-value: sampletype MOV10_overexpression vs control 
DataFrame with 6 rows and 6 columns
               baseMean log2FoldChange      lfcSE       stat    pvalue       padj
              <numeric>      <numeric>  <numeric>  <numeric> <numeric>  <numeric>
1/2-SBSRNA4  45.6520399     0.26976764 0.18775752  1.4367874 0.1507784 0.25242910
A1BG         61.0931017     0.20999700 0.17315013  1.2128030 0.2252051 0.34444163
A1BG-AS1    175.6658069    -0.05197768 0.12366259 -0.4203185 0.6742528 0.77216278
A1CF          0.2376919     0.02237286 0.04577046  0.4888056 0.6249793         NA
A2LD1        89.6179845     0.34598540 0.15901426  2.1758136 0.0295692 0.06725157
A2M           5.8600841    -0.27850841 0.18051805 -1.5428286 0.1228724 0.21489067
```
> *NOTE:* The results table looks very much like a data frame and in many ways it can be treated like one (i.e when accessing/subsetting data). However, it is important to recognize that it is actually stored in a `DESeqResults` object. When we start visualizing our data, this information will be helpful. 


Let's go through some of the columns in the results table to get a better idea of what we are looking at. To extract information regarding the meaning of each column we can use `mcols()`:

	mcols(res_tableOE, use.names=T)

* `baseMean`: mean of normalized counts for all samples
* `log2FoldChange`: log2 fold change
* `lfcSE`: standard error
* `stat`: Wald statistic
* `pvalue`: Wald test p-value
* `padj`: BH adjusted p-values
 

***

**Exercise**

1. Create a contrasts vector called `contrast_kd` for the Mov10_knockdown comparison to control.
2. Use that contrasts vector to extract a results table and store that to a variable called `res_tableKD`.  
3. Create a contrasts vector for the Mov10_overexpression comparison to *all other samples*.

*** 


### Summarizing results and identifying DEGs

To summarize the results table, a handy function in DESeq2 is `summary()`. Confusingly it has the same name as the function used to inspect data frames. This function when called with a DESeq results table as input, will summarize the results at a given FDR threshold. 

	## Summarize results
	summary(res_tableOE)
	

```  
out of 19748 with nonzero total read count
adjusted p-value < 0.1
LFC > 0 (up)     : 3657, 19% 
LFC < 0 (down)   : 3897, 20% 
outliers [1]     : 0, 0% 
low counts [2]   : 3912, 20% 
(mean count < 4)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results
```

In addition to the number of genes up- and down-regulated at and FDR < 0.1, the function also reports the number of genes that were tested (genes with non-zero total read count), and the number of genes not included in multiple test correction due to a low mean count (which in our case is < 4).

The default FDR threshold is set to `alpha = 0.1`, which is quite liberal. Let's try changing that to `0.05` -- *how many genes are we left with*?

The FDR threshold on it's own doesn't appear to be reducing the number of significant genes. With large significant gene lists it can be hard to extract meaningful biological relevance. To help increase stringency, one can also add a fold change threshold. The `summary()` function doesn't have an argument for fold change threshold, but instead we can use the base R function `subset()`.

Let's first create variables that contain our threshold criteria:

	### Set thresholds
	padj.cutoff <- 0.05
	lfc.cutoff <- 1

The `lfc.cutoff` is set to 1; remember that we are working with log2 fold changes so this translates to an actual fold change of 2 which is pretty reasonable. Now let's setup our **`subset()` function**. Start building from the inside out:

	subset(res_tableOE)

We need to add our selection criteria. The first is our FDR threshold:

	subset(res_tableOE, padj < padj.cutoff)

Now let's add in the log2 fold change criteria. Because we want both up- and down-regulated genes we will use the absolute value of the fold change using the `abs(log2FoldChange)` function:

	subset(res_tableOE, padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff)

Now, finally we will put all of that inside the `summary()` function. This is a fast way of getting overall statistics and deciding whether our threshold is still too liberal or perhaps overly stringent.

	summary(subset(res_tableOE, padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff), alpha =0.05)


Does this reduce our results? How many genes are up-regulated and down-regulated at this new threshold?

We should have a total of 884 genes (682 up-regulated and 202 down-regulated) that are significantly differentially expressed. To denote these genes as significant we can add a column in our results table. The column will be a logical vector, where `TRUE` means the gene passes our threshold and `FALSE` means it fails.

	# Add a threshold vector
	threshold <- res_tableOE$padj < padj.cutoff & 
                   abs(res_tableOE$log2FoldChange) > lfc.cutoff
                   
To add this vector to our results table we can use the `$` notation to create the column on the left hand side of the assignment operator, and the assign the vector to it:

	res_tableOE$threshold <- threshold                

Now we can easily check how many genes are significant by using the `which()` function:

	length(which(res_tableOE$threshold))

***

**Exercise**

1. Explore the results table summary for the **Mov10_knockdown comparison to control**. How many genes are differentially expressed using the default thresholds?
2. Using the same thresholds as above (`padj.cutoff < 0.05` and `lfc.cutoff = 1`), report the number of genes that are up- and down-regulated in Mov10_knockdown compared to control.
3. Add a new column called `threshold` to the `res_tableKD` which contains a logical vector denoting genes as being differentially expressed or not.

*** 

> **NOTE: on p-values set to NA**
> > 
> 1. If within a row, all samples have zero counts, the baseMean column will be zero, and the log2 fold change estimates, p-value and adjusted p-value will all be set to NA.
> 2. If a row contains a sample with an extreme count outlier then the p-value and adjusted p-value will be set to NA. These outlier counts are detected by Cook’s distance. 
> 3. If a row is filtered by automatic independent filtering, for having a low mean normalized count, then only the adjusted p-value will be set to NA. 
>


## Visualizing the results

One way to visualize results would be to simply plot the expression data for a handful of our top genes. We could do that by picking out specific genes of interest, for example Mov10:

	# Plot expression for single gene
	plotCounts(dds, gene="MOV10", intgroup="sampletype")
	
![topgene](../img/topgen_plot.png)

### Volcano plot

This would be great to validate a few genes, but for more of a global view there are other plots we can draw. A commonly used one is a volcano plot; in which you have the log transformed adjusted p-values plotted on the y-axis and log2 fold change values on the x-axis. There is no built-in function for the volcano plot in DESeq2, but we can easily draw it using `ggplot2`. First, we will need to create a `data.frame` object from the results, which is currently stored in a `DESeqResults`  object:

	# Create dataframe for plotting
	df <- data.frame(res_tableOE)

Now we can start plotting. The `geom_point` object is most applicable, as this is essentially a scatter plot:

```
	# Volcano plot
	ggplot(df) +
  		geom_point(aes(x=log2FoldChange, y=-log10(padj), colour=threshold)) +
  		xlim(c(-2,2)) +
  		ggtitle('Mov10 overexpression') +
  		xlab("log2 fold change") + 
 		ylab("-log10 adjusted p-value") +
  		theme(legend.position = "none",
        	plot.title = element_text(size = rel(1.5)),
        	axis.title = element_text(size = rel(1.5)),
        	axis.text = element_text(size = rel(1.25)))  
```

![volcano](../img/volcanoplot-1.png)


### Heatmap

Alternatively, we could extract only the genes that are identifed as significant and the plot the expression of those genes using a heatmap.


First, let's sort the results file by adjusted p-value:
	
	### Sort the results tables
	res_tableOE_sorted <- res_tableOE[order(res_tableOE$padj), ]
	res_tableKD_sorted <- res_tableKD[order(res_tableKD$padj), ]
	
Now let's get the gene names for those significant genes:

	### Get significant genes
	sigOE <- row.names(res_tableOE_sorted)[which(res_tableOE_sorted$threshold)]
	sigKD <- row.names(res_tableKD_sorted)[which(res_tableKD_sorted$threshold)]
	
We can then use those genes to select the corresponding rows from the normalized data matrix:

	### Extract normalized expression for significant genes
	norm_OEsig <- normalized_counts[sigOE,]

Now let's draw the heatmap using `pheatmap`:

	### Annotate our heatmap (optional)
	annotation <- data.frame(sampletype=meta[,'sampletype'], 
                         row.names=row.names(meta))

	### Set a color palette
	heat.colors <- brewer.pal(6, "YlOrRd")
	
	### Run pheatmap
	pheatmap(norm_OEsig, color = heat.colors, cluster_rows = T, show_rownames=F,
	annotation= annotation, border_color=NA, fontsize = 10, scale="row",
         fontsize_row = 10, height=20)
         
![sigOE_heatmap](../img/sigOE_heatmap.png)       


***

**Exercise**

Generate two figures for the KD-control comparison: a volcano plot and a heatmap. Save both images to file.

***



## Hypothesis testing: Likelihood ratio test (LRT)

An alternative to pair-wise comparisons is to **analyze all levels of a factor at once**. By default the Wald test is used to generate the results table, but DESeq2 also offers the LRT which is used to identify any genes that show change in expression across the three levels. This type of test can be especially useful in analyzing time course experiments. 

To use the LRT, we use the `DESeq()` function but this time adding two arguments: 1) to specify that we want to use the LRT `test` and 2) the `reduced` model:

	
	### Likelihood ratio test
	dds_lrt <- DESeq(dds, test="LRT", reduced = ~ 1)

Since our model only has one factor (`sampletype`), the reduced model is just the intercept. The LRT is comparing the full model to the reduced model to identify significant genes. The p-values are determined solely by the difference in deviance between the full and reduced model formula (not fold changes). Generally, this test will result in a larger number of genes than the individual pair-wise comparisons. While the LRT is a test of significance for differences of any level of the factor, one should not expect it to be exactly equal to the union of sets of genes using Wald tests (alhtough there will be substantial overlap).

Let's take a look at the results table:

	# Extract results
	res_LRT <- results(dds_lrt, test="LRT")
	
You will find that similar columns are reported for the LRT test. One thing to note is, even though there are fold changes present they are not directly associated with the actual hypothesis test. Thus, when filtering significant genes from the LRT we use only the FDR as our threshold. *How many genes are significant at `padj < 0.05`?*

	length(which(res_LRT$padj < padj.cutoff))
	
Similar to our other result tables, let's add in a column to denote which genes are significant:

	res_LRT$threshold <- res_LRT$padj < padj.cutoff


Having this colum will allow us to make some quick comparisons as to whether we see an overlap with our pair-wise Wald test results.

	# Get sig gene lists
	LRTgenes <- row.names(res_LRT)[which(res_LRT$threshold)]
	OEgenes <- row.names(res_tableOE)[which(res_tableOE$threshold)]
	KDgenes <- row.names(res_tableKD)[which(res_tableKD$threshold)]

How many genes from the Mov10 overexpression Wald test are contained in the LRT gene set? And for the Mov10 knockdown? 

The number of significant genes observed from the LRT is quite high. We are **unable to set a fold change criteria here since the statistic is not generated from any one pairwise comparison.** This list includes genes that can be changing in any number of combinations across the three factor levels. It is advisable to instead increase the stringency on our criteria and lower the FDR threshold.

***

**Exercise**

1. Using a more stringent cutoff of `padj < 0.001`, count how many genes are significant using the LRT method.
2. Set the variables `OEgenes` and `KDgenes`to contain the genes that meet the  threshold `padj < 0.001`.
3. Find the overlapping number of genes between these gene sets and the genes from LRT at `padj < 0.0001`.

***

## Exporting significant gene lists

The next step in our workflow is interpretation of gene lists using various tools for functional analysis. Depending on the tool you choose to use downstream, you will require different information from the results table as input. To be safe it is wise to keep atleast one copy of the full results table with relevant information. 
	
Let's use the `write.table()` function to write the ordered results to file:

	### Write sorted results to file
	write.table(res_tableOE_sorted, file="results/results_OE_sortedPval.txt", sep="\t", quote=F, col.names=NA)
	
	write.table(res_tableKD_sorted, file="results/results_KD_sortedPval.txt", sep="\t", quote=F, col.names=NA)

One of the tools we will be using for functional analysis (gProfiler) will require only the gene names of the significant genes, but ordered by adjusted p-value. These lists we had created above for visualization.

To write these lists to file we will use the `write()` function which will write the contents to file on single line, or if `ncol` is specified, into a certain number of columns:

	### Write genes to file
	write(sigOE, file="results/Mov10_oe_logFC_1_pVal_0.05.txt", ncol=1)
	write(sigKD, file="results/Mov10_kd_logFC_1_pVal_0.05.txt", ncol=1)
	

	
## Saving the project

Now we are set up for functional analysis of our gene lists. Make sure you save your R session as you quit RStudio to your DEanalysis project, so you don't lose all your work from this DE analysis module!

---
*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*

* *Materials and hands-on activities were adapted from [RNA-seq workflow](http://www.bioconductor.org/help/workflows/rnaseqGene/#de) on the Bioconductor website*
