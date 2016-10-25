---
title: "Data manipulation"
authors: Meeta Mistry and Mary Piper
date: "Tuesday, June 28, 2016"
---
Approximate time: 30 min

## Learning Objectives
* Understand and implement functions nested within other functions


## Nested functions

Thus far, to perform any specific task, we have executed every function separately; if we wanted to use the results of a function for downstream purposes, we saved the results to a variable. As you become more comfortable with R, you will find that it is more efficient to code using nested functions, or functions within other functions, which will allow you to execute multiple commands at the same time.

Even if you decide to avoid writing nested functions for the time being, you should still have experience reading and understanding them. The key to understanding nested functions is to **read from the inside out**.

### Nested functions practice #1

You realize that you forgot to include important metadata regarding sex of your samples in your `metadata` file. You would like to add this data to your `metadata` dataframe, and using functions separately would require us to execute three separate steps:  

**Step 1:** Create the `sex` vector: 
	
	sex <- c("M","F","M","M","F","M","M","F","M","M","F","M")
	
**Step 2:** Turn the `sex` vector into a factor variable:
	 
	 sex_fr <- factor(sex)
	 
**Step 3:** Use the `cbind` function to add the column to the **end** of the `metadata` dataframe: 

	metadata2 <- cbind(metadata, sex=sex_fr)

Instead of performing all three steps, we would like to create a nested function. **To create a nested function, simply replace the variable name with its contents**. We could combine steps 1 and 2 by replacing `sex` in **Step 2** with it's contents (`c("M","F","M","M","F","M","M","F","M","M","F","M")`):

	sex_fr <- factor(c("M","F","M","M","F","M","M","F","M","M","F","M"))
	metadata3 <- cbind(metadata, sex=sex_fr)
	
It is possible to combine all steps, but your code would be difficult to read, so we don't recommend doing this:

	metadata4 <- cbind(metadata,
			sex=factor(c("M","F","M","M","F","M","M","F","M","M","F","M")))

### Nested functions practice #2			
Now, let's say that you are interested in counting the number samples in your dataset that have  "Wt" genotype within our `metadata` file. Obviously, for our small file, we could just look at the file, but if we had too many samples to count, we could do the following:

**Step 1:** Determine the **location** of samples with `genotype` equal to "Wt":
	
	wt_loc <- which(metadata$genotype == "Wt")
	

**Step 2:** Determine the number of samples with `genotype` "Wt":
	
	length(wt_loc)
	
Alternatively, we could combine the steps:

	length(which(metadata$genotype == "Wt"))
	

Learning to understand nested functions is a critical part of your mastery of R. Not only will their use improve your efficiency, but nested functions are frequently encountered in help forums and R package documentation, so understanding them is critical to your learning process. 


---

*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*

* *The materials used in this lesson is adapted from work that is Copyright © Data Carpentry (http://datacarpentry.org/). 
All Data Carpentry instructional material is made available under the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0).*
