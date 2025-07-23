# Example Python Notebook for Analyzing Data from Skyline

I've I use two Skyline reports for the generation of the input data.
- I use the standard Replicates report in the document grid to export meta data annotations from within the Skyline document
- I use the included `MJM Protein Total Areas.skyr` to generate the Protein Quant matrix from the document grid

I've divided my analysis into three notebooks to make it easier to find things.
- 1-import-skyline-output.ipynb
  - imports the data from the two matrixes, parses the output, and cleans the replicate labels, puts samples into groups, etc...
  - I do a check for missing data. With the Skyline imputation there should be none.
- 2-normalize-data.ipynb
  - plots the distribution of intensities of the raw data
  - performs normalization using both median and VSN normalization
  - uses PCA to change for batch effects and grouping
  - looks at the correlation of the controls
  - performs clustering on the different controls to confirm that the different controls are different as expected.
- 3-differential-analysis.ipynb
  - T-test analysis with multiple testing correction
    - I make a number of plots and output a table.
    - I also perform a Mann-Whitney analysis
    - I perform a binary classification using an SVM with cross validation. I also perform some feature selection using Mann Whitney and also by how well the features are replicated in the cross validation steps.
    - I also have implemented a general linear model for unpaired data and a mixed effects model for paired data. I also have implemented the ability to use additional metadata fields as covariates

These analyses are nothing fancy and this is just for use as an example of how I begin looking at our differential analysis data.
