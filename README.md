# HhsFunctionalRankings

This workflow is a companion to the to the paper "Integrating Network-based Functional 
Enrichment with Genetic Association to Prioritize Causal Genes in Quantitative Trait 
Loci." It contains all code required to reproduce the results and figures in the manuscript. Data are available at Figshare [https://figshare.com/articles/Data_required_to_run_HhsFunctionalRankings_workflow/8205356]

The workflow is in an R markdown file (Hhs.Rmd)

## Directory Organization
To replicate the analysis in the manuscript, run the code in Hhs.Rmd.
This code expects the following directories within the project directory.

* *code*: A directory containing all code subdirectories. These are available on this GitHub page.
* *data*: A directory containing all data files used in the analysis. This is available at Figshare [https://figshare.com/articles/Data_required_to_run_HhsFunctionalRankings_workflow/8205356]
* *results*: A directory containing all files produced by the analysis. The output is also available at Figshare [https://figshare.com/articles/Data_required_to_run_HhsFunctionalRankings_workflow/8205356]. However, downloading these results is not necessary as running this workflow will reproduce those results. 

## Package Requirements
This workflow was built in R/3.5.1

An interntet connection is required to run this workflow, as it queries Biomart to retrieve information about genes.

The following is a list of packages required by this workflow. Please install these prior to running the workflow

* e1071
* RColorBrewer
* doParallel
* foreach
* gProfileR
* pheatmap
* igraph
* biomaRt
* knitr
* DescTools
* here
* emma
* rPref
* evd

