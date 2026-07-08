# Exploratory notebooks

This folder contains exploratory notebooks used to develop the workflow. 

* `01-cnmf-results-exploration.Rmd`: This notebook looks at the results from running `cNMF` on a range of k values for a subset of samples.

* `02-metaprogram-evaluation.Rmd`: This notebook looks at the results from generating metaprograms for a subset of samples. 

* `03-metaprogram-comparison.Rmd`: This notebook compares metaprograms generated across different values of k for the full cohort as of January 2026. 

* `04-sample-diversity.Rmd`: This notebook looks at approaches and metrics for measuring sample diversity within metaprograms. 

* `05-metaprogram-correlation.Rmd`: This notebook explores the distribution of correlation-related metrics in an attempt to define outliers and set thresholds for identiying good metaprograms.

* `06-gene-set-uniqueness.Rmd`: This notebook explores metrics used to measure gene set uniqueness and biological specificity for identifying good metaprograms. 

* `07-top-gene-overlap.Rmd`: This notebook looks at metrics that can be used to measure the overlap in the top genes for each metaprogram. 

* `08-cell-scoring-metrics.Rmd`: This notebook looks at metrics used to measure variations in cell scores across metaprograms. 

* `09-metaprogram-cnv-correlation.Rmd`: This notebook looks at the correlation of CNVs with metaprogram scores for the full Sean Karl cohort. 

* `10-shared-metaprogram-metrics.html`: This notebook explores metrics used to compare metaprograms generated across different sample cohorts in an effort to identify shared metaprograms. 

* `11-all-shared-metrics.html`: This notebook evaluates the output from comparing all subsets to the full Sean Karl cohort to generate saturation plots. 
