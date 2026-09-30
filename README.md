# Replication data and code for "Dual-Capacity Signaling: How Police Show Regime Strength on Douyin"

Authors: Yiqiang Wang, Yingqi Huang, Haohan Chen

From this folder, run `rmarkdown::render("Replication_Code_FigTab.Rmd")`.
The R Markdown file recalculates tables and figures included in the manuscript and associated appendix.

`data` contains four processed inputs:

- `analysis.rds`: the signaling outcomes, engagement variables, event
  indicators, and city-level covariates used in the models. 
- `table_01_labels.csv`: six binary coder/model labels for Table 1.
- `description_validation_ai.rds`: the eight evaluator ratings needed for
  Appendix Figure 2.
- `description_validation_ai_human.rds`: the 16 evaluator ratings and one
  eligibility flag needed for Appendix Figure 3.

- Machine information
  - System: Windows Server 2022 Datacenter
  - RAM: 32 GB
- Software information
  - R version 4.6.1 (2026-06-24)
  - Packages used (please see the `sessionInfo()` output below for the full list of packages and their versions)

```r
## R version 4.6.1 (2026-06-24 ucrt)
## Platform: x86_64-w64-mingw32/x64
## Running under: Windows Server 2022 x64 (build 26100)
## 
## Matrix products: default
##   LAPACK version 3.12.1
## 
## locale:
## [1] LC_COLLATE=English_United States.utf8 
## [2] LC_CTYPE=English_United States.utf8   
## [3] LC_MONETARY=English_United States.utf8
## [4] LC_NUMERIC=C                          
## [5] LC_TIME=English_United States.utf8    
## 
## time zone: Asia/Shanghai
## tzcode source: internal
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
## [1] knitr_1.51
## 
## loaded via a namespace (and not attached):
##  [1] digest_0.6.39     R6_2.6.1          fastmap_1.2.0     xfun_0.59        
##  [5] cachem_1.1.0      htmltools_0.5.9   rmarkdown_2.31    lifecycle_1.0.5  
##  [9] cli_3.6.6         sass_0.4.10       jquerylib_0.1.4   compiler_4.6.1   
## [13] rstudioapi_0.19.0 tools_4.6.1       evaluate_1.0.5    bslib_0.11.0     
## [17] yaml_2.3.12       otel_0.2.0        jsonlite_2.0.0    rlang_1.3.0
```
