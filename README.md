# Genetic approaches for studying microsporidia

Data and R code used to generate the figures in:

> Reinke AW. Genetic approaches for studying microsporidia. \*Philosophical Transactions of the Royal Society B\* (in revision).

This review surveys the three main genetic approaches that have been applied to microsporidia: RNA interference (RNAi), genome transformation/modification, and experimental evolution. For each approach, published studies were compiled into curated tables (Tables S1–S3), and these tables were analysed and plotted with the scripts in this repository.

## Repository structure

```
microsporidia-genetic-approaches/
├── README.md
├── data/
│   ├── Table\_S1\_RNAi.xlsx
│   ├── Table\_S2\_transformation.xlsx
│   ├── fig3B\_transformation\_summary.csv
│   └── Table\_S3\_experimental\_evolution.xlsx
├── scripts/
│   ├── fig2\_rnai\_analysis.R
│   ├── fig3\_genome\_transformation.R
│   ├── fig4\_experimental\_evolution.R
│   └── figS1\_extraction\_validation.R
└── figures/
    ├── fig2\_rnai\_analysis.pdf / .png
    ├── fig3\_genome\_transformation.pdf / .png
    ├── fig4\_experimental\_evolution.pdf / .png
    └── figS1\_extraction\_validation.pdf / .png
```

## Figures, scripts and data

|Figure|Script|Input data|Description|
|-|-|-|-|
|Figure 1|—|—|Schematic overview of genetic approaches (drawn, not generated with code)|
|Figure 2|`fig2\_rnai\_analysis.R`|Table S1|RNAi studies: publications over time, targets by protein category, knockdown phenotype magnitude, mRNA knockdown vs phenotype, off-target expression, and phenotypes examined|
|Figure 3|`fig3\_genome\_transformation.R`|Table S2; `fig3B\_transformation\_summary.csv`|Genome-transformation studies over time and summary table of transformation studies|
|Figure 4|`fig4\_experimental\_evolution.R`|Table S3|Experimental evolution studies over time and selection conditions and outcomes|
|Figure S1|`figS1\_extraction\_validation.R`|Table S1|Comparison of values extracted from figures by a generative AI model with raw data provided by the original authors|

## Data

* **Table S1 (`Table\_S1\_RNAi.xlsx`)**: Papers applying RNAi against microsporidia genes.

  * `Master\_Bibliography`: one row per paper, with parasite species, host, targets and assays.
  * `Categories of RNAi targets`: the likely protein category assigned to each RNAi target.
  * `Quant\_Timepoint\_Data`: quantitative values extracted from published figures for each paper, gene, timepoint and measurement type (target mRNA, off-target mRNA and growth phenotypes), expressed relative to the negative control.
  * `Author data comparison`: values extracted from figures compared with raw data published by the authors.
* **Table S2 (`Table\_S2\_transformation.xlsx`)**: Papers applying transformation and genomic modification to microsporidia (sheet `Master\_Bibliography`).
* **`fig3B\_transformation\_summary.csv`**: The DNA delivery approach, fluorescent marker, selection and reported result for each of the three transformation studies in Table S2, shown as the table in Figure 3B.
* **Table S3 (`Table\_S3\_experimental\_evolution.xlsx`)**: Papers applying experimental evolution to microsporidia (sheet `Experimental\_Evolution\_Biblio`), including the selection condition, the assay performed after passage and whether a change was observed.

Quantitative values in Table S1 were extracted from published figures using a generative AI model, and at least one value from every study (and several values from studies contributing multiple data points) was checked manually against the source figure. Figure S1 shows how these extracted values compare with raw data published by the original authors.

## Running the code

All scripts are run from the repository root, so that the relative paths to `data/` and `figures/` resolve correctly. For example:

```bash
Rscript scripts/fig2\_rnai\_analysis.R
```

Or, in R/RStudio, set the working directory to the repository root and run:

```r
source("scripts/fig2\_rnai\_analysis.R")
```

Each script reads its table from `data/` and writes a vector PDF and a 600-dpi PNG to `figures/`. The regression statistics for Figure S1 are also printed to the console.

## Requirements

Figures were generated in R version 4.6.1 using ggplot2 and patchwork. The following packages are required:

```r
install.packages(c("readxl", "dplyr", "tidyr", "ggplot2", "scales", "stringr",
                   "forcats", "broom", "patchwork", "gridExtra", "gtable",
                   "cellranger", "ggtext", "janitor"))
```

## Statistics

Values from RNAi experiments were expressed as log2-transformed ratios of treatment to control. For each paper, gene and phenotype, only the timepoint with the strongest reduction in phenotype was used. Relationships between target mRNA knockdown and either phenotype or off-target mRNA levels were assessed for each group by ordinary least-squares linear regression (`lm`). For each regression, the coefficient of determination (R²) is reported, and the slope was tested against zero by a two-sided t-test at p < 0.05, without correction for multiple comparisons.

## Use of AI

Generative AI models (Claude, ChatGPT and Gemini) were used to help write the code in this repository and to extract values from published figures. All code, extracted values and outputs were checked by the author.

## Citation

If you use these data or code, please cite the paper above.

## Contact

Aaron Reinke — aaron.reinke@utoronto.ca
Department of Molecular Genetics, University of Toronto

