# Source code for "Dysregulated RNA Splicing Targets Cancer via Proteotoxic Stress"

## Abstract

Aberrant RNA splicing is a hallmark of aggressive cancers. We show that oncogenic transcriptional and metabolic cues, including c-MYC activation and/or lactate accumulation, reprogram splicing by upregulating the spliceosome factor USP39 to sustain malignant proteostasis. Disruption of USP39 compromises splicing fidelity, drives accumulation of aberrant proteins, and induces proteotoxic stress that selectively eliminates cancer cells and promotes tumor regression in vivo, while sparing non-transformed tissues. Leveraging this vulnerability, we develop first-generation therapies that induce proteotoxic RNA splicing by targeting USP39. These findings establish a synthetic lethal interaction between oncogenic splicing and proteostasis collapse and identify proteotoxic RNA splicing as a therapeutic strategy for aggressive solid tumors, such as hepatocellular carcinoma.

## Source code

The source code can be found hosted on a webpage [here](https://felixhaidle.github.io/prieto_garcia_et_al_2026/), or alternatively in a GitHub repository [here](https://github.com/felixhaidle/prieto_garcia_et_al_2026) and is split into three parts:

* [RNA-seq MAJIQ Data Processing](https://felixhaidle.github.io/prieto_garcia_et_al_2026/rnaseq_raw_data_processing.html)
* [Junction count based Alternative Splicing Analysis in TCGA LIHC Cohort](https://felixhaidle.github.io/prieto_garcia_et_al_2026/comparing_rnaseq_and_tcga_results.html); Fig. 2C, S2B
* [Read Coverage based Alternative Splicing Analysis in TCGA LIHC Cohort](https://felixhaidle.github.io/prieto_garcia_et_al_2026/recount_paper.html); Fig. 1I, S1M

## Tools and packages

* The TCGA cohort analysis was performed with the help of [**recount3**](https://doi.org/10.1186/s13059-021-02533-6)
