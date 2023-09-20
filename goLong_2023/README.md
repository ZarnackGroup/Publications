# SF3B1 mutations primarily affect proximal alternative splice site in CLL and MDS patients.

## Abstract

Alternative splicing plays a critical role in generating transcriptome diversity, and aberrant splicing is frequently observed in cancer. Mutations in the splicing factor SF3B1 are particularly common in patients with chronic lymphocytic leukaemia (CLL) and myelodysplastic syndromes (MDS), with an opposite predictive model. In order to get insight into the effect of SF3B1 mutations on the transcriptome we used long-read sequencing data derived from MDS and CLL patient cohorts, as well as cell lines. Our results revealed that SF3B1 mutations specifically alter the usage of 3’ alternative splice sites especially within short proximity to the canonical splice sites. Moreover, disease-specific differences in the SF3B1 mutations effect seem to emerge from different transcriptomic profiles rather than different mechanism of SF3B1 mutation itself. Full isoform information enabled to predict the functional consequences of the aberrant splicing and gain mechanistic insights into the role of mutated SF3B1 in splicing.

## Source code

The source code is split into three parts:

* The [analysis of the isoseq data](https://github.com/ZarnackGroup/go_long2023/tree/main/1%20Isoseq%20Analysis) ; Figures 1-4
* The [analysis of the iCLIP data](https://github.com/ZarnackGroup/go_long2023/tree/main/2%20Binding%20Site%20Analysis); Figure 5
* The [integration of isoseq and iCLIP data](https://github.com/ZarnackGroup/go_long2023/tree/main/3%20Integration%20Analysis); Figure 6

## Tools and packages

* The isoseq analysis for long read transcriptome sequencing was performed with the *IsoTools python module* ([GitHub](https://github.com/MatthiasLienhard/isotools))
* The iCLIP analysis was performed with the *BindingSiteFinder* R package ([GitHub](https://github.com/ZarnackGroup/BindingSiteFinder), [Bioconductor](https://www.bioconductor.org/packages/release/bioc/html/BindingSiteFinder.html)).

## Reference
TBA
