# Source code for "Long-read transcriptome sequencing of CLL and MDS patients uncovers molecular effects of SF3B1 mutations"

## Abstract

Mutations in splicing factor 3B subunit 1 (SF3B1) frequently occur in patients with chronic lymphocytic leukemia (CLL) and myelodysplastic syndromes (MDSs). These mutations have different effects on the disease prognosis with beneficial effect in MDS and worse prognosis in CLL patients. A full-length transcriptome approach can expand our knowledge on SF3B1 mutation effects on RNA splicing and its contribution to patient survival and treatment options. We applied long-read transcriptome sequencing (LRTS) to 44 MDS and CLL patients, as well as two pairs of isogenic cell lines with and without SF3B1 mutations, and found >60% of novel isoforms. Splicing alterations were largely shared between cancer types and specifically affected the usage of introns and 3' splice sites. Our data highlighted a constrained window at canonical 3' splice sites in which dynamic splice-site switches occurred in SF3B1-mutated patients. Using transcriptome-wide RNA-binding maps and molecular dynamics simulations, we showed multimodal SF3B1 binding at 3' splice sites and predicted reduced RNA binding at the second binding pocket of SF3B1K700E Our work presents the hitherto most-complete LRTS study of the SF3B1 mutation in CLL and MDS and provides a resource to study aberrant splicing in cancer. Moreover, we showed that different disease prognosises result most likely from the different cell types expanded during carcinogenesis rather than different mechanisms of action of the mutated SF3B1. These results have important implications for understanding the role of SF3B1 mutations in hematological malignancies and other related diseases. 

## Source code

The source code can be found [here](https://github.com/ZarnackGroup/go_long2023/tree/main) and is split into three parts:

* The [analysis of the IsoSeq data](https://github.com/ZarnackGroup/go_long2023/tree/main/1%20Isoseq%20Analysis); Figures 1-4
* The [analysis of the iCLIP data](https://github.com/ZarnackGroup/go_long2023/tree/main/2%20Binding%20Site%20Analysis); Figure 5
* The [integration of IsoSeq and iCLIP data](https://github.com/ZarnackGroup/go_long2023/tree/main/3%20Integration%20Analysis); Figure 6

## Tools and packages

* The IsoSeq analysis for long read transcriptome sequencing was performed with the *IsoTools python module* ([GitHub](https://github.com/MatthiasLienhard/isotools))
* The iCLIP analysis was performed with the *BindingSiteFinder* R package ([GitHub](https://github.com/ZarnackGroup/BindingSiteFinder), [Bioconductor](https://www.bioconductor.org/packages/release/bioc/html/BindingSiteFinder.html)).

## Reference
[1] Pacholewska A*, Lienhard M*, Brüggemann M*, Hänel H, Bilalli L, Königs A, Heß F, Becker K, Köhrer K, Kaiser J, Gohlke H, Gattermann N, Hallek M, Herling CD, König J, Grimm C, Herwig R#, Zarnack K#, Schweiger MR#. Long-read transcriptome sequencing of CLL and MDS patients uncovers molecular effects of SF3B1 mutations. Genome Res. 2024 Nov 20;34(11):1832-1848. https://doi.org/10.1101/gr.279327.124
