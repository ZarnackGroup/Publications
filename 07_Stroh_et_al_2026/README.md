# Source code for "m6A sites in the coding region trigger translation-dependent mRNA decay"

## Abstract

During S phase, the genome is replicated in a tightly regulated spatiotemporal order described as DNA replication timing (RT). Discontinuous lagging-strand synthesis produces Okazaki fragments whose strand-specific distribution reflects replication dynamics. Here, we present RepliCNN, a deep learning framework based on one- dimensional convolutional neural networks to predict RT from Okazaki fragment distributions obtained from strand-specific 3′ DNA end sequencing methods such as GLOE-Seq, TrAEL-seq, or OK-Seq. RepliCNN also automatically annotates replication origins, termination zones, replication fork directionality, and origin efficiency genome-wide from a single dataset. Benchmarking on public and in-house human and yeast datasets using leave-one-chromosome-out cross-validation demonstrates high predictive accuracy in both wild-type and perturbation experiments, enabling comprehensive analyses of replication dynamics from strand-specific DNA 3′ end sequencing data.

## Source code

The source code is deposited at https://github.com/ZarnackGroup/RepliCNN_paper.

## Reference

[1] Stroh D, Zilio N, Pabba MK, Roukos V, Cardoso MC, Ulrich HD, Zarnack K. RepliCNN: High-resolution inference of the DNA replication program from strand-specific 3′ DNA end sequencing. bioRxiv. 2026. https://doi.org/10.64898/2026.03.12.710907
