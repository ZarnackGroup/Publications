# Souce code for "Deep and accurate detection of m6A RNA modifications using miCLIP2 and m6Aboost machine learning"

## Abstract

N6-methyladenosine (m6A) is the most abundant internal RNA modification in eukaryotic mRNAs and influences many aspects of RNA processing, such as RNA stability and translation. miCLIP (m6A individual-nucleotide resolution UV crosslinking and immunoprecipitation) is an antibody-based approach to map m6A sites in the transcriptome with single-nucleotide resolution. However, due to broad antibody reactivity, reliable identification of m6A sites from miCLIP data remains challenging. Here, we present several experimental and computational innovations, that significantly improve transcriptome-wide detection of m6A sites. Based on the recently developed iCLIP2 protocol, the optimised miCLIP2 results in high-complexity libraries from less input material, which yields a more comprehensive representation of m6A sites. Next, we established a robust computational pipeline to identify true m6A sites from our miCLIP2 data. The analyses are calibrated with data from Mettl3 knockout cells to learn the characteristics of m6A deposition, including a significant number of m6A sites outside of DRACH motifs. In order to make these results universally applicable, we trained a machine learning model, m6Aboost, based on the experimental and RNA sequence features. Importantly, m6Aboost allows prediction of genuine m6A sites in miCLIP data without filtering for DRACH motifs or the need for Mettl3 depletion. Using m6Aboost, we identify thousands of high-confidence m6A sites in different murine and human cell lines, which provide a rich resource for future analysis. Collectively, our combined experimental and computational methodology greatly improves m6A identification [1].

## Links for source code and packages

The source code of "Deep and accurate detection of m6A RNA modifications using miCLIP2 and m6Aboost machine learning" is split in three parts:    

* The source code for the *bin-based differential methylation analysis* can be found 
at https://github.com/Codezy99/miCLIP2-DMA-source-code
* For running the m6Aboost model to predict m6A sites from miCLIP2 data, 
we made a package `m6Aboost` which can be found 
at https://bioconductor.org/packages/devel/bioc/html/m6Aboost.html
* The transcript metaprofiles were generated with the package `cliProfiler` which 
can be found at https://github.com/Codezy99/cliProfiler

## Reference
[1] Körtel N*, Rücklé C*, Zhou Y*, Busch A, Hoch-Kraft P, Sutandy FXR, Haase J, Pradhan M, Musheev M, Ostareck D, Ostareck-Lederer A, Dieterich C, Hüttelmaier S, Niehrs C, Rausch O, Dominissini D, König J#, Zarnack K#. Deep and accurate detection of m6A RNA modifications using miCLIP2 and m6Aboost machine learning. Nucleic Acids Res. 2021 Sep 20;49(16):e92. https://doi.org/10.1093/nar/gkab485.
