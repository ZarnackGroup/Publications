# Souce codes of "Deep and accurate detection of m6A RNA modifications using miCLIP2 and m6Aboost machine learning"

## Abstract of the paper

N6-methyladenosine (m6A) is the most abundant internal RNA modification in eukaryotic mRNAs and influences many aspects of RNA processing, such as RNA stability and translation. miCLIP (m6A individual-nucleotide resolution UV crosslinking and immunoprecipitation) is an antibody-based approach to map m6A sites in the transcriptome with single-nucleotide resolution. However, due to broad antibody reactivity, reliable identification of m6A sites from miCLIP data remains challenging. Here, we present several experimental and computational innovations, that significantly improve transcriptome-wide detection of m6A sites. Based on the recently developed iCLIP2 protocol, the optimised miCLIP2 results in high-complexity libraries from less input material, which yields a more comprehensive representation of m6A sites. Next, we established a robust computational pipeline to identify true m6A sites from our miCLIP2 data. The analyses are calibrated with data from Mettl3 knockout cells to learn the characteristics of m6A deposition, including a significant number of m6A sites outside of DRACH motifs. In order to make these results universally applicable, we trained a machine learning model, m6Aboost, based on the experimental and RNA sequence features. Importantly, m6Aboost allows prediction of genuine m6A sites in miCLIP data without filtering for DRACH motifs or the need for Mettl3 depletion. Using m6Aboost, we identify thousands of high-confidence m6A sites in different murine and human cell lines, which provide a rich resource for future analysis. Collectively, our combined experimental and computational methodology greatly improves m6A identification.

## The links for the source code and packages

The source code of "Deep and accurate detection of m6A RNA modifications using miCLIP2 and m6Aboost machine learning" are splited in three parts:    

* The source code of bin-based differential methylation analysis could be found 
in [https://github.com/Codezy99/miCLIP2-DMA-source-code](https://github.com/Codezy99/miCLIP2-DMA-source-code)
* For running the m6Aboost to predict the 

