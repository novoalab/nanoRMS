# NanoRMS: predicting NANOpore Rna Modification Stoichiometry
Prediction of RNA modification stoichiometry from per-read current intensity information in direct RNA sequencing datasets

## General description
This code uses Nanopolish eventalign output files and performs the following steps:

* Step 1. Convert Nanopolish eventalign outputs into processed output for each 15-mer region 
* Step 2. Visualization of the per-read results (PCA, per-read current intensities)
* Step 3. Stoichiometry prediction, using either KMEANS or KNN.

Note, the KNN should only be used if you have a training set with 0% modified (e.g. IVT, knockout condition) and 100% (or quasi-100%) modified sites (e.g. rRNAs, IVT products). If this is not the case, KMEANS is recommended.


## Requirements/dependencies

* Nanopolish (tested with version xxx)
* R (tested version 3.6.3)

## Running the code:


## Citation: 

Begik O*, Lucas MC*, Ramirez JM, Milenkovic I, Cruciani S, Vieira HGS, Medina R, Liu H, Sas-Chen A, Mattick JS, Schwartz S and Novoa EM. Decoding ribosomal RNA modification dynamics at single molecule resolution. bioRxiv 2020. doi: https://doi.org/10.1101/2020.07.06.189969

## Questions/doubts?
Please open an issue in the GitHub repo if you have any questions/doubts/suggestions about how to use this software. Thanks!
