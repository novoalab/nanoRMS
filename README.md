# NanoRMS: predicting NANOpore Rna Modification Stoichiometry
Prediction and visualization of RNA modification stoichiometry in direct RNA sequencing datasets from per-read information 

![alt text](./img/init_fig.png "init_fig")


## Table of Contents  
- [General Description](#General-description)
- [Prediction of RNA modified sites](#Prediction-of-RNA-modified-sites)
- [RNA modification stoichiometry estimation using Nanopolish resquiggling (not recommended)](#RNA-modification-stoichiometry-estimation-using-nanopolish-resquiggling)
- [RNA modification stoichiometry estimation using Tombo resquiggling (recommended)](#RNA-modification-stoichiometry-estimation-using-tombo-resquiggling)
- [Visualization of per-read current intensities at individual sites](Visualization-of-per-read-current-intensities-at-individual-sites)
- [Dependencies and versions](#Dependencies-and-versions)
- [Citation](#Citation) 
- [Contact](#Contact) 
 

## General description
* NanoRMS predicts modification stoichiometries by classifying reads into modified/unmodified, based on their per-read features (current intensity, dwell time and/or trace).
* NanoRMS can be run in single mode (1sample) or paired mode (2 samples).
* NanoRMS can run both unsupervised (e.g. KMEANS, Aggregative Clustering, GMM) and supervised machine learning algorithms (e.g. KNN, Random Forest). The later will require pairwise samples where one of the conditions is a knockout.
* NanoRMS can predict stoichiometry from Nanopolish resquiggled reads or from Tombo resquiggled reads. The latter is the recommended option.

## Prediction of RNA modified sites

### 1. Extract base-calling features using Epinano-RMS 

#### Create a dictionary file from reference fasta
Please note that we are using  a modified version of [EpiNano](https://github.com/enovoa/EpiNano/blob/master/README.md) that was specifically developed for nanoRMS. 

Firstly, we need to create a dictionary file from the reference fasta file
```
java -jar epinano_RMS/picard.jar CreateSequenceDictionary REFERENCE=reference.fasta OUTPUT= reference.fasta.dict
```

#### Run Epinano-RMS
Requires python3

```
python3 epinano_RMS/epinano_rms.py -R <reference_file> -b <bam_file> -s epinano_RMS/sam2tsv
```
Example using test data: 

```
python3 epinano_RMS/epinano_rms.py -R test_data/yeast_rRNA_ref -b test_data/wt_sorted.bam -s epinano_RMS/sam2tsv
```

### 2. Predict RNA modifications

#### a) Single sample RNA modification prediction (i.e. "de novo" prediction)

Single sample '*de novo*' RNA modification prediction has been tested for predicting pseudouridine RNA modifications in mitochondrial rRNAs, and the novel predicted sites were afterwards validated using CMC-based probing followed by sequencing), validating 2 out of the 2 sites that were predicted in all 3 biological replicates. 

This code relies on the identification of pseudouridine base-calling error 'signatures', which allows us to predict RNA modifications de novo in individual samples, as long as the stoichiometry of modification is sufficiently high (i.e. to be distinguished from background base-calling error of direct RNA sequencing).


General usage: 
```
Rscript predict_singleSample.R <epinanofile_rep1> <epinanofile_rep2> <epinanofile_rep3> 
```

Example using test data (prediction of pseudouridine sites on mitochondrial ribosomal RNAs): 
```
Rscript predict_singleSample.R wt_epinano.csv sn3_epinano.csv sn36_epinano.csv
```

#### b) Paired sample RNA modification prediction (i.e. "differential-error"-based prediction)

Pseudouridine is not always present in high stoichiometries (e.g. rRNAs), but can also be present in low stoichiometries (e.g. in mRNAs). Please note that pseudouridines in mRNAs cannot be accurately predicted using "de novo" mode, because the background nanopore 'error' is too similar to the 'error' caused by the presence of pseudouridine. 
For such cases, we can predict differentially pseudouridylated sites by identifying which sites show pseudouridine differential error signatures between two conditions, as shown below. This type of pairwsie comparison can be done for WT-KO, or between two conditions (e.g. normal-heat stress). 


General usage: 
```
Rscript predict_twoSample_transcript.R <epinanofile_rep1_normal> <epinanofile_rep1_heatshock> <epinanofile_rep2_normal> <epinanofile_rep2_heatshock>
```

Example using test data: 
```
Rscript predict_twoSample_transcript.R ncRNA_normal_rep1_epinano.csv ncRNA_heatshock_rep1_epinano.csv ncRNA_normal_rep2_epinano.csv ncRNA_heatshock_rep2_epinano.csv
```

## RNA modification stoichiometry estimation using Nanopolish resquiggling 

This version is deprecated. If you still wish to use it, you can find the details and code [here](https://github.com/novoalab/nanoRMS/blob/master/README_nanoRMS_nanopolish.md) 


## RNA modification stoichiometry estimation using Tombo resquiggling 

To use this version, you can find the details [here](https://github.com/novoalab/nanoRMS/blob/master/per_read/README.md)

0. First download test data
This dataset is described in more details in per_read directory. 
```bash
(cd per_read && wget https://public-docs.crg.es/enovoa/public/lpryszcz/src/nanoRMS/per_read/guppy3.0.3.hac -q --show-progress -r -c -nc -np -nH --cut-dirs=6 --reject="index.html*")
```

1. Retrieve per-read features from all samples. 
```bash
ref=per_read/guppy3.0.3.hac/Saccharomyces_cerevisiae.R64-1-1_firstcolumn.ncrna.fa
per_read/get_features.py --rna -f $ref -t 6 -i per_read/guppy3.0.3.hac/*WT??C/workspace/*.fast5
```

2. Estimate modification frequency difference between two samples  
Note, you'll need to provide candidate positions that are likely modified. Those were identified earlier (see ,
so here we'll just generate BED file from existing candidate file. 
```bash
# prepare BED
f=per_read/results/predictions_ncRNA_WT30C_WT45C.tsv.gz
zgrep -v X.Ref $f |awk -F'\t' 'BEGIN {OFS = FS} {print $1,$2-1,$2,".",100,"+"}' > $f.bed

# calculate modification frequency difference between control and sample of interest
per_read/get_freq.py -f $ref -b $f.bed -o $f.bed.tsv.gz -1 per_read/guppy3.0.3.hac/*WT30C/workspace/*.fast5.bam -2 per_read/guppy3.0.3.hac/*WT45C/workspace/*.fast5.bam
```

![alt text](img/get_freq.png "Stoichometry output")

Note, the candidate position file (`-b $f.bed`) has to refer only to reference positions that have some reads aligned -
otherwise the `get_freq.py` will fail with an error `ValueError: Found array with 0 sample(s) (shape=(0,6)) while a minimu of 1 is required`.  
Please note that KMEANS does not accurately assign directionality of the stoichiometry change, whereas KNN does (because KMEANS randomly assigns one cluster as "modified" and another as "unmodified". Thus, to know the directionality of the change for KMEANS stoichiometry predictions, you will need to infer that from the directionality of mismatch error in that given position. If you don't care about the directionality of the change, but just about the effect size of the change, you can just take the absolute values of the predicted stoichiometry changes.

## Visualization of per-read current intensities at individual sites

### 1. Data pre-processing

Firstly, generate a collapsed Nanopolish event align output, by collapsing all the multiple observations for a given position from a same read.

```
python3 per_read_mean.py <event_align_file>
```

Example using test data:

```
python3 per_read_mean.py test_data/data1_eventalign_output.txt
```
Secondly, create 15-mer windows of per-read current intensities centered in positions of interest

The output of Nanopolish event align generated in the previous step is used as input in this script.

```
Rscript --vanilla nanopolish_window.R positions_file <input_table> <label>
```

Example using test data:

```
Rscript --vanilla nanopolish_window.R test_data/positions test_data/data1_eventalign_output.txt_processed_perpos_mean.csv data1
```


### 2. Visualization of current intensity information:

#### 2.1 Density plots

```
Rscript --vanilla density_nanopolish.R <window_file1> <window_file2> <window_file3(optional)> <window_file4(optional)>
```

Example using test data:

```
Rscript --vanilla nanopolish_density_plot.R test_data/sn34_window_file.tsv test_data/wt_window_file.tsv
```

![alt text](./img/density.png "Density")


#### 2.2. Mean current intensity plots centered in the modified sites
```
Rscript --vanilla nanopolish_meanlineplot.R <window_file1> <window_file2> <window_file3(optional)> <window_file4(optional)>
```
Example using test data:

```
Rscript --vanilla nanopolish_meanlineplot.R test_data/sn34_window_file.tsv test_data/wt_window_file.tsv
```


![alt text](./img/mean_current.png "Mean_current")


#### 2.3 Per-read current intensity plots centered in the modified sites
```
Rscript --vanilla nanopolish_perreadlineplot.R <window_file1> <window_file2> <window_file3(optional)> <window_file4(optional)>
```
Example using test data:

```
Rscript --vanilla nanopolish_perreadlineplot.R test_data/sn34_window_file.tsv test_data/wt_window_file.tsv
```


![alt text](./img/per_read_current.png "Per_read")


#### 2.4 PCA plots from the per-read 15-mer current intensity data
```
Rscript --vanilla nanopolish_pca.R <window_file1.tsv> <window_file2.tsv> <window_file3.tsv(optional)> <window_file4.tsv(optional)>
```

Example using test data:

```
Rscript --vanilla nanopolish_pca.R test_data/sn34_window_file.tsv test_data/wt_window_file.tsv
```

![alt text](./img/pca.png "PCA")

## Dependencies and versions

Software | Version 
--- | ---
matplotlib | 3.1.2
numba | 0.51.1
numpy | 1.19.4
ont-fast5-api | 3.1.6
ont-tombo | 1.5.1
openpyxl | 3.0.5
pandas | 0.24.2
scipy | 1.5.2
seaborn | 0.11.0
sklearn | 0.23.1


## Citation: 

Begik O*, Lucas MC*, Pryszcz LP, Ramirez JM, Medina R, Milenkovic I, Cruciani S, Vieira HGS, Liu H, Sas-Chen A, Mattick JS, Schwartz S and Novoa EM. Quantitative profiling of native RNA modifications and their dynamics using nanopore sequencing. bioRxiv 2020. doi: https://doi.org/10.1101/2020.07.06.189969v3

## Contact
Please open an issue in the GitHub repo if you have any questions/doubts/suggestions about how to use this software. Thanks!
