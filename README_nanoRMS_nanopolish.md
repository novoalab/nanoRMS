
## Predicting RNA modification stoichiometry using Nanopolish  (not recommended)

To predict RNA modification stoichiometry using Nanopolish, nanoRMS will need: i) Nanopolish eventalign output files, and ii) a list of predicted candidate RNA modification sites. It then performs the following steps:

* 1. Collapse Nanopolish eventalign output
* 2. Convert Nanopolish eventalign outputs into processed output for each 15-mer region 
* 3. Visualization of the per-read results (PCA, per-read current intensities) -- optional step, but highly recommended to see how your samples look like in terms of modified/unmodified reads
* 4. Stoichiometry prediction, using either KMEANS or KNN.


### STEP 1. Run Nanopolish on your FAST5: getting per-read current intensities
Before you start, you first need to run **[Nanopolish](https://github.com/jts/nanopolish) index** and **Nanopolish eventalign** on the raw FAST5 reads using the following command line: 

```bash
nanopolish index -d fast5_file fastq_file -s sequencing_summary_file
#Index file should be in the same path as the fastq file

nanopolish eventalign \
    --reads path/to/fastq_file \
    --bam  path/to/bam_file \
    --genome path/to/reference_file.fasta \
    --scale-events > output.txt
```

### STEP 2. Run EpiNano on your BAM: getting differential base-calling 'errors'

#### 2.1 Obtain *EpiNano* base-calling error information from mapped BAM files:
```
## Convert BAM to TSV
samtools view -h bamfile.bam | java -jar sam2tsv.jar -r reference_file  > bamfile.bam.tsv 

## Run EpiNano
python3.7 TSV_to_Variants_Freq.py3 -f bamfile.bam.tsv -t 10 

```
Example using test data:
```
samtools view -h test_data/test.bam | java -jar sam2tsv.jar -r test_data/yeast_rRNA_ref  > test.bam.tsv 

python3.7 TSV_to_Variants_Freq.py3 -f test.bam.tsv -t 10 
```

#### 2.2. Then, convert *EpiNano* outputs into Summed_Errors, using the code below: 
```
Rscript summed_errors.R <epinano_file> <summed_errors_output_file>
```

Example using test data:

```
Rscript summed_errors.R test_data/wt_epinano.csv wt
```
#### 2.3. You can visualize your EpiNano results using the following code (optional):

* a) Scatterplot representations
```
Rscript epinano_scatterplot.R <input1> <label1> <input2> <label2> <feature>

## Feature can be: "mis", "ins", "del", "q_mean", "q_median", "q_std" and "sum" if you are using the input with sum errors

```

Example using test data:

```
Rscript epinano_scatterplot.R test_data/wt_epinano.csv wt test_data/sn34ko_epinano.csv sn34ko mis
```
![alt text](./img/mis_scatter.png "Mismatch Scatter Plot")


* b) Per-transcript representations

```
Rscript epinano_barplot.R input1 label1 input2 label2 feature
```
Example with test data:
```
Rscript epinano_barplot.R test_data/wt_epinano.csv wt test_data/sn34ko_epinano.csv sn34ko mis
```


![alt text](./img/delta_mis_barplot.png "Delta Mismatch Barplot Plot")


Result: There are two regions that show distinct mismatch profiles when comparing WT and snR34-KO, one centered in 25s_rRNA:2880 and one centered in 25s_rRNA:2826. These are actually the two exact locations that are expected to be affected by snR34 depletion. You can check predicted target sites of snR34 as well as of other snoRNAs in yeast [here](https://people.biochem.umass.edu/sfournier/fournierlab/snornadb/snrs/snr34_ta.php).


### STEP 3. Run nanoRMS:

#### 3.1. Pre-processing the Nanopolish event align output 
Generate a collapsed Nanopolish event align output, by collapsing all the multiple observations for a given position from a same read.

```
python3 per_read_mean.py <event_align_file>
```

Example using test data:

```
python3 per_read_mean.py test_data/data1_eventalign_output.txt
```


#### 3.2. Create 15-mer windows of per-read current intensities centered in positions of interest

The output of Nanopolish event align generated in the previous step is used as input in this script.

```
Rscript --vanilla nanopolish_window.R positions_file <input_table> <label>
```


Example using test data:

```
Rscript --vanilla nanopolish_window.R test_data/positions test_data/data1_eventalign_output.txt_processed_perpos_mean.csv data1
```


#### 3.3. Visualize current intensity information of modified sites (optional)

#### Distribution of current intensities at the modified site (position 0)

```
Rscript --vanilla density_nanopolish.R <window_file1> <window_file2> <window_file3(optional)> <window_file4(optional)>
```

Example using test data:

```
Rscript --vanilla nanopolish_density_plot.R test_data/sn34_window_file.tsv test_data/wt_window_file.tsv
```

![alt text](./img/density.png "Density")


#### Mean current intensity plots centered in the modified sites
```
Rscript --vanilla nanopolish_meanlineplot.R <window_file1> <window_file2> <window_file3(optional)> <window_file4(optional)>
```
Example using test data:

```
Rscript --vanilla nanopolish_meanlineplot.R test_data/sn34_window_file.tsv test_data/wt_window_file.tsv
```


![alt text](./img/mean_current.png "Mean_current")


#### Per-read current intensity plots centered in the modified sites
```
Rscript --vanilla nanopolish_perreadlineplot.R <window_file1> <window_file2> <window_file3(optional)> <window_file4(optional)>
```
Example using test data:

```
Rscript --vanilla nanopolish_perreadlineplot.R test_data/sn34_window_file.tsv test_data/wt_window_file.tsv
```


![alt text](./img/per_read_current.png "Per_read")


#### PCA plots from the per-read 15-mer current intensity data
```
Rscript --vanilla nanopolish_pca.R <window_file1.tsv> <window_file2.tsv> <window_file3.tsv(optional)> <window_file4.tsv(optional)>
```

Example using test data:

```
Rscript --vanilla nanopolish_pca.R test_data/sn34_window_file.tsv test_data/wt_window_file.tsv
```


#### Export each individual position for clustering
```
Rscript --vanilla nanopolish_export_each_position.R <window_file1.tsv> <window_file2.tsv> <window_file3.tsv(optional)> <window_file4.tsv(optional)>
```

Example using test data:

```
Rscript --vanilla nanopolish_export_each_position.R test_data/bc2_window_file.tsv
```



![alt text](./img/pca.png "PCA")


### STEP 4. Estimation of RNA modification stoichiometry 

#### a) Using KMEANS clustering

```
R --vanilla < read_clustering.R --args <file1.tsv> <file2.tsv> kmeans
```

Example using test data:
```
R --vanilla < read_clustering.R --args test_data/25s_2880.wt.15mer.perread.h.tsv test_data/25s_2880.sn34.15mer.perread.h.tsv kmeans
```


#### b) Using KMEANS on PCAed data
```
R --vanilla < read_clustering.R --args <file1.ts>v <file2.tsv> kmeans_pca
```
Example using test data:
```
R --vanilla < read_clustering.R --args test_data/25s_2880.wt.15mer.perread.h.tsv test_data/25s_2880.sn34.15mer.perread.h.tsv kmeans_pca
```

![alt text](./img/KMEANS_PCA_plots.png "KMEANS_PCA_plots")


#### c) Using K-nearest neighbour (KNN) 
Important note: this option should be only used when 0% and 100% (or similar) modified datasets are available for training, e.g. rRNAs coupled with KO condition. 

```
R --vanilla < read_clustering.R --args <file1.tsv> <file2_UNM.tsv> <validation_set.tsv(optional)> knn
```

Example using test data:
```
R --vanilla < read_clustering.R --args test_data/25s_2880.wt.15mer.perread.h.tsv test_data/25s_2880.sn34.15mer.perread.h.tsv knn
```

Example using test data that includes predictions in independent validation data:

```
R --vanilla < read_clustering.R --args test_data/25s_2880.wt.15mer.perread.h.tsv test_data/25s_2880.sn34.15mer.perread.h.tsv test_data/25s_2880.sn3.15mer.perread.h.tsv test_data/25s_2880.sn36.15mer.perread.h.tsv  knn_validation
```


![alt text](./img/KNN_plots.png "KNN_plots")
(Note: PCA is not used for stoichiometry prediction when using KNN - the PCA has only been used for illustration purposes)

