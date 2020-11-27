# NanoRMS: predicting NANOpore Rna Modification Stoichiometry

Below you'll find information how to retrieve pre-read features and estimate modification stoichiometry from Nanopore raw data. 

0. First install all dependencies
- from conda: `conda install hdf5 minimap2`
- from pip: `pip install jupyterlab matplotlib numba numpy ont-fast5-api ont-tombo[full] openpyxl pandas seaborn scipy sklearn`
- manually: [nanopolish](https://github.com/jts/nanopolish)

1. Then download the data
```bash
# downsampled Fast5 for CC, rRNA & ncRNA (7GB)
wget https://public-docs.crg.es/enovoa/public/lpryszcz/src/nanoRMS/per_read/guppy3.0.3.hac -q --show-progress -r -c -nc -np -nH --cut-dirs=6 --reject="index.html*"

# pre-computed BAMs for mRNA (42GB) - this is only needed if you wish to analyse yeast mRNA for pseudoU
wget https://public-docs.crg.es/enovoa/public/lpryszcz/src/nanoRMS/per_read.bam/bams -q --show-progress -r -c -nc -np -nH --cut-dirs=6 --reject="index.html*"
```

CC, rRNA and ncRNA are downsampled to up to 1,000 reads per each reference
in order to facillicate reasonable size and speed-up analysis.
It wasn't feasible to downsample mRNA dataset, thus you'll need to download it from SRA (see paper for accession). 

2. Retrieve per-read features from all samples.
```bash
# for rRNA and CC # ~10 minutes
./get_features.py --rna -f cc_yeast_rrna.fa -t 6 -i guppy3.0.3.hac/{CC,*snR,*wt}*/workspace/*.fast5

# for yeast ncRNA # ~10 minutes
./get_features.py --rna -f guppy3.0.3.hac/Saccharomyces_cerevisiae.R64-1-1_firstcolumn.ncrna.fa -t 6 -i guppy3.0.3.hac/*WT??C/workspace/*.fast5
# and link output files
cd guppy3.0.3.hac && for f in *_WT??C/workspace/batch0.fast5.bam*; do d=`echo $f|cut -f1 -d'/'`; ln -s $f $d.`basename $f|cut -f3- -d"."`; done && ../
```

BAM files from yeast mRNAs Pus1/4 KO and heat-shock are already preprocessed. For the raw Fast5 files check the paper. 

3. Run nanopolish if you wish to compare it with tombo. You can skip this step if you don't want to compare the tools. 
```bash
ref=cc_yeast_rrna.fa 
for d in guppy3.0.3.hac/*_{snR*,wt}/; do
 if [ ! -s $d/fq.gz.bam.events.gz ]; then
  echo `date` $d;
  cat $d/*.fastq.gz > $d/fq.gz;
  nanopolish index -d $d $d/fq.gz 2> /dev/null;
  minimap2 -ax map-ont $ref $d/fq.gz 2> /dev/null | samtools sort -o $d/fq.gz.bam;
  samtools index $d/fq.gz.bam;
  nanopolish eventalign -n -t 6 -q 10 --progress --signal-index --scale-events --reads $d/fq.gz --bam $d/fq.gz.bam --genome $ref | gzip > $d/fq.gz.bam.events.gz;
 fi;
done; date
```

4. Finally you can analyse the data using one of the Jupyter Notebooks:
- rRNA_density.ipynb: plot feature densities for rRNA (Fig S6)
- rRNA_stoichometry.ipynb: benchmarking of methods and features for varying stoichiometries (Fig S5C)  
and estimation of modification frequency across all samples (Fig 3J)
- rRNA_tombo_vs_nanopolish.ipynb: comparison of Nanopolish and Tombo (Fig S5 A-B)
- mRNA.ipynb: 
