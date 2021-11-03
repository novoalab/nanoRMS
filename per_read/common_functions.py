"""
Here we store all functions that are used across Jupyter notebooks
"""

import csv, gzip, os, matplotlib.pyplot as plt, numpy as np, pandas as pd, pysam, sys
import seaborn as sns#; sns.set()
import eif_new as iso_new
from sklearn.cluster import AgglomerativeClustering, KMeans
from sklearn.svm import OneClassSVM
from sklearn.ensemble import IsolationForest, RandomForestClassifier
from sklearn.mixture import GaussianMixture, BayesianGaussianMixture
from sklearn.neighbors import KNeighborsClassifier
from datetime import datetime
from collections import Counter
from multiprocessing import Pool

# it's only DNA as in SAM U should be A
base2complement = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}

# nanopolish parser
def mer_generator(handle, k=15):
    """Report consecutive k-mers from nanopolish output"""
    # data handle
    pcontig, pread_name, mer_data = 0, 0, []
    rd = csv.reader(handle, delimiter="\t", quotechar='"')
    header = rd.__next__() #; print(header)
    for i, r in enumerate(rd):
        #if i>100000: break
        if not i%10000:
            sys.stderr.write(" %s \r"%i)
        contig, position, reference_kmer, read_name, strand, event_index, event_level_mean, event_stdv, event_length, model_kmer, model_mean, model_stdv, standardized_level, start_idx, end_idx = r[:15]
        # skip undetermined model_kmers
        #if "N" in model_kmer: continue
        # get int and float
        position, event_level_mean, event_length = int(position), float(event_level_mean), float(event_length)
        # start over
        if pcontig!=contig or pread_name!=read_name:
            #if len(mer_data)==k: yield pcontig, ppos, mer_data[:-1]
            pread_name, pcontig, ppos, mer_data = read_name, contig, position, [[]]
        # define data to store
        data = (event_level_mean, event_length)
        # add to previous mer
        if position == ppos: mer_data[-1].append(data)
        # add new mer position
        elif position == ppos+1: 
            mer_data.append([data, ])
            ppos = position
        # start new mer
        else: 
            if len(mer_data)==k: yield pcontig, ppos, mer_data[:-1]
            ppos, mer_data = position, [[data, ]]
        # report middle position only if full mer
        if len(mer_data)==k+1: # skip last pos, since it's still not complete
            yield pcontig, ppos-1, mer_data[:-1]
            mer_data = mer_data[1:]
            
def nanopolish2regions(fn, regions, nn=1, maxReads=2000):
    """Create dictionary of nanopolish regions"""
    k = 2*nn+1
    pos2data = {(ref, pos): [] for ref, pos, mt in regions}
    for ref, pos, data in mer_generator(gzip.open(fn, "rt"), k):
        if (ref, pos-nn) not in pos2data: continue
        # get weithted average of events at every position
        si = [np.average([e[0] for e in d], weights=[e[1] for e in d]) for d in data]
        pos2data[(ref, pos-nn)].append(si)
    return pos2data

# get coverage in reads per each reference position
def pass_filters(a, mapq=10):
    if a.mapq<mapq or a.is_secondary or a.is_supplementary or a.is_qcfail or a.is_duplicate: return False
    return True

def get_coverage(regions1, fnames1, sample2nanopolish1):
    """Return coverage from Nanopolish"""
    pos2count1 = {(ref, pos): [sum([1 for a in pysam.AlignmentFile(fn[:-10]).fetch(ref, pos-1, pos) if pass_filters(a)]) for fn in fnames1] 
                                   for ref, pos, mt in regions1}
    # get number of resquiggled reads from tombo
    tombo1 = ["guppy3.0.3.hac/%s/workspace/batch0.fast5.bam"%fn.split("/")[-2] for fn in fnames1]
    tombo_p2c1 = {(ref, pos): [sum([1 for a in pysam.AlignmentFile(fn).fetch(ref, pos-1, pos) if pass_filters(a)]) for fn in tombo1] 
                  for ref, pos, mt in regions1}
    # combine
    names1 = ["%s %s"%(n, fn.split("/")[-2]) for n in ("coverage", "nanopolish", "tombo") for fn in fnames1]
    df4c1 = pd.DataFrame([[ref, pos, *cov, *[len(sample2nanopolish1[i][(ref, pos)]) for i, fn in enumerate(fnames1)], *tombo_p2c1[(ref, pos)]] 
                          for (ref, pos), cov in pos2count1.items()], columns=["chrom", "pos", *names1])
    return df4c1

def get_coverage2(regions1, fnames1, sample2nanopolish1): 
    pos2count1 = {(ref, pos): [sum([1 for a in pysam.AlignmentFile(fn[:-10]).fetch(ref, pos-1, pos) if pass_filters(a)]) for fn in fnames1] 
                                   for ref, pos, mt in regions1}
    # get number of resquiggled reads from tombo
    tombo1 = ["guppy3.0.3.hac/%s/workspace/batch0.fast5.bam"%fn.split("/")[-2] for fn in fnames1]
    tombo_p2c1 = {(ref, pos): [sum([1 for a in pysam.AlignmentFile(fn).fetch(ref, pos-1, pos) if pass_filters(a)]) for fn in tombo1] 
                  for ref, pos, mt in regions1}
    nanopolish_p2c = {(ref, pos): [len(sample2nanopolish1[i][(ref, pos)]) for i, fn in enumerate(fnames1)] for ref, pos, mt in regions1}
    # combine
    dframes = []
    names = [fn.split("/")[-2].split("_")[-1] for fn in fnames1]
    for n, d in zip(("minimap2", "nanopolish", "tombo"), 
                    (pos2count1, nanopolish_p2c, tombo_p2c1)):
        df = pd.DataFrame([[ref, pos, n, *d[(ref, pos)]] for ref, pos, mt in regions1], columns=["chrom", "pos", "name", *names])
        dframes.append(df)
    df4c1 = pd.concat(dframes).reset_index()
    return df4c1

def get_coverage3(regions1, fnames1, sample2nanopolish1, mod="pU"):
    pos2count1 = {(ref, pos): [sum([1 for a in pysam.AlignmentFile(fn[:-10]).fetch(ref, pos-1, pos) if pass_filters(a)]) for fn in fnames1] 
                                   for ref, pos, mt in regions1}
    # get number of resquiggled reads from tombo
    tombo1 = ["guppy3.0.3.hac/%s/workspace/batch0.fast5.bam"%fn.split("/")[-2] for fn in fnames1]
    tombo_p2c1 = {(ref, pos): [sum([1 for a in pysam.AlignmentFile(fn).fetch(ref, pos-1, pos) if pass_filters(a)]) for fn in tombo1] 
                  for ref, pos, mt in regions1}
    nanopolish_p2c = {(ref, pos): [len(sample2nanopolish1[i][(ref, pos)]) for i, fn in enumerate(fnames1)] for ref, pos, mt in regions1}
    # combine
    dframes = []
    names = [fn.split("/")[-2].split("_")[-1] for fn in fnames1]
    strain2idx = {n: i for i, n in enumerate(names)}
    for n, d in zip(("nanopolish", "tombo"), (nanopolish_p2c, tombo_p2c1)):
        df = pd.DataFrame([[ref, pos, n, m, d[(ref, pos)][strain2idx[s]]/1000] for ref, pos, mt in regions1 
                           for s, m in zip(("wt", mt), (mod, "unmod"))], 
                          columns=["chrom", "pos", "name", "base", "resquiggled"])
        dframes.append(df)
    df = pd.concat(dframes).reset_index()
    return df

# FastA/BAM parsers
def get_revcomp(bases):
    """Return reverse comlement"""
    return "".join(base2complement[b] for b in bases[::-1])

def fasta2bases(fastafn, ref, start, end, strands="+-", n=3):
    """Generator of individual bases from FastA file.

    The output consists of: 
    - position in reference (1-based)
    - strand integer (0 for plus, 1 for minus)
    - strand as +/-
    - base (complement for -)
    - 7-mer (+/- n bases surrounding given position)
    """
    fasta = pysam.FastaFile(fastafn)
    ref2len = {r: l for r, l in zip(fasta.references, fasta.lengths)}
    if ref not in ref2len: #fasta.references:
        raise StopIteration
    for pos, refbase in enumerate(fasta.fetch(ref, start, end), start+1):
        refbase = refbase.upper()
        # combine before start NNN (if needed) sequence from ref and after start NNN (if needed)
        mer = "N"*(n-pos+1) + "".join(fasta.fetch(ref, pos-n-1 if pos>n+1 else 0, pos+n)) + "N"*(pos-ref2len[ref]+n)
        mer = mer.upper() # need to be upper case
        for si, strand in enumerate(strands):
            if si:
                refbase = base2complement[refbase]
                mer = get_revcomp(mer)
            yield pos, si, strand, refbase, mer

def moving_average(a, n=5):
    """Calculate moving average including first n-1 objects"""
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    ret[n-1:] *= 1 / n
    ret[:n-1] *= 1 / np.arange(1, n)
    return ret

def bam2data(bam, ref, start, end, rna=True, nn=1, features=["si", "tr"],
             maxDepth=100000, mapq=20, dtype="float16", verbose=1, logger=sys.stderr.write):
    """Generator of data for consecutive positions from ref:start-end region"""
    sam = pysam.AlignmentFile(bam)#; print(ref, start, end)
    # get dt_shift
    f2idx = {f: i for i, f in enumerate(features)}
    dt_shift_keys = list(filter(lambda k: k.startswith("dt") and k!="dt0", f2idx.keys()))
    dt_shift = 0 if not len(dt_shift_keys) else int(dt_shift_keys[0][2:]) # dt10 > 10
    # update end position with shift
    end += dt_shift # here for DNA it's a bit complicated as we'd need to do start-=dt_shift
    # this is needed later
    requested_tags = list(filter(lambda f: not f.startswith("dt"), features))
    if dt_shift or "dt0" in features: requested_tags.append("dt")
    # here only SI & MP # here dt_shift should be read from the feature

    id_tags = np.empty(maxDepth, dtype=object) # store id per read 
    id_tags[:] = ""

    data = np.empty((len(features), maxDepth, end-start), dtype=dtype)
    # solve missing trace at deletions in the read
    data[:] = -1 # store -1 instead of 0 (trace with -1 will be skipped)
    strands = np.zeros((maxDepth, end-start), dtype="int8") # 1: +/FOR; -1: -/REV; 0: no alignment/read
    readi = 0
    for a in sam.fetch(ref, start, end):
        # filter by mapq
        if a.mapq<mapq: continue
        # make sure first position of read always corresponds to first pos in data
        while a.pos>start: # consider skipping first/last 5-15 bases
            # report data for reads from + & - strand separately
            for strand in (1, -1):
                flt = strands[:readi, nn] == strand
                yield (id_tags[:readi].tolist(), data[:, :readi][:, flt, :2*nn+1])
            # strip position 0 from arrays
            data = data[:, :, 1:]
            strands = strands[:, 1:]
            start += 1
        # define read start & end for current data view
        s, e = start-a.pos, a.aend-a.pos if a.aend<=end else end-a.pos
        # and region end
        re = e-s
        # get data from tags
        tags = {k: v for k, v in a.tags}
        # turn exonic blocks back into genomic coordinates
        if "bs" in tags and len(tags["bs"])>2:
            # get blocks as 2D array (start-end) with exonic intervals of the read
            blocks = np.array(tags["bs"]).reshape(-1, 2)-tags["bs"][0]
            # take care only about requested features
            _tags = {}
            for f in requested_tags:
                # storing 1s is importand as dt is log2 obs/exp, thus it can't be 0s
                _tags[f] = np.ones(a.reference_length, dtype=dtype)
            # mark exonic block in strands
            read_strands = np.zeros(a.reference_length, dtype="int8")
            # store block info
            pe = 0
            for bs, be in blocks:
                # mark exonic block in read_strands
                read_strands[bs:be] = -1 if a.is_reverse else 1
                # store exon block into new tags
                blen = be-bs
                for f in requested_tags:
                    #print(f, bs, be, be-bs, pe, be-pe)
                    _tags[f][bs:be] = tags[f][pe:pe+blen]
                pe += blen
            # replace tags & udpate exonic strands
            tags = _tags
            strands[readi, :re] = read_strands[s:e]
        else:
            # mark entire read as stand
            strands[readi, :re] = -1 if a.is_reverse else 1
        # here we need to add special treatment for dt!
        if "dt0" in f2idx or dt_shift:
            # normalise dwell time by moving average and store as log2
            dt = np.array(tags["dt"])
            dt = np.log2(dt / moving_average(dt)) #dt.mean())
        # store
        for j, k in enumerate(features, 0): #for k, j in f2idx.items(): #
            # dwell-time for position 0
            if k=="dt0": data[j, readi, :re] = dt[s:e]
            # shifted dwell time
            elif k.startswith("dt"):
                if rna and not a.is_reverse or not rna and a.is_reverse: 
                    if e>s+dt_shift: # make sure enough alignment here # and len(dt[s+dt_shift:e]):
                        data[j, readi, :re-dt_shift] = dt[s+dt_shift:e]
                elif e-dt_shift>s: # and here as well len(dt[s:e-dt_shift]): 
                    data[j, readi, :re-dt_shift] = dt[s:e-dt_shift]
            # normalise trace - this isn't needed cause we'll do min-max anyway
            elif k.startswith("t"): data[j, readi, :re] = np.array(tags[k][s:e], dtype=dtype)/255
            # and remaining features
            else:
                data[j, readi, :re] = tags[k][s:e]
        
        id_tags[readi] = tags["ID"] if 'ID' in tags.keys() else ''

        readi += 1
        # clean-up only if maxDepth reached
        if readi>=maxDepth:
            if verbose: logger("[INFO] maxDepth reached for %s:%s-%s @ %s\n"%(ref, start, end, bam))
            # get algs that still carry data
            ## here read has strand over from 0 to end (not excluding introns)
            ne = strands[:, 0]!=0 # np.all(strands!=0, axis=0)#?
            readi = ne.sum() # update readi
            if readi>=maxDepth: # if all reads are still aligned, remove random 25% of algs
                ne[np.random.randint(0, len(ne), int(0.25*maxDepth))] = False
                readi = ne.sum() # update readi
            # update strands & data
            _strands, _data, _id_tags = np.zeros_like(strands), np.zeros_like(data),  np.zeros_like(id_tags)
            _strands[:readi] = strands[ne]
            _data[:, :readi] = data[:, ne]
            _id_tags[:readi] = id_tags[ne]
            strands, data, id_tags = _strands, _data, _id_tags
    # report last bit from region
    for pos in range(strands.shape[1]-nn):
        # report data for reads from + & - strand separately
        for strand in (1, -1):
            flt = strands[:readi, pos+nn] == strand

            yield (id_tags[:readi].tolist(), data[:, :readi][:, flt, pos:pos+2*nn+1])

# functions we'll need to load the data 
def load_data(fasta, bams, regions, features, max_reads=1000, strands="+-", nn=1):
    """Return features for positions of interest"""
    # get storage
    k = 2*nn+1
    fi = 0
    sam = pysam.AlignmentFile(bams[0])
    region2data = {}
    for ri, (ref, pos, _) in enumerate(regions, 1):
        sys.stderr.write(" %s / %s %s:%s \r"%(ri, len(regions), ref, pos))
        start, end = pos-1, pos
        # extend start/end by nn and end by dt_shift
        ##this is for RNA, for DNA start start needs to be -dt_shift
        parsers = [bam2data(bam, ref, start-nn if start>=nn else 0, end+2*nn, True, 
                            nn, features, max_reads) for bam in bams]
        refparser = fasta2bases(fasta, ref, start, end, strands)
        for ((pos, _, strand, refbase, mer), *calls) in zip(refparser, *parsers):
            if strand=="+":
                region2data[(ref, pos)] = (mer, [np.hstack(c) for c in calls])
    return region2data

def load_data_reps(fasta, bams, regions, features, strains, strains_unique, maxReads=10000, strands="+-", nn=1):
    """Return features for positions of interest"""
    # get storage
    k = 2*nn+1
    fi = 0
    strain2idx = {s: idx for idx, s in enumerate(strains_unique)}
    region2data = {}
    for ri, (ref, pos, strand) in enumerate(regions, 1):
        if type(strand)==float: strand="+" # sometimes strand is missing, assume +
        start, end = pos-1, pos
        sys.stderr.write(" %s / %s %s:%s-%s \r"%(ri, len(regions), ref, start, end))
        # extend start/end by nn and end by dt_shift
        ##this is for RNA, for DNA start start needs to be -dt_shift
        parsers = [bam2data(bam, ref, start-nn if start>=nn else 0, end+2*nn, True, 
                            nn, features, maxReads) for bam in bams]
        refparser = fasta2bases(fasta, ref, start, end, strands)
        for ((pos, _, _strand, refbase, mer), *calls) in zip(refparser, *parsers):
            if _strand==strand:
                sdata = [[], []] #np.hstack(c) for c in calls]
                sid = [[], []]
                for c, s in zip(calls, strains): 
                    sdata[strain2idx[s]].append(np.hstack(c[1])) # feature data
                    sid[strain2idx[s]].extend(c[0]) # read ids
                # merge replicates
                region2data[(pos, strand)] = (mer, sid, [np.vstack(sd) for sd in sdata])

    return region2data

def get_data_mix(unmod, mod, frac, max_reads):
    """Return sample containing mod[:frac*max_reads] and unmod[:(1-frac)*max_reads]"""
    mod_n = int(round(frac*max_reads))
    unmod_n = max_reads-mod_n #int(round((1-frac)*max_reads))
    return np.vstack([unmod[:unmod_n], mod[:mod_n]])

def load_data_stoichometry(fasta, bams, regions, features, samples, fracs, 
                           maxReads=1000, strands="+-", nn=1):
    """Return features for positions of interest"""
    # get storage
    k = 2*nn+1
    fi = 0
    sam = pysam.AlignmentFile(bams[0])
    region2data = {}
    sample2idx = {s: i for i, s in enumerate(samples)}; print(sample2idx)
    for ri, (ref, pos, mt) in enumerate(regions, 1):
        sys.stderr.write(" %s / %s %s:%s \r"%(ri, len(regions), ref, pos))
        start, end = pos-1, pos
        # extend start/end by nn and end by dt_shift
        ##this is for RNA, for DNA start start needs to be -dt_shift
        parsers = [bam2data(bam, ref, start-nn if start>=nn else 0, end+2*nn, True, 
                            nn, features, maxReads) for bam in bams]
        refparser = fasta2bases(fasta, ref, start, end, strands)
        for ((pos, _, strand, refbase, mer), *calls) in zip(refparser, *parsers):
            if strand=="+":
                sample2data = [np.hstack(c) for c in calls]
                # get min number of reads
                max_reads = int(min(map(len, sample2data))/3)#; print(ref, pos, mt, max_reads, [s.shape for s in sample2data])
                # first get 2 fully unmodified and 1 fully modified sample - those reads won't be used later on
                data_frac = [sample2data[sample2idx[mt]][max_reads:2*max_reads], # this will be used as 0 sample
                             sample2data[sample2idx[mt]][-max_reads:], sample2data[sample2idx["wt"]][-max_reads:], # those two will be training set
                            ] 
                # the get samples with given fractions of modified reads
                data_frac += [get_data_mix(sample2data[sample2idx[mt]], 
                                           sample2data[sample2idx["wt"]], frac, max_reads) 
                              for frac in fracs]
                region2data[(ref, pos)] = (mer, data_frac)
    return region2data

def load_data_train_test_val(fasta, bams, regions, features, samples, maxReads=1000, strands="+-", nn=1):
    """Return features for positions of interest"""
    # get storage
    k = 2*nn+1
    fi = 0
    sam = pysam.AlignmentFile(bams[0])
    region2data = {}
    sample2idx = {s: i for i, s in enumerate(samples)}; print(sample2idx)
    for ri, (ref, pos, mt) in enumerate(regions, 1):
        sys.stderr.write(" %s / %s %s:%s \r"%(ri, len(regions), ref, pos))
        start, end = pos-1, pos
        # extend start/end by nn and end by dt_shift
        ##this is for RNA, for DNA start start needs to be -dt_shift
        parsers = [bam2data(bam, ref, start-nn if start>=nn else 0, end+2*nn, True, 
                            nn, features, maxReads) for bam in bams]
        refparser = fasta2bases(fasta, ref, start, end, strands)
        for ((pos, _, strand, refbase, mer), *calls) in zip(refparser, *parsers):
            if strand=="+":
                sample2data = [np.hstack(c) for c in calls]
                # get min number of reads
                max_reads = int(min(map(len, sample2data))/3)#; print(ref, pos, mt, max_reads, [s.shape for s in sample2data])
                # first get 2 fully unmodified and 1 fully modified sample - those reads won't be used later on
                data_frac = [sample2data[sample2idx[mt]][max_reads:2*max_reads], # this will be used as 0 sample
                             sample2data[sample2idx[mt]][-max_reads:], sample2data[sample2idx["wt"]][-max_reads:], # those two will be training set
                            ] 
                # get a bit of every sample
                data_frac += [sd[:max_reads] for sd in sample2data]
                region2data[(ref, pos)] = (mer, data_frac)
    return region2data

# functions we'll need to plot
def get_modfreq_from_quantiles_many_samples(scores_per_sample, q=0.1):
    """Return modification frequency calculated using quantiles method"""
    freqs = np.zeros(len(scores_per_sample))
    minc = min(map(len, scores_per_sample))
    q1, q2 = np.quantile(np.concatenate([s[:minc] for s in scores_per_sample]), [q, 1-q])
    for i, _scores in enumerate(scores_per_sample): 
        confs = [(_scores<q1).sum(), (_scores>q2).sum()]
        if not sum(confs): continue
        mod_freq = confs[1]/sum(confs)
        freqs[i] = mod_freq
    return freqs

def get_mod_freq_clf(df, cols, chr_pos, strains, clf, method="GMM"):
    """Predict modification frequency using single classifier"""
    results = []
    for cp in chr_pos:
        # min-max normalisation
        _df = df.loc[(df["chr_pos"]==cp)&(df.Strain.isin(strains)), cols+["Strain"]]
        _X = min_max_norm(_df[cols].to_numpy().astype("float"))
        # get fit and clusters
        clusters = clf.fit_predict(_X)
        # for outlier method, store outliers (-1) as cluster_1 and normal (1) as cluster_0
        if max(clusters)>1: 
            clusters[clusters!=0] = 1
        elif -1 in clusters and 1 in clusters: # outlier method
            clusters[clusters==1] = 0
            clusters[clusters<0] = 1
        # get modification freuqency - simply number of 1s over all for each sample
        freqs = [clusters[_df["Strain"]==s].mean() for s in strains]
        results.append((cp, method, *freqs, ", ".join(map(str, strains[1:]))))
    return results
            
def min_max_norm(X):
    """Return (X-min(X))/(max(X)-min(X))"""
    #return X # no min_max_norm ;)
    Xmax, Xmin = X.max(axis=0), X.min(axis=0)
    sel = Xmin!=Xmax
    if sel.sum():
        X[:, sel] = (X[:, sel] - Xmin[sel])/(Xmax[sel] - Xmin[sel]) # here if min==max div by 
    return X

def get_mod_freq_two_step(df, cols, chr_pos, strains, method="GMM+eIF", clf_name="GMM", 
                          clf=GaussianMixture(n_components=4, random_state=0), 
                          clf2_name="eIF", clf2=iso_new.iForest(random_state=0), 
                          OFFSET=None):
    """Predict modification frequency using """
    results = []
    for cp in chr_pos:
        _df = df.loc[(df["chr_pos"]==cp)&(df.Strain.isin(strains)), cols+["Strain"]]
        _X = min_max_norm(_df[cols].to_numpy().astype("float"))
        # get clusters from GMM using only SIGNAL INTENSITY
        clusters = clf.fit_predict(_X) #[:,:3]
        c2i = Counter(clusters)#; print(c2i)
        # get outliers using every cluster as training sset
        mod_freqs = np.zeros((len(c2i), len(strains)))
        mod_freqs1 = np.zeros_like(mod_freqs)
        for cl in list(c2i.keys())[:3]:
            Xtrain = _X[clusters==cl]
            if len(Xtrain)<3: continue # this is arbitrary value
            scores = clf2.fit(Xtrain).score_samples(_X)
            offset = (max(scores)-min(scores))/2 if not OFFSET else OFFSET
            y_pred = scores>offset
            # get mod_freq from outlier score cut-off
            mod_freqs1[cl] = [y_pred[_df["Strain"]==s].mean() for s in strains]
            # and using quantile method
            mod_freqs[cl] = get_modfreq_from_quantiles_many_samples([scores[_df["Strain"]==s] for s in strains])

        # pick cluster that gave the largest difference in mod_freq between any two samples
        extremes = np.vstack([np.nanmin(mod_freqs, axis=1), np.nanmax(mod_freqs, axis=1)])
        mod_freq_idx = np.abs(np.diff(extremes, axis=0)).argmax()#; print(mod_freq_idx)
        # and report
        #results.append((cp, "%s+%s_c"%(clf_name, clf2_name), *mod_freqs1[mod_freq_idx], 
        #                ", ".join(map(str, strains[1:]))))
        results.append((cp, method, *mod_freqs[mod_freq_idx], ", ".join(map(str, strains[1:]))))
    return results

def get_mod_freq_clf_train_test(df, cols, chr_pos, strains, train_samples, 
                                clf=KNeighborsClassifier(), method="KNN"):
    """Predict modification frequency using single classifier"""
    results = []
    for cp in chr_pos:
        # train classifier using train sampels: unmod and mod
        _df = df.loc[(df["chr_pos"]==cp)&(df.Strain.isin(train_samples)), cols+["Strain"]]
        X_train = min_max_norm(_df[cols].to_numpy().astype("float"))
        y_train = _df.Strain==train_samples[-1]
        clf.fit(X_train, y_train)
        # min-max normalisation
        _df = df.loc[(df["chr_pos"]==cp)&(df.Strain.isin(strains)), cols+["Strain"]]
        _X = min_max_norm(_df[cols].to_numpy().astype("float"))
        # get fit and clusters
        clusters = clf.predict(_X) # this will return 0 (unmodified) and 1 (modified)
        # get modification freuqency - simply number of 1s over all for each sample
        freqs = [clusters[_df["Strain"]==s].mean() for s in strains]
        results.append((cp, method, *freqs, ", ".join(map(str, strains[1:]))))
    return results

def generate_figures_and_xls(outdir, cols_starts, region2data, ext, xls, group2pos, feature_names, samples):
    """Generate figures and tables"""
    all_freqs = []
    # concatenate all pos and samples into one dataframe
    dframes = []
    for ri, (ref, pos) in enumerate(region2data.keys()): #regions): #[3]#; print(ref, pos, mt)
        mer, calls = region2data[(ref, pos)]
        for c, s in zip(calls, samples): 
            df = pd.DataFrame(c, columns=feature_names)
            df["Strain"] = s
            df["chr_pos"] = "%s:%s"%(ref, pos)
            dframes.append(df)
    # read all tsv files
    df = pd.concat(dframes).dropna().reset_index()
    chr_pos, strains = df["chr_pos"].unique(),  df["Strain"].unique()    
    
    # compare individual methods
    for clf, method in (
                        (iso_new.iForest(ntrees=100, random_state=0), "GMM+eIF"), 
                        (GaussianMixture(random_state=0, n_components=2), "GMM"), 
                        (AgglomerativeClustering(n_clusters=2), "AggClust"), 
                        (KMeans(n_clusters=2), "KMeans"), 
                        (OneClassSVM(), "OCSVM"), 
                        (IsolationForest(random_state=0), "IF"), 
                        (iso_new.iForest(ntrees=100, random_state=0), "eIF"), 
                        (KNeighborsClassifier(), "KNN"), 
                        (RandomForestClassifier(), "RF"), 
                        ):
        fname = method
        print(fname)
        outfn = os.path.join(outdir, "%s.%s"%(fname, ext))        
        results = []
        for i, cols_start in enumerate(cols_starts, 1):
            # narrow down the features to only signal intensity & trace
            cols = list(filter(lambda n: n.startswith(cols_start), feature_names)); cols #, "DT"
            # compare all samples to 0%
            s0 = samples[0]
            for s in samples[3:]: 
                with np.errstate(under='ignore'):
                    if "+" in method:
                        clf2_name = method.split("+")[-1]
                        results += get_mod_freq_two_step(df, cols, chr_pos, [s0, s], "_".join(cols_start),  
                                                         OFFSET=0.5, clf2_name=clf2_name, clf2=clf)
                    elif method in ("KNN", "RF"):
                        results += get_mod_freq_clf_train_test(df, cols, chr_pos, [s0, s], samples[1:3], clf, "_".join(cols_start))
                    else:
                        results += get_mod_freq_clf(df, cols, chr_pos, [s0, s], clf, "_".join(cols_start))
                    
        # and store mod_freq predicted by various methods
        freqs = pd.DataFrame(results, columns=["chr_pos", "features", "mod_freq wt", "mod_freq strain", "strain"])
        freqs["diff"] = freqs.max(axis=1)-freqs.min(axis=1); freqs
        for name, pos in group2pos.items(): #(("negative", negatives), ("pU", pU_pos), ("Nm", Nm_pos)):
            freqs.loc[freqs["chr_pos"].isin(pos), "group"] = name
        #freqs.to_csv(outfn, sep="\t"); freqs.head()
        freqs.to_excel(xls, fname, index=False)
        # plot differences between methods
        for group, pos in group2pos.items():
            freqs.loc[freqs["chr_pos"].isin(pos), "modification"] = group
        #g = sns.catplot(x="strain", y="diff", hue="features", col="modification", data=freqs, kind="box")#, palette="Blues")
        g = sns.catplot(x="strain", y="diff", hue="features", col="modification", data=freqs, kind="point", ci=None)#, palette="Blues")
        fig = g.fig
        fig.suptitle(method)
        for ax in fig.axes:
            ax.set_xlabel("Expected mod_freq")
            ax.set_ylabel("Observed mod_freq [absolute difference between wt & mt]")
            ax.set_ylim(0, 1)
        fig.savefig(outfn)
        plt.close() # clear axis
        freqs["name"] = fname
        all_freqs.append(freqs)
    return all_freqs

def generate_figures_and_xls_all_strains(outdir, cols_starts, region2data, ext, xls, group2pos, feature_names, samples):
    """Generate figures and tables"""
    all_freqs = []
    # concatenate all pos and samples into one dataframe
    dframes = []
    for ri, (ref, pos) in enumerate(region2data.keys()): #regions): #[3]#; print(ref, pos, mt)
        mer, calls = region2data[(ref, pos)]
        for c, s in zip(calls, samples): 
            df = pd.DataFrame(c, columns=feature_names)
            df["Strain"] = s
            df["chr_pos"] = "%s:%s"%(ref, pos)
            dframes.append(df)
    # read all tsv files
    df = pd.concat(dframes).dropna().reset_index()
    chr_pos, strains = df["chr_pos"].unique(),  df["Strain"].unique()    
    # compare individual methods
    for clf, method in (
                        (KMeans(n_clusters=2), "KMeans"), 
                        (KNeighborsClassifier(), "KNN"), 
                        #(iso_new.iForest(ntrees=100, random_state=0), "GMM+eIF"), 
                        (GaussianMixture(random_state=0, n_components=2), "GMM"), 
                        (AgglomerativeClustering(n_clusters=2), "AggClust"), 
                        #(OneClassSVM(), "OCSVM"), 
                        (IsolationForest(random_state=0), "IF"), 
                        #(iso_new.iForest(ntrees=100, random_state=0), "eIF"), 
                        (RandomForestClassifier(), "RF"), 
                        ):
        fname = method
        for i, cols_start in enumerate(cols_starts, 1):
            results = []
            feat_name = "_".join(cols_start)
            fname = "%s.%s"%(method, feat_name); print(fname)
            outfn = os.path.join(outdir, "%s.%s"%(fname, ext))
            # narrow down the features to only signal intensity & trace
            cols = list(filter(lambda n: n.startswith(cols_start), feature_names))#; print(cols) #, "DT"
            # compare all samples to 0%
            s0 = samples[0]
            for s in samples[3:]: 
                with np.errstate(under='ignore'):
                    if "+" in method:
                        clf2_name = method.split("+")[-1]
                        results += get_mod_freq_two_step(df, cols, chr_pos, [s0, s], feat_name,  
                                                         OFFSET=0.5, clf2_name=clf2_name, clf2=clf)
                    elif method in ("KNN", "RF"):
                        results += get_mod_freq_clf_train_test(df, cols, chr_pos, [s0, s], samples[1:3], clf, feat_name)
                    else:
                        results += get_mod_freq_clf(df, cols, chr_pos, [s0, s], clf, feat_name)
                    
            # and store mod_freq predicted by various methods
            freqs = pd.DataFrame(results, columns=["chr_pos", "features", "mod_freq wt", "mod_freq strain", "strain"])
            freqs["diff"] = freqs.max(axis=1)-freqs.min(axis=1); freqs
            for name, pos in group2pos.items(): #(("negative", negatives), ("pU", pU_pos), ("Nm", Nm_pos)):
                freqs.loc[freqs["chr_pos"].isin(pos), "group"] = name
            #freqs.to_csv(outfn, sep="\t"); freqs.head()
            freqs.to_excel(xls, fname, index=False)
            # plot differences between methods
            for group, pos in group2pos.items():
                freqs.loc[freqs["chr_pos"].isin(pos), "modification"] = group
            #return freqs
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 5))#, sharey="all")
            sns.barplot(x="chr_pos", y="mod_freq strain", hue="strain", edgecolor="white", palette=["#f8786fff", "#7aae02ff", "#00bfc2ff", "#c67afeff"], 
                        data=freqs[(freqs["features"]==feat_name)&(freqs["group"]=="pU")], ax=ax1)
            sns.barplot(x="chr_pos", y="mod_freq strain", hue="strain", edgecolor="white", palette=["#ed823aff", "#1c6ca9ff", "#35d1bbff", "#c978fdff"], 
                        data=freqs[(freqs["features"]==feat_name)&(freqs["group"]=="Nm")], ax=ax2)
            ax1.set_ylabel("Per-site stoichiometry"); ax2.set_ylabel("")
            ax1.get_legend().remove(); ax2.get_legend().remove()#ax1.legend([]); ax2.legend([])
            ax1.set_ylim(0, 1); ax2.set_ylim(0, 1); #ax2.set(aspect=1.7)
            ax1.set_title("pU modifications"); ax2.set_title("Nm modifications")
            fig.suptitle(fname)
            fig.savefig(outfn)
            plt.close() # clear axis
            freqs["name"] = fname
            all_freqs.append(freqs)
    return all_freqs

def plot_figures(outdir, df, mt, strains_unique, hue=[], ext="pdf"):
    # join with predictionsxs
    fnames = df["features"].unique()
    fig, axes = plt.subplots(len(fnames), 1, figsize=(12, 5*len(fnames)))
    fig.suptitle(mt)
    #df_predicted = df[df["Prediction"]=="Predicted"]
    #sns.boxplot(x="method", y="diff", hue="features", data=df_predicted, ax=ax) # 
    # plot boxplot with stipplot
    for ai, (ax, fname) in enumerate(zip(axes, fnames)): 
        sns.boxplot(x="method", y="diff", hue="group", data=df[df["features"]==fname], ax=ax, color=".8", showfliers=False)#, width=0.8)
        sns.stripplot(x="method", y="diff", hue="group", data=df[df["features"]==fname], ax=ax, dodge=True)
        ax.set_ylabel("Absolute difference between %s & WT"%mt)
        ax.set_title(fname); ax.set_xlabel("")
        if not ai: ax.legend(bbox_to_anchor=(0, 1.1, 1, 0), loc="lower left", mode="expand", ncol=2) #bbox_to_anchor=(1.01, 1), borderaxespad=0)
        else: ax.get_legend().remove() # get rid of legend for subsequent plots
    fig.savefig(os.path.join(outdir, "%s.boxplot.%s"%(mt, ext)))
    
    # plot scatterplot
    methods = df.method.unique() #fnames = df.features.unique()
    groups = df.group.unique(); groups
    colors =  sns.color_palette(n_colors=len(groups))#"flare"
    markers = [".", "1", "2", "o", "o"]
    f = "SI_TR"
    fig, axes = plt.subplots(1, len(methods), figsize=(5*len(methods), 5), sharex="all", sharey="all")
    for ai, (ax, m) in enumerate(zip(axes, methods)):
        #g = sns.scatterplot(*strains_unique[::-1], hue="group", data=df[(df["features"]==f)&(df["method"]==m)], ax=ax); ax.get_legend().remove()
        for c, g, r in zip(colors, groups, markers): 
            ax.scatter(*strains_unique[::-1], color=c, alpha=0.75, marker=r, label=g, 
                        data=df[(df["features"]==f)&(df["method"]==m)&(df["group"]==g)])
        ax.set_title(m)
        ax.set_xlabel(strains_unique[1]); ax.set_ylabel(strains_unique[0])
        ax.plot(np.linspace(0, 1, 50), np.linspace(0, 1, 50), "grey")
    lgd = ax.legend(bbox_to_anchor=(1.01, 1), borderaxespad=0)
    #ax.legend(bbox_to_anchor=(-2.5, 1.1, 2.5, 0), loc="lower left", ncol=3)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    fig.suptitle("{} {}".format(mt, f))
    fig.savefig(os.path.join(outdir, "%s.scatter.%s"%(mt, ext)), bbox_extra_artists=(lgd,), bbox_inches='tight')
        
def plot_boxplot(outdir, df, mt, method, ext="pdf"):
    fnames = df["features"].unique()
    fig, axes = plt.subplots(len(fnames), 1, figsize=(7, 5*len(fnames)))
    fig.suptitle(mt)
    df = df.sort_values("New_Status")
    # plot boxplot with stipplot
    for ai, (ax, fname) in enumerate(zip(axes, fnames)): 
        sns.boxplot(x="New_Status", y="diff", hue="Prediction", data=df[df["features"]==fname], ax=ax, color=".8", showfliers=False)#, width=0.8)
        sns.stripplot(x="New_Status", y="diff", hue="Prediction", data=df[df["features"]==fname], ax=ax, dodge=True)
        ax.set_ylabel("Absolute difference between %s & WT"%mt)
        ax.set_title(fname); ax.set_xlabel("")
        if not ai: ax.legend(bbox_to_anchor=(0, 1.1, 1, 0), loc="lower left", mode="expand", ncol=2) #bbox_to_anchor=(1.01, 1), borderaxespad=0)
        else: ax.get_legend().remove() # get rid of legend for subsequent plots
    fig.savefig(os.path.join(outdir, "%s.boxplot.%s.%s"%(mt, method, ext)))
    
def plot_density(outdir, sdata, mt, group, ref, pos, strand, mer, feature_names, colors, ext="pdf"):
    """Plot and save density plot for given position"""
    fig, axes = plt.subplots(1, len(feature_names), figsize=(4*len(feature_names), 4))
    fig.suptitle("{} {} {}:{}{} {}".format(mt, group, ref, pos, strand, mer))
    for fi, (ax, f) in enumerate(zip(axes, feature_names)):
        for si, (s, c) in enumerate(zip((mt, "wt"), colors)):
            sns.kdeplot(sdata[si][:, fi], color=c, linewidth=2, shade=True, alpha=.5, legend=False, ax=ax)
            ax.set_xlabel(f); ax.set_ylabel("")
        axes[0].set_ylabel("Density")
    fig.savefig(os.path.join(outdir, "{}:{}{}.{}".format(ref, pos, strand, ext)))
    plt.close()

# classifiers and mod_freq estimators
def get_freq(y_pred, cov):
    freq = []
    ps = 0
    for c in cov:
        freq.append(y_pred[ps:ps+c].mean())
        ps+=c
    return freq

def get_freq_clf(region2data, strains_unique, cols_starts, feature_names, clf=KNeighborsClassifier(), clf_name="KNN"):
    """Return data frame"""
    rows = []
    for cols_start in cols_starts:
        features = "_".join(cols_start)
        cidx = [i for i, n in enumerate(feature_names) if n.startswith(cols_start)]
        sys.stderr.write(" %s          \r"%(features, ))
        for (ref, pos, strand), (mer, data) in region2data.items():
            pos_info = "{}:{}{}".format(ref, pos, strand)
            cov = list(map(len, data))
            X = np.vstack(data)[:, cidx] # get only columns corresponding to features of interests
            #X = min_max_norm(X) # minmax_norm
            y = np.zeros(len(X)) # KO
            y[len(data[0]):] = 1 # WT - here many may be unmodified
            clf.fit(X, y) # here we train and predict on the same dataset
            y_pred = clf.predict(X)
            freq = get_freq(y_pred, cov)
            #print(ref, pos, strand, cov, freq)
            rows.append((pos_info, clf_name, features, *cov, *freq))
    # get df with all predicitons
    df = pd.DataFrame(rows, columns=["chr_pos", "method", "features", *["%s cov"%s for s in strains_unique], *strains_unique])
    df["diff"] = abs(df[strains_unique[1]]-df[strains_unique[0]])
    return df
