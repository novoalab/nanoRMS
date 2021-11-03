#!/usr/bin/env python3
desc="""Requiggle basecalled FastQ files and features in BAM file. 

For all reference bases we store (as BAM comments):
- normalised signal intensity mean [tag si:B,f]
- reference base probability [tag tr:B:C] retrieved from guppy (trace scaled 0-255)
- dwell time [tag dt:B:C] in signal step capped at 255

All features are matched versus padded reference sequnce blocks 
ie excluding introns and large (padded) deletions from reference. 
Those blocks (2-D array of start & ends) are stored as flattened 1-D array [tag bs:B:i]
ie. exons [(8114, 8244), (8645, 8797)] will be stored as array('I', [8114, 8244, 8645, 8797]). 

--rna will automatically enable spliced alignments. 
"""
epilog="""Author: l.p.pryszcz+git@gmail.com
Cologne/Barcelona/Mizerów, 17/06/2020
"""

import itertools, json, os, resource, scipy, subprocess, sys, numpy as np, pysam, tempfile
from tombo import tombo_stats, resquiggle, tombo_helper
from tombo._default_parameters import OUTLIER_THRESH, SHIFT_CHANGE_THRESH, SCALE_CHANGE_THRESH, RNA_SAMP_TYPE, DNA_SAMP_TYPE, COLLAPSE_RNA_STALLS, COLLAPSE_DNA_STALLS, STALL_PARAMS#, FM_OFFSET_DEFAULT
from ont_fast5_api.fast5_interface import get_fast5_file
from datetime import datetime
from multiprocessing import Pool
from array import array
from copy import deepcopy
# add PATH - needed by fast5_to_fastq.py
os.environ["PATH"] = "%s:%s"%(':'.join(sys.path), os.environ["PATH"]) 

VERSION = '0.11b'
DEFAULT_STALL_PARAMS = tombo_helper.stallParams(**STALL_PARAMS)
USE_START_CLIP_BASES = resquiggle.USE_START_CLIP_BASES

# only DNA bases as in SAM U is always referred as T
bases = "ACGT"
base2idx = {b: i for i, b in enumerate(bases)}
base2complement = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
# add lower-case for get_aligned_pairs as it reports substitutions as lower-case
for b, i in list(base2idx.items()): base2idx[b.lower()] = i
for b, c in list(base2complement.items()): base2complement[b.lower()] = c

def minimap2_proc(ref, fast5, threads=1, spliced=0, sensitive=1): 
    """Run minimap2 and return its stdout"""
    mode = ["-axmap-ont", ]
    if spliced:
        mode = ["-axsplice", "-uf"]
    args1 = ["minimap2", "--MD", "-Y", "-t%s"%threads] + mode
    if sensitive:
        args1 += ["-k7", "-w5", "-m20", "-A6", "-B4"]
    args1 += [ref, "-"]
    # fast5_to_fastq
    args0 = ["fast5_to_fastq.py", "-i%s"%fast5]
    proc0 = subprocess.Popen(args0, stdout=subprocess.PIPE)
    # minimap2
    proc1 = subprocess.Popen(args1, stdin=proc0.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return proc1

def adjust_map_res(map_res, seq_samp_type, rsqgl_params, TRIM_RNA_ADAPTER=False):
    if seq_samp_type.name == RNA_SAMP_TYPE:
        if TRIM_RNA_ADAPTER:
            # trim DNA adapter off of RNA signal
            adapter_end = tombo_stats.trim_rna(map_res.raw_signal, rsqgl_params)
            # trim off adapter
            map_res = map_res._replace(raw_signal=map_res.raw_signal[adapter_end:])

        # flip raw signal for re-squiggling
        map_res = map_res._replace(raw_signal=map_res.raw_signal[::-1])

    elif seq_samp_type.name == DNA_SAMP_TYPE and USE_START_CLIP_BASES:
        # flip raw signal, genome and start clip seqs for re-squiggling
        map_res = map_res._replace(
            raw_signal=map_res.raw_signal[::-1],
            genome_seq=map_res.genome_seq[::-1])

    if ((COLLAPSE_RNA_STALLS and seq_samp_type.name == RNA_SAMP_TYPE) or
        (COLLAPSE_DNA_STALLS and seq_samp_type.name == DNA_SAMP_TYPE)):
        map_res = map_res._replace(stall_ints=tombo_stats.identify_stalls(map_res.raw_signal, DEFAULT_STALL_PARAMS))

    return map_res

def adjust_rsqgl_res(rsqgl_res, all_raw_signal, seq_samp_type, USE_START_CLIP_BASES=False):
    if seq_samp_type.name == DNA_SAMP_TYPE and USE_START_CLIP_BASES:
        # flip raw signal and events back for storage in genome direction
        rev_rsrtr = (all_raw_signal.shape[0] -
                     rsqgl_res.read_start_rel_to_raw -
                     rsqgl_res.segs[-1])
        rev_segs = -1 * (rsqgl_res.segs[::-1] - rsqgl_res.segs[-1])
        rsqgl_res = rsqgl_res._replace(
            read_start_rel_to_raw=rev_rsrtr, segs=rev_segs,
            genome_seq=rsqgl_res.genome_seq[::-1],
            raw_signal=rsqgl_res.raw_signal[::-1])

    return rsqgl_res

def map_read(a, faidx, seq_samp_type, std_ref, ref2len):
    """Get resquiggle result with read alignement info"""
    seq_data = tombo_helper.sequenceData(seq=a.seq, id=a.qname, mean_q_score=np.mean(a.query_qualities))
    # get chrom, start and end
    chrm, ref_start, ref_end = a.reference_name, a.reference_start, a.reference_end
    # store strand & number of clipped bases relative to read sequence
    if a.is_reverse:
        strand = "-"
        num_start_clipped_bases = len(seq_data.seq) - a.qend
        num_end_clipped_bases = a.qstart
    else:
        strand = "+"
        num_start_clipped_bases = a.qstart
        num_end_clipped_bases = len(seq_data.seq) - a.qend
    
    # 'ID', 'Subgroup', 'ClipStart', 'ClipEnd', 'Insertions', 'Deletions', 'Matches', 'Mismatches'
    align_info = tombo_helper.alignInfo(seq_data.id, "", num_start_clipped_bases, num_end_clipped_bases,
                                        0, 0, a.alen, 0) # this isn't used anywhere, so just don't bother computing it!
    
    # extract genome sequence from mappy aligner
    # expand sequence to get model levels for all sites (need to handle new
    # sequence coordinates downstream)
    start_skip = 0
    # get exonic blocks
    blocks = get_exonic_blocks(a)
    align_info.blocks = deepcopy(blocks)
    dnstrm_bases = std_ref.kmer_width - std_ref.central_pos - 1 
    if ((seq_samp_type.name == RNA_SAMP_TYPE and strand == '+') or
        (seq_samp_type.name == DNA_SAMP_TYPE and strand == '-' and USE_START_CLIP_BASES) or
        (seq_samp_type.name == DNA_SAMP_TYPE and strand == '+' and not USE_START_CLIP_BASES)):
        if ref_start < std_ref.central_pos: 
            start_skip = std_ref.central_pos-ref_start
            ref_start = std_ref.central_pos
        ref_seq_start = ref_start - std_ref.central_pos
        ref_seq_end = ref_end + dnstrm_bases
    else:
        if ref_start < dnstrm_bases: 
            start_skip = dnstrm_bases-ref_start
            ref_start = dnstrm_bases
        ref_seq_start = ref_start - dnstrm_bases
        ref_seq_end = ref_end + std_ref.central_pos
    # update blocks start & end with kmer specific shifts - this sequence won't be saved! 
    blocks[0][0] = ref_seq_start
    blocks[-1][1] = ref_seq_end
    # get exonic sequence
    genome_seq = "".join([faidx.fetch(chrm, s, e) for s, e in blocks])
    # get missing bases in the end
    end_skip = 0 if blocks[-1][1]<=ref2len[chrm] else blocks[-1][1]-ref2len[chrm]
    # enlarge genome seq by missing bits from ends with (random!) bases - As for now
    if start_skip or end_skip:
        genome_seq = "A"*start_skip + genome_seq + "A"*end_skip
    if strand == '-':
        genome_seq = tombo_helper.rev_comp(genome_seq)
    # store enlarged genome for P-value calculation, so no trimming needed later :)
    genome_seq = genome_seq.upper() #.upper() is important to correctly process soft-masked sequences
    align_info.refseq = genome_seq.upper() # res.genome_seq is altered during find_adaptive_assignment
    genome_loc = tombo_helper.genomeLocation(ref_start, strand, chrm)
    return tombo_helper.resquiggleResults(align_info, genome_loc, genome_seq, seq_data.mean_q_score)
   
def get_exonic_blocks(a):
    """Return exonic blocks this is start-end reference-based 
    for consecutive exons covered by given read.
 
    Note, those are not necesarily exact exons, just exons infered from read alignment. 
    """
    blocks = []
    s = e = a.pos
    # iter read blocks
    for code, bases in a.cigar:
        # count blocks that alter reference positions (ignore ie insertions [1])
        if code in (0, 2, 7, 8): e += bases
        # exclude introns - those are reported as reference-padded alignment part
        elif code == 3:
            blocks.append([s, e])
            s = e + bases
            e = s
    # store exon after last intron (or entire transcript if no introns)
    blocks.append([s, e])
    return blocks

def resquiggle_reads(multifast5_fn, aligner, ref, seq_samp_type, std_ref, rsqgl_params, 
                     outlier_thresh=OUTLIER_THRESH, max_scaling_iters=3, max_per_ref=0, 
                     valid_bases=set(list('ACGT'))):
    ref2c = {}
    # process reads from multi fast5
    faidx = pysam.FastaFile(ref)
    ref2len = {r: l for r, l in zip(faidx.references, faidx.lengths)}#; ref2len
    f5file = get_fast5_file(multifast5_fn, mode="r")
    for a in aligner:
        # process only given number of reads per reference
        if max_per_ref:
            contig = a.reference_name #map_results.genome_loc.Chrom
            if contig in ref2c:
                if ref2c[contig]>=max_per_ref: continue
            else: ref2c[contig] = 0
        # skip reads without alignment or secondary/qcfails
        if a.is_unmapped or a.is_secondary or a.is_qcfail:
            yield None, "No alignment" if a.is_unmapped else "Secondary alignment"
            continue
        # get alignment data
        map_results = map_read(a, faidx, seq_samp_type, std_ref, ref2len)
        # make sure only ACGT chars in reference!
        if set(map_results.genome_seq).difference(valid_bases):
            yield None, "Non-ACGT sequence" # instead maybe just replace by random char?
            continue
        # extract data from FAST5
        read = f5file.get_read(a.qname) #read_id)
        all_raw_signal = read.get_raw_data(scale=False)
        map_results = map_results._replace(raw_signal=all_raw_signal)
        try:
            # this causes sometimes TomboError: Read event to sequence alignment extends beyond bandwidth
            map_results = adjust_map_res(map_results, seq_samp_type, rsqgl_params)
            rsqgl_res = resquiggle.resquiggle_read(map_results, std_ref, rsqgl_params, outlier_thresh)
            n_iters = 1
            while n_iters < max_scaling_iters and rsqgl_res.norm_params_changed:
                rsqgl_res = resquiggle.resquiggle_read(map_results._replace(scale_values=rsqgl_res.scale_values),
                                                       std_ref, rsqgl_params, outlier_thresh)
                n_iters += 1
        except Exception as inst:
            yield None, str(inst)
            continue
        rsqgl_res = adjust_rsqgl_res(rsqgl_res, all_raw_signal, seq_samp_type)
        # add alignment and read as those are needed later
        rsqgl_res.a, rsqgl_res.read = a, read
        # update ref counter 
        if ref2c: ref2c[contig] += 1
        yield rsqgl_res, ""

def get_norm_mean(raw, segs): 
    """Return raw signal means for given segments."""
    return np.array([raw[segs[i]:segs[i+1]].mean() for i in range(len(segs)-1)])

def get_trace_for_reference_bases(a, read, rna, func=np.mean):
    """Return reference-aligned trace for tr (ref base), tA, tC, tG, tT"""
    def get_bidx_fwd(b): return base2idx[b] 
    def get_bidx_rev(b): return base2idx[base2complement[b]] 
    # trace for reference bases
    tr = np.zeros(a.reference_length, dtype="uint8")
    # trace and move data from read
    bcgrp = read.get_latest_analysis("Basecall_1D")
    trace = read.get_analysis_dataset(bcgrp, "BaseCalled_template/Trace")
    if trace is None:
        logger("[ERROR] Trace table is missing in Fast5 file! Basecall Fast5 files again using --fast5_out option. ")
        return tr
    move = read.get_analysis_dataset(bcgrp, "BaseCalled_template/Move")
    move_pos = np.append(np.argwhere(move==1).flatten(), len(trace)) # add end of trace
    # combine flip & flop probabilities
    ## here we get sum of flip & flop. maybe get just one? but flop is usually lower...
    trace[:, :len(bases)] += trace[:, len(bases):]
    trace = trace[:, :len(bases)]
    # here we need to remember that DNA 5'>3', but RNA 3'>5'
    # plus the strand matters
    if a.is_reverse: # for REV alg
        get_bidx = get_bidx_rev # take complement base
        if not rna: move_pos = move_pos[::-1] # reverse move_pos for DNA
    else: # for FWD alg
        get_bidx = get_bidx_fwd # take base
        if rna: move_pos = move_pos[::-1] # reverse move_pos for RNA
    # process aligned bases - that's quite elegant, right? :P
    ## with_seq require MD tags: in minimap2 use --MD and -Y (soft-clip supplementary)
    for qi, ri, b in a.get_aligned_pairs(with_seq=True, matches_only=True): 
        # get start & end in trace-space
        s, e = move_pos[qi:qi+2]
        if s>e: s, e = e, s # fix s, e for reversed move_pos
        tr[ri-a.reference_start] = func(trace[s:e, get_bidx(b)], axis=0)
    return tr

def get_trace_for_all_bases(a, read, rna, func=np.mean):
    """Return reference-aligned trace for tr (ref base), tA, tC, tG, tT"""
    def get_bidx_fwd(b): return base2idx[b] 
    def get_bidx_rev(b): return base2idx[base2complement[b]] 
    # trace for reference bases
    tr = np.zeros((a.reference_length,5), dtype="uint8") # one column per base + canonical col
    # trace and move data from read
    bcgrp = read.get_latest_analysis("Basecall_1D")
    trace = read.get_analysis_dataset(bcgrp, "BaseCalled_template/Trace")
    if trace is None:
        logger("[ERROR] Trace table is missing in Fast5 file! Basecall Fast5 files again using --fast5_out option. ")
        return tr
    move = read.get_analysis_dataset(bcgrp, "BaseCalled_template/Move")
    move_pos = np.append(np.argwhere(move==1).flatten(), len(trace)) # add end of trace
    # combine flip & flop probabilities
    ## here we get sum of flip & flop. maybe get just one? but flop is usually lower...
    trace[:, :len(bases)] += trace[:, len(bases):]
    trace = trace[:, :len(bases)]
    # here we need to remember that DNA 5'>3', but RNA 3'>5'
    # plus the strand matters
    if a.is_reverse: # for REV alg
        get_bidx = get_bidx_rev # take complement base
        if not rna: move_pos = move_pos[::-1] # reverse move_pos for DNA
    else: # for FWD alg
        get_bidx = get_bidx_fwd # take base
        if rna: move_pos = move_pos[::-1] # reverse move_pos for RNA
    # process aligned bases - that's quite elegant, right? :P
    ## with_seq require MD tags: in minimap2 use --MD and -Y (soft-clip supplementary)
    for qi, ri, b in a.get_aligned_pairs(with_seq=True, matches_only=True): 
        # get start & end in trace-space
        s, e = move_pos[qi:qi+2]
        if s>e: s, e = e, s # fix s, e for reversed move_pos
        tr[ri-a.reference_start,0] = func(trace[s:e, 0], axis=0)
        tr[ri-a.reference_start,1] = func(trace[s:e, 1], axis=0)
        tr[ri-a.reference_start,2] = func(trace[s:e, 2], axis=0)
        tr[ri-a.reference_start,3] = func(trace[s:e, 3], axis=0)
        tr[ri-a.reference_start,4] = func(trace[s:e, get_bidx(b)], axis=0)
    return tr

def process_fast5(fast5, ref, rna=True, sensitive=False):
    """Process individual Fast5 files"""
    outfn = "%s.bam"%fast5 #.d2r

    # uncomment if you don't wish to recompute previously computed bam files
    # if os.path.isfile(outfn): return outfn
    faidx = pysam.FastaFile(ref)
    ref2len = {r: l for r, l in zip(faidx.references, faidx.lengths)}
    # load model & its parameters
    if rna:
        seq_samp_type = tombo_helper.seqSampleType('RNA', True)
        rsqgl_params = tombo_stats.load_resquiggle_parameters(seq_samp_type)
        std_ref = tombo_stats.TomboModel(seq_samp_type=seq_samp_type)
        spliced = True
    else:
        seq_samp_type = tombo_helper.seqSampleType('DNA', False)
        rsqgl_params = tombo_stats.load_resquiggle_parameters(seq_samp_type)
        spliced = False
        std_ref = tombo_stats.TomboModel(seq_samp_type=seq_samp_type)
    # get resquiggle parameters
    i, errors = 0, {} 
    # prep aligner, signal model and parameters
    aligner = minimap2_proc(ref, fast5, sensitive=sensitive, spliced=spliced)
    sam = pysam.AlignmentFile(aligner.stdout)
    # open unsorted bam for saving alignements with features
    tmp = tempfile.NamedTemporaryFile(delete=False); tmp.close()
    bam_unsorted = pysam.AlignmentFile(tmp.name, "wb", header=sam.header)

    for i, (res, err) in enumerate(resquiggle_reads(fast5, sam, ref, seq_samp_type, std_ref, rsqgl_params), 1):
        #if i>200: break
        if not i%100: sys.stderr.write(" %s - %s reads skipped: %s \r"%(i, sum(errors.values()), str(errors)))
        if not res:
            if err not in errors: errors[err] = 1
            else: errors[err] += 1
            continue
        # get pysam alignment object & exonic blocks
        a, blocks = res.a, res.align_info.blocks

        # get signal intensity means
        si = get_norm_mean(res.raw_signal, res.segs)
        # catch problems - here exonic seq will have different length
        if len(si)!=sum([e-s for s, e in blocks]): #a.reference_length:
            region = "%s:%s-%s"%(a.reference_name, a.reference_start, a.reference_end)
            print(a.qname, region, sam.lengths[a.reference_id], a.reference_length, len(si), blocks)
        # get dwell times capped at 255 to fit uint8 (1 byte per base)
        dt = res.segs[1:]-res.segs[:-1]
        dt[dt>255] = 255
        # get reference-aligned base probabilities: tr (ref base)
        tr = get_trace_for_all_bases(a, res.read, rna) # trA, trC, trG, trT, (canonical) tr
        if a.is_reverse: si, dt = si[::-1], dt[::-1]
        # and finally set tags matching refseq
        ## but if alignment reaches seq end the end signal/probs will be wrong!
        ## same at exon-intron boundaries
        a.set_tag("bs", array("i", np.array(blocks).flatten()))
        a.set_tag("si", array("f", si))
        a.set_tag("dt", array("B", dt))
        # tr correspond to reference base
        # get exonic tr
        exonic_pos = np.concatenate([np.arange(s, e) for s, e in blocks])
        tr = tr[exonic_pos-a.pos]

        a.set_tag("tA", array("B", tr[:,0]))
        a.set_tag("tC", array("B", tr[:,1]))
        a.set_tag("tG", array("B", tr[:,2]))
        a.set_tag("tT", array("B", tr[:,3]))
        a.set_tag("tr", array("B", tr[:,4]))

        # add quality scores
        a.set_tag("QQ", array("B", a.query_qualities))

        # read id
        a.set_tag('ID', a.qname)

        # store read alignment with additional info
        bam_unsorted.write(a)

    # close tmp, sort, index & clean-up
    bam_unsorted.close()
    pysam.sort("-o", outfn, tmp.name)
    pysam.index(outfn)
    os.unlink(tmp.name)
    # write error report
    with open('%s.json'%outfn, 'w') as f:
        errors["Alignements"] = i # store number of alignements
        f.write(json.dumps(errors)) #
    return outfn

def mod_encode(fnames, fasta, threads=1, rna=True, sensitive=False, mem=1):
    """Process multiple directories from Fast5 files"""
    # no need to have more threads than input directories ;) 
    if threads > len(fnames):
        threads = len(fnames)
    # use pool if more than 1 thread, otherwise just itertools
    if threads>1: p = Pool(threads, maxtasksperchild=1)
    else: p = itertools
    # get arguments for func
    args = [(fn, fasta, rna, sensitive) for fn in fnames]# if not os.path.isfile("%s.bam"%fn)]
    # return list of outputs
    return list(p.starmap(process_fast5, args))    

def memory_usage(childrenmem=True, div=1024.):
    """Return memory usage in MB including children processes"""
    mem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / div
    if childrenmem:
        mem += resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss / div
    return mem
        
def logger(info, add_timestamp=1, add_memory=1, out=sys.stderr):
    """Report nicely formatted stream to stderr"""
    info = info.rstrip('\n')
    memory = timestamp = ""
    if add_timestamp:
        timestamp = "[%s]"%str(datetime.now()).split(".")[0] 
    if add_memory:
        memory = " [mem: %5.0f MB]"%memory_usage()
    out.write("%s %s%s\n"%(timestamp, info, memory))

def main():
    import argparse
    usage   = "%(prog)s -v" #usage=usage, 
    parser  = argparse.ArgumentParser(description=desc, epilog=epilog, \
                                      formatter_class=argparse.RawTextHelpFormatter)
  
    parser.add_argument('--version', action='version', version=VERSION)   
    parser.add_argument("-v", "--verbose", action="store_true", help="verbose")    
    parser.add_argument("-i", "--input", nargs="+", help="input Fast5 file(s)")
    parser.add_argument("--rna", action='store_true', help="project is RNA sequencing [DNA]")
    parser.add_argument("-f", "--fasta", required=1, help="reference FASTA file")
    parser.add_argument("-t", "--threads", default=1, type=int, help="number of cores to use [%(default)s]")
    parser.add_argument("-s", "--sensitive", action='store_true', help="use sensitive alignment")

    o = parser.parse_args()
    if o.verbose: 
        sys.stderr.write("Options: %s\n"%str(o))

    # encode tombo output into BAM files
    logger("Processing %s file(s)..."%len(o.input))
    bamfiles = mod_encode(o.input, o.fasta, o.threads, o.rna, o.sensitive)
        
if __name__=='__main__': 
    t0 = datetime.now()
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\nCtrl-C pressed!      \n")
    #except IOError as e:
    #    sys.stderr.write("I/O error({0}): {1}\n".format(e.errno, e.strerror))
    dt = datetime.now()-t0
    sys.stderr.write("#Time elapsed: %s    \n"%dt)

