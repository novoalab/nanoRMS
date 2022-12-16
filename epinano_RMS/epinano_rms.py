#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys,os,re,io,pysam
import shutil, fileinput
import glob, itertools
import subprocess
import argparse
import multiprocessing as mp
from multiprocessing import Process, Manager
from functools import partial
from sys import __stdout__
import dask
import dask.dataframe as dd
import pandas as pd
from collections import defaultdict
from collections import OrderedDict
import numpy as np

#~~~~~~~~~~~~~~~~~~~~ private function ~~~~~~~~
# func1 subprocess call linux cmmands
def touch(fname):
    if os.path.exists(fname):
        os.utime(fname, None)
    else:
        open(fname, 'a').close()

def openfile(f):
	if f.endswith ('.gz'):
		fh = gzip.open (f,'rt')
	elif f.endswith ('bz') or f.endswith ('bz2'):
		fh = bz2.open(f,'rt')
	else:
		fh = open(f,'rt')
	return fh

def spot_empty_tsv (tsv):
        ary = []
        cnt = 0
        with open (tsv,'r')  as fh:
                for l in fh:
                        if cnt <2:
                                ary.append (l)
                        else:
                                break
                        cnt += 1
        return True if len (ary)>1 else False


def split_tsv_for_per_site_var_freq(tsv,folder, q, number_threads,  num_reads_per_chunk=4000):
	'''
	'''
	head = next(tsv)
	firstline = next (tsv)
	current_rd = firstline.split()[0]
	rd_cnt = 1
	idx = 0
	out_fn = "{}/CHUNK_{}.txt".format(folder, idx)
	out_fh = open (out_fn, 'w')
	#chunk_out = [] # open ("CHUNK_{}.txt".format(idx),'w')
	#chunk_out.append(firstline)
	print (firstline.rstrip(), file=out_fh)
	try:
		for line in tsv:
			rd = line.split()[0]
			if current_rd != rd:
				rd_cnt += 1
				current_rd = rd
				if ((rd_cnt-1) % num_reads_per_chunk == 0 and rd_cnt >= num_reads_per_chunk):
					q.put ((idx, out_fn)) #.close()
					idx += 1
					out_fn = "{}/CHUNK_{}.txt".format(folder,idx)
					out_fh = open (out_fn, 'w')
			print (line.rstrip(), file=out_fh)
		out_fh.close()
		q.put((idx, out_fn))
	except:
		raise
		sys.stderr.write("split tsv file on reads failed\n")
	finally:
		for _ in range(number_threads):
			q.put(None)

def  proc_small_freq (small_freq_fn):
    df = pd.read_csv (small_freq_fn)
    df['pos'] = df['pos'].astype(str)
    df['index'] = df[['#Ref','pos','base','strand']].apply (lambda x: "-EPIJN-".join(x),axis=1)
    df.drop(['#Ref','pos','base','strand'], axis=1, inplace=True)
    df.set_index(['index'], inplace=True)
    df['qual'] = df['qual'].replace(r':{2,}',':',regex=True)
    df['qual'] = df['qual'].replace(r':$','',regex=True)
    df['bases'] = df['bases'].replace(r':{2,}',':', regex=True)
    df['bases'] = df['bases'].replace(r':$','', regex=True)
    df[['_A_', '_C_', '_G_', '_T_']] = df['bases'].str.split(pat=':',  expand=True)
    df.drop (['bases'],axis=1, inplace=True)
    return df

def file_exist (file):
	return os.path.exists (file)

def _rm (file):
	os.remove (file)

def stdin_stdout_gen (stdin_stdout):
	'''
	generator for subprocess popen stdout
	'''
	for l in stdin_stdout:
		if isinstance (l,bytes):
			yield (l.decode('utf-8'))
		else:
			yield l

def print_from_stdout (stdout_lst, outputfh):
	for i, o in enumerate (stdout_lst):
		for l in o:
			if l.decode().startswith ('#'):
				if i >1 :
					continue
			outputfh.write(l.decode())
#~~~~~~~

def java_bam_to_tsv (bam_file,  reference_file, sam2tsv):
	'''
	type: reference types,i.e., trans or genome
	'''

	awk_forward_strand = """ awk '{if (/^#/) print $0"\tSTARAND"; else print $0"\t+"}' """
	awk_reverse_strand = """ awk '{if (/^#/) print $0"\tSTARAND"; else print $0"\t-"}' """
	cmds = []
	cmd1 = (f"samtools view -h -F 3860 {bam_file} | java -jar  {sam2tsv} -r {reference_file} "
					f"| {awk_forward_strand} ")
	cmd2 = (f"samtools view -h -f 16 -F 3844 {bam_file} | java -jar  {sam2tsv} -r {reference_file} "
					f" | {awk_reverse_strand}")
	cmds = [cmd1,cmd2]
	return cmds
# data frame

def tsv_to_freq_multiprocessing_with_manager (tsv_reads_chunk_q, out_dir):
	'''
	mutliprocessing
	produced with sam2tsv.jar with strand information added
	read read-flags reference       read-pos        read-base       read-qual       ref-pos ref-base                cigar-op                strand
	a3194184-d809-42dc-9fa1-dfb497d2ed6a    0       cc6m_2244_T7_ecorv      0       C       #       438     G       S       +
	'''
	for idx, tsv_small_chunk_fn in iter (tsv_reads_chunk_q.get, None):
		filename = "{}/small_{}.freq".format(out_dir, idx)
		outh = open (filename,'w')
		mis = defaultdict(int) # mismatches
		mat = defaultdict (int) #matches
		ins = defaultdict(int) # insertions
		dele = defaultdict(int) # deletions
		cov = OrderedDict ()  # coverage
		ins_q = defaultdict(list)
		aln_mem = []  #read, ref, refpos; only store last entry not matching insertion
		pos = defaultdict(list) # reference positions
		base = {} # ref base
		qual = defaultdict(list)
		read_bases = defaultdict (dict)
        #READ_NAME     FLAG    CHROM   READ_POS  BASE   QUAL  REF_POS REF  OP   STRAND
        #read read-flags        reference       read-pos        read-base       read-qual       ref-pos ref-base                cigar-op                strand
		tsv_small_chunk = open (tsv_small_chunk_fn,'r')
		for line in tsv_small_chunk:
			if line.startswith ('#'):
				continue
			ary = line.rstrip().split()
			if ary[-2] in ['M','m']:
				k = (ary[2], int (ary[-4]), ary[-1]) #
				cov[k] = cov.get(k,0) + 1
				aln_mem = []
				aln_mem.append((ary[0],ary[2],int(ary[-4]), ary[-1]))
				qual[k].append (ord(ary[-5])-33)
				base[k] = ary[-3].upper()
				read_bases[k][ary[4]] = read_bases[k].get(ary[4], 0) + 1
				if (ary[-3] != ary[4]):
					mis[k] += 1
				else:
					mat[k] += 1
			if ary[-2] == 'D':
				k = (ary[2], int(ary[-4]), ary[-1])
				cov[k] = cov.get(k,0) + 1
				aln_mem = []
				aln_mem.append((ary[0],ary[2],int(ary[-4]), ary[-1]))
				base[k] = ary[-3].upper()
				dele[k] = dele.get(k,0) + 1
			if ary[-2] == 'I':
				last_k = aln_mem[-1][1],aln_mem[-1][2],aln_mem[-1][3] # last alignment with match/mismatch/del
				next_k = (ary[2], last_k[1] + 1,last_k[2])
				if last_k[0] != ary[2]:
					pass
				ins_k_up = (ary[0], ary[2], last_k[1],last_k[2])
				ins_k_down = (ary[0], ary[2], last_k[1] + 1,last_k[2])
				if (ins_k_down) not in ins_q:
					ins[next_k] = ins.get(next_k,0) + 1
					ins_q[ins_k_down].append(ord(ary[-5])-33)
				if (ins_k_up) not in ins_q:
					ins[last_k] = ins.get(last_k,0) + 1
					ins_q[ins_k_up].append(ord(ary[-5])-33)
		header = '#Ref,pos,base,cov,mat,mis,ins,del,qual,strand,bases\n'
		outh.write(header)
		os.remove(tsv_small_chunk_fn)

		for k in cov.keys():
			depth = cov.get (k,0)
			Mis = mis.get (k,0)
			Mat = mat.get (k,0)
			Del = dele.get (k,0)
			q_lst = qual.get (k,[0])
			try:
				q_lst = ':'.join (map (str, q_lst))+':'  # dataframe sum
				num_ins = ins.get (k,0)
				bases_counts = "0:0:0:0:"
				if k in read_bases:
					bases_counts = ":".join ([str(read_bases[k].get(l,0)) for l in 'ACGT'])
				inf = "{},{},{},{},{},{},{},{},{},{},{}:\n".format (k[0], k[1], base[k], depth, Mat, Mis, num_ins, Del, q_lst, k[2], bases_counts)
				outh.write (inf)
			except:
				sys.stderr.write ("file {} {} does not work\n".format (tsv,k))


def df_is_not_empty(df):
	'''
	 input df is a df filtred on reference id
	 if is is empty: next (df.iterrows()) does not work
	 otherwise it returns a row of df
	'''
	try:
		next (df.iterrows())
		return True
	except:
		return False


def _tsv_gen_ (bam_fn, ref_fn, sam2tsv_jar):
	cmds = java_bam_to_tsv (bam_fn, ref_fn, sam2tsv_jar) #, args.type)
	cmd1 = subprocess.Popen ((cmds[0]), stdout=subprocess.PIPE, stderr = subprocess.PIPE,shell=True)
	cmd2 = subprocess.Popen ((cmds[1]), stdout=subprocess.PIPE, stderr = subprocess.PIPE,shell=True)
	returncode1 = cmd1.returncode
	returncode2 = cmd2.returncode
	if any ([returncode1, returncode2] ):
		res1 = cmd1.communicate()
		res2 = cmd2.communicate()
		print (res1[1], res2[1], file=sys.stderr)
		exit()
	return itertools.chain (stdin_stdout_gen (cmd1.stdout), stdin_stdout_gen (cmd2.stdout))




#~~~~~~~~~~~~~~~~~~~~~~~ main () ~~~~~~~~~~~~~~~~~~~~~~~
def main ():
	parser = argparse.ArgumentParser()
	parser.add_argument ('-R','--reference', type=str, required=True, help='samtools faidx indexed reference file')
	parser.add_argument ('-b', '--bam', type=str, required=True, help='bam file; if given; no need to offer reads file; mapping will be skipped')
	parser.add_argument ('-s', '--sam2tsv',type=str, required=True, default='',help='/path/to/sam2tsv.jar; needed unless a sam2tsv.jar produced file is already given')
	parser.add_argument ('-n', '--number_cpus', type=int, default=4,  help='number of CPUs')
	parser.add_argument ('-d', '--delete', action='store_true', help = 'delete intermediate files')
	args=parser.parse_args()

#~~~~~~~~~~~~~~~~~~~~~~~ prepare for analysis ~~~~~~~~~~~~~~
	prefix = ''
	if args.reference:
		if not file_exist (args.reference):
			sys.stderr.write (args.reference, 'does not exist')
			exit()
		dict_fn = args.reference + '.dict'
		if not file_exist (dict_fn):
			sys.stderr.write (dict_fn, 'needs to be created using picard.jar CreateSequenceDictionary')
			exit()
		ref_faidx = args.reference +'.fai'
		if not file_exist (ref_faidx):
			sys.stderr.write (ref_faidx, 'needs to be created with samtools faidx')
			exit()
	if args.bam:
		bam_file = args.bam
		if not file_exist (bam_file):
			sys.stderr.write (bam_file+' does not exist; please double check!\n')
			exit()
		else:
			if not file_exist (args.sam2tsv):
				sys.stderr.write ("Please offer correctly path to sam2tsv.jar\n".format(args.sam2tsv))
				exit()
			if not os.path.exists (bam_file+'.bai'):
				print (bam_file)
				sys.stderr.write ('bam file not indexed!\nstarting indexing it ...')
				pysam.index (bam_file)
			if not args.reference :
				sys.stderr.write('requires reference file that was used for reads mapping\n')
			prefix = re.sub (r'.bam$', '', bam_file) #  bam_file.replace('.bam','')

#~~~~~~~~~~~~~~~~  SAM2TSV ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
################# funciton run commands ###########################
#~~~~~~~~~~~~~~~~ split tsv  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	tsv_gen  = _tsv_gen_(args.bam, args.reference, args.sam2tsv)
	tmp_dir = prefix + '.tmp_splitted_base_freq'
	
	progress_fn = ".{}.done_splitting".format(tmp_dir)
		
	if not os.path.exists(progress_fn):
		if  os.path.exists(tmp_dir):
			shutil.rmtree (tmp_dir)
			sys.stderr.write ("{} already exists, will overwrite it\n".format(tmp_dir))
		os.mkdir (tmp_dir)
		number_threads = args.number_cpus
		manager = Manager()
		q = manager.Queue(args.number_cpus)
#~~~~~~~~~~~~~~~~ compute per site variants frequecies ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#1 calculate variants frequency for each small splitted file
		processes = []
		ps = Process (target = split_tsv_for_per_site_var_freq, args = (tsv_gen, tmp_dir, q, number_threads, 2500))
		processes.append (ps)
		for _ in range(number_threads):
			ps = Process (target= tsv_to_freq_multiprocessing_with_manager, args = (q, tmp_dir))
			processes.append (ps)
		for ps in processes:
			ps.daemon = True
			ps.start()
		for ps in processes:
			ps.join()
		touch (".{}.done_splitting".format(tmp_dir))
#2 combine small files and produce varinats frequencies per ref-position
	small_freq_fns = [os.path.join (tmp_dir, f) for f in os.listdir(tmp_dir) if f.startswith('small_')]
	out = open (prefix + '.per.site.baseFreq.csv', 'w')
	print ('#Ref,pos,base,strand,cov,mean_q,median_q,std_q,mis,ins,del,ACGT', file=out)
	ddf_lst = [] 
	for f in small_freq_fns:
        	df = proc_small_freq (f)
	        ddf = dd.from_pandas(df, npartitions=2)
        	ddf_lst.append(ddf)
	ddf_cat = dd.concat (ddf_lst, axis=1)
	
	for r in  ddf_cat.iterrows ():
		index, var = r[0],r[1]
		var_df = pd.DataFrame (np.split(np.array(var), len(r[1])/10)) #10: cov, mat, mis, ins, del, qual, _A_, _C_, _G_, _T_
		var_df.columns = ['cov', 'mat', 'mis', 'ins', 'del', 'qual', '_A_', '_C_','_G_', '_T_']
		cov=var_df['cov'].sum()
		mat="{:.6f}".format(var_df['mat'].sum()/cov)
		mis="{:.6f}".format(var_df['mis'].sum()/cov)
		ins="{:.6f}".format(var_df['ins'].sum()/cov)
		dele = "{:.6f}".format(var_df['del'].sum()/cov)
		qual = var_df['qual'].dropna().sum() #apply(str).replace(np.nan,'',regex=True).sum()
		qual = np.array(qual.split(':')).astype(int)
		qmn, qme, qst = "{:.6f}".format(np.mean(qual)), "{:.6f}".format(np.median(qual)), "{:.6f}".format(np.std(qual))
		ACGTs = [var_df['_A_'].dropna().astype(int).sum(),  var_df['_C_'].dropna().astype(int).sum() ,		
		var_df['_G_'].dropna().astype(int).sum() , var_df['_T_'].dropna().astype(int).sum()]
		index = index.replace('-EPIJN-',',')
		out.write ("{},{},{},{},{},{},{},{},{}\n".format(index,cov,qmn,qme,qst,mis,ins,dele,":".join(map (str,ACGTs))))
# ~~~~~~~~~~~~~~~~~ delete intermediate files
	if args.delete:
		pool = mp.Pool(args.number_cpus)
		tmp_files = glob.glob("{}/small*".format(tmp_dir))
		pool.map(_rm,  tmp_files)
		shutil.rmtree(tmp_dir)
if __name__ == "__main__":
	main()
