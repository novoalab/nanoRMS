#!/usr/bin/env python 
# -*- coding: utf-8 -*- 

import sys,os,re,io
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



def pysam_bam_to_tsv (bam):
	bamfh = pysam.AlignmentFile(bam,'rb')
	out_tsv_fh = open (bam + '.tsv', 'w')
	header = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format("#READ_NAME","FLAG","CHROM","READ_POS","BASE", "QUAL","REF_POS","REF","OP",'STRAND')
	out_tsv_fh.write(header)
	for read in bamfh.fetch():
		o1, o2, o3  = read.query_name, read.flag, read.reference_name
		query_seq = read.query_sequence
		pairs = read.get_aligned_pairs(with_seq=True)
		pairs = clean_soft_hard_clippings(pairs)
		pairs = clean_soft_hard_clippings(pairs[::-1])
		pairs = pairs[::-1]
		strand = '-' if read.is_reverse else '+'
		op =''
		for p in pairs:
			try:
				o9 = variant_typing(p)
				op = o9
			except:
				sys.stderr.write ("{}\t{}\t{} is problematic\n".format (read.reference_name, read.query_name, p) )
				exit()
			if op in ['D']:
				o4, o5, o6 = '.', '.', '.'
				o7, o8 = p[1] + 1, p[2]
			elif op in ['I'] :
				o4,o5,o6 = p[0],query_seq[int(p[0])],read.query_qualities[p[0]]
				o7,o8 = '.','.'
			else:
				o4,o5,o6, o7, o8= p[0], query_seq[int(p[0])].upper(), read.query_qualities [p[0]], int (p[1]) + 1, p[2].upper()
			out_tsv_fh.write ("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(o1,o2,o3,o4,o5,o6,o7,o8,op,strand))
	out_tsv_fh.close()
	return (bam+'.tsv')


def split_tsv_for_per_site_var_freq(tsv, q, number_threads,  num_reads_per_chunk=4000):
	'''

	'''
	head = next(tsv)
	firstline = next (tsv)
	current_rd = firstline.split()[0]
	rd_cnt = 1
	idx = 0
	chunk_out = [] # open ("CHUNK_{}.txt".format(idx),'w')
	chunk_out.append(firstline)
	try:
		for line in tsv:
			rd = line.split()[0]
			if current_rd != rd:
				rd_cnt += 1
				current_rd = rd
				if ((rd_cnt-1) % num_reads_per_chunk == 0 and rd_cnt >= num_reads_per_chunk):
					q.put ((idx, chunk_out)) #.close()
					idx += 1
					chunk_out = [] #open ("CHUNK_{}.txt".format(idx),'w')
			chunk_out.append(line)
		q.put ((idx, chunk_out))
	except:
		raise
		sys.stderr.write("split tsv file on reads failed\n")
	finally:
		for _ in range(number_threads):
			q.put(None)

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

def java_bam_to_tsv (bam_file,  reference_file, sam2tsv, type):
	'''
	type: reference types,i.e., trans or genome 
	'''
	
	awk_forward_strand = """ awk '{if (/^#/) print $0"\tSTARAND"; else print $0"\t+"}' """
	awk_reverse_strand = """ awk '{if (/^#/) print $0"\tSTARAND"; else print $0"\t-"}' """
	cmds = []

	if type.lower().startswith ("t"):	
		cmd =  f"samtools view -h -F 3860 {bam_file} | java -jar  {sam2tsv} -r {reference_file} "\
			f" | {awk_forward_strand}"		
		cmds = [cmd]
	else:
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
	for idx, tsv_small_chunk in iter (tsv_reads_chunk_q.get, None):
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
	

def df_proc (df, outfn):
	'''
	input is a dataframe for either forward or reverse strand
	'''
	if not df_is_not_empty (df):
		print ("empty dataframe for {}".format(outfn), file=sys.stderr)
		return None
	outfh = open (outfn, 'w')
	header = "#Ref,pos,base,strand,cov,q_mean,q_median,q_std,mis,ins,del,ACGT_freq"
	print (header, file=outfh)
	gb = df.groupby(['#Ref','pos','base','strand']).agg({
               'cov':['sum'],
               'mis':['sum'],
               'ins':['sum'],
               'del':['sum'],
               'qual':['sum'],
		'bases':['sum']})
	for row in gb.itertuples():
		index = ",".join (map (str,row[0]))
		cov, mis, ins, _del, qual,bases = row[1:]
		mis = '%0.5f' % (mis/cov)
		ins = '%0.5f' % (ins/cov)
		_de = '%0.5f' % (_del/cov)
		q = np.array ([x for x in qual.split(':') if x ]).astype(int)  
		qmn,qme,qst = '%0.5f' % np.mean(q), '%0.5f' % np.median(q), '%0.5f' % np.std(q)	
		ACGT_freq = []
		acgt = ''
		try:
			bases = np.array(  [ ele for ele in  bases.split(':') if ele] ).astype(int)	
			acgt = np.array([0,0,0,0])
			for x in range (0, len(bases), 4):
				acgt = acgt + np.array (bases[x:x+4])
			outfh.write ("{},{},{},{},{},{},{},{},{}\n".format(index,cov,qmn,qme,qst,mis,ins,_de,":".join(map (str,acgt))))
		except:
			print ('warning:',row,file=sys.stderr)
			raise 





#~~~~~~~~~~~~~~~~~~~~~~~ main () ~~~~~~~~~~~~~~~~~~~~~~~
def main ():
	parser = argparse.ArgumentParser()
	required_args = parser.add_argument_group ('Required Arguments')
	required_args.add_argument ('-R','--reference', help='samtools faidx indexed reference file')
	required_args.add_argument ('-b', '--bam', type=str, help='bam file; if given; no need to offer reads file; mapping will be skipped')
	required_args.add_argument ('-s', '--sam2tsv',type=str, default='',help='/path/to/sam2tsv.jar; needed unless a sam2tsv.jar produced file is already given')
	parser.add_argument ('-f','--file', type=str, help='tsv file generated by sam2tsv.jar; if given, reads mapping and sam2tsv conversion will be skipped')
	parser.add_argument ('-n', '--number_cpus', type=int, default=4,  help='number of CPUs') 
	parser.add_argument ('-T', '--type', type=str, default="t" ,help="reference types, which is either g(enome) or t(ranscriptome);")
	args=parser.parse_args()

#~~~~~~~~~~~~~~~~~~~~~~~ prepare for analysis ~~~~~~~~~~~~~~ 
	tsv_gen = None  # generator 
	prefix = '' 
 
	def _tsv_gen ():
		if not args.file:
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
						sys.stderr.write ('bam file not indexed!\nstarting indexing it ...')
						os.system ('samtools index ' + bam_file + '.bai')
					if not args.reference :
						sys.stderr.write('requires reference file that was used for reads mapping\n')
					prefix = bam_file.replace('.bam','')
					cmds = java_bam_to_tsv (bam_file, args.reference, args.sam2tsv, args.type)
					if args.type[0].lower() == 't': #mapping to transcriptome; only one sam2tsv.jar command 
						cmd = subprocess.Popen ((cmds[0]), stdout=subprocess.PIPE, stderr=subprocess.PIPE,shell=True )
						returncode = cmd.returncode 
						if returncode:
							print (res[1], file=sys.stderr)
							exit()
						tsv_gen = stdin_stdout_gen (cmd.stdout)
					elif args.type[0].lower() == 'g': #mapping to genome; sam2tsv.jar caled twice for + and - strands 
						cmd1 = subprocess.Popen ((cmds[0]), stdout=subprocess.PIPE, stderr = subprocess.PIPE,shell=True)
						cmd2 = subprocess.Popen ((cmds[1]), stdout=subprocess.PIPE, stderr = subprocess.PIPE,shell=True)
						returncode1 = cmd1.returncode
						returncode2 = cmd2.returncode
						if any ([returncode1, returncode2] ):
							res1 = cmd1.communicate()
							res2 = cmd2.communicate() 
							print (res1[1], res2[1], file=sys.stderr)
							exit()
						tsv_gen = itertools.chain (stdin_stdout_gen (cmd1.stdout), stdin_stdout_gen (cmd2.stdout))
		else:
			if  args.file:
				tsv_file = args.file 
				prefix = tsv_file.replace ('.tsv','')
				if os.path.exists (args.file):
					fh = openfile (tsv_file)
					firstline = fh.readline()
					fh.close()
					if len (firstline.rstrip().split()) != 10:
						sys.stderr.write('tsv file is not in right format!')
						sys.stderr.write('tsv files should contain these columns {}\n'.format("#READ_NAME     FLAG    CHROM   READ_POS        BASE    QUAL    REF_POS REF     OP      STARAND"))
					sys.stderr.write (tsv_file + ' already exists; will skip reads mapping and sam2tsv conversion \n')			
					tsv_gen = openfile (tsv_file)
				else:
					sys.stderr.write (tsv_file + ' does not exist; please double check \n')
					exit()
		return tsv_gen, prefix 
#~~~~~~~~~~~~~~~~  SAM2TSV ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
################# funciton run commands ########################### 
#~~~~~~~~~~~~~~~~ split tsv  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	tsv_gen, prefix = _tsv_gen()

	tmp_dir = prefix + '.tmp_splitted'
	if  os.path.exists(tmp_dir):
		shutil.rmtree (tmp_dir)
		sys.stderr.write ("{} already exists, will overwrite it\n".format(tmp_dir))

	if not os.path.exists (tmp_dir):
		os.mkdir (tmp_dir)

		number_threads = args.number_cpus 

		manager = Manager()
		q = manager.Queue(args.number_cpus)
#~~~~~~~~~~~~~~~~ compute per site variants frequecies ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#1 calculate variants frequency for each small splitted file 
		processes = []

		ps = Process (target = split_tsv_for_per_site_var_freq, args = (tsv_gen, q, number_threads, 2500))
		processes.append (ps)

		for _ in range(number_threads):
			ps = Process (target= tsv_to_freq_multiprocessing_with_manager, args = (q, tmp_dir))
			processes.append (ps)
		for ps in processes:
			ps.daemon = True
			ps.start()
		for ps in processes:
			ps.join()

#2 combine small files and produce varinats frequencies per ref-position
#persite_var = prefix +'.per_site.var.csv'
	df = dd.read_csv ("{}/small_*freq".format(tmp_dir))
	df = df.compute()
	df_fr = df[df['strand']=="+"]
	out = prefix + '.per.site.fwd.csv'

	if args.type.lower() == 't' :
		df_proc (df_fr, out)
	elif args.type.lower() == 'g':
		df_rev = df[df['strand'] == "-"]
		out_rev = prefix + '.per.site.rev.csv'
		if args.number_cpus > 1:
			processes = [None, None]
			processes[0] = Process (target = df_proc, args = (df_fr,out))
			processes[1] = Process (target = df_proc, args = (df_rev,out_rev))
			for ps in processes:
				ps.aemon = True
				ps.start()		
			for ps in processes:
				ps.join()
		else:
			df_proc (df_fr, out)
			df_proc (df_rev, out_rev)
		
	#var_files = df_proc (tmp_dir, prefix, 2)


	if  os.path.exists(tmp_dir):
		pool = mp.Pool(args.number_cpus)
		tmp_files = glob.glob("{}/small*".format(tmp_dir))
		pool.map(_rm,  tmp_files)
		shutil.rmtree(tmp_dir)
#print ("per site variants frequencies table has been generated", file = sys.stderr)
#3 sliding window per site variants --> for making predicitons 
if __name__ == "__main__":
	main()
