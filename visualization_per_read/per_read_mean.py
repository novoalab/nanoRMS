import sys
import pandas as pd
infile = sys.argv[1]
inp=pd.read_csv(infile,sep='\t') #this step will take ages if the file is huge!
inp = inp [inp['model_kmer'] != 'NNNNN'] #Remove NNNNN values
grouped_multiple_mean_inp = inp.groupby(['contig', 'position','reference_kmer', "read_index"]).agg({'event_level_mean':['mean']}) #Collapse multiple observations from the same read
grouped_multiple_mean_inp = grouped_multiple_mean_inp.reset_index()
grouped_multiple_mean_inp.columns =  grouped_multiple_mean_inp.columns.droplevel(-1)
grouped_multiple_mean_inp.to_csv('{}_processed_perpos_mean.csv'.format(infile), sep='\t', index = False) #Export the file