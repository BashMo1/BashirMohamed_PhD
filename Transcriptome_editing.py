#!/usr/bin/env python

#Imports
import re
from Bio import SeqIO

#Functions
def read_in_fasta(afasta):
	'''Reads in a fasta file to a dictionary'''
	fasta_dict = {}
	fasta_sequences = SeqIO.parse(open(afasta),'fasta')
	for fasta in fasta_sequences:
		fasta_dict[fasta.id] = str(fasta.seq).upper()
	return fasta_dict
    
def filter_fasta(fasta_dict):
    '''filters fasta and returns a dict for each region that passes the filtering'''
        
    filtered_dict = {}
    fp_dict = {}
    cds_dict = {}
    tp_dict = {}
    
    for k, v in fasta_dict.items():
        total_length = int(re.findall(r"\|\d+\|",k)[0].strip('|'))
        transcript = k.split('|')[0]
        fp_search = re.findall(r"\|\UTR5:\d+\-\d+\|",k)
        if fp_search != []:
            fp_len = int(fp_search[0].strip('|').split('-')[1])
            
            tp_search = re.findall(r"\|\UTR3:\d+\-\d+\|",k)
            if tp_search != []:
                tp_range = tp_search[0].split(':')[1].strip('|')
                tp_len = int(tp_range.split('-')[1]) - int(tp_range.split('-')[0]) + 1
                
                cds_search = re.findall(r"\|\CDS:\d+\-\d+\|",k)
                cds_range = cds_search[0].split(':')[1].strip('|')
                cds_len = int(cds_range.split('-')[1]) - int(cds_range.split('-')[0]) + 1
                
                if total_length == fp_len + cds_len + tp_len:
                    if cds_len % 3 == 0:
                        fp_seq = v[:fp_len]
                        cds_seq = v[fp_len:-tp_len]
                        tp_seq = v[-tp_len:]
                        if len(fp_seq) + len(cds_seq) + len(tp_seq) == total_length:
                            if cds_seq[:3] == "ATG":
                                if cds_seq[-3:] == "TAA" or cds_seq[-3:] == "TGA" or cds_seq[-3:] == "TAG":
                                    filtered_dict[k] = v
                                    fp_dict[transcript] = fp_seq
                                    cds_dict[transcript] = cds_seq
                                    tp_dict[transcript] = tp_seq
    return [filtered_dict,fp_dict,cds_dict,tp_dict]

def write_fasta(dictionary, outfyle, LW=80):
	'''takes a dictionary and writes fasta'''
	with open(outfyle, 'w') as g:
		for k,v in dictionary.items():
			g.write('>' + k + '\n')
			for i in range(0, len(v), LW):
				g.write(v[i:i+LW] + '\n')

#Main Function
def main():
    master_fasta = read_in_fasta("gencode.v34.pc_transcripts.fa")
    filtered_dict,fp_dict,cds_dict,tp_dict = filter_fasta(master_fasta)
    write_fasta(filtered_dict,'gencode.v34.pc_transcripts_filtered.fa')
    write_fasta(fp_dict,'gencode.v34.pc_transcripts_5UTRs.fa')
    write_fasta(cds_dict,'gencode.v34.pc_transcripts_CDSs.fa')
    write_fasta(tp_dict,'gencode.v34.pc_transcripts_3UTRs.fa')

if __name__ == '__main__': 
    main()
 
