#!/usr/bin/env python3

import os
import sys
import argparse
import numpy as np
import subprocess as sp
from Bio import SeqIO

def alleleFormatCallers(allelesBN,outDir):
	seq={}
	with open(allelesBN) as record:
		for record in SeqIO.parse(record,"fasta"):
			sequence = str(record.seq)
			sequenceUpdate = sequence.replace("-","")
			seq[record.id]=sequenceUpdate

	out_filename = os.path.join(outDir,'reference_alleles.fasta')
	with open(out_filename, "w") as output_file:
		for key in seq.keys():
			name = key.split("_")
			ids = name[0].replace("_","",1)
			refCheck = name[-1]
			if refCheck == '1':
				output_file.write(">"+ids+name[1]+"_1"+"\n"+str(seq[key])+"\n")




def main():

	parser = argparse.ArgumentParser()
	parser.add_argument("-bn",required=True, help="REQUIRED: BioNumerics alleleDB")
	parser.add_argument("-output",required=True, help="REQUIRED: Output directory for reformatted BioNumerics allele DB")


	args = parser.parse_args()
	allelesDB = args.bn
	outDir = args.output
	alleleFormatCallers(allelesDB,outDir)

if __name__ == '__main__':
	main()
	