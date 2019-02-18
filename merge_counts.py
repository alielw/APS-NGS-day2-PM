#!/usr/bin/python2.6
# -*- coding: utf-8 -*-
#Take htseq-count output files, extracts read counts, prints into a new file. Input files are called
#sample.htseq
#==============================================================================
import argparse
import sys
import os
from collections import defaultdict
#==============================================================================
#Command line options==========================================================
#==============================================================================
parser = argparse.ArgumentParser()
parser.add_argument("htseq", type=str,
					help="Folder of htseq-count files")
# This checks if the user supplied any arguments. If not, help is printed.
if len(sys.argv) == 1:
	parser.print_help()
	sys.exit(1)
# Shorten the access to command line arguments.
args = parser.parse_args()
#==============================================================================
#Functions=====================================================================
#==============================================================================

def get_expr(file,dict_FPKM, dict_count):
	with open(file, "r") as infile:
		for line in infile:
			line = line.rstrip()
			if line.startswith("gene_id"):
				pass
			else:
				gene = line.split("\t")[0]
				FPKM = line.split("\t")[6]
				count = line.split("\t")[4]
				dict_FPKM[gene].append(float(FPKM))
				dict_count[gene].append(count)
	return dict_FPKM, dict_count

def list_folder(infolder):
    '''Returns a list of all files in a folder with the full path'''
    return [os.path.join(infolder, f) for f in os.listdir(infolder) if f.endswith(".htseq")]

#==============================================================================
#Main==========================================================================
#==============================================================================
def main():

	folder = list_folder(args.htseq)
	print folder

	#extract read counts

	header = []
	count_dict = defaultdict(list)
	for file in folder:
		sample = os.path.basename(file).split(".htseq")[0]
		header.append(sample)
		print sample
		with open(file, "r") as infile:
			for line in infile:
				line = line.rstrip()
				if not line.startswith("__"):
					gene = line.split("\t")[0]
					count = line.split("\t")[1]
					count_dict[gene].append(count)
	print "Number of transcripts =", len(count_dict)

	#check expression data for each sample

	for gene in count_dict:
		if len(count_dict[gene]) == len(header):
			pass
		else:
			print "ERROR - missing expression data"
			sys.exit()

	#print read counts into one file

	outfile = os.path.dirname(args.htseq)+"/"+"htseq.merged.txt"
	print "Printing to ", outfile

	with open(outfile, "w") as out:
		out.write("gene_id")
		for i in header:
			out.write("\t")
			out.write(i)
		out.write("\n")
		for gene in count_dict:
			out.write(gene)
			for i in count_dict[gene]:
				out.write("\t")
				out.write(i)
			out.write("\n")

if __name__ == '__main__':
	main()