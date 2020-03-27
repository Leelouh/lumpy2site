#!/usr/bin/env python3

import sys,re,os,cgi,csv, getopt, argparse, time, math, shutil
from copy import deepcopy
from random import randrange, uniform
from numpy import *
import scipy as sp
from pandas import *
import argparse

parser = argparse.ArgumentParser(add_help=True)

parser.add_argument('-in', action='store', type=str, help="file from extract_bnd_with_evid.py", dest="in")
parser.add_argument('-out', action='store', type=str, help="output file with results", dest="out")
parser.add_argument('-size', action='store', type=int, help="size of the transposon ex: -size 2133", dest="size")

args = vars(parser.parse_args())

if (args["in"] is None or args["out"] is None or args["size"] is None): #check whether the arguments have been passed
	parser.print_help()
	sys.exit()
else :
	pass

input = args["in"]

output = args["out"]

te_size = args["size"]


def lumpy_2_final (input, te_size, output):
	lines = input.split("\n")
	i = 0

	dict_chr_line = {}
	line_chr = ""

	for line in lines:
		search_site = re.search('^chr', line)
		if (search_site):

			if (line not in dict_chr_line):
				dict_chr_line[line] = ""
				line_chr = line
			first_evidence = True
			count_7_col = False
			count_8_col = False
			i=i+1
			#recupInfoGFF
			Chr_splitted = re.split('\s+',line_chr)
			nb_chr = Chr_splitted[0] #chr nb
			bkpt = Chr_splitted[1] #exact breakpoint/site position insertion on genomics
			details = Chr_splitted[7] #information about the nb of reads supporting the event, etc.

			tab_details= re.split(';', details)
			strand = tab_details[1]
			SU = tab_details[9]

			while (i<len(lines)):
				search_evidence = re.search('^\s+Evidence', lines[i])
				search_next_site = re.search('^chr', lines[i])

				if (search_next_site):
					del best_value
					break
				elif (search_evidence):
					line_splitted = re.split('\s+',lines[i])

					startTpson = int(line_splitted[7])
					endTpson = int(line_splitted[8])

					#info gff
					startChr = line_splitted[4] #start position on genomics
					endChr = line_splitted[5] #end pos on genomics
					fq = line_splitted[2] #read name that supports the information

					if first_evidence == True : #First evidence
						first_evidence = False

						if (startTpson < (int(te_size) - endTpson)) : #we set up the algo: we start with a count on the 7th column...
							count_7_col = True
						else : #ou de 8
							count_8_col = True
					if count_7_col == True : #if we count on 7
						if 'best_value' not in locals() : #we check that best_value is not yet defined and we store the first value of col 7 of the evidence as well as the line
							best_value = startTpson
							best_end = endTpson
							bkptTpson = startTpson
							myL = lines[i] + " " + str(bkptTpson)
							dict_chr_line[line_chr] = myL

						else :
							if (startTpson < best_value):
								best_value = startTpson
								best_end = endTpson
								bkptTpson = startTpson
								myL = lines[i] + " " + str(bkptTpson)
								dict_chr_line[line_chr] = myL

							elif (startTpson <= best_value and endTpson >= best_end): #If ever, in the following lines, the value of col 7 that we read is < the stored value, we overwrite the said col 7.
								best_value = startTpson
								best_end = endTpson
								bkptTpson = startTpson
								myL = lines[i] + " " + str(bkptTpson)
								dict_chr_line[line_chr] = myL


					elif count_8_col == True : #if we count on 8
						if 'best_value' not in locals() : #we check that best_value is not yet defined and we store the first value of col 7 of the evidence as well as the line
							best_value = endTpson
							best_start = startTpson
							bkptTpson = endTpson
							myL = lines[i] + " " + str(bkptTpson)
							dict_chr_line[line_chr] = myL

						else :
							if (endTpson > best_value): #If ever, in the following lines, the value of collar 8 that we read is > to the one stored, we overwrite the said value
								best_value = endTpson
								best_start = startTpson
								bkptTpson = endTpson
								myL = lines[i] + " " + str(bkptTpson)
								dict_chr_line[line_chr] = myL

							elif (endTpson >= best_value and startTpson <= best_start):
								best_value = endTpson
								best_start = startTpson
								bkptTpson = endTpson
								myL = lines[i] + " " + str(bkptTpson)
								dict_chr_line[line_chr] = myL

				i=i+1

	for key in dict_chr_line :
		Chr_splitted = re.split('\s+',key)
		nb_chr = Chr_splitted[0] #chr nb
		bkpt = Chr_splitted[1] #breakpoint/site position "exact" insertion on genomics
		details = Chr_splitted[7] #information about the nb of reads supporting the event, etc.
		bkpt_tpson_main = Chr_splitted[4]

		tab_details= re.split(';', details)
		strand = tab_details[1]
		SU = tab_details[9]

		evidence_splitted = re.split('\s+',dict_chr_line[key])

		startTpson = int(evidence_splitted[7])
		endTpson = int(evidence_splitted[8])

		#info gff
		startChr = evidence_splitted[4] #start position on genomics
		endChr = evidence_splitted[5] #genomic pos fin
		fq = evidence_splitted[2] #name of the read that supports the info

		bkptTp=evidence_splitted[15]


		gff_infos = [nb_chr,'Lumpy','Insertion',str(startChr),str(endChr),'.','+','.','BREAKPOINT='+bkpt+';'+strand+';READ='+fq+';TPSON='+str(startTpson)+'_'+str(endTpson)+';BKPTpson='+bkptTp+';BKPTpsonRef='+bkpt_tpson_main]

		with open(output, 'a') as my_file_results:
			my_file_results.write('\t'.join(gff_infos)+ "\n")


with open(input, 'r') as input_file_data : 
	input_file = input_file_data.read()


lumpy_2_final(input_file, te_size, output)
