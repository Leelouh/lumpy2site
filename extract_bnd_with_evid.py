#!/usr/bin/env python3

import sys,re,os,cgi,csv, getopt, argparse, time, math, shutil
from copy import deepcopy
from random import randrange, uniform
from numpy import *
import scipy as sp
from pandas import *
import argparse

parser = argparse.ArgumentParser(add_help=True)

parser.add_argument('-in', action='store', type=str, help="VCF file output from lumpy", dest="in")
parser.add_argument('-out', action='store', type=str, help="output file with results", dest="out")

parser.add_argument('-name', action='store', type=str, help="name of the transposon for boundary to filter ex: -name IFP2_neo", dest="name")

args = vars(parser.parse_args())

if (args["in"] is None or args["out"] is None or args["name"] is None): #check whether the arguments have been passed
	parser.print_help()
	sys.exit()
else :
	pass

input = args["in"]

output = args["out"]

te_name = args["name"]

def do_BND_evidence_search (input, te_name, output) :
	lines = input_file.split("\n")
	i = 0

	for line in lines :

		search_BND = re.search("SVTYPE\=BND", line) #my research on BND

		if search_BND :

			line_splitted = line.split("\t")
			line_splitted = list(filter(None, line_splitted))
			if len(line_splitted) >= 5 :
				if te_name in line_splitted[0] and te_name not in line_splitted[4] or te_name in line_splitted[4] and te_name not in line_splitted[0] :
					line_BND = line
					list_BND_evidences = [] #The list will contain all the BND
					j = i #an incrementer to go back

					while j > 0 : #We check that we're not in a weird case (go to the end of the list if your first evidences are at the very beginning of the file)

						j = j - 1 #We're disincrementing to go backwards

						search_evidence = re.search("Evidence\:", lines[j]) #We're looking for the evidence

						if search_evidence : #If there's evidence, we fill out the evidence list.

							list_BND_evidences.append(lines[j])
						else : #If there's no evidence, we finish the search for evidence.
							break

					with open(output, 'a') as my_file_results:
						if len(list_BND_evidences) != 0 :
							my_file_results.write(str(line_BND) + "\n")

							for lines_evidence in list_BND_evidences :
								my_file_results.write(str(lines_evidence) + "\n")

		i += 1

with open(input, 'r') as input_file_data : #The input file is opened for read
	input_file = input_file_data.read()


do_BND_evidence_search(input_file, te_name, output) #The main function is launched
