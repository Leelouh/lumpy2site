#!/usr/bin/env python3

import sys,re,os,cgi,csv, getopt, argparse, time, math, shutil
from copy import deepcopy
from random import randrange, uniform
from numpy import *
import scipy as sp
from pandas import *
import argparse

parser = argparse.ArgumentParser(add_help=True)
#parser = argparse.ArgumentParser()

# parser.add_argument('-in', action='store', type=argparse.FileType('r'), help="VCF file output from lumpy", dest="in")
# parser.add_argument('-out', action='store', type=argparse.FileType('a'), help="output file with results", dest="out")

parser.add_argument('-in', action='store', type=str, help="VCF file output from lumpy", dest="in")
parser.add_argument('-out', action='store', type=str, help="output file with results", dest="out")

parser.add_argument('-name', action='store', type=str, help="name of the transposon for boundary to filter ex: -name IFP2_neo", dest="name")

# #parser.add_argument('-a', action="store_true", default=False)
# parser.add_argument('-b', action="store", dest="b")
# parser.add_argument('-c', action="store", dest="c", type=int)

# print('print_usage output:')
# parser.print_usage()
# print()
#
# print('print_help output:')
# parser.print_help()
#
# print(parser.parse_args())
args = vars(parser.parse_args())

if (args["in"] is None or args["out"] is None or args["name"] is None): #vérifier si les arguments ont été passé
	parser.print_help()
	sys.exit()
else :
	pass

#input=args["in"]
#print(input)

#print("fichier entrée {}".format(args["in"]))

# parser.add_argument('--in', action='store',required=True, type=argparse.FileType('r'), help="VCF file output from lumpy", dest="in")
# parser.add_argument('--out', action='store',required=True, type=argparse.FileType('a'), help="output file with results", dest="out")
# parser.add_argument('--name', action='store',required=True, type=str, help="name of the transposon for boundary to filter ex: --name IFP2_neo", dest="name")
#
#
# args = parser.parse_args('--input --output --name'.split())
#
# args = parser.parse_args()

input = args["in"]

output = args["out"]

te_name = args["name"]

# def create_output_dir (output_dir) : #Fonction pour creer un dir
#
# 	directory = output_dir
# 	if not os.path.isdir(directory):
# 		os.makedirs(directory)

def do_BND_evidence_search (input, te_name, output) : #Fonction principale
	lines = input_file.split("\n")

	i = 0

	for line in lines :

		search_BND = re.search("SVTYPE\=BND", line) #ton grep sur BND

		if search_BND : #Si ca marche :

			line_splitted = line.split("\t")
			line_splitted = list(filter(None, line_splitted))
			# if len(line_splitted) >= 5 :
			# 	if te_name in line_splitted[0] and te_name not in line_splitted[4] or te_name in line_splitted[4] and te_name not in line_splitted[0] :
			if len(line_splitted) >= 5 :
				if te_name in line_splitted[0] and te_name not in line_splitted[4] or te_name in line_splitted[4] and te_name not in line_splitted[0] :

					line_BND = line

					list_BND_evidences = [] #La liste qui contiendra tous les BND

					j = i #Un incrémenteur pour revenir en arrière

					while j > 0 : #On vérifie quand même qu'on soit pas dans un cas bizarre (aller à la fin de la liste si tes premières évidences sont au tout début du fichier)

						j = j - 1 #On désincrémente pour aller en arrière

						search_evidence = re.search("Evidence\:", lines[j]) #On cherche les évidences

						if search_evidence : #Si y'a évidence on remplit la liste d'évidences

							list_BND_evidences.append(lines[j])
						else : #Si y'a pas évidence on termine la recherche d'évidences.
							break

					with open(output, 'a') as my_file_results:
						if len(list_BND_evidences) != 0 :
						#if len(list_BND_evidences) == 0 :
							#print ("nothing here !")
						#else :
							my_file_results.write(str(line_BND) + "\n")

							for lines_evidence in list_BND_evidences :
								my_file_results.write(str(lines_evidence) + "\n")
			#else :
				#print ("this is strange")

		i += 1

with open(input, 'r') as input_file_data : #On ouvre le fichier d'entrée pour le lire
	input_file = input_file_data.read()

#create_output_dir(sys.argv[2]) #On crée un dictionnaire pour mettre les résultats dedans

#output_file = args.out


do_BND_evidence_search(input_file, te_name, output) #On lance la fonction principale
