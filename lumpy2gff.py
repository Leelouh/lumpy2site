#!/usr/bin/env python3

#usage : lumpy2site.py results.txt taille_tpson/plasmide
#lancer sur results.txt
#~/scripts/python/lumpy2site/lumpy2site_inprocess.py res.txt 2156
#/home/lhelou/scripts/python/lumpy2site/lumpy2site_inprocess.py /home/lhelou/ancestralSpecies/Analyses/profil_insertion/analyse_pgep/second_run/IFP2_IFP2/trimo/extract_20/results.txt 2156 > res030120

import sys,re,os,cgi,csv, getopt, argparse, time, math, shutil
from copy import deepcopy
from random import randrange, uniform
from numpy import *
import scipy as sp
from pandas import *

def create_output_dir (output_dir): #Fonction pour creer un dir
	directory = output_dir
	if not os.path.isdir(directory):
		os.makedirs(directory)

def lumpy_2_final (input_file, te_size):
	lines = input_file.split("\n")
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
			nb_chr = Chr_splitted[0] #numero chr
			bkpt = Chr_splitted[1] #position du breakpoint/site insertion "exact" sur le génomique
			details = Chr_splitted[7] #informations relatives aux nb de reads supportant l'événement, etc.

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
					startChr = line_splitted[4] #position de début sur génomique
					endChr = line_splitted[5] #pos fin génomique
					fq = line_splitted[2] #nom du read qui supporte l'info

					if first_evidence == True : #Première évidence
						first_evidence = False

						if (startTpson < (int(te_size) - endTpson)) : #on met en place l'algo : on part sur le comptage de 7
							count_7_col = True
						else : #ou de 8
							count_8_col = True
					if count_7_col == True : #si on compte sur 7
						if 'best_value' not in locals() : #on verifie que best_value n'est pas encore défini et on stock la première valeur de la col 7 de l'évidence ainsi que la ligne
							best_value = startTpson
							best_end = endTpson
							dict_chr_line[line_chr] = lines[i]

						else :
							if (startTpson < best_value): #nontestée
								best_value = startTpson
								best_end = endTpson
								dict_chr_line[line_chr] = lines[i]
							elif (startTpson <= best_value and endTpson >= best_end): # Si jamais, dans les lignes suivantes, la valeur de la col 7 qu'on lit est < a celle qui est stockee, on ecrase ladite
								best_value = startTpson
								best_end = endTpson
								dict_chr_line[line_chr] = lines[i]

					elif count_8_col == True : #si on compte sur 8
						if 'best_value' not in locals() : #on verifie que best_value n'est pas encore défini et on stock la première valeur de la col 7 de l'évidence ainsi que la ligne
							best_value = endTpson
							best_start = startTpson
							dict_chr_line[line_chr] = lines[i]
						else :
							if (endTpson > best_value):# and startTpson <= best_start) : # Si jamais, dans les lignes suivantes, la valeur de la col 8 qu'on lit est > a celle qui est stockee, on ecrase ladite valeur
								best_value = endTpson
								best_start = startTpson
								dict_chr_line[line_chr] = lines[i]
							elif (endTpson >= best_value and startTpson <= best_start):
								best_value = endTpson
								best_start = startTpson
								dict_chr_line[line_chr] = lines[i]
				i=i+1

	for key in dict_chr_line :
		Chr_splitted = re.split('\s+',key)
		nb_chr = Chr_splitted[0] #numero chr
		bkpt = Chr_splitted[1] #position du breakpoint/site insertion "exact" sur le génomique
		details = Chr_splitted[7] #informations relatives aux nb de reads supportant l'événement, etc.

		tab_details= re.split(';', details)
		strand = tab_details[1]
		SU = tab_details[9]

		evidence_splitted = re.split('\s+',dict_chr_line[key])

		startTpson = int(evidence_splitted[7])
		endTpson = int(evidence_splitted[8])

		#info gff
		startChr = evidence_splitted[4] #position de début sur génomique
		endChr = evidence_splitted[5] #pos fin génomique
		fq = evidence_splitted[2] #nom du read qui supporte l'info

		gff_infos = [nb_chr,'Lumpy','Insertion',str(startChr),str(endChr),'.','+','.','BREAKPOINT='+bkpt+';'+strand+';READ='+fq+';TPSON='+str(startTpson)+'_'+str(endTpson)]
		print ('\t'.join(gff_infos))


with open(sys.argv[1], 'r') as input_file_data : #On ouvre le fichier d'entrée pour le lire
	input_file = input_file_data.read()

te_size = sys.argv[2]

lumpy_2_final(input_file, te_size)
