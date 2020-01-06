#!/usr/bin/env python3

#usage : lumpy2site.py results.txt taille_tpson/plasmide
#lancer sur results.txt
#~/scripts/python/lumpy2site/lumpy2site_inprocess.py res.txt 2156

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
			#best_value = False
			#print (line)
			i=i+1
			while (i<len(lines)):
				search_evidence = re.search('^\s+Evidence', lines[i])
				search_next_site = re.search('^chr', lines[i])

				if (search_next_site):
					del best_value
					break
				elif (search_evidence):
					line_splitted = re.split('\s+',lines[i])

					start = int(line_splitted[7])
					end = int(line_splitted[8])
					#print ("start =",start)
					if first_evidence == True : #Première évidence
						#best_value = False
						first_evidence = False

						if (start < (int(te_size) - end)) : #on met en place l'algo : on part sur le comptage de 7
						#if int(line_splitted[7]) < (int(te_size) - int(line_splitted[8])) :
							#print("start ",start,"\ntpson-end ",int(te_size) - end,"\n")
							count_7_col = True
							#best_value = False
						else : #ou de 8
							count_8_col = True
							#best_value = False
					if count_7_col == True : #si on compte sur 7
						if 'best_value' not in locals() :
						#if best_value == False : # On lance et on stocke la premiere valeur de la col 7 de l'evidence (et la ligne
							best_value = start
							dict_chr_line[line_chr] = lines[i]
							best_end = end
							#print ("bestvalue false deb : ",best_value)
						else :
							#best_value
							if (start <= best_value and end >= best_end): # Si jamais, dans les lignes suivantes, la valeur de la col 7 qu'on lit est < a celle qui est stockee, on ecrase ladite
								#if (end>best_end):
								best_value = start
								dict_chr_line[line_chr] = lines[i]
								best_end = end
							#print ("la ligne", lines[i])
						# elif (start > best_value):
						# 	pass


					elif count_8_col == True : #si on compte sur 8
						if 'best_value' not in locals() :
						#if best_value == False : # On lance et on stocke la premiere valeur de la col 8 de l'evidence (et la ligne)
							best_value = end
							dict_chr_line[line_chr] = lines[i]
							best_start = start
							#print ("bestvalue false deb : ",best_value)
						else :
							if (end >= best_value and best_start >= start) : # Si jamais, dans les lignes suivantes, la valeur de la col 8 qu'on lit est > a celle qui est stockee, on ecrase ladite valeur
								#if (best_start<=start):
								best_value = end
								dict_chr_line[line_chr] = lines[i]
					#print (start,"\t",end,"\n")

				i=i+1

#	print (dict_chr_line)

	for key in dict_chr_line :
		print(key)
		print(dict_chr_line[key])

with open(sys.argv[1], 'r') as input_file_data : #On ouvre le fichier d'entrée pour le lire
	input_file = input_file_data.read()

te_size = sys.argv[2]

lumpy_2_final(input_file, te_size)
