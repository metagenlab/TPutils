#!/usr/bin/python

import sys
import os
from Bio import Entrez, SeqIO
import re
from Bio.Seq import Seq

''' this script generate output for EBI submission
take a fasta file with scaffolds and generate :

- an agp file describing the order of contigs
- a fasta file with contigs (instead of scaffolds)

'''



from optparse import OptionParser
parser = OptionParser()

parser.add_option("-i", "--input",dest="input_file",action="store", type="string", help="file with scaffolds as fasta format", metavar="FILE")
parser.add_option("-o", "--file",dest="output",action="store", type="string", help="output name. Write file with contigs (scaffolds split  in contigs is gaps > 1 N)")
parser.add_option("-n", "--organism",dest="organism",action="store", type="string",default="unknown", help="organisme name")
parser.add_option("-s", "--strain",dest="strain",action="store", type="string",default="unknown", help="strain name")
parser.add_option("-d", "--project",dest="project",action="store", type="string",default="unknown", help="project id")



(options, args) = parser.parse_args()







##############

input_file = 		options.input_file				# file with scaffolds as fasta format

output_contigs = 	str(options.output)+".fa"		# file with contigs (scaffolds split  in contigs is gaps > 10 N)
output_agp =		str(options.output)+".agp"	# agp file for submission

output_embl =		str(options.output)+".embl"	# embl file for submission


##optionnal -- leave blank otherwise

organism = options.organism
strain = options.strain
project_id = options.project
experimental_details = ""

#############


def contigs_to_list (input_file) :
	
	
	'''create 2 list with scaffolds, contigs and gaps
	input must be a fasta with scaffolds'''
	
	L = []	# 2D-list of contigs - [scaffolds] [contigs sequence]
	L0 = []	### 2D-list of gaps - [scaffolds] [N sequence]
	
	scaf_iterator = 1

	# for each scaffold
	
	for record in SeqIO.parse(open(input_file, "r"), "fasta"):
		
		#for contigs
		
		seq=str(record.seq)
		cont = re.split("N{1,}", seq)
		L .append(cont)
		
		#for gaps
		
		l_tmp = re.findall("N{1,}", seq)
		for i in range (len(l_tmp)) :
			L0.append(len(l_tmp[i]) + 1)	#+1 to ajust position
		
		print ("Scaffold "+str(scaf_iterator)+ " : split in "+ str(len(cont)) + " contigs")

		scaf_iterator += 1
		
	print ("------------\nNumber of contigs (>1xN):", len(L0)+1, "+", scaf_iterator-2,"contigs-scaffolds")
		
	return L, L0
	
	

def list_to_fasta_contigs (L, output_contigs):
	
	records = []
	nb_contig = 1
	
	# for each scaffold

	for i in range(len(L)) :
		
		# for each contigs

		for j in range(len(L[i])) :
			
			# create the id of the contig
			z = str(nb_contig)
			while len(z) < 5:
				z ="0"+z
			id_contig = "contig"+z
			
			# create the id of the scaffold
			z = str(i+1)
			while len(z) < 5:
				z ="0"+z
			id_scaffold = "scaffold"+z
			
			
			# create the description of the contig
			desc = "[Scaffold name="+ id_scaffold+"]"
			desc +=	"[Contig name="+id_contig+"]"
			if (project_id != "") : desc +=	"[Registered Genome Project ID="+project_id+"]"
			if (organism != "") : desc += "[Organism="+organism+"]"
			if (strain != "") : desc += "[Strain name="+strain+"]"
			if (experimental_details != "") : desc += "[Experimental details="+experimental_details+"]"
			
			# add the sequence to a list
			seq = Seq(L[i][j])
			tmp = SeqIO.SeqRecord(seq, id=id_contig, description=desc)
			records.append(tmp)
	
			nb_contig +=1
			
	SeqIO.write(records, output_contigs, "fasta")
	
	
	
def list_to_AGP (L, L0, output_agp) :
	
	tmp = ""		#string to print
	nb_contig = 1
	
	# for each scaffold

	for i in range(len(L)):
		
		z = str(i+1)
		while len(z) < 5:
			z ="0"+z
		id_scaff =	"scaffold"+z
		
		# for each contigs

		c = 1			# line iterator
		nb_base = 1		#position of bases
		
		for j in range(len(L[i])) :
			
			#for contig
			
			z = str(nb_contig)
			while len(z) < 5: z ="0"+z
			id_contig = "contig"+z
			
			tmp += id_scaff						#name of the scaffold
			tmp += "\t"+ str(nb_base) 				#position 1 contig
			tmp += "\t"+ str(nb_base +len(L[i][j])-1 )	#position 2 contig
			tmp += "\t" + str(c)					#number of the contig/gap
			tmp += "\tW"							# W ?
			tmp += "\t"+id_contig+"\t1"			#contig name + (1) ?
			tmp += "\t"+ str( len(L[i][j]) )			#length of contig
			tmp += "\t+\n"						#end string
			
			nb_base += len(L[i][j])
			#print nb_base
			#if nb_base<101 : print ("Warning, contig"+id_scaff+"-"+id_contig+" under 100bp")
			c += 1
			nb_contig+=1
			
			#for gaps
			
			if j==len(L[i])-1 : continue

			tmp += id_scaff						#name of the scaffold
			tmp += "\t"+ str(nb_base) 				#position 1 gap
			tmp += "\t"+ str(nb_base +L0[j]-1 )		#position 2 gap
			tmp += "\t" + str(c)					#number of contig/gap
			tmp += "\tN"							# N ?
			tmp += "\t"+str(L0[j])+"\tscaffold"		# gap length + type 
			tmp += "\tyes\tpaired-ends\n"			# end line
			
			nb_base += L0[j]
			c+=1
			
	
	f=open(output_agp, "w")
	f.write(tmp)
	f.close()
	
	#print tmp


def fasta_contigs_to_embl_format (input_file, output_embl):
	
	return
	
	### a faire .......
	
	for record in SeqIO.parse(open(input_file), "fasta") :
		
		 print (record)
		 
		 
	'''	 
	output_handle = open(input_file, "r")
	SeqIO.write(sequences, output_handle, "fasta")
	output_handle.close()
	'''
	
	
def gb_to_embl (gb_file, embl_file):
	
	count = SeqIO.convert(gb_file, "genbank", embl_file, "embl")
	print ("Converted %i records" % count)

#gb_to_embl(path+"/Klebsiella.gbk", path+"/Klebsiella.embl")

###############

L, L0= contigs_to_list (input_file)
list_to_fasta_contigs (L, output_contigs)
list_to_AGP(L, L0, output_agp)

fasta_contigs_to_embl_format(output_contigs, output_embl)

