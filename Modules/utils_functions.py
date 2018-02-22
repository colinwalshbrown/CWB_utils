import sys
import getopt
import os
import time
from Bio import SeqIO

################
## From a sps name, finds the position of the sps in a genome aln
# Input: sps, map folder
# Output: the position index= 0; 1; 2...
#################
def FindsSpsPosInAln (sps, genomes_folder):
	try:
		fgenomes=open(genomes_folder+"/genomes", 'r')
	except IOError, e:
		print "Unknown genome file: ",genomes_folder+"/genomes"
		sys.exit()
	for l in fgenomes: # finds the position of the sps in the aln.
		for i in range(len(l.split())):
			if (sps==l.split()[i]):
				sps_pos=i
				break
	fgenomes.close()
	return (sps_pos)


################
## From a sps position in an aln and a fragment (actually, a particular line in the map), retrieves the info of the fragment for the sps of interest
# Input: sps position in the aln (0,1,2...), list containing the map (separator=space), considered line i, fragment number
# Output: beginning of the fragment in the aln, end of the fragment in the aln, orientation of the fragment 
#################
def retrieveFragmentInfo (sps_pos, lmap, i, fragment):
	debut=4*sps_pos
	if fragment == lmap[i][0]: # chooses the good alignment file fragment to open
		beg_fragment=int(lmap[i][debut+2])
		end_fragment=int(lmap[i][debut+3])
		orient = lmap[i][debut+4] #strand of the fragment in the alignment
		return(beg_fragment, end_fragment, orient)
			
################
## adds a record in a dictionary that contains records of aln fragments
# input = dictionary to use, fragment number to be added, map_folder, format of the aln
# output= completed dictionary
#################
def CompleteFragmentDict (dict_fragment, fragment, map_folder, format):
	if format=="fasta":
		ext=".mfa"
	elif format=="clustal":
		ext=".aln"
	frag_handle = open(map_folder+"/"+fragment+ext, "rU")
	dict_fragment[fragment]=list()
	for record in SeqIO.parse(frag_handle, format) :
		dict_fragment[fragment].append(record)
	frag_handle.close()
#	print fragment
	return(dict_fragment)
			
##################
## From a seq position index, outputs a position index in the aln
##################
def seqpos2alnpos (dict_fragment, fragment, sps, beg_frag, end_frag, orientation, pos_genome): # returns the position in the aln corresponding to the position in the genome.
	for record in [0,1,2]:
		if sps in dict_fragment[fragment][record].id:
			n=0
			for i in range(len(dict_fragment[fragment][record].seq)):
				if dict_fragment[fragment][record].seq[i] != "-":
					n=n+1
				if (orientation == "+"):
					if beg_frag + n == pos_genome:
						return (i)
						break
				else:
					if end_frag - n == pos_genome:
						return (i)
						break

##################
## From a position index in an aln fragment, outputs a seq position in the genome (the opposite to the previous function)
##################
def alnpos2seqpos (dict_fragment, fragment, sps, beg_frag, end_frag, orientation, pos_aln): # returns the position in the aln corresponding to the position in the genome.
	for record in [0,1,2]:
		if sps in dict_fragment[fragment][record].id:
			n=0
			for i in range(len(dict_fragment[fragment][record].seq)):
				if dict_fragment[fragment][record].seq[i] != "-":
					n=n+1
				if (orientation == "+"):
					if i==pos_aln:
						pos_genome = beg_frag + n
						return (pos_genome)
						break
				else:
					if i==pos_aln:
						pos_genome = end_frag - n
						return (pos_genome)
						break



##################
## From 1 position in a genome, gives the corresponding position in an aln.
# map_folder: contains the file ("map") where one can find for each aln file the corresponding genome positions for all genomes
# genomes_folder: contains the genome file ("genomes") that gives the order of the genomes in the map file and the alns
# Returns: fragment, pos_a, orient: output fragment name, aln position, orientation
##################
def genomePosition2alnPosition (sps, chr, pos_g, lmap, map_folder, genomes_folder, dict_fragment): 

	fragment = "NA"
	pos_a = 0
	orient = "NA"

## Finds position of the sps in the map file
	sps_pos=FindsSpsPosInAln (sps, genomes_folder) # finds the position of the sps in the aln.
	debut=4*sps_pos

## Finds good aln fragment and coordinates of the fragment in the aln.
	for i in range(len(lmap)):	
		if (chr==lmap[i][debut+1]) or (chr=="chr"+lmap[i][debut+1]): 
			if ((int(lmap[i][debut+2]) <= pos_g) and (int(lmap[i][debut+3]) >= pos_g)):
				fragment = lmap[i][0] # chooses the good alignment file fragment to open
				beg_fragment, end_fragment, orient= retrieveFragmentInfo (sps_pos, lmap, i, fragment)
				break
	# if no fragment could be found, return "NA" values.
	if fragment=="NA":
		return("NA","NA","NA")
	# completes the fragment dictionary if it does not already contains the aln #fragment
	if not dict_fragment.__contains__(fragment):
		dict_fragment=CompleteFragmentDict(dict_fragment, fragment, map_folder, "fasta")
	# opens the aln fragment (from the dictionary completed or not with the previous command) and screens it to find the good coordinate pos_a.
	pos_a = seqpos2alnpos(dict_fragment, fragment, sps, beg_fragment, end_fragment, orient, pos_g)
	return(fragment, pos_a, orient)



##################
# From 1 position in an aln, gives the corresponding position in a genome.
# map_folder: contains the file ("map") where one can find for each aln file the corresponding genome positions for all genomes
# genomes_folder: contains the genome file ("genomes") that gives the order of the genomes in the map file and the alns
# Input: fragment, pos_a, orient: output fragment name, aln position, 
# Returns 
##################
def alnPosition2genomePosition (sps, fragment, pos_a, orient, lmap, map_folder, genomes_folder, dict_fragment): 

## Finds position of the sps in the map file
	sps_pos=FindsSpsPosInAln (sps, genomes_folder) # finds the position of the sps in the aln.
	debut=4*sps_pos
	

# finds the infos for the aln fragment
	for i in range(len(lmap)):	
		if ((int(lmap[i][0]) == int(fragment))):
			beg_fragment, end_fragment, orient= retrieveFragmentInfo (sps_pos, lmap, i, fragment)
			chr = lmap[i][debut+1] # chooses the good chr
			break
										
	# opens the aln fragment and screens it to find the good coordinate pos_a.
	handle = open(map_folder+"/"+fragment+".mfa", "rU")
	pos_g = alnpos2seqpos (dict_fragment, fragment, sps, beg_fragment, end_fragment, orient, pos_a)
	handle.close()

	return(chr, pos_g)


#################
# opens a gff file and enters the info into a list of list (more precisely a list of lines of the gff).
#################
def gff2list(gff_file, liste_numeric): #liste_numeric=liste of column indexes that should be considered numeric
	try:
		flist=open(gff_file,'r')
	except IOError, e:
		print "Unknown gff file: ", gff_file
		sys.exit()
	lgff=list()
	print "...entering the gff file "+gff_file+" in a list..."
	for l in flist:
		lsplit=l.split()
		if "##FASTA" in lsplit[0]:
			break
		elif "#" not in lsplit[0]:
			if "Chromosome" not in lsplit[0] and "tart" not in lsplit[4]: #removes comment and header lines
				lgff.append(l.strip().split("\t"))
				for i in liste_numeric: # transforms the designated columns into numeric values.
					lgff[len(lgff)-1][i]=float(lgff[len(lgff)-1][i])
	print "the length of the gff list is "+str(len(lgff))
	flist.close()
	return (lgff)
	
#################
# opens a pfopo (peak fitting orthologous peaks output) and puts it in a list of lists
#################
def pfopo2list(pfopo_file, n, liste_numeric): # n= the number of sps in the table, liste_numeric=liste of columns to be considered numeric
	try:
		flist=open(pfopo_file,'r')
	except IOError, e:
		print "Unknown pfopo file: ", pfopo_file
		sys.exit()
	lpfopo=list() # first list level=lines in the table
	for l in flist: 
		lpfopo.append(l.split()) #second list level= data
		if "#" not in l.split()[0]:
			if "Sps" not in l.split()[0] and "tart" not in l.split()[4]: #will transform certain cata into numeric values, except for the header line, if present.
				for i in liste_numeric: # transforms the designated columns into numeric values.
					lpfopo[len(lpfopo)-1][i]=float(lpfopo[len(lpfopo)-1][i])
		
	flist.close()
	return (lpfopo)


#################
## finds the position of the sps that is looked at in the pfopo file.
#################
def findsPositionSpsPfopo(linput,sps):
	print "Finding the sps position in the input file..."
	sps_position="NA"
	index=0
	for i in range(len(linput[1])):
		if linput[0][i]=="Sps" or "Sps_" in linput[0][i]:
			if linput[1][i]==sps:
				sps_position=index
				break
			index=index+1
	if sps_position=="NA":
		print "The requested sps "+sps+" has not been found in the pfopo table. Script avorted."
		sys.exit()
	return (sps_position)

#################
## finds the column numbers of the different categories of the pfopo list.
# the parameter list_annotated indicates whether the annotation parameters are present in the input table. Default is False
# Annotation parameters are given if present
#################
def entersPfopoParametersInList(linput, sps_position):
	list_annotated=False
	list_parameters=list()
#	list_parameters=[0,0,0,0,0,0,0,0,0,0,0]
	spsIndex = 0
	chrIndex = 0
	heightIndex = 0
	genPosIndex = 0
	fragmentIndex = 0
	alnPosIndex = 0
	orientIndex = 0
	five_prime_name = 0
	five_prime_dist = 0
	three_prime_name = 0
	three_prime_dist = 0

	index=0
	Sps_index_read=False
	while index <= len(list_parameters) :
		if "Sps" in linput[0][index]:
			if not Sps_index_read:
				spsIndex = index
				list_parameters.append(index)
#				list_parameters[0] = spsIndex
				Sps_index_read=True
			else:
				break # on a fait la boucle
			index=index + 1
		if "Chr" in linput[0][index]:
			chrIndex = index
			list_parameters.append(index)
			index=index + 1
		if "Peak_height" in linput[0][index]:
			heightIndex = index
			list_parameters.append(index)
			index=index + 1
		if "Genome_position" in linput[0][index]:
			genPosIndex = index
			list_parameters.append(index)
			index=index + 1
		if "Fragment" in linput[0][index]:
			fragmentIndex = index
			list_parameters.append(index)
			index=index + 1
		if "Alignment_position" in linput[0][index]:
			alnPosIndex = index
			list_parameters.append(index)
			index=index + 1
		if "orientation" in linput[0][index]:
			orientIndex = index
			list_parameters.append(index)
			index=index + 1
		if "five_prime_gene" in linput[0][index]:
			five_prime_name = index
			list_parameters.append(index)
			list_annotated=True # indicates that the annotations are present in the table
			index=index + 1
		if "five_prime_dist" in linput[0][index]:
			five_prime_dist = index
			list_parameters.append(index)
			index=index + 1
		if "three_prime_gene" in linput[0][index]:
			three_prime_name = index
			list_parameters.append(index)
			index=index + 1
		if "three_prime_dist" in linput[0][index]:
			three_prime_dist = index
			list_parameters.append(index)
			index=index + 1
		index=index+1
	for i in range(len(list_parameters)):
		list_parameters[i]=list_parameters[i]+sps_position*len(list_parameters)
	return(list_parameters, list_annotated)

#################
## enters the genome in a dictionnaire: 
# for each chromosome=key, puts a string=sequence of the chromosome, chr_prst is a parameter that monitors whether the chromosome names start with "chr" (that will be removed by cka)
#################
def enters_seq_genome_dans_dict(genome, dict_genome, chr_prst): 
	test_premier=True #tests whether we are reading the first chromosome (for function enters_seq_genome_dans_dict)
	chr=""
	sequence_chr=""
	try:
		fgenome=open(genome, 'r')
	except IOError, e:
		print " Could not enter the genome sequence in a dictionary because UNKNOWN GENOME FILE: ",genome
		sys.exit()

	print("Reading genome file...")
	for l in fgenome:
		if (">" in l):
			if (test_premier):
				test_premier=False
			else:
				dict_genome[chr]=sequence_chr
				sequence_chr=""
				chr=""
			chr=l.split(">")[1].strip()
			if "chr" in chr[0:3]:
				chr_prst=1
		else:
			sequence_chr=sequence_chr+l.strip()
	fgenome.close()
	dict_genome[chr]=sequence_chr
	return (dict_genome, chr_prst)

#################
#This function extracts the best BS from a patser output, and returns its characteristics
#The best BS is either the closest to the center (criterion="center") or the best p-value (criterion="pvalue")
#################
def extractBestPeakFromPatserOutput(file, criterion):
	try:
		f=open(file, 'r')
	except IOError, e:
		print "Unknown patser file: ",file
		sys.exit()
	starting=False
	sequenceToBS = dict() # a dictionary that links the sequence to all BS found in this sequence
	for l in f:
		if ("average score above numerically calculated cutoff" in l):
			starting=True
		elif starting:
			if ":" in l: #line to analyse
				liste = l.split()
				if "ln(p-value)" in l:
					if liste[0] in sequenceToBS:
						sequenceToBS[liste[0] ].append([liste[2], liste[4], liste[6]]) # one triplet of info per BS
					else:
						sequenceToBS[liste[0] ] = list()
						sequenceToBS[liste[0] ].append([liste[2], liste[4], liste[6]]) # one triplet of info per BS
				else:
					if liste[0] in sequenceToBS:
						sequenceToBS[liste[0] ].append([liste[2], liste[4], "10"]) # one triplet of info per BS
					else:
						sequenceToBS[liste[0] ] = list()
						sequenceToBS[liste[0] ].append([liste[2], liste[4], "10"]) # one triplet of info per BS
					
			
	f.close()
	#done with file parsing, now for each sequence we extract the best BS.
	sequenceToBestBS = dict()
	for k,v in sequenceToBS.iteritems():
		sequenceToBestBS[k] = list()
		best = list()
		if (criterion == "center"):
			begend = k.split(":")[1].split("-")
			beg = int(begend[0])
			end = int(begend[1])
			mid = (end - beg) /2
			mindist = 100000000000000
			for l in v: #we go through all the BS triplets
				posStr = l[0]
				pos = int()
				comp = False
				if ("C" in posStr):
					comp = True
					pos = int(posStr.replace("C",""))
				else:
					pos = int(posStr)
				if (abs (pos - mid) < mindist): # better candidate BS than the current one
					mindist = abs (pos - mid)
					best = list()
					best.append(str(pos))
					if comp:
						best.append("C")
					else :
						best.append("F")
					best.append(l[1])
					best.append(l[2])
		elif (criterion == "pvalue"):
			minPValue = 100
			for l in v: #we go through all the BS triplets
				posStr = l[0]
				pos = int()
				comp = False
				if ("C" in posStr):
					comp = True
					pos = int(posStr.replace("C",""))
				else:
					pos = int(posStr)
				if (float(l[2]) < minPValue): # better candidate BS than the current one
					minPValue = float(l[2])
					best = list()
					best.append(str(pos))
					if comp:
						best.append("C")
					else :
						best.append("F")
					best.append(l[1])
					best.append(l[2])
		sequenceToBestBS[k] = best
	return (sequenceToBestBS)
 

#######
# Parcourt les sgrs de tout le genome et les met dans 2 tableaux.
#######
def sgr2dict(sgr):
	try:
		fsgr=open(sgr,'r')
	except IOError, e:
		print "Unknown file sgr: ", sgr
		sys.exit()
	time1=time.clock()
	num_pos=1	
	dsgr=dict()
	print ("Reading sgr: "+sgr)
	sgr_value_current_chr=list()
	previous_chr=""
	current_chr=""
	for l in fsgr:
		current_chr=l.split()[0]
		position_in_chr=int(l.split()[1])
		if current_chr!= previous_chr or previous_chr=="":
			if current_chr != previous_chr and previous_chr!="":
				dsgr[previous_chr]=sgr_value_current_chr
				sgr_value_current_chr=list()
			previous_chr=current_chr
			if (not dsgr.__contains__(current_chr)):
				dsgr[current_chr]=list()
				sgr_value_current_chr.append("0")
				#dsgr[current_chr].append("NA")
				num_pos=1
				previous_chr=current_chr
		while(num_pos < position_in_chr):
			sgr_value_current_chr.append("0")	
			#dsgr[l.split()[0]].append("NA")
			num_pos=num_pos+1		

		sgr_value_current_chr.append(l.split()[2])
#		dsgr[l.split()[0]].append(l.split()[2])
		num_pos=num_pos+1
	time2=time.clock()
	fsgr.close()
	print ("time "+str(time2-time1))
	print ("Sgr "+sgr+" file contains "+str(len(dsgr))+" entries.")
	return(dsgr)


#################
# from a position, a sgr dictionary and a window, gives a list containing sgr values of a window centered on the given position.
#################
def listSgrWindow(dsgr, chr, position, window):
	lsgr_out=list()
	lsgr_out=dsgr[chr][int(position)-int(int(window)/2):int(position)+int(int(window)/2)+1]
	return(lsgr_out)
