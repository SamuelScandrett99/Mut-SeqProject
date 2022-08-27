#!/usr/bin/env python

import sys
import getopt		# remove
from tkinter import *
import random
import Bio.SeqIO
from collections import OrderedDict
import numpy as np
import pandas as pd	
import matplotlib.pyplot as plt
import logomaker as lm
import PIL
from PIL import Image
from decimal import Decimal
import os


def main(): 
	
	# Declaring global CONSTANTS throughout simulation.
	# Since these are stated by the user at the start of the program they do not
	# change, but are referenced in various functions throughout.

	global total_cycles
	global runs
	global error_rate
	global SNVc
	global INSc
	global DELc
	global log_output
	global logFile
	global output
	global output_name
	global figures
	global outputCovGraph
	global bases
	global default
	global mutdict
	global no_logos

	# Lists/dictionaries used throughout

	bases = ['A', 'T', 'C', 'G']
	mutdict = {}   		# dict of no. of point mutations
	mutpositions = {} 	# count of how many muts at each nuc position (mutPlot)
	nopositions = {}	# count of positions unmutated (mutPlot)
	filteredMuts = [] 	# master list of all filtered mutations (1 mutation)
	mutcoverage = []	# % coverage per run (1: 1%)
	coverage = {}		# mutations per base location (1: ['A','C'])
	masterMutDict = {}	# Dictionary of all sequences, rather than lists
	filterMutDict = {}


	# Error-prone polymerase probabilities

	default = [25,25,25,25]		# no bias



	##################
	### PARAMETERS ###
	################## 

	total_cycles = None			# Number of cycles
	runs = None					# Number of parallel runs
	error_rate = None			# Polymerase error rate 
	SNVc = Decimal("0.82")		# Probability of point mutation
	INSc = Decimal("0.06")		# Probability of insertion	
	DELc = Decimal("0.12")		# Probability of deletion	
	log_output = False		
	output = False
	output_name = None
	figures = True
	outputCovGraph = True
	


	##############################
	### COMMAND LINE ARGUMENTS ###
	##############################

	argv = sys.argv[1:]

	if len(sys.argv) <= 7:
		titleScreen()
		sys.exit(1)
	
	try:
		opts, args = getopt.getopt(argv, "i:o:c:r:e:lxys:d:n:", ["input=",
									 	 "cycles=",
										 "runs=",
										 "error-rate=",
										 "log-output-file",
										 "output=",
										 "no-figures",
										 "no-coverage"
										 "snv="
										 "deletion="
										 "insertion="])
								 
	except getopt.GetoptError as err:
		print(err)
		titleScreen()
		sys.exit(2)	
	
	for opt, arg in opts:
		if opt in ('-i', '--input'):
			pbp1a = open(arg)
			for gene in Bio.SeqIO.parse(pbp1a, format = "fasta"):
				identifier, sequence = gene.id, str(gene.seq)
			# No. stacked logos in plot
			if len(sequence) > 500:
				no_logos = 5
			else:
				extract = str(len(sequence))
				no_logos = int(extract[0])
			
		elif opt in ('-c', '--cycles'):
			total_cycles = int(arg)
		elif opt in ('-r', '--runs'):
			runs = int(arg)
		elif opt in ('-e', '--error-rate'):	
			error_rate = float(arg)
		elif opt in ('-l', '--log-output-file'):
			log_output = True
		elif opt in ('-o', '--output'):
			output = True
			output_name = arg
		elif opt in ('-x', '--no-figures'):
			figures = False	
			outputCovGraph = False
		elif opt in ('-y', '--no-coverage'):
			outputCovGraph = False
		elif opt in ('-s', '--snv'):
			SNVc = Decimal(arg) 
		elif opt in ('-d', '--deletion'):
			DELc = Decimal(arg) 
		elif opt in ('-n', '--insertion'):
			INSc = Decimal(arg) 
	
	
	#Write mutlog file#
	logFileName = ('mutLog%s%s%s%s.txt' % (
			identifier,total_cycles,runs,error_rate))
	if log_output == True:
		logFile = open(logFileName, 'w')

	#Write fasta file
	if output == True:
		mutFasta = open(output_name, 'w')



	########################
	### BEGIN SIMULATION ###
	########################

	# Error if mut probabilities not = 1
	total = [SNVc, INSc, DELc]
	if sum(total) != 1:
		print("Probability does not equal 1 (100%), adjust accordingly")
		return

	runningScreen()
	run_iterator = 1
	while run_iterator <= runs:
		
		logAppend("Run: %s" % (run_iterator))
		final_DNA, final_muts = cycle(sequence)	

		# Store all final mutant sequences (final_muts) in dictionary w/ counter
		# as needed for the ALL logoplot - masterMutDict
		for item in final_muts:
			dictCounter(item, masterMutDict)

		# Filter out all mutants != 1 mutation after each run & store in dictionary
		# w/ counter again as needed for FILTERED logoplot - filteredMuts
		if outputCovGraph == True:
			filterMutDict, coverage = filtration(masterMutDict, filterMutDict, coverage, sequence)
			mutcoverage = calculatingCoverage(coverage, mutcoverage, sequence)

		# If output is enabled, write the DNA at the end of each run to FASTA file
		mutant = 1
		if output == True:
			fastaSave(final_DNA, identifier, mutFasta, mutant)
			

		run_iterator = run_iterator + 1
	
	
	#This gives end coverage rather than per run - faster
	if outputCovGraph == False:
		filterMutDict, coverage = filtration(masterMutDict, filterMutDict, coverage, sequence)
		mutcoverage = calculatingCoverage(coverage, mutcoverage, sequence)
		covName = None
	
	if outputCovGraph == True:
		covName = runs_coverage(mutcoverage, identifier)

	
	#Figures (default is true)
	if figures == True:
		mpsName = mutPerSeq(final_muts, identifier)	
		pmcName = point_mutcount(identifier)
		lpALLName = logoplot(masterMutDict, identifier, "ALL")
		lpFILTName = logoplot(filterMutDict, identifier, "FILTERED")
		mpmutsName = mutplot(filterMutDict, mutpositions, None, identifier, "muts")
		mpbothName = mutplot(filterMutDict, mutpositions, nopositions, identifier, "both")
		if mutcoverage[-1] != 100:
			uncoveredList = invertingLogo(coverage, sequence, bases)
			lpINVName = logoplot(uncoveredList, identifier, "FILTERED_INVERTED")
	else:
		mpsName = None
		pmcName = None
		lpALLName = None
		lpFILTName = None
		mpmutsName = None
		mpbothName = None
		lpINVName = None

	#if outputCovGraph == True:
	#	covName = runs_coverage(mutcoverage, identifier)	


	#Add % coverage to optional log file
	lengthMaster = dictLens(masterMutDict)
	logAppend("\nTotal mutants: %s" % (lengthMaster))
	lengthFilter = dictLens(filterMutDict)
	logAppend("Total filtered mutants (with one mutation): %s" % ((lengthFilter)))
	logAppend("%s cycles + %s runs gives %.2f%% coverage" % (
		total_cycles, runs, mutcoverage[-1])) 
	
	outputScreen(lengthMaster, lengthFilter, mutcoverage, mpsName, pmcName, lpALLName, lpFILTName, 
		mpmutsName, mpbothName, lpINVName, covName, logFileName)
	

def cycle(seqDNA):
	
	all_seqs = []		# all products from each cycle (NOT INCLUDING START SEQ)
	all_seqs_mut = []	# all products only showing mutations

	temp_list1 = []		
	temp_list2 = []	
	mut_temp_list1 = []
	mut_temp_list2 = []
	
	output_DNA = []		# Final sequences from each run (after total_cycles) 
	output_muts = []	# Final sequences with only mutant bases and fullstops

	cycle_no = 1	# iterator for cycle no.
	
	temp_list1.append(seqDNA)
	
	new = "." * len(seqDNA)
	mut_temp_list1.append(new)
	
	while cycle_no <= total_cycles:
		logAppend("Cycle: %s" % (cycle_no))
		
		for each, mut_each in zip(temp_list1, mut_temp_list1):	
			
			copy = each
			mutCopy = mut_each
			
			all_seqs.append(copy)
			all_seqs_mut.append(mutCopy)
			
			temp_list2.append(all_seqs[-1])
			mut_temp_list2.append(all_seqs_mut[-1])
			
			mutation(each, mut_each, default, all_seqs, all_seqs_mut)
			
			temp_list2.append(all_seqs[-1])
			mut_temp_list2.append(all_seqs_mut[-1])
		
		# Clear large lists from memory
		temp_list1.clear()
		temp_list1.extend(temp_list2)
		temp_list2.clear()	
		mut_temp_list1.clear()
		mut_temp_list1.extend(mut_temp_list2)
		mut_temp_list2.clear()
		all_seqs.clear()	

		cycle_no = cycle_no + 1	
	
	output_DNA.extend(temp_list1) 
	output_muts.extend(mut_temp_list1)
	
	return output_DNA, output_muts	
	

def mutation(insert, mutSeq, polBiases, allSeqs, allSeqsMuts):
		
	gene_length = len(insert)
	n = 1
	while n <= gene_length:
		x = round(random.random(), 5)
		
		# If random number = error rate reroll
		while x == error_rate:
			x = round(random.random(), 5)
		
		if x < error_rate:
			
			nucleotide = insert[n-1]	
			
			# Randomiser for selecting which mutation type	
			m = round(random.random(), 2)
			
			while True:
				if m == 0:
					m = round(random.random(), 2)
				else:
					break 
		
			INSlower = 1 - INSc		
			DELlower = INSlower - DELc	
			SNVlower = DELlower - SNVc	
		
			SNVupper = 0 + SNVc
			DELupper = SNVupper + DELc
			INSupper = DELupper + INSc
			
			# Point mutation
			if SNVlower <= m <= SNVupper:
				point = random.choices(bases, weights=polBiases)
				point = str(point[0])
				while True:			
					if point == nucleotide:
						point = random.choices(bases, weights=polBiases)
						point = str(point[0])
					else:
						break

				changed = nucleotide + " to " + point
				dictCounter(changed, mutdict)	
				logAppend("Point mutation %s at %s" % (changed, n-1))	
				insert = insert[:n-1] + point + insert[n:] 
				mutSeq = mutSeq[:n-1] + point + mutSeq[n:] 
				
			# Deletion mutation	
			elif DELlower <= m <= DELupper:
				logAppend("Deletion at %s" % (n-1))
				insert = insert[:n-1] + insert[n:]
				mutSeq = mutSeq[:n-1] + mutSeq[n:]
				gene_length -= 1
				
			# Insertion mutation	
			elif INSlower <= m <= INSupper:
				ins = random.choices(bases, weights=default)
				ins = str(ins[0])
				befaft = random.randint(1,2)
				if befaft == 1:		# insert before n
					insert = insert[:n-1] + ins + insert[n-1:]
					mutSeq = mutSeq[:n-1] + ins + mutSeq[n-1:]
					logAppend("Insertion of %s before %s" % (ins, n-1))
				elif befaft == 2:	# insert after n
					insert = insert[:n] + ins + insert[n:]
					mutSeq = mutSeq[:n] + ins + mutSeq[n:]
					logAppend("Insertion of %s after %s" % (ins, n-1))	 
				gene_length += 1
	
		n += 1			
	allSeqs.append(insert)
	allSeqsMuts.append(mutSeq)			


def fastaSave(selection, mutName, fileName, mutNo):
	
	wrapping = 60
	
	if output == True:
		for each in selection:
	
			name = (">%s_M%s" % (mutName, mutNo))
			fileName.write(name + "\n")
		
			for x in range(0, len(each), wrapping):
				fileName.write(each[x:x + wrapping] + "\n")
			mutNo += 1

def logAppend(text):
	
	if log_output == True:
		logFile.write(text + "\n")
			

def dictLens(dict_to_use):
	
	tempList = []
	for item, count in dict_to_use.items():
		tempList.extend([item for i in range(count)])

	length = len(tempList)

	return length
	

def filtration(currentDict, futureDict, locDict, insert):

	# Going through masterDict and appending to filtered dictionary
	futureDict = {}
	locDict = {}
	for item, count in currentDict.items():
		mutcount = 0
		falsemutcount = False
		for coord, (base, original) in enumerate(zip(item, insert)):
			if base.isupper():
				mutcount += 1
				pos = coord
				mut = base
				if base == original:
					falsemutcount = True
			
		if mutcount == 1 and len(item) == len(insert) and falsemutcount == False:
			futureDict[item] = count

			if 1+pos not in locDict:
				locDict[1+pos] = list()
			if base not in locDict[1+pos]:
				locDict[1+pos].append(mut)

	return futureDict, locDict
	
	
def calculatingCoverage(locDict, runningCov, insert):

	# Calculate cov each run
	count = 0
	for each, each2 in locDict.items():
		count += len(each2)
	
	cov_all = len(insert) * 3
	cov = (count/cov_all) * 100
	runningCov.append(cov)

	return runningCov


def invertingLogo(dictOfBaseCov, insert, nucleotides):

    # calculating every possible base change
    totalCov = {}
    for loc, b in enumerate(insert):
        totalCov[1+loc] = list()
        for each in nucleotides:
            if each != b:
                totalCov[1+loc].append(each)

    # Dictionary of bases not changed to    
    uncovered = {}
    for (a, b) in dictOfBaseCov.items():
        for (x, y) in totalCov.items():
            if a == x:
                difference = (list(set(y) - set(b)))
                uncovered[a] = difference
            elif x in uncovered:
                continue
            elif a != x:
                uncovered[x] = y

    #generate 4 strings of each of the missing bases
    uncovlist = []
    for atgc in nucleotides:
        template = ""
        for loc, letters in uncovered.items():
            if atgc in letters:
                template += atgc
            else: 
                template += "'"
        uncovlist.append(template)
    
    return uncovlist


def titleScreen():

	print("\n\n\t\t"+"-"*41)
	print("\n\t\t\tError-Prone PCR Simulator")
	print("\t\t\t    Samuel Scandrett")
	print("\n\t\t"+"-"*41)
	print("\n\n"+"-"*4+" USAGE "+"-"*4)
	print("\nepPCR.py -i (input_file) -c (cycles) -r (runs) -e (error_rate) -o (output_file)")
	print("\n"+"-"*4+" PARAMETERS "+"-"*4)
	print("""\n-c\t--cycles\t\t<int>\t\tThe number of cycles
-d\t--deletion\t\t<float>\t\tThe probability for a deletion mutation (default: 0.12)
-e\t--error-rate\t\t<float>\t\tThe error rate as decimal
-i\t--input\t\t\t<str>\t\tThe file to have epPCR performed on
-l\t--log-output-file\t\t\tEnable log file production containing each mutation and percentage coverage
-n\t--insertion\t\t<float>\t\tThe probability for an insertion mutation (default: 0.06)
-o\t--output\t\t<str>\t\tThe name of the output file in fasta format
-r\t--runs\t\t\t<int>\t\tThe number of runs (template molecules)
-s\t--snv\t\t\t<float>\t\tThe probability for an SNV mutation (default: 0.82)
-x\t--no-figures\t\t\t\tDisable figure production
-y\t--no-coverage\t\t\t\tDisable production of coverage graph (FASTER)\n""")


def runningScreen():
	
	print("\n\n\t\t"+"-"*41)
	print("\n\t\t\tError-Prone PCR Simulator")
	print("\t\t\t    Samuel Scandrett")
	print("\n\t\t"+"-"*41)
	print("\nError-Prone PCR Parameters:")
	print("\n\tCycles:\t\t\t%s" % (total_cycles))
	print("\tRuns:\t\t\t%s" % (runs))
	print("\tError Rate:\t\t%s" % (error_rate))
	print("\tSNV probability:\t%s" % (SNVc))
	print("\tInsertion probability:\t%s" % (INSc))
	print("\tDeletion probability:\t%s" % (DELc))
	print("\tOutput log file:\t%s" % (log_output))
	print("\tProduce figures:\t%s" % (figures))
	print("\tProduce coverage graph:\t%s" % (outputCovGraph))
	print("\tOutput fasta:\t\t%s" % (output))
	

def outputScreen(totalMuts, oneMutMuts, perCov, mutperseqName, pointcountName, logoALL, logoFILT, 
	mutplotmuts, mutplotboth, logoINV, covg, logF):
	
	print("\nCoverage:")
	print("\n\tTotal mutants: %s" % ((totalMuts)))
	print("\tTotal filtered mutants (with one mutation): %s" % ((oneMutMuts)))
	print("\t%s cycles + %s runs gives %.2f%% coverage" % (
		total_cycles, runs, perCov[-1])) 
	
	print("\nOutput Files:\n")
	if output == True:
		print("\tMutant fasta:\t\t%s" % (output_name))
	if figures == True:
		print("\tMutations per seq:\t%s" % (mutperseqName))
		print("\tSNV count:\t\t%s" % (pointcountName))
		print("\tMutation plots:\t\t%s" % (mutplotmuts))
		print("\t\t\t\t%s" % (mutplotboth))
		print("\tLogo plots:\t\t%s" % (logoALL))
		print("\t\t\t\t%s" % (logoFILT))
	if perCov[-1] != 100 and logoINV != None:
		print("\t\t\t\t%s" % (logoINV))
	if outputCovGraph == True:
		print("\tCoverage:\t\t%s" % (covg))
	if log_output == True:
		print("\tLog file:\t\t%s" % (logF))
	if (output == False and figures == False and 
		outputCovGraph == False and log_output == False):
		print("\tNone")
	print("")	

	
def dictCounter(dictItem, dictName):
	
	# if not in dict then add and make value = 0
	# if is in dict then value will increase by 1
	
	if dictItem not in dictName.keys():
		dictName[dictItem] = 0
	dictName[dictItem] += 1


def logoplot(input_mutants, geneIdentifier, name):

	# If the input is a dict, if not treat as a list
	list_of_mutants = []
	if type(input_mutants) == dict:
		for item, count in input_mutants.items():
			list_of_mutants.extend([item for i in range(count)])
	else:
		list_of_mutants = input_mutants


	# Plotting Logoplot in segments of 100 bases
	im_list = []	
	length = 0

	cscheme = {
        "'": 'white',
        'ATG': [0, 0.5, 0],
        'T': [1, 0, 0],
        'C': [0, 0, 1],
        'G': [1, 0.65, 0]}

	while length != no_logos*100:
		
		im_name = ('%s_logo%s%s-%s_%sc%sr%s.jpg' % (
			geneIdentifier, name, length, length+100, 
			total_cycles, runs, error_rate))
		im_list.append(im_name)
		
		segment = (list(map(lambda i:i[length:length+100], list_of_mutants)))
		counts_mat = lm.alignment_to_matrix(segment)

		seqs_logo = lm.Logo(counts_mat, 
				font_name = 'DejaVu Sans', 
				color_scheme = cscheme,
				figsize = (30,4))
	
		seqs_logo.ax.set_xlabel("Nucleotide number")
		seqs_logo.ax.set_title("Logo plot of %s bases %s to %s" % (
					geneIdentifier, length, (length+100)))
		
		plt.savefig(im_name)	
		#plt.show()
		plt.close()
		
		length = length + 100
	
	# Compile all logos of 100 bases
	o = [PIL.Image.open(x) for x in im_list]
	min_shape = sorted([(np.sum(x.size), x.size) for x in o])[0][1]
	ims_all = np.vstack([np.asarray(x.resize(min_shape)) for x in o])
	ims_all = PIL.Image.fromarray(ims_all)
	figureName = ('logoCompiled%s_%sc%sr%s_%s.jpg'
		% (name, total_cycles, runs, error_rate, geneIdentifier))
	ims_all.save(figureName)
	
	# Delete individual logos after compiling
	for individual_logo in im_list:
		os.remove(individual_logo)

	return figureName

def point_mutcount(geneIdentifier):
	
	mutdict_percent = {}
	total = sum(mutdict.values())
	
	for entry, p in mutdict.items():	# raw no. of mutations 
		percent = p * 100 / total	# converted to percentages
		mutdict_percent[entry] = percent
		
	sorted_mutdict = OrderedDict(sorted(mutdict_percent.items())) # alphabetically
	
	cols = ['green','green','green','blue','blue','blue',
		'orange','orange','orange','red','red','red']
		
	plt.bar(sorted_mutdict.keys(), sorted_mutdict.values(), color = cols)
	plt.title('Percentage of each point mutation')
	plt.xlabel('Point mutation')
	plt.ylabel('Percentage of mutations (%)')
	plt.xticks(rotation = 'vertical')	
	figureName = ('SNVcount_%sc%sr%s_%s.svg' % (total_cycles, runs, 
		error_rate, geneIdentifier))
	plt.savefig(figureName)
	plt.close()
	
	return figureName
	
def mutPerSeq(mps_name, geneIdentifier):
	
	mutperseq = {}
	
	for per_seq in mps_name:
		count = 0
		for base in per_seq:
			if base.isupper():
				count += 1
		if count in mutperseq:
			mutperseq[count] += 1
		else:
			mutperseq[count] = 1
	
	perpercent = {}
	total = sum(mutperseq.values())
	for entry, p in mutperseq.items():
		percent = p * 100 / total
		perpercent[entry] = percent
	
	perpercent = OrderedDict(sorted(perpercent.items()))
		
	plt.plot(perpercent.keys(), perpercent.values(), '-x')
	plt.title("Number of mutations per mutant DNA product")
	plt.xlabel("Number of mutations")
	plt.ylabel("Percentage of sequences")
	plt.grid()
	figureName = ('mutPerSeq_%sc%sr%s_%s.svg' % (total_cycles, runs, 
		error_rate, geneIdentifier))
	plt.savefig(figureName)
	plt.close()

	return figureName
	
		
def mutplot(baseDict, posDict, noposDict, geneIdentifier, name):
	
	temp = []
	for item, count in baseDict.items():
		temp.extend([item for i in range(count)])

	# count abundance of mutations at positions (mutPlot)
	for each in temp:
		for pos, base in enumerate(each):
			if base.isupper():
				dictCounter(1+pos, posDict)
			elif noposDict != None:
				dictCounter(1+pos, noposDict)

	plt.figure(figsize=(16,8))	
	plt.plot(posDict.keys(), posDict.values(), 'kD', ms = 1.5)
	if noposDict != None:
		plt.plot(noposDict.keys(), noposDict.values(), 'rD', ms = 1.5)
	plt.xlabel("Nucleotide position")
	plt.ylabel("Frequency of mutations")
	plt.title("Total mutations for each nucleotide position")
	figureName = ('mutPlotFILTERED%s_%sc%sr%s_%s.svg' % (name,
		total_cycles, runs, error_rate, geneIdentifier))
	plt.savefig(figureName)
	plt.close()

	return figureName
	

def runs_coverage(covDict, geneIdentifier):
	
	x = range(runs)
	count = 0
	plt.plot(x, covDict)
	
	for each in covDict:
		count += 1
		if each == 100:
			maxcov = each
			plt.plot(count, maxcov, 'kx')
			plt.axvline(x = count, color = 'k', ls = '--',
				label = 'Runs for 100% coverage')
			break
			
	plt.xlabel("Number of runs")
	plt.ylabel("Percentage coverage (%)")
	plt.title("Percentage coverage of mutants with one mutation")
	plt.grid()
	plt.axhline(y=100, color = 'red', ls = '--', alpha = 0.5,
		label = '100% coverage')
	plt.axhline(y=99, color = 'orange', ls = '--', alpha = 0.5,
		label = '99% coverage')
	plt.axhline(y=95, color = 'gold', ls = '--', alpha = 0.5,
		label = '95% coverage')
	figureName = ('coverage_%sc%sr%s_%s.svg' % (total_cycles, runs, 
		error_rate, geneIdentifier))
	plt.savefig(figureName)
	plt.close()
	
	return figureName


main() 
