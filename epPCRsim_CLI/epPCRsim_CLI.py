# Imports
import sys
import getopt
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
	global log_file
	global output
	global output_name
	global figures
	global output_cov_graph
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
	master_mut_dict = {}		# Dictionary of all sequences, rather than lists
	filter_mut_dict = {}


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
	output_cov_graph = True
	


	##############################
	### COMMAND LINE ARGUMENTS ###
	##############################

	argv = sys.argv[1:]

	if len(sys.argv) <= 7:
		title_screen()
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
		title_screen()
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
			output_cov_graph = False
		elif opt in ('-y', '--no-coverage'):
			output_cov_graph = False
		elif opt in ('-s', '--snv'):
			SNVc = Decimal(arg) 
		elif opt in ('-d', '--deletion'):
			DELc = Decimal(arg) 
		elif opt in ('-n', '--insertion'):
			INSc = Decimal(arg) 
	
	
	#Write mutlog file#
	log_file_name = ('mutLog%s%s%s%s.txt' % (
			identifier,total_cycles,runs,error_rate))
	if log_output == True:
		log_file = open(log_file_name, 'w')

	#Write fasta file
	if output == True:
		mut_fasta = open(output_name, 'w')



	########################
	### BEGIN SIMULATION ###
	########################

	# Error if mut probabilities not = 1
	total = [SNVc, INSc, DELc]
	if sum(total) != 1:
		print("Probability does not equal 1 (100%), adjust accordingly")
		return

	running_screen()
	run_iterator = 1
	while run_iterator <= runs:
		
		log_append("Run: %s" % (run_iterator))
		final_DNA, final_muts = cycle(sequence)	

		# Store all final mutant sequences (final_muts) in dictionary w/ counter
		# as needed for the ALL logoplot - master_mut_dict
		for item in final_muts:
			dict_counter(item, master_mut_dict)

		# Filter out all mutants != 1 mutation after each run & store in dictionary
		# w/ counter again as needed for FILTERED logoplot - filteredMuts
		if output_cov_graph == True:
			filter_mut_dict, coverage = filtration(master_mut_dict, filter_mut_dict, coverage, sequence)
			mutcoverage = calculating_coverage(coverage, mutcoverage, sequence)

		# If output is enabled, write the DNA at the end of each run to FASTA file
		mutant = 1
		if output == True:
			fasta_save(final_DNA, identifier, mut_fasta, mutant)
			

		run_iterator = run_iterator + 1
	
	
	#This gives end coverage rather than per run - faster
	if output_cov_graph == False:
		filter_mut_dict, coverage = filtration(master_mut_dict, filter_mut_dict, coverage, sequence)
		mutcoverage = calculating_coverage(coverage, mutcoverage, sequence)
		cov_name = None
	
	if output_cov_graph == True:
		cov_name = runs_coverage(mutcoverage, identifier)

	
	#Figures (default is true)
	if figures == True:
		mps_name = mut_per_seq(final_muts, identifier)	
		pmc_name = point_mutcount(identifier)
		lp_all_name = logoplot(master_mut_dict, identifier, "ALL")
		lp_filt_name = logoplot(filter_mut_dict, identifier, "FILTERED")
		mp_muts_name = mutplot(filter_mut_dict, mutpositions, None, identifier, "muts")
		mp_both_name = mutplot(filter_mut_dict, mutpositions, nopositions, identifier, "both")
		if mutcoverage[-1] != 100:
			uncovered_list = inverting_logo(coverage, sequence, bases)
			lp_inv_name = logoplot(uncovered_list, identifier, "FILTERED_INVERTED")
	else:
		mps_name = None
		pmc_name = None
		lp_all_name = None
		lp_filt_name = None
		mp_muts_name = None
		mp_both_name = None
		lp_inv_name = None


	#Add % coverage to optional log file
	length_master = dict_lens(master_mut_dict)
	log_append("\nTotal mutants: %s" % (length_master))
	length_filter = dict_lens(filter_mut_dict)
	log_append("Total filtered mutants (with one mutation): %s" % ((length_filter)))
	log_append("%s cycles + %s runs gives %.2f%% coverage" % (
		total_cycles, runs, mutcoverage[-1])) 
	
	output_screen(length_master, length_filter, mutcoverage, mps_name, pmc_name, lp_all_name, lp_filt_name, 
		mp_muts_name, mp_both_name, lp_inv_name, cov_name, log_file_name)
	

def cycle(seq_DNA):
	
	all_seqs = []		# all products from each cycle (NOT INCLUDING START SEQ)
	all_seqs_mut = []	# all products only showing mutations

	temp_list1 = []		
	temp_list2 = []	
	mut_temp_list1 = []
	mut_temp_list2 = []
	
	output_DNA = []		# Final sequences from each run (after total_cycles) 
	output_muts = []	# Final sequences with only mutant bases and fullstops

	cycle_no = 1	# iterator for cycle no.
	
	temp_list1.append(seq_DNA)
	
	new = "." * len(seq_DNA)
	mut_temp_list1.append(new)
	
	while cycle_no <= total_cycles:
		log_append("Cycle: %s" % (cycle_no))
		
		for each, mut_each in zip(temp_list1, mut_temp_list1):	
			
			copy = each
			mut_copy = mut_each
			
			all_seqs.append(copy)
			all_seqs_mut.append(mut_copy)
			
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
	

def mutation(insert, mut_seq, pol_biases, all_seqs, all_seqs_muts):
		
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
				point = random.choices(bases, weights=pol_biases)
				point = str(point[0])
				while True:			
					if point == nucleotide:
						point = random.choices(bases, weights=pol_biases)
						point = str(point[0])
					else:
						break

				changed = nucleotide + " to " + point
				dict_counter(changed, mutdict)	
				log_append("Point mutation %s at %s" % (changed, n-1))	
				insert = insert[:n-1] + point + insert[n:] 
				mut_seq = mut_seq[:n-1] + point + mut_seq[n:] 
				
			# Deletion mutation	
			elif DELlower <= m <= DELupper:
				log_append("Deletion at %s" % (n-1))
				insert = insert[:n-1] + insert[n:]
				mut_seq = mut_seq[:n-1] + mut_seq[n:]
				gene_length -= 1
				
			# Insertion mutation	
			elif INSlower <= m <= INSupper:
				ins = random.choices(bases, weights=default)
				ins = str(ins[0])
				befaft = random.randint(1,2)
				if befaft == 1:		# insert before n
					insert = insert[:n-1] + ins + insert[n-1:]
					mut_seq = mut_seq[:n-1] + ins + mut_seq[n-1:]
					log_append("Insertion of %s before %s" % (ins, n-1))
				elif befaft == 2:	# insert after n
					insert = insert[:n] + ins + insert[n:]
					mut_seq = mut_seq[:n] + ins + mut_seq[n:]
					log_append("Insertion of %s after %s" % (ins, n-1))	 
				gene_length += 1
	
		n += 1			
	all_seqs.append(insert)
	all_seqs_muts.append(mut_seq)			


def fasta_save(selection, mut_name, file_name, mut_no):
	
	wrapping = 60
	
	if output == True:
		for each in selection:
	
			name = (">%s_M%s" % (mut_name, mut_no))
			file_name.write(name + "\n")
		
			for x in range(0, len(each), wrapping):
				file_name.write(each[x:x + wrapping] + "\n")
			mut_no += 1

def log_append(text):
	
	if log_output == True:
		log_file.write(text + "\n")
			

def dict_lens(dict_to_use):
	
	temp_list = []
	for item, count in dict_to_use.items():
		temp_list.extend([item for i in range(count)])

	length = len(temp_list)

	return length
	

def filtration(current_dict, future_dict, loc_dict, insert):

	# Going through masterDict and appending to filtered dictionary
	future_dict = {}
	loc_dict = {}
	for item, count in current_dict.items():
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
			future_dict[item] = count

			if 1+pos not in loc_dict:
				loc_dict[1+pos] = list()
			if base not in loc_dict[1+pos]:
				loc_dict[1+pos].append(mut)

	return future_dict, loc_dict
	
	
def calculating_coverage(loc_dict, running_cov, insert):

	# Calculate cov each run
	count = 0
	for each, each2 in loc_dict.items():
		count += len(each2)
	
	cov_all = len(insert) * 3
	cov = (count/cov_all) * 100
	running_cov.append(cov)

	return running_cov


def inverting_logo(dict_of_base_cov, insert, nucleotides):

    # calculating every possible base change
    total_cov = {}
    for loc, b in enumerate(insert):
        total_cov[1+loc] = list()
        for each in nucleotides:
            if each != b:
                total_cov[1+loc].append(each)

    # Dictionary of bases not changed to    
    uncovered = {}
    for (a, b) in dict_of_base_cov.items():
        for (x, y) in total_cov.items():
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


def title_screen():

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


def running_screen():
	
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
	print("\tProduce coverage graph:\t%s" % (output_cov_graph))
	print("\tOutput fasta:\t\t%s" % (output))
	

def output_screen(totalMuts, oneMutMuts, perCov, mut_per_seqName, pointcountName, logoALL, logoFILT, 
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
		print("\tMutations per seq:\t%s" % (mut_per_seqName))
		print("\tSNV count:\t\t%s" % (pointcountName))
		print("\tMutation plots:\t\t%s" % (mutplotmuts))
		print("\t\t\t\t%s" % (mutplotboth))
		print("\tLogo plots:\t\t%s" % (logoALL))
		print("\t\t\t\t%s" % (logoFILT))
	if perCov[-1] != 100 and logoINV != None:
		print("\t\t\t\t%s" % (logoINV))
	if output_cov_graph == True:
		print("\tCoverage:\t\t%s" % (covg))
	if log_output == True:
		print("\tLog file:\t\t%s" % (logF))
	if (output == False and figures == False and 
		output_cov_graph == False and log_output == False):
		print("\tNone")
	print("")	

	
def dict_counter(dict_item, dict_name):
	
	# if not in dict then add and make value = 0
	# if is in dict then value will increase by 1
	
	if dict_item not in dict_name.keys():
		dict_name[dict_item] = 0
	dict_name[dict_item] += 1


def logoplot(input_mutants, gene_identifier, name):

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
			gene_identifier, name, length, length+100, 
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
					gene_identifier, length, (length+100)))
		
		plt.savefig(im_name)	
		#plt.show()
		plt.close()
		
		length = length + 100
	
	# Compile all logos of 100 bases
	o = [PIL.Image.open(x) for x in im_list]
	min_shape = sorted([(np.sum(x.size), x.size) for x in o])[0][1]
	ims_all = np.vstack([np.asarray(x.resize(min_shape)) for x in o])
	ims_all = PIL.Image.fromarray(ims_all)
	figure_name = ('logoCompiled%s_%sc%sr%s_%s.jpg'
		% (name, total_cycles, runs, error_rate, gene_identifier))
	ims_all.save(figure_name)
	
	# Delete individual logos after compiling
	for individual_logo in im_list:
		os.remove(individual_logo)

	return figure_name

def point_mutcount(gene_identifier):
	
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
	figure_name = ('SNVcount_%sc%sr%s_%s.svg' % (total_cycles, runs, 
		error_rate, gene_identifier))
	plt.savefig(figure_name)
	plt.close()
	
	return figure_name
	
def mut_per_seq(mps_name, gene_identifier):
	
	mut_per_seq = {}
	
	for per_seq in mps_name:
		count = 0
		for base in per_seq:
			if base.isupper():
				count += 1
		if count in mut_per_seq:
			mut_per_seq[count] += 1
		else:
			mut_per_seq[count] = 1
	
	perpercent = {}
	total = sum(mut_per_seq.values())
	for entry, p in mut_per_seq.items():
		percent = p * 100 / total
		perpercent[entry] = percent
	
	perpercent = OrderedDict(sorted(perpercent.items()))
		
	plt.plot(perpercent.keys(), perpercent.values(), '-x')
	plt.title("Number of mutations per mutant DNA product")
	plt.xlabel("Number of mutations")
	plt.ylabel("Percentage of sequences")
	plt.grid()
	figure_name = ('mut_per_seq_%sc%sr%s_%s.svg' % (total_cycles, runs, 
		error_rate, gene_identifier))
	plt.savefig(figure_name)
	plt.close()

	return figure_name
	
		
def mutplot(base_dict, pos_dict, nopos_dict, gene_identifier, name):
	
	temp = []
	for item, count in base_dict.items():
		temp.extend([item for i in range(count)])

	# count abundance of mutations at positions (mutPlot)
	for each in temp:
		for pos, base in enumerate(each):
			if base.isupper():
				dict_counter(1+pos, pos_dict)
			elif nopos_dict != None:
				dict_counter(1+pos, nopos_dict)

	plt.figure(figsize=(16,8))	
	plt.plot(pos_dict.keys(), pos_dict.values(), 'kD', ms = 1.5)
	if nopos_dict != None:
		plt.plot(nopos_dict.keys(), nopos_dict.values(), 'rD', ms = 1.5)
	plt.xlabel("Nucleotide position")
	plt.ylabel("Frequency of mutations")
	plt.title("Total mutations for each nucleotide position")
	figure_name = ('mutPlotFILTERED%s_%sc%sr%s_%s.svg' % (name,
		total_cycles, runs, error_rate, gene_identifier))
	plt.savefig(figure_name)
	plt.close()

	return figure_name
	

def runs_coverage(cov_dict, gene_identifier):
	
	x = range(runs)
	count = 0
	plt.plot(x, cov_dict)
	
	for each in cov_dict:
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
	figure_name = ('coverage_%sc%sr%s_%s.svg' % (total_cycles, runs, 
		error_rate, gene_identifier))
	plt.savefig(figure_name)
	plt.close()
	
	return figure_name


main() 
