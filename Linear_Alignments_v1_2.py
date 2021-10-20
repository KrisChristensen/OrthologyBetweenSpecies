##########################################################
### Import Necessary Modules

import argparse		               #provides options at the command line
import sys		               #take command line arguments and uses it in the script
import gzip		               #allows gzipped files to be read
import re		               #allows regular expressions to be used
import textwrap		               #allows the use of textwrapping for long sequences
from collections import defaultdict
import operator

"""
This script is meant to be called from another script.  It retains alignments or anything with the proper format that are linear in their progression and having a common orientation.  This version is made for long distance alignments.
"""
sort = {}
lengths = {}

def linearAlignments(input_string, maximum_scaff_diff, maximum_chrom_diff, minimum_length, minimum_total_length):      
	##################
	##################
	###First put into a hash to sort by the first distances and then by the second distances
	sort.clear()
        numAligns = 0
	for pos in input_string.split(","):
		psstart, psend, pcstart, pcend, scaf_len, orientation = pos.split("\t")[0:6]
                numAligns += 1
		if float(psstart) in sort:
		        if float(pcstart) in sort[float(psstart)]:
				sort[float(psstart)][float(pcstart)] += ",{}\t{}\t{}\t{}\t{}\t{}".format(psstart,psend,pcstart,pcend,scaf_len,orientation)
		        else:
				sort[float(psstart)][float(pcstart)] = "{}\t{}\t{}\t{}\t{}\t{}".format(psstart,psend,pcstart,pcend,scaf_len,orientation)
		else:
			sort[float(psstart)] = {}
			sort[float(psstart)][float(pcstart)] = "{}\t{}\t{}\t{}\t{}\t{}".format(psstart,psend,pcstart,pcend,scaf_len,orientation)
	sys.stderr.write("\tFound #Alignments: {}\n\n".format(numAligns))

	##################
	##################
	###Sort by first positions and then second position, then start creating linear alignments within a max distance set above
	###If an alignment could go to multiple linear alignments, it goes to the best, but if there are two "best," then it doesn't go to any linear alignments.
	len_count = 1
	final_list = "NA"
	lengths.clear()
	len_count = linear(False, len_count, maximum_scaff_diff, maximum_chrom_diff, "+")
	sys.stderr.write("\tFound #Pos Alignments: {}\n".format(len(lengths.keys())))
	final_list = finallist(final_list, minimum_length, minimum_total_length)
	lengths.clear()
	len_count = linear(True, len_count, maximum_scaff_diff, maximum_chrom_diff, "-")
	sys.stderr.write("\tFound #Neg Alignments: {}\n".format(len(lengths.keys())))
	final_list = finallist(final_list, minimum_length, minimum_total_length)
	###Finally the positions of the filtered alignments are returned
	final_positions = []
	for final in final_list.split(","):
		final_position = 0
		for t1 in input_string.split(","):
			##################################
			psstart, psend, pcstart, pcend, scaf_len, orientation = t1.split("\t")[0:6]
			t1 = "{}\t{}\t{}\t{}\t{}\t{}".format(psstart,psend,pcstart,pcend,scaf_len,orientation)
			##################################
			if str(final) == str(t1):
				final_positions.append(final_position)
			final_position += 1
	if len(final_positions) > 0:
		return(final_positions)
	else:
		return("")

def linear(d, l, msd, mcd, ori):
	direction = d
	len_count = l
	maximum_scaff_diff = msd
	maximum_chrom_diff = mcd
        main_orientation = ori
	for scaffoldStart in sorted(sort, key=float):
		for chromStart in sorted(sort[scaffoldStart], key=float, reverse=direction):
			for mul in sort[scaffoldStart][chromStart].split(","):
				(tpsstart,tpsend,tpcstart,tpcend,tscaf_len,torient) = mul.split("\t")
				#sys.stderr.write("\tposition: ss:{} se:{} cs:{} ce:{} sl:{}\n".format(tpsstart, tpsend, tpcstart, tpcend, tscaf_len))
				best_distance = "NA"
				best_count = 0
				best_position = "NA"
				if main_orientation == torient:
					for ls in lengths:
						all_lengths = lengths[ls].split(",")
						(lpsstart,lpsend,lpcstart,lpcend,lpscaff_len,lporient) = all_lengths[-1].split("\t")
						distance_between_scaf = float(tpsstart) - float(lpsend)
						distance_between_chr = float(tpcstart) - float(lpcend)
						if direction: ###Reverse direction/orientation
							 distance_between_chr = float(lpcstart) - float(tpcend)
						if float(distance_between_scaf) >= 0 and float(distance_between_chr) >= 0 and float(distance_between_scaf) <= float(maximum_scaff_diff) and float(distance_between_chr) <= float(maximum_chrom_diff) and lporient == torient:
							if best_distance == "NA":
								best_distance = float(distance_between_scaf) + float(distance_between_chr) 
								best_count = 1
								best_position = ls
							else:
								if float(best_distance) > float(distance_between_scaf) + float(distance_between_chr):
									best_distance = float(distance_between_scaf) + float(distance_between_chr)
									best_count = 1
									best_position = ls
								elif float(best_distance) == float(distance_between_scaf) + float(distance_between_chr):
									best_count += 1
					if int(best_count) == 1:
						lengths[best_position] += ",{}\t{}\t{}\t{}\t{}\t{}".format(tpsstart,tpsend,tpcstart,tpcend,tscaf_len,torient)
					
					else:
						lengths[len_count] = "{}\t{}\t{}\t{}\t{}\t{}".format(tpsstart,tpsend,tpcstart,tpcend,tscaf_len,torient)
						len_count += 1
	return(len_count)

def finallist(f, m, mtl):
	final_list = f
	minimum_length = m
	minimum_total_length = mtl
	total_kept = 0
	###All of the linear alignments are checked for a minimum length (based on the first sequence, in my case the scaffold vs the genome, the scaffold)
	for ls in lengths:
		total_length = 0
		total_nucleotides = 0
		length_start = "NA"
		length_end = "NA"
		number_of_alignments = 0
		for tl in lengths[ls].split(","):
			(tlpsstart,tlpsend,tlpcstart,tlpcend,tlscaf_len,tlorient) = tl.split("\t")
			number_of_alignments += 1
			if length_start == "NA":
				length_start = tlpsstart
				length_end = tlpsend
			else:
				if float(length_start) > float(tlpsstart):
					length_start = tlpsstart
				if float(length_end) < float(tlpsend):
					length_end = tlpsend
			total_nucleotides += int(tlscaf_len)
		total_length = float(length_end) - float(length_start)
		#sys.stderr.write("\t\t\t#Alignments: {}\tStart: {}\tFinish: {}\tNucleotides: {}\n".format(number_of_alignments, length_start, length_end, total_nucleotides))
		if float(total_length) >= float(minimum_length) and int(total_nucleotides) >= int(minimum_total_length):
			total_kept += 1
			if final_list == "NA":
				final_list = lengths[ls]
			else:
				final_list += ",{}".format(lengths[ls])
	sys.stderr.write("\t\tKept: {}\n".format(total_kept))
	return(final_list)


if __name__ == '__main__':
	test = "0.1	0.2	0.1	0.2	12	-,0.1	0.2	0.2	0.3	8	+,0.2	0.3	0.3	0.4	9	+,0.2	0.2	0.5	0.5	11	-,0.3	0.4	0.4	0.5	22	+"
	smax_difference = 0.1
	cmax_difference = 0.1
	min_length = 0.2
	min_tot_len = 0

	output = linearAlignments(test, smax_difference, cmax_difference, min_length, min_tot_len)	
	print "{}".format(output)	
