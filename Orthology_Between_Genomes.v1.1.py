##########################################################
### Import Necessary Modules

import argparse		               #provides options at the command line
import sys		               #take command line arguments and uses it in the script
import gzip		               #allows gzipped files to be read
import re		               #allows regular expressions to be used

##########################################################
### Command-line Arguments
parser = argparse.ArgumentParser(description="A script to find orthologous genes between genomes")
parser.add_argument("-gffs", help = "The location of the gff file for the subject (database) species", default=sys.stdin, required=True)
parser.add_argument("-gffq", help = "The location of the gff file for the query species", default=sys.stdin, required=True)
parser.add_argument("-alignn", help = "The location of the alignment file between the species (nucleotide)--should be filtered already", default=sys.stdin, required=True)
parser.add_argument("-alignp", help = "The location of the alignment file between the species (protein)--should be filtered already", default=sys.stdin, required=True)
parser.add_argument("-fais", help = "The location of the index for the subject fasta file", default=sys.stdin, required=True)
parser.add_argument("-faiq", help = "The location of the index for the query fasta file", default=sys.stdin, required=True)
args = parser.parse_args()

#########################################################
### Open file (object-oriented programming)

class GlobalVar():
	###Makes variables accessible to all functions###
	gene_position = {}
	gene_position["subject"] = {}
	gene_position["query"] = {}
	protein_position = {}
	protein_position["subject"] = {}
	protein_position["query"] = {}
	protein_name_2_gene_name = {}
	nucleotide_alignments1_start = {}
	nucleotide_alignments1_end = {}
	nucleotide_alignments2_start = {}
	nucleotide_alignments2_end = {}
	chromosome_lengths = {}
	protein_alignments = {}
	orthologs = {}
	ortholog_conflict = {}

class OpenFile():
    	### Opens the file and directs it to a line-reader and stores information into variables used later
	def __init__ (self, filename, ty, ty2):
		"""Opens the file (accepts gzipped files)"""
		if re.search(".gz$", filename):
			self.file = gzip.open(filename, 'rb')
		else:
			self.file = open(filename, 'r')
		if ty == "gff":
			sys.stderr.write("Reading gff {} file\n".format(ty2))
			self.readLinesGFF(self.file, ty2)
		elif ty == "aln":
			sys.stderr.write("Reading alignment ({}) file\n".format(ty2))          
			self.readLinesAln(self.file, ty2)
		elif ty == "fai":
			sys.stderr.write("Reading fai ({}) file\n".format(ty2))          
			self.readLinesFai(self.file, ty2)			

	def readLinesGFF(self, f, ty2):
		"""Reads the lines from gff file and outputs them to global variable"""
		self.rna_name_2_gene_name = {}
		self.filename = f
		self.type = ty2
		self.gene_count = 0
		self.linked_proteins = 0
		for line in self.filename:
			line = line.rstrip('\n')
			if not re.search("^#", line):
				self.list_of_columns = line.split("\t")
				###Finds the link between gene id and protein name (it can be tricky as unique gene names may be the same between species)
				if re.search("^gene", self.list_of_columns[2]) and re.search("protein_coding", self.list_of_columns[8]):
					self.gene_count += 1
					self.chrom = self.list_of_columns[0]
					self.start = self.list_of_columns[3]
					self.end = self.list_of_columns[4]
					self.gene_name = "NA"
					self.gene_id = "NA"
					self.gene_unique = "NA"
					###The GFF file does not contain uniform column 8's########
					###########################################################
					for self.info in self.list_of_columns[8].split(";"):
						if re.search("^ID=gene", self.info):
							self.gene_unique = self.info[3:]
						elif re.search("^Dbxref=GeneID:", self.info):
							self.gene_id = self.info[14:]
						elif re.search("^Name=", self.info):
							self.gene_name = self.info[5:]
					###########################################################
					if self.gene_name == "NA" or self.gene_id == "NA" or self.gene_unique == "NA":
						sys.stderr.write("Warning, gene at Chr {}, Pos {} {} not named correctly\n".format(self.chrom, self.start, self.end))
					if self.gene_unique in GlobalVar.gene_position[self.type]:
						sys.stderr.write("Warning: gu:{}, gi:{}, gene:{} found more than once\n".format(self.gene_unique, self.gene_id, self.gene_name))
					else:
						GlobalVar.gene_position[self.type][self.gene_unique] = "{}\t{}\t{}\t{}\t{}".format(self.gene_id, self.gene_name, self.chrom, self.start, self.end)
				###The mRNA id is needed to link the protein id to the gene id
				elif re.search("^mRNA", self.list_of_columns[2]):
					self.rna_id = "NA"
					self.gene_unique = "NA"
					###########################################################
					for self.info in self.list_of_columns[8].split(";"):
						if re.search("^ID=rna", self.info):
							self.rna_id = self.info[3:]
						elif re.search("^Parent=gene", self.info):
							self.gene_unique = self.info[7:]
					###########################################################
					if self.rna_id == "NA" or self.gene_unique == "NA":
						sys.stderr.write("Warning, RNA at Chr {}, Pos {} {} not named correctly\n".format(self.list_of_columns[0], self.list_of_columns[3], self.list_of_columns[4]))
					if self.rna_id in self.rna_name_2_gene_name:
						sys.stderr.write("Warning: RNA_id: {} found more than once\n".format(self.rna_id))
						sys.stderr.write("Line: {}\n".format(line))
					else:
						self.rna_name_2_gene_name[self.rna_id] = self.gene_unique
				###The CDS information contains the protein id
				elif re.search("^CDS", self.list_of_columns[2]):
					self.parent_rna = "NA"
					self.gene_id = "NA"
					self.protein_name = "NA"
					self.gene_name = "NA"
					self.product_name = "NA"
					self.gene_unique = "NA"
					###########################################################
					for self.info in self.list_of_columns[8].split(";"):
						if re.search("^Parent=rna", self.info):
							self.parent_rna = self.info[7:]
						elif re.search("^Dbxref=GeneID:", self.info):
							self.gene_id = self.info.split(",")[0][14:]
						elif re.search("^Name=", self.info):
							self.protein_name = self.info[5:]
						elif re.search("^gene=", self.info):
							self.gene_name = self.info[5:]
						elif re.search("^product=", self.info):
							self.product_name = self.info[8:]					
					###########################################################
					if self.parent_rna in self.rna_name_2_gene_name:
						self.gene_unique = self.rna_name_2_gene_name[self.parent_rna]
					if (self.parent_rna == "NA" or self.gene_id == "NA" or self.protein_name == "NA" or
						self.gene_name == "NA" or self.product_name == "NA" or self.gene_unique == "NA"):
						sys.stderr.write("Warning, Protein at Chr {}, Pos {} {} not named correctly\n".format(self.list_of_columns[0], self.list_of_columns[3], self.list_of_columns[4]))
					if self.protein_name not in GlobalVar.protein_name_2_gene_name and self.gene_unique in GlobalVar.gene_position[self.type]:
						self.linked_proteins += 1
						self.gi, self.gn, self.chr, self.str, self.en = GlobalVar.gene_position[self.type][self.gene_unique].split("\t")
						GlobalVar.protein_name_2_gene_name[self.protein_name] = "{}\t{}\t{}\t{}\t{}\t{}\t{}".format(self.gene_unique, self.gene_id, self.gene_name, self.product_name, self.chr, self.str, self.en)
						###Puts the protein in order on chromosome
						if self.chr in GlobalVar.protein_position[self.type]:
							if self.str in GlobalVar.protein_position[self.type][self.chr]:
								GlobalVar.protein_position[self.type][self.chr][self.str] += ",{}".format(self.protein_name)
							else:
								GlobalVar.protein_position[self.type][self.chr][self.str] = "{}".format(self.protein_name)
						else:
							GlobalVar.protein_position[self.type][self.chr] = {}
							GlobalVar.protein_position[self.type][self.chr][self.str] = "{}".format(self.protein_name)


		sys.stderr.write("\tFound {} genes in gff file\n\tFound {} linked proteins\n\n".format(self.gene_count, self.linked_proteins))
		self.filename.close() 

	def readLinesAln(self, f, ty2):
		"""Reads alignment files and outputs them to a global variable"""
		self.filename = f
		self.type = ty2
		self.subject_proteins = 0
		self.protein_aligns = 0
		self.nucleotide_aligns = 0
		if self.type == "nucleotide":
			#GlobalVar.nucleotide_alignments_start/end = {}
			####Finds the minimum and maximum positions for orthologous chromosomes for both chromosomes
			self.header = self.filename.readline()
			for line in self.filename:
				line = line.rstrip('\n')
				self.nucleotide_aligns += 1
				self.query, self.subject, self.qstart, self.qend, self.sstart, self.send, self.length, self.orientation = line.split("\t")
				if str(self.query) != str(self.subject): ###The same chromosome is not allowed to align to itself (if doing homeologous analysis)
					####Chromosome 1###
					if self.query in GlobalVar.nucleotide_alignments1_start:
						if self.subject in GlobalVar.nucleotide_alignments1_start[self.query]:
							if float(GlobalVar.nucleotide_alignments1_start[self.query][self.subject]) > float(self.qstart):
								GlobalVar.nucleotide_alignments1_start[self.query][self.subject] = float(self.qstart)
						else:
							GlobalVar.nucleotide_alignments1_start[self.query][self.subject] = float(self.qstart)
						if self.subject in GlobalVar.nucleotide_alignments1_end[self.query]:
							if float(GlobalVar.nucleotide_alignments1_end[self.query][self.subject]) < float(self.qend):
								GlobalVar.nucleotide_alignments1_end[self.query][self.subject] = float(self.qend)
						else:
							GlobalVar.nucleotide_alignments1_end[self.query][self.subject] = float(self.qend)
					else:
						GlobalVar.nucleotide_alignments1_start[self.query] = {}
						GlobalVar.nucleotide_alignments1_end[self.query] = {}
						GlobalVar.nucleotide_alignments1_start[self.query][self.subject] = float(self.qstart)
						GlobalVar.nucleotide_alignments1_end[self.query][self.subject] = float(self.qend)
					####Chromosome 2###
					if self.subject in GlobalVar.nucleotide_alignments2_start:
						if self.query in GlobalVar.nucleotide_alignments2_start[self.subject]:
							if float(GlobalVar.nucleotide_alignments2_start[self.subject][self.query]) > float(self.sstart):
								GlobalVar.nucleotide_alignments2_start[self.subject][self.query] = float(self.sstart)
						else:
							GlobalVar.nucleotide_alignments2_start[self.subject][self.query] = float(self.sstart)
						if self.query in GlobalVar.nucleotide_alignments2_end[self.subject]:
							if float(GlobalVar.nucleotide_alignments2_end[self.subject][self.query]) < float(self.send):
								GlobalVar.nucleotide_alignments2_end[self.subject][self.query] = float(self.send)
						else:
							GlobalVar.nucleotide_alignments2_end[self.subject][self.query] = float(self.send)
					else:
						GlobalVar.nucleotide_alignments2_start[self.subject] = {}
						GlobalVar.nucleotide_alignments2_end[self.subject] = {}
						GlobalVar.nucleotide_alignments2_start[self.subject][self.query] = float(self.sstart)
						GlobalVar.nucleotide_alignments2_end[self.subject][self.query] = float(self.send)
			sys.stderr.write("\tFound {} alignments\n\n".format(self.nucleotide_aligns))
		elif self.type == "protein":
			#GlobalVar.protein_alignments = {}
			####Finds protein matches passing filtering
			for line in self.filename:
				line = line.rstrip('\n')
				self.query,self.subject,self.per_id,self.aln_len,self.mismatch,self.gap,self.qstart,self.qend,self.sstart,self.send,self.evalue,self.bit = line.split("\t")
				#if str(self.query) != str(self.subject): ###The same protein is not allowed to align to itself (if doing homeologous analysis)
				self.protein_aligns += 1
				if self.subject in GlobalVar.protein_alignments:
					GlobalVar.protein_alignments[self.subject][self.query] = "NA"
				else:
					self.subject_proteins += 1
					GlobalVar.protein_alignments[self.subject] = {}
					GlobalVar.protein_alignments[self.subject][self.query] = "NA"
			sys.stderr.write("\tFound {} subject proteins and {} alignments\n\n".format(self.subject_proteins, self.protein_aligns))
		self.filename.close()

	def readLinesFai(self, f, ty2):
		"""Reads the lines from the fai files and outputs them to a global variable"""
		#GlobalVar.chromosome_lengths = {}
		self.filename = f
		self.type = ty2
		self.seq_count = 0
		self.total_length = 0
		for line in self.filename:
			line = line.rstrip('\n')
			self.chr, self.length, self.un1, self.un2, self.un3 = line.split("\t")
			GlobalVar.chromosome_lengths[self.chr] = self.length
			self.seq_count += 1
			self.total_length += int(self.length)
		sys.stderr.write("\tFound {} sequences, with a total length of {} nucleotides\n\n".format(self.seq_count, self.total_length))

class Analysis():
    	### Begins the analysis ###
	def __init__ (self):
		"""This analysis identifies the orthologous protein alignments"""
		###First orthologous chromosomes are identified
		self.already_used = {}
		for self.chromosome1 in GlobalVar.nucleotide_alignments2_start:
			for self.chromosome2 in GlobalVar.nucleotide_alignments2_start[self.chromosome1]:
				self.chr1le = GlobalVar.chromosome_lengths[self.chromosome1]
				self.chr1st = int(int(self.chr1le)*float(GlobalVar.nucleotide_alignments2_start[self.chromosome1][self.chromosome2]))
				self.chr1en = int(int(self.chr1le)*float(GlobalVar.nucleotide_alignments2_end[self.chromosome1][self.chromosome2]))
				self.chr2le = GlobalVar.chromosome_lengths[self.chromosome2]
				self.chr2st = int(int(self.chr2le)*float(GlobalVar.nucleotide_alignments1_start[self.chromosome2][self.chromosome1]))
				self.chr2en = int(int(self.chr2le)*float(GlobalVar.nucleotide_alignments1_end[self.chromosome2][self.chromosome1]))
				self.window = self.chr1st
				self.list_of_proteins1 = {}
				###Finds all of the proteins in the orthologous region for the subject chromosome
				for self.position in GlobalVar.protein_position["subject"][self.chromosome1]:
					if int(self.position) >= int(self.chr1st) and int(self.position) <= int(self.chr1en):
						for self.protein in GlobalVar.protein_position["subject"][self.chromosome1][self.position].split(","):
							self.list_of_proteins1[self.protein] = 0
				###Finds all of the proteins in the orthologous region for the query chromosome
				self.list_of_proteins2 = {}
				for self.position in GlobalVar.protein_position["query"][self.chromosome2]:
					if int(self.position) >= int(self.chr2st) and int(self.position) <= int(self.chr2en):
						for self.protein in GlobalVar.protein_position["query"][self.chromosome2][self.position].split(","):
							self.list_of_proteins2[self.protein] = 0
				###Finds the potential othologous protein matches (and change them to genes)
				#GlobalVar.protein_name_2_gene_name[self.protein_name] = "{}\t{}\t{}\t{}\t{}\t{}\t{}".
				#format(self.gene_unique, self.gene_id, self.gene_name, self.product_name, self.chr, self.str, self.en)
				self.already_used.clear()
				for self.protein1 in self.list_of_proteins1:
					for self.protein2 in self.list_of_proteins2:
						###These are all of the proteins in the orthologous regions and we check to see if they align as well###
						###Check to see if proteins align and the alignment passed the thresholds###
						if self.protein1 in GlobalVar.protein_alignments and self.protein2 in GlobalVar.protein_alignments[self.protein1]:
							self.gene_u1, self.gene_id1, self.gene_name1, self.product_name1, self.chr1, self.str1, self.en1 =  GlobalVar.protein_name_2_gene_name[self.protein1].split("\t")
							self.gene_u2, self.gene_id2, self.gene_name2, self.product_name2, self.chr2, self.str2, self.en2 =  GlobalVar.protein_name_2_gene_name[self.protein2].split("\t")
							####The gff file has multiple gene names that are not unique so I used the ID, but this is shared between gff ###
							####s was added for the subject and q for query to make them unique between###
							self.gene_name_used1 = "s" + self.gene_u1
							self.gene_name_used2 = "q" + self.gene_u2
							if self.gene_name_used1 in self.already_used and str(self.gene_name_used2) != str(self.already_used[self.gene_name_used1]):
								GlobalVar.ortholog_conflict[self.gene_name_used1] = 1
							elif self.gene_name_used2 in self.already_used and str(self.gene_name_used1) != str(self.already_used[self.gene_name_used2]):
								GlobalVar.ortholog_conflict[self.gene_name_used1] = 1
							else:
								GlobalVar.orthologs[self.gene_name_used1] = self.gene_name_used2
								self.already_used[self.gene_name_used1] = self.gene_name_used2
								self.already_used[self.gene_name_used2] = self.gene_name_used1

class Output():
    	### Outputs the results ###
	def __init__ (self):
		"""This outputs the proteins and corresponding orthologs in order"""
		print "Subject.Chromosome\tSubject.Start\tSubject.End\tSubject.Gene\tQuery.Chromosome\tQuery.Start\tQuery.End\tQuery.Gene"
		self.already_printed = {}
		for self.chromosome in sorted(GlobalVar.protein_position["subject"]):
			for self.position in sorted(GlobalVar.protein_position["subject"][self.chromosome].iterkeys(), key=int):
				for self.protein in GlobalVar.protein_position["subject"][self.chromosome][self.position].split(","):
					###################################################################################
					self.gene_u1, self.gene_id1, self.gene_name1, self.product_name1, self.chr1, self.str1, self.en1 =  GlobalVar.protein_name_2_gene_name[self.protein].split("\t")
					self.gene_name_used1 = "s" + self.gene_u1
					if self.gene_name_used1 not in self.already_printed:
						if self.gene_name_used1 in GlobalVar.orthologs and self.gene_name_used1 not in GlobalVar.ortholog_conflict: 
							self.gene_name_used2 = GlobalVar.orthologs[self.gene_name_used1]
							self.gene_u2 = self.gene_name_used2[1:]
							self.gid2, self.gn2, self.chr2, self.str2, self.end2 = GlobalVar.gene_position["query"][self.gene_u2].split("\t")
							print "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(self.chr1,self.str1,self.en1,self.gene_name1,self.chr2,self.str2,self.end2,self.gn2)
						else:
							print "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(self.chr1,self.str1,self.en1,self.gene_name1,"NA","NA","NA","NA")
						self.already_printed[self.gene_name_used1] = 1


if __name__ == '__main__':
	GlobalVar()
	OpenFile(args.gffs, "gff", "subject")
	OpenFile(args.gffq, "gff", "query")
	OpenFile(args.alignn, "aln", "nucleotide")
	OpenFile(args.alignp, "aln", "protein")
	OpenFile(args.fais, "fai", "subject")
	OpenFile(args.faiq, "fai", "query")
	Analysis()
	Output()
