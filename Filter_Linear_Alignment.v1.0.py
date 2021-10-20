##########################################################
### Import Necessary Modules

import argparse		               #provides options at the command line
import sys		               #take command line arguments and uses it in the script
import gzip		               #allows gzipped files to be read
import re		               #allows regular expressions to be used
import GeneralOverlap	      	       #allows overlaps to be identified

##########################################################
### Command-line Arguments
parser = argparse.ArgumentParser(description="A script to filter overlapping linear alignments, retaining only the best (and sections that do not overlap)")
parser.add_argument("-align", help = "The location of the alignment file that has already been linearly filtered", default=sys.stdin, required=True)
args = parser.parse_args()

#########################################################
### Open file (object-oriented programming)

class GlobalVar():
	###Makes variables accessible to all functions###
	nucleotide_alignments1_start = {} ###Chrom1 => Chrom2 => start
	nucleotide_alignments1_end = {}   ###Chrom1 => Chrom2 => end
	nucleotide_alignments2_start = {} ###Chrom2 => Chrom1 => start
	nucleotide_alignments2_end = {}   ###Chrom2 => Chrom1 => end
	lengths = {}    ###Chrom1 => Chrom2 => added lengths
			###Chrom2 => Chrom1 => added lengths

class OpenFile():
    	### Opens the file and directs it to a line-reader and stores information into variables used later
	def __init__ (self, filename, ty, ty2):
		"""Opens the file (accepts gzipped files)"""
		if re.search(".gz$", filename):
			self.file = gzip.open(filename, 'rb')
		else:
			self.file = open(filename, 'r')
		if ty == "aln":
			sys.stderr.write("Reading alignment ({}) file\n".format(ty2))          
			self.readLinesAln(self.file, ty2)			


	def readLinesAln(self, f, ty2):
		"""Reads alignment files and outputs them to a global variable"""
		self.filename = f
		self.type = ty2
		self.nucleotide_aligns = 0
		####Finds the minimum and maximum positions for orthologous chromosomes for both chromosomes
		self.header = self.filename.readline()
		for line in self.filename:
			line = line.rstrip('\n')
			self.nucleotide_aligns += 1
			self.query, self.subject, self.qstart, self.qend, self.sstart, self.send, self.length, self.orientation = line.split("\t")
			if str(self.query) != str(self.subject)	and not re.search ("^NW_", self.query): 
			###The same chromosome is not allowed to align to itself (if doing homeologous analysis) and scaffolds are not supported
				####Chromosome 1###
				if self.query in GlobalVar.nucleotide_alignments1_start:
					if self.subject in GlobalVar.nucleotide_alignments1_start[self.query]:
						GlobalVar.lengths[self.query][self.subject] += int(self.length)
						if float(GlobalVar.nucleotide_alignments1_start[self.query][self.subject]) > float(self.qstart):
							GlobalVar.nucleotide_alignments1_start[self.query][self.subject] = float(self.qstart)
					else:
						GlobalVar.lengths[self.query][self.subject] = int(self.length)
						GlobalVar.nucleotide_alignments1_start[self.query][self.subject] = float(self.qstart)
					if self.subject in GlobalVar.nucleotide_alignments1_end[self.query]:
						if float(GlobalVar.nucleotide_alignments1_end[self.query][self.subject]) < float(self.qend):
							GlobalVar.nucleotide_alignments1_end[self.query][self.subject] = float(self.qend)
					else:
						GlobalVar.nucleotide_alignments1_end[self.query][self.subject] = float(self.qend)
				else:
					GlobalVar.nucleotide_alignments1_start[self.query] = {}
					GlobalVar.nucleotide_alignments1_end[self.query] = {}
					GlobalVar.lengths[self.query] = {}
					GlobalVar.nucleotide_alignments1_start[self.query][self.subject] = float(self.qstart)
					GlobalVar.nucleotide_alignments1_end[self.query][self.subject] = float(self.qend)
					GlobalVar.lengths[self.query][self.subject] = int(self.length)
				####Chromosome 2###
				if self.subject in GlobalVar.nucleotide_alignments2_start:
					if self.query in GlobalVar.nucleotide_alignments2_start[self.subject]:
						GlobalVar.lengths[self.subject][self.query] += int(self.length)
						if float(GlobalVar.nucleotide_alignments2_start[self.subject][self.query]) > float(self.sstart):
							GlobalVar.nucleotide_alignments2_start[self.subject][self.query] = float(self.sstart)
					else:
						GlobalVar.lengths[self.subject][self.query] = int(self.length)
						GlobalVar.nucleotide_alignments2_start[self.subject][self.query] = float(self.sstart)
					if self.query in GlobalVar.nucleotide_alignments2_end[self.subject]:
						if float(GlobalVar.nucleotide_alignments2_end[self.subject][self.query]) < float(self.send):
							GlobalVar.nucleotide_alignments2_end[self.subject][self.query] = float(self.send)
					else:
						GlobalVar.nucleotide_alignments2_end[self.subject][self.query] = float(self.send)
				else:
					GlobalVar.nucleotide_alignments2_start[self.subject] = {}
					GlobalVar.nucleotide_alignments2_end[self.subject] = {}
					GlobalVar.lengths[self.subject] = {}
					GlobalVar.nucleotide_alignments2_start[self.subject][self.query] = float(self.sstart)
					GlobalVar.nucleotide_alignments2_end[self.subject][self.query] = float(self.send)
					GlobalVar.lengths[self.subject][self.query] = int(self.length)
		sys.stderr.write("\tFound {} alignments\n\n".format(self.nucleotide_aligns))
		self.find_overlap()
		self.filename.close()

	def find_overlap(self):
		"""Overlaps are detected and removed"""
		self.remove = {}
		for self.subject in GlobalVar.nucleotide_alignments2_start:
			for self.query1 in GlobalVar.nucleotide_alignments2_start[self.subject]:
				self.length1 = GlobalVar.lengths[self.subject][self.query1]
				self.start1 = GlobalVar.nucleotide_alignments2_start[self.subject][self.query1]
				self.finish1 = GlobalVar.nucleotide_alignments2_end[self.subject][self.query1]
				for self.query2 in GlobalVar.nucleotide_alignments2_start[self.subject]:
					if (self.query1 != self.query2 and "{}:{}".format(self.subject, self.query2) not in self.remove and 
						"{}:{}".format(self.subject, self.query1) not in self.remove):
						self.length2 = GlobalVar.lengths[self.subject][self.query2]
						self.start2 = GlobalVar.nucleotide_alignments2_start[self.subject][self.query2]
						self.finish2 = GlobalVar.nucleotide_alignments2_end[self.subject][self.query2]
						self.overlap, self.lenOver = GeneralOverlap.Overlap(self.start1, self.finish1, self.start2, self.finish2)
						sys.stderr.write("Sub: {}, Quer1: {}, Quer2: {}\n\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(self.subject, self.query1, self.query2, self.length1, self.start1, self.finish1, self.length2, self.start2, self.finish2, self.overlap))
						if self.overlap == "overlap1" or self.overlap == "overlap2":
							if self.length1 > self.length2:
								self.remove["{}:{}".format(self.subject, self.query2)] = 1
							else:
								self.remove["{}:{}".format(self.subject, self.query1)] = 1
						elif self.overlap == "overlap3":
							if self.length1 > self.length2:
								GlobalVar.nucleotide_alignments2_end[self.subject][self.query2] = float(GlobalVar.nucleotide_alignments2_start[self.subject][self.query1])
							else:
								GlobalVar.nucleotide_alignments2_start[self.subject][self.query1] = float(GlobalVar.nucleotide_alignments2_end[self.subject][self.query2])
						elif self.overlap == "overlap4":
							if self.length1 > self.length2:
								GlobalVar.nucleotide_alignments2_start[self.subject][self.query2] = float(GlobalVar.nucleotide_alignments2_end[self.subject][self.query1])
							else:
								GlobalVar.nucleotide_alignments2_end[self.subject][self.query1] = float(GlobalVar.nucleotide_alignments2_start[self.subject][self.query2])

		for self.pair in self.remove:
			self.s1, self.q1 = self.pair.split(":")
			sys.stderr.write("\tRemoved1: {}\t{}\n".format(self.s1, self.q1))
			del GlobalVar.nucleotide_alignments2_start[self.s1][self.q1]
			del GlobalVar.nucleotide_alignments1_start[self.q1][self.s1]

		self.remove.clear()
		for self.subject in GlobalVar.nucleotide_alignments1_start:
			for self.query1 in GlobalVar.nucleotide_alignments1_start[self.subject]:
				self.length1 = GlobalVar.lengths[self.query1][self.subject]
				self.start1 = GlobalVar.nucleotide_alignments1_start[self.subject][self.query1]
				self.finish1 = GlobalVar.nucleotide_alignments1_end[self.subject][self.query1]
				for self.query2 in GlobalVar.nucleotide_alignments1_start[self.subject]:
					if (self.query1 != self.query2 and "{}:{}".format(self.subject, self.query2) not in self.remove and 
						"{}:{}".format(self.subject, self.query1) not in self.remove):
						self.length2 = GlobalVar.lengths[self.query2][self.subject]
						self.start2 = GlobalVar.nucleotide_alignments1_start[self.subject][self.query2]
						self.finish2 = GlobalVar.nucleotide_alignments1_end[self.subject][self.query2]
						self.overlap, self.lenOver = GeneralOverlap.Overlap(self.start1, self.finish1, self.start2, self.finish2)
						if self.overlap == "overlap1" or self.overlap == "overlap2":
							if self.length1 > self.length2:
								self.remove["{}:{}".format(self.subject, self.query2)] = 1
							else:
								self.remove["{}:{}".format(self.subject, self.query1)] = 1
						elif self.overlap == "overlap3":
							if self.length1 > self.length2:
								GlobalVar.nucleotide_alignments1_end[self.subject][self.query2] = float(GlobalVar.nucleotide_alignments1_start[self.subject][self.query1])
							else:
								GlobalVar.nucleotide_alignments1_start[self.subject][self.query1] = float(GlobalVar.nucleotide_alignments1_end[self.subject][self.query2])
						elif self.overlap == "overlap4":
							if self.length1 > self.length2:
								GlobalVar.nucleotide_alignments1_start[self.subject][self.query2] = float(GlobalVar.nucleotide_alignments1_end[self.subject][self.query1])
							else:
								GlobalVar.nucleotide_alignments1_end[self.subject][self.query1] = float(GlobalVar.nucleotide_alignments1_start[self.subject][self.query2])
						
		for self.pair in self.remove:
			self.s1, self.q1 = self.pair.split(":")
			sys.stderr.write("\tRemoved2: {}\t{}\n".format(self.s1, self.q1))
			del GlobalVar.nucleotide_alignments1_start[self.s1][self.q1]

class OpenFileII():
    	### Opens the file and directs it to a line-reader and stores information into variables used later
	def __init__ (self, filename, ty, ty2):
		"""Opens the file (accepts gzipped files)"""
		if re.search(".gz$", filename):
			self.file = gzip.open(filename, 'rb')
		else:
			self.file = open(filename, 'r')
		if ty == "aln":
			sys.stderr.write("Reading alignment ({}) file (again)\n".format(ty2))          
			self.readLinesAln(self.file, ty2)			


	def readLinesAln(self, f, ty2):
		"""Reads alignment files and outputs them to a global variable"""
		self.filename = f
		self.type = ty2
		self.nucleotide_aligns = 0
		####Finds the minimum and maximum positions for orthologous chromosomes for both chromosomes
		self.header = self.filename.readline().rstrip('\n')
		print "{}".format(self.header)
		for line in self.filename:
			line = line.rstrip('\n')
			self.query, self.subject, self.qstart, self.qend, self.sstart, self.send, self.length, self.orientation = line.split("\t")
			if str(self.query) != str(self.subject): ###The same chromosome is not allowed to align to itself (if doing homeologous analysis)
				if (self.query in GlobalVar.nucleotide_alignments1_start and 
					self.subject in GlobalVar.nucleotide_alignments1_start[self.query] and
					self.subject in GlobalVar.nucleotide_alignments2_start and 
					self.query in GlobalVar.nucleotide_alignments2_start[self.subject]):
					self.q1start = GlobalVar.nucleotide_alignments1_start[self.query][self.subject]
					self.q1end = GlobalVar.nucleotide_alignments1_end[self.query][self.subject]
					self.s1start = GlobalVar.nucleotide_alignments2_start[self.subject][self.query]
					self.s1end = GlobalVar.nucleotide_alignments2_end[self.subject][self.query]
					if (float(self.qstart) >= float(self.q1start) and float(self.qstart) <= float(self.q1end) and
						float(self.qend) >= float(self.q1start) and float(self.qend) <= float(self.q1end) and
						float(self.sstart) >= float(self.s1start) and float(self.sstart) <= float(self.s1end) and
						float(self.send) >= float(self.s1start) and float(self.send) <= float(self.s1end)): 
						print "{}".format(line)
						self.nucleotide_aligns += 1
		sys.stderr.write("\tKept {} alignments\n\n".format(self.nucleotide_aligns))
		self.filename.close()



if __name__ == '__main__':
	GlobalVar()
	OpenFile(args.align, "aln", "nucleotide")
	OpenFileII(args.align, "aln", "nucleotide")
