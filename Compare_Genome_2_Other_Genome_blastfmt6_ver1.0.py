##########################################################
### Import Necessary Modules

import argparse		               #provides options at the command line
import sys		               #take command line arguments and uses it in the script
import gzip		               #allows gzipped files to be read
import re		               #allows regular expressions to be used
import textwrap		               #allows the use of textwrapping for long sequences
import Linear_Alignments_v1_2	       #allows linear alignments to be found and retained

##########################################################
### Command-line Arguments
parser = argparse.ArgumentParser(description="A script to filter non-linear global alignments from blast alignments (outfmt 6)")
parser.add_argument("-aln", help = "The location of an alignment file with a genome aligned to another genome (blast outfmt 6)", default=sys.stdin, required=True)
parser.add_argument("-qfasta", help = "The location of a fasta file with the query scaffolds for sizes", default=sys.stdin, required=True)
parser.add_argument("-sfasta", help = "The location of a fasta file with the query scaffolds for sizes", default=sys.stdin, required=True)
parser.add_argument("-smax", help = "The maximum distance a linear alignment can be apart (in fraction of query), default:0.01", default=0.01)
parser.add_argument("-cmax", help = "The maximum distance a linear alignment can be apart (in fraction of subject), default:0.01", default=0.01)
parser.add_argument("-minl", help = "The minimum length (in fraction of query) that a linear alignment can be to be retained, default:0.1", default=0.1)
parser.add_argument("-minal", help = "The minimum actual length that a linear alignment can be to be retained, default:20000", default=20000)
args = parser.parse_args()


#########################################################
### Open file (object-oriented programming)


class OpenFasta():
    query_lengths = {}
    subject_lengths = {}
    ### Opens the file either using a regular mechanism or opens it after uncompressing the data
    def __init__ (self, f1, f2):
        """Initiates the class by opening a fasta file either compressed or uncompressed"""
        if re.search(".gz$", f1):
            self.fasq = gzip.open(f1, 'rb')
        else:
            self.fasq = open(f1, 'r')
        if re.search(".gz$", f2):
            self.fass = gzip.open(f2, 'rb')
        else:
            self.fass = open(f2, 'r')  
        sys.stderr.write("\nOpened query fasta file\nOpened subject fasta file\n")          
        self.readLinesFasta(self.fasq, "query")
        self.readLinesFasta(self.fass, "subject")

    ### Reads the fasta file and returns length of sequences
    def readLinesFasta(self, fa, ty):
        """Reads fasta file line by line and returns the length of each sequence to a dictionary/hash and other info to user"""
        self.fasta = fa
        self.type = ty
        self.header = "NA"
        self.nucleotide_length = 0
        for line in self.fasta:
            line = line.rstrip('\n')
            if re.search("^\>", line):
                if self.header == "NA":
                    self.header = line[1:].split(" ")[0]
                else:
                    if self.type == "query":
                        OpenFasta.query_lengths[self.header] = self.nucleotide_length
                    else:
                        OpenFasta.subject_lengths[self.header] = self.nucleotide_length
                    self.nucleotide_length = 0
                    self.header = line[1:].split(" ")[0]
            elif re.search("\w", line):
		self.nucleotide_length += int(len(line))
        self.chrom_count = 0
        self.chrom_length = 0
        if self.type == "query":
            OpenFasta.query_lengths[self.header] = self.nucleotide_length
            for self.head in OpenFasta.query_lengths:
                self.chrom_count += 1
                self.chrom_length += int(OpenFasta.query_lengths[self.head])
            sys.stderr.write("Finished reading query fasta file: Found {} sequence(s), {} nucleotide(s)\n".format(self.chrom_count, self.chrom_length))
        else:
            OpenFasta.subject_lengths[self.header] = self.nucleotide_length
            for self.head in OpenFasta.subject_lengths:
                self.chrom_count += 1
                self.chrom_length += int(OpenFasta.subject_lengths[self.head])
            sys.stderr.write("Finished reading subject fasta file: Found {} sequence(s), {} nucleotide(s)\n".format(self.chrom_count, self.chrom_length))      
        self.fasta.close()


class OpenAln():
    MainAlignments = {} ### query => subject => "queryPercentStart\tqueryPercentEnd\tsubjectPercentStart\tsubjectPercentEnd\tlength\torientation"
    ### Opens the file either using a regular mechanism or opens it after uncompressing the data
    def __init__ (self, aln):
        """Opens an alignment file output by the blast program outfmt 6 (compressed or uncompressed)"""
        if re.search(".gz$", aln):
            self.aln = gzip.open(aln, 'rb')
        else:
            self.aln = open(aln, 'r')
        sys.stderr.write("\nOpened blast alignment (format outfmt 6)\n")
        print "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format("Query", "Subject", "Start.PercentQueryLength",
            "End.PercentQueryLength", "Start.PercentSubjectLength", "End.PercentSubjectLenth", "LengthAlignment", "Orientation")            
        self.readLinesAln(self.aln)

    def readLinesAln(self, a):
        """Alignment file is read in line by line, but only alignments passing the filtering are returned"""
        self.open_aln = a
        for line in self.open_aln:
            line = line.rstrip('\n')
            (self.queryID, self.subjectID, self.percentID, self.length, self.mismatch, self.gapopen, self.queryStart, self.queryEnd, 
                self.subjectStart, self.subjectEnd, self.evalue, self.bitscore) =  line.split("\t")  
            self.orientation = "+"
            if int(self.subjectEnd) < int(self.subjectStart):
                self.orientation = "-"
                self.temp = self.subjectStart
                self.subjectStart = self.subjectEnd
                self.subjectEnd = self.temp
            (self.queryStart, self.queryEnd, self.subjectStart, self.subjectEnd) = (float(self.queryStart)/int(OpenFasta.query_lengths[self.queryID]), 
		float(self.queryEnd)/int(OpenFasta.query_lengths[self.queryID]), float(self.subjectStart)/int(OpenFasta.subject_lengths[self.subjectID]), 
                float(self.subjectEnd)/int(OpenFasta.subject_lengths[self.subjectID])) ###Puts the distances in percent of chromosome
            if not self.queryID in OpenAln.MainAlignments:
                OpenAln.MainAlignments[self.queryID] = {}
                OpenAln.MainAlignments[self.queryID][self.subjectID] = "{}\t{}\t{}\t{}\t{}\t{}".format(self.queryStart, self.queryEnd, self.subjectStart, self.subjectEnd,
                    self.length, self.orientation)
            else:
                if self.subjectID in OpenAln.MainAlignments[self.queryID]:
                    OpenAln.MainAlignments[self.queryID][self.subjectID] += ",{}\t{}\t{}\t{}\t{}\t{}".format(self.queryStart, self.queryEnd, self.subjectStart, 
                        self.subjectEnd, self.length, self.orientation)
                else:
                    OpenAln.MainAlignments[self.queryID][self.subjectID] = "{}\t{}\t{}\t{}\t{}\t{}".format(self.queryStart, self.queryEnd, self.subjectStart, 
                        self.subjectEnd, self.length, self.orientation)
        self.open_aln.close()
        sys.stderr.write("Finished reading alignment file\n\n")
        self.checkLinearity()

    def checkLinearity(self):
        for self.query in OpenAln.MainAlignments:
            for self.subject in OpenAln.MainAlignments[self.query]:
                self.long_list = OpenAln.MainAlignments[self.query][self.subject].split(",")
                sys.stderr.write("\nStarting Analysis for {} vs {}\n".format(self.query, self.subject))
                self.final = Linear_Alignments_v1_2.linearAlignments(OpenAln.MainAlignments[self.query][self.subject], args.smax, args.cmax, args.minl, args.minal)
                for self.aln in self.final:
                    (self.qstart, self.qend, self.sstart, self.send, self.len, self.orientation) = self.long_list[self.aln].split("\t")
                    print "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(self.query, self.subject, self.qstart, self.qend, self.sstart, self.send, self.len, self.orientation)

if __name__ == '__main__':
    open_fasta = OpenFasta(args.qfasta, args.sfasta)
    open_aln = OpenAln(args.aln)
