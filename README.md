# OrthologyBetweenSpecies
A method to identify orthologous and syntenic genes between closely related species.  This pipeline identifies orthologs based on synteny and sequence similarity.  This pipeline removes potential ortholog assignments if there are multiple locations in which it could be placed (i.e. multiple syntenic regions).  Future versions will hopefully address this issue in a more refined manner.

<!-- TABLE OF CONTENTS -->
<details open="open">
  <summary>Table of Contents</summary>
  <ol>
    <li><a href="#requirements">Requirements</a></li>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#license">License</a></li>
  </ol>
</details>

<!-- requirements -->
## Requirements

These scripts have been tested with Python 2.7 and will not work with later versions unless they are modified.

These scripts work with alignments from blast (blast format 6 is required).  We assume the user has already generated blast databases below.  Please see blast documentation for further details and for an explanation of command line options used below.

Files: Both genome fasta files, both set of protein sequences for entire genome in fasta format, and both .gff files for the genomes, samtools faidx indexes of both genomes (.fai)

<!-- usage -->
## Usage

1) Align genomes of closely related species.  This was done using blast for this pipeline (hopefully future versions will include a quicker software solution).

       Example (for closely related species): 
       blastn -task megablast -db genome1 -query genome2.fasta -evalue 0.000001 -max_target_seqs 3 -max_hsps 20000 -outfmt 6 -word_size 40\
       -perc_identity 96 -lcase_masking -soft_masking false -num_threads 10 -out genome1_vs_genome2.aln
       
2) Filter non-linear global alignments

       Note: Requires the python script Linear_Alignments_v1_2.py to be in the same directory
       
       Example:
       python Compare_Genome_2_Other_Genome_blastfmt6_ver1.0.py -aln genome1_vs_genome2.aln -qfasta genome2.fasta -sfasta genome1.fasta -minl 0.01 -minal 30000 \
       > Filtered_genome1_vs_genome2.aln
       
       Help:
       python Compare_Genome_2_Other_Genome_blastfmt6_ver1.0.py -h
  
3) Keep the best overlapping linear global alignment

       Note: Requires the python script GeneralOverlap.py to be in the same directory
       
       Example:
       python Filter_Linear_Alignment.v1.0.py -align Filtered_genome1_vs_genome2.aln > Filtered2_genome1_vs_genome2.aln
       
       Help:
       python Filter_Linear_Alignment.v1.0.py -h
       
4) Align all protein sequences from genome1 to genome2

       Example:
       blastp -db genome1 -query genome2.faa -max_target_seqs 3 -max_hsps 20 -evalue 0.01 -outfmt 6 -num_threads 10 -out genome1_vs_genome2.protein.aln
     
5) Filter protein alignments

       Example:
       python Filter_Alignments_Blast_Fmt6_Protein_ver1.0.py -aln_file genome1_vs_genome2.protein.aln -fastas genome1.faa -fastaq genome2.faa -min_per 80 \
       -min_aln_per 80 > Filtered.genome1_vs_genome2.protein.aln
   
       Help:
       python Filter_Alignments_Blast_Fmt6_Protein_ver1.0.py -h

6) Generate a list of orthologous genes between the two genomes

       Example:
       python Orthology_Between_Genomes.v1.1.py -gffs genome1.gff -gffq genome2.gff -alignn Filtered2_genome1_vs_genome2.aln \
       -alignp Filtered.genome1_vs_genome2.protein.aln -fais genome1.fai -faiq genome2.fai > orthologs.txt
       
       Help:
       python Orthology_Between_Genomes.v1.1.py -h



<!-- license -->
## License 

Distributed under the MIT License.
