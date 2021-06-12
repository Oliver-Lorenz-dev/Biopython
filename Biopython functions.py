# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

#biopython
import Bio
from Bio.Seq import Seq

a = Seq('ATGATCTCGTAA')

import matplotlib.pyplot as plt
from collections import Counter

'''function plots graph of base frequency in a sequence'''
def base_frequency(dna):
    base_frequency = Counter(dna)
    plt.bar(base_frequency.keys(),base_frequency.values())

base_frequency(a)

#comp rev comp
a.complement()
a.reverse_complement()

#translation
mRNA = a.transcribe()
mRNA.translate()
#direct from DNA
a.translate()

# mrna to DNA
mRNA.back_transcribe()

#GC
from Bio.SeqUtils import GC
GC(a)

#AT
def at_content(seq):
    result = float(seq.count('A') + seq.count('T'))/len(seq) * 100
    return result

#nucleotide search
from Bio.SeqUtils import nt_search

main_seq = Seq('ACTATTGATT')
subseq = Seq('ATT')

nt_search(str(main_seq),str(subseq))

#alignments
from Bio import pairwise2
from Bio.pairwise2 import format_alignment

seq1 = Seq('ACTCGT')
seq2 = Seq('ATTCG')
'''function to perform global alignment on 2 sequences'''
def global_alignment(seq1,seq2):
    alignments = pairwise2.align.globalxx(seq1,seq2)
    return (format_alignment(*alignments[0]))

print(global_alignment(seq1,seq2))
'''function to perform local alignment on 2 sequences'''
def local_alignment(seq1,seq2):
    alignments = pairwise2.align.localxx(seq1,seq2)
    return (format_alignment(*alignments[0]))

print(local_alignment(seq1,seq2))

# Get the alignment by only the score
alignment2 = pairwise2.align.globalxx(seq1,seq2,one_alignment_only=True,score_only=True)
print(alignment2)

#reading in files
from Bio import SeqIO
#fasta
def read_fasta(file):
    for record in SeqIO.parse(file,"fasta"):
        print(record)

read_fasta('sequence.fasta')

#just sequence

dna_record = SeqIO.read("sequence.fasta","fasta")
dna_seq = dna_record.seq
#print(dna_seq)

#genbank
# Reading GenBank 
# Load A GenBank File
def read_genbank():
    for record in SeqIO.parse("sequence.gb","genbank"):
        print(record)

read_genbank()

# Reading PDB Files
'''from Bio.PDB import PDBParser,MMCIFParser
parser = PDBParser()
structure = parser.get_structure("6LU7","6lu7.pdb")
model  = structure[0]'''
# Check for atoms
'''for chain in model:
    print(f'Chain {chain},Chain_ID{chain.id}')
    for residue in chain:
        for atom in residue:
            print(atom)'''

#BLAST
from Bio.Blast import NCBIWWW
#Read in seq
covid_record = SeqIO.read("covid_sequence_MT385497.fasta","fasta") 
covid_dna = covid_record.seq
#BLAST



with NCBIWWW.qblast("blastn","nt",covid_dna) as result_handle:
    with open("result_blast_covid.xml","w") as xml_file:
        xml_file.write(result_handle)
        
    
    






