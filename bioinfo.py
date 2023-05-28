#import Seq module from Bio.seq library
from Bio.Seq import Seq
a = Seq("ATGGGCCAAACTAATT")
print("The Sequence is", a)

#Calculating GC Content
from Bio.SeqUtils import gc_fraction
GC = 100 * (a.count("G") + a.count("C"))/ (len(a))
GC_1 = gc_fraction(a)
print("The GC content is", GC)
print("The GC content is", GC_1)

#finding complement and reverse complement 
b = Seq("ATCCGGGAACCTACGACGAATC")
complement = a.complement()
reverse_complement = a.reverse_complement()
print("Sequence =", b)
print("The complement of the sequence is", a)
print("The reverse complement is", a)

#transcription and translation of the sequence
#trancription
coding_seq = Seq("ATTGCCATGCCATGAAGCTTAGATGCAAGATACAGGACA")
print("Coding sequence =", coding_seq)
mRNA = coding_seq.transcribe()
print("biopython_mRNA =", mRNA)
true_mRNA = coding_seq.reverse_complement().transcribe()
print("Transcribed mRNA =", true_mRNA)

#translation
translate = true_mRNA.translate()
print("Translated mRNA =",translate)

#termination at stop codons
true_mRNA = coding_seq.reverse_complement()
stop_true_mRNA = true_mRNA.translate(to_stop=True)
print("Stop codon peptide =", stop_true_mRNA)

#translation tables
from Bio.Data import CodonTable
translation_table = CodonTable.ambiguous_dna_by_id[2]
print("translation table", translation_table)