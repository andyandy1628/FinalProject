#Question 1 DNA: Count "A","T","C","G" respectively in a DNA sequence
def count(dnaSeq)
    print(dnaSeq.count("A"),dnaSeq.count("C"),dnaSeq.count("G"),dnaSeq.count("T"))
    
#Question 2 RNA: Transcribe the complementary chain of a DNA sequence into RNA
def transcription(dnaSeq):
    print(dnaSeq.replace("T","U"))
    
#Question 3 REVC: Print out the reverse complement of a DNA sequence
def complement_dna(dnaSeq):
    """helper function that sets the rules of complementing DNA sequences."""
    #The complement rule is "A" to "T", "T" to "A", "C" to "G", "G" to "C". 
    if dnaSeq in "A":
         return "T"
    elif dnaSeq in "T":
         return "A"
    elif dnaSeq in "C":
         return "G"
    else:
         return "C"
         
def reverse_complement(dnaSeq):
    reverse =  dnaSeq[:: -1]
    complement = ''
    #Get the complementary sequence using the rules.
    for each in reverse:
        complement += complement_dna(each)
    print(complement)

#Question 4 HAMM: Count point mutations between two DNA sequences
def count_point_mutation(dnaSeq1,dnaSeq2):
    count = 0
    for i in range(len(dnaSeq1)):
        if dnaSeq1[i] != dnaSeq2[i]:
            count += 1
    print(count)
    
#Question 5 PROT: Translate rna into protein
def translation(rnaSeq):
    #RNA translation tab;e
    RNA_codon_table = {
        'UUU':'F', 'CUU':'L', 'AUU':'I', 'GUU':'V', 
        'UUC':'F', 'CUC':'L', 'AUC':'I', 'GUC':'V', 
        'UUA':'L', 'CUA':'L', 'AUA':'I', 'GUA':'V', 
        'UUG':'L', 'CUG':'L', 'AUG':'M', 'GUG':'V',                  
        'UCU':'S', 'CCU':'P', 'ACU':'T', 'GCU':'A', 
        'UCC':'S', 'CCC':'P', 'ACC':'T', 'GCC':'A', 
        'UCA':'S', 'CCA':'P', 'ACA':'T', 'GCA':'A', 
        'UCG':'S', 'CCG':'P', 'ACG':'T', 'GCG':'A', 
        'UAU':'Y', 'CAU':'H', 'AAU':'N', 'GAU':'D', 
        'UAC':'Y', 'CAC':'H', 'AAC':'N', 'GAC':'D', 
        'UAA':'Stop', 'CAA':'Q', 'AAA':'K', 'GAA':'E', 
        'UAG':'Stop', 'CAG':'Q', 'AAG':'K', 'GAG':'E', 
        'UGU':'C', 'CGU':'R', 'AGU':'S', 'GGU':'G', 
        'UGC':'C', 'CGC':'R', 'AGC':'S', 'GGC':'G', 
        'UGA':'Stop', 'CGA':'R', 'AGA':'R', 'GGA':'G', 
        'UGG':'W', 'CGG':'R', 'AGG':'R', 'GGG':'G',
    }
    protein =""  
    #Translate rna in groups of three. Ignore the stop codon at last
    for i in range(0, len(rnaSeq)-3, 3): 
        codon = rnaSeq[i:i + 3] 
        protein+= RNA_codon_table[codon] 
    print(protein)
    

    
    

