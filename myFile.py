#Import math for some calculations in this file
import math 

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
    
   
#Question 6 IPRB: Calculate the probability of producing offsprings with a dominant allele
def pr_dominant(k,m,n):
    #total population of mating organisms
    population = k + m + n
    #total number of offsprings produced if any two organisms can mate
    All_Offspring = math.comb(population, 2)
    #total number of offsprings with dominant phenotype
    #Any two organisms of k, any groups of one k and one m, any groups of one k and one n, half of the groups of one m and one n, 3/4 of any two organisms of m can produce offsprings with a dominant allele
    Dominant_Offspring = math.comb(k,2) + k*m + k*n + 0.5*m*n + 0.75*(math.comb(m,2))
    #Calculate the probability of producing offsprings with dominant phenotype
    probability = Dominant_Offspring/All_Offspring
    print(probability)

#Question 7 GC: Select the name of the dna sequence which has the biggest gc content from a FASTA file. Calculate the biggest gc content. 
def GC_content(data):
    #Open the data file. Data refers to "rosalind_gc.txt", a file provided by the question
    dnafile = open(data,"r")
    #Read the file and get the raw data
    rawdata = dnafile.read()
    #Creates a dictionary
    dic = {}
    #The file is a fasta file, so identify  ">" to split different dna sequences from each other. Ignore every title line starting with ">"
    for dna in rawdata.split(">")[1:]:
        parts = dna.split("\n")
        #The title line starting with ">"
        dna = parts[0]
        #Get the dna sequence connecting different lines together
        seq = ''.join(parts[1:])
        #Calculate the gc content
        gc = 100 * (seq.count("G") + seq.count("C")) / float(len(seq))
        #Use the dictionary. Key is dna name and value is the gc content of dna
        dic[gc] = dna
    #print out the biggest gc content and the corresponding dna name
    print(dic[max(dic)],max(dic))
    

    
    

