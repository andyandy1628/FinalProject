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
    

#Question 8 SUBS: Find all the locations of motifs in a DNA sequence
def find_motif(dnaSeq, motif):
    #define the locations as a list
    locations = []
    #For loop goes through the DNA sequence
    for i in range(0,len(dnaSeq)):
        #Check if parts of the DNA sequence equals to the motif
        if dnaSeq[i:i+len(motif)] == motif:
            #Rules in the question: the index starts at "1" rather than "0"
            locations.append(i+1)
    print(locations)
    
#Question 9 FIB: Calculate the total number of rabbits pairs after given months and every rabbit pair can produce given amount of new rabbits
def rabbits(months, rabbit_pairs_produced):
    #In month 1 and month 2, there are only one pairs of rabbits
    #Define the number of rabbits in following months as intergers
    F1 = 1
    F2 = 1
    F3 = 0
    #months + 1 because need to consider the rabbits produced by the last generation
    #Start to calculate the rabbits produced
    #[3:] because small rabbits in month 3 grows big enough to reproduce in month 4
    for i in range(months+1)[3:]:
        #Rabbits grow bigger to reproduction age
        F3 = F2
        #F2 equals rabbits produced by first generation plus current rabbit numbers in F2
        #Based on Fibonacci sequence
        F2 = F1*rabbit_pairs_produced + F2
        F1 = F3
    print(F2)
    
    
#Question 10 PERM: Permutations of positive intergers smaller than or equal to a given number
def factorial(number):
    #calculate the factorial of number
    #factorial of 1 is 1
    if number == 1:
        return 1
    #e.g. factorial of 5 is 5*4*3*2*1
    else:
        return number*factorial(number-1)

def all_number(number):
    #Append all positive integers smaller than the given number to a list
    all = []
    for n in range(number):
        all.append(n+1)
    return all

def permutation(all, position, last):
    #Position is 0 and last is len(all)
    #Give all the permutations of the number
    #If current position(starts from 0) is the end of the number list
    if position == last:
        #If n is the last number in the list, start a new line
        for n in all:
            if n != all[-1]:
                print(n, end= "")
            else:
                print(n, end= "\n")

    #If current position is not the end of the number list
    else:
        for index in range(position, last):
            #adjust to the right position
            all[index], all[position] = all[position], all[index]
            #permutation starting with next number in the list
            permutation(all, position + 1, last)
            all[index], all[position] = all[position], all[index]
            
#Question 11 LCSM: Find the longest common substring among multiple DNA sequences
def sortfile(filename):
    with open(filename)as file:
        lines = file.read().splitlines()
    #Figuring out a way to give the shortest DNA sequence because it is the longest possible common substring
    

    

    
    

