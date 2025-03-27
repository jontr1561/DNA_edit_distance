import random
import matplotlib.pyplot as plt
from Bio import Entrez, SeqIO


# function generating random DNA sequences of length 400
def rand_dna():
    dna_choice = ['A', 'C', 'T', 'G']  # Possible entries in the DNA sequence
    sequence = random.choices(dna_choice, k=400)  # creates the sequence of length 400
    return sequence


# function that given string1 and string2 calculates the minimum edit distance
def opt(string1, string2):
    x = len(string1)  # gets length of string1
    y = len(string2)  # gets length of string2
    table = [[0] * (y + 1) for _ in range(x + 1)]  # initializes a table to track calculations
    for i in range(x+1):  # inputs entries for the first column of the table
        table[i][0] = i
    for j in range(y+1):  # inputs entries for the first row of the table
        table[0][j] = j
    for i in range(1, x+1):  # fills out the rest of the entries in the table
        for j in range(1, y+1):
            if string1[i-1] == string2[j-1]:
                table[i][j] = table[i-1][j-1]
            else:  # calculates the minimum edit distances by reference to the table
                table[i][j] = 1 + min(table[i-1][j], table[i][j-1], table[i-1][j-1])
    return table[x][y]


# Generates 20 random DNA sequences and calculates the edit distance for each pair of sequences
def random_trial():
    edit_results = []
    for iteration in range(20):  # runs this sequence 20 times
        string1 = rand_dna()  # generates a random DNA sequence of length 400
        string2 = rand_dna()  # generates a random DNA sequence of length 400
        result = opt(string1, string2)  # finds minimum edit distance for the two sequences
        edit_results.append(result)  # adds the result to a running list of results from past iterations
    print(edit_results)
    plt.hist(edit_results)  # creates a histogram modeling the edit distance result array
    plt.show()  # displays the histogram


# Dictionary with the keys being the name of the species and values being the accession code for data bank
speciesdata = {
    'German Neanderthal': 'AF011222',
    'Russian Neanderthal': 'AF254446',
    'European Human': 'X90314',
    'Mountain Gorilla Rwanda': 'AF089820',
    'Chimp Troglodytes': 'AF176766',
    'Puti Orangutan': 'AF451972',
    'Jari Orangutan': 'AF451964',
    'Western Lowland Gorilla': 'AY079510',
    'Eastern Lowland Gorilla': 'AF050738',
    'Chimp Schweinfurthii': 'AF176722',
    'Chimp Vellerosus': 'AF315498',
    'Chimp Verus': 'AF176731'
}


# function that given a list of accession codes runs the edit distance program pair-wise
def real_data(sequence_list):
    counter = 0
    species1 = 'blank'
    species2 = 'blank'
    result_list = []
    for species in sequence_list:  # iterates through the entire sequence list
        if counter == 0:
            species1 = species  # if it is the first species in the pair-wise comparison, assigns species1
            counter += 1  # increments count
        elif counter == 1:
            species2 = species  # if it is the second species in the pair-wise comparison assigns species2
            counter = 0  # resets count to0
        result = opt(species1, species2)  # runs the edit distance program on current iteration of sequence pair
        result_list.append(result)  # appends result to list of results
    return result_list


# Provide your email address to NCBI
Entrez.email = "jihtran@ucdavis.edu"


def fetch_sequence(accession):
    # Fetch the sequence from NCBI using the accession number
    handle = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta")
    record = SeqIO.read(handle, "fasta")
    handle.close()
    return record


def main():
    random_trial()  # edit distance algorithm initiated on random trials
    accessions = []  # initiates an empty list
    for species in speciesdata.values():  # for each species, grabs the accession code
        accessions.append(species)  # enters in accessions codes into a list
    sequence_list = []  # initiates an empty list
    for numbers in accessions:  # for loop that extracts DNA sequence for each species and appends it to the list
        sequence_record = fetch_sequence(numbers)
        sequence = sequence_record.seq
        sequence_list.append(sequence)
    real_trial = real_data(sequence_list)
    plt.hist(real_trial)  # creates histogram representing edit distance for real DNA sequences
    plt.show()  # displays histogram
    print(real_trial) # displays the list of edit distance values


if __name__ == "__main__":
    main()
