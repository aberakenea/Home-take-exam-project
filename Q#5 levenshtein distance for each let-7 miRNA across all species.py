import os

from Levenshtein import distance

def calculate_levenshtein_distance(file_path):

    if not os.path.isfile(file_path):

        print("The specified file does not exist.")

        return

    let7_sequences = []

    total_let7_miRNA = 0

    total_distance = 0

    total_pairs = 0

    with open(file_path, 'r') as file:

        species = ""

        sequence = ""

        for line in file:

            if line.startswith('>'):

                header = line[1:].strip()

                species = extract_species(header)

                sequence = ""

            else:

                sequence = line.strip()

            if 'let-7e' in header and species and sequence:


                let7_sequences.append((species, sequence))


                total_let7_miRNA += 1

    for i in range(len(let7_sequences) - 1):

        for j in range(i + 1, len(let7_sequences)):


            species1, seq1 = let7_sequences[i]


            species2, seq2 = let7_sequences[j]


            levenshtein_distance = distance(seq1, seq2)

            total_distance += levenshtein_distance

            total_pairs += 1

            print(f"Species: {species1} - {species2} | Levenshtein Distance: {levenshtein_distance}")


    if total_pairs > 0:

        average_distance = total_distance / total_pairs

        print(f"Total 'let-7' miRNAs: {total_let7_miRNA}")

        print(f"Average Levenshtein Distance of 'let-7' miRNAs: {average_distance:.2f}")

def extract_species(header):

    species = header.split(' ')[0]

    return species

# File path

file_path = r"C:\Users\wku\Advancedprogram\project\GC_calc-complexity-project\mature.fa"

calculate_levenshtein_distance(file_path)