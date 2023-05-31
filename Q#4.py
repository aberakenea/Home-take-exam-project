import Levenshtein
import matplotlib.pyplot as plt

let7a_mirnas = []

file_path = r"C:\Users\wku\Advancedprogram\project\GC_calc-complexity-project\mature.fa"

with open(file_path, "r") as file:
    miRNA_name = ""
    sequence = ""
    for line in file:
        if line.startswith(">"):
            if miRNA_name != "":
                if miRNA_name.startswith("hsa-let-7"):
                    let7a_mirnas.append((miRNA_name, sequence))
            miRNA_name = line.strip()[1:]
            sequence = ""
        else:
            sequence += line.strip()

    if miRNA_name.startswith("-let-7"): 
        let7a_mirnas.append((miRNA_name, sequence))

print("let-7 miRNA codes:")
for mirna, _ in let7a_mirnas:
    print(mirna)

# Calculate pairwise Levenshtein distance
distances = []
for mirna, seq in let7a_mirnas:
    print("Pairwise Levenshtein distances for", mirna)
    for other_mirna, other_seq in let7a_mirnas:
        if mirna != other_mirna:
            distance = Levenshtein.distance(seq, other_seq)
            distances.append(distance)
            print(f"Levenshtein distance between {mirna} and {other_mirna}: {distance}")

# Generate histogram plot
plt.hist(distances, bins=15)  # Adjust the number of bins as needed
plt.xlabel("Levenshtein Distance")
plt.ylabel("Frequency")
plt.title("Distribution of Levenshtein Distances")
plt.show()
