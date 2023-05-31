import os
import Levenshtein
import matplotlib.pyplot as plt

def extract_let7_code(header):
    if 'let-7' in header:
        code = header.split('let-7')[1].strip()[0]
        return f"let-7{code}"
    return ""

file_path = r"C:\Users\wku\Advancedprogram\project\GC_calc-complexity-project\mature.fa"
if not os.path.isfile(file_path):
    print("The specified file does not exist.")
    exit()

let7_sequences = {}
let7_frequency = {}

with open(file_path, 'r') as file:
    let7_code = ""
    sequence = ""
    for line in file:
        if line.startswith('>'):
            header = line[1:].strip()
            let7_code = extract_let7_code(header)
            sequence = ""
        else:
            sequence = line.strip()

        if let7_code and sequence:
            let7_sequences.setdefault(let7_code, []).append(sequence)
            let7_frequency[let7_code] = let7_frequency.get(let7_code, 0) + 1

total_sequences_count = 0
avg_distances = []
frequencies = []
miRNA_families = []

for let7_code, sequences in let7_sequences.items():
    total_distance = 0
    total_pairs = 0

    if len(sequences) < 2:
        continue

    for i in range(len(sequences) - 1):
        for j in range(i + 1, len(sequences)):
            total_distance += Levenshtein.distance(sequences[i], sequences[j])
            total_pairs += 1

    average_distance = total_distance / total_pairs
    frequency = let7_frequency.get(let7_code, 0)
    avg_distances.append(average_distance)
    frequencies.append(frequency)
    miRNA_families.append(let7_code)
    print(f"The Average Levenshtein distance among all pairs for miRNA family {let7_code}: {average_distance:.2f}")
    print(f"The frequency of miRNA family {let7_code}: {frequency}")
    total_sequences_count += frequency

print(f"Total sequences count: {total_sequences_count}")

# Plotting
x_pos = range(len(miRNA_families))
fig, ax1 = plt.subplots()
ax1.bar(x_pos, avg_distances, align='center', alpha=0.5)
ax1.set_ylabel('Average Levenshtein Distance')
ax1.set_title('Average Levenshtein Distance and Frequency for miRNA Families')

ax2 = ax1.twinx()
ax2.plot(x_pos, frequencies, 'r')
ax2.set_ylabel('Frequency')

plt.xticks(x_pos, miRNA_families)
plt.xlabel('miRNA Family')

plt.show()