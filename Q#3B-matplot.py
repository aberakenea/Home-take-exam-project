'''
Created on May 28, 2023

@author: wku
'''
import matplotlib.pyplot as plt
import numpy as np
import os

def extract_let7_code(header):
    if 'let-7' in header:
        code = header.split('let-7')[1].split()[0]
        return f"let-7{code[0]}"
    return ""

def extract_species_code(header):
    start_index = header.find('>') + 1
    end_index = header.find('-', start_index)
    if start_index < end_index:
        return header[start_index:end_index]
    return ""

def let7_family_presence(file_path):
    if not os.path.isfile(file_path):
        print("The specified file does not exist.")
        return

    let7_species = {}

    with open(file_path, 'r') as file:
        current_species = ""
        for line in file:
            if line.startswith('>'):
                header = line[1:].strip()
                current_species = extract_species_code(header)
            else:
                let7_code = extract_let7_code(header)
                if let7_code:
                    let7_species.setdefault(current_species, {}).setdefault(let7_code, 0)
                    let7_species[current_species][let7_code] += 1

    # Filter species with frequency count less than 10
    let7_species_filtered = {species: counts for species, counts in let7_species.items() if sum(counts.values()) >= 10}
    species_list = list(let7_species_filtered.keys())
    let7_codes = list(set().union(*[d.keys() for d in let7_species_filtered.values()]))

    # Create a matrix to store the presence counts
    presence_matrix = np.zeros((len(species_list), len(let7_codes)))

    for i, species in enumerate(species_list):
        for j, let7_code in enumerate(let7_codes):
            presence_matrix[i, j] = let7_species_filtered[species].get(let7_code, 0)

    # Set the positions of the bars on the x-axis
    x = np.arange(len(species_list))
    # Set the width of the bars
    bar_width = 0.35
    # Create the figure and axes
    fig, ax = plt.subplots(figsize=(10, 6))

    # Plot the stacked bars
    bottom = np.zeros(len(species_list))
    for i, let7_code in enumerate(let7_codes):
        ax.bar(x, presence_matrix[:, i], bottom=bottom, label=let7_code)
        bottom += presence_matrix[:, i]

        # Add count labels for individual let-7 miRNAs
        for j, species in enumerate(species_list):
            count = presence_matrix[j, i]
            if count > 0:
                ax.text(x[j], bottom[j] - count / 2, str(int(count)), ha='center', va='center')

    # Add labels and title
    ax.set_xlabel('Species')
    ax.set_ylabel('Presence Count')
    ax.set_title('Presence of let-7 Family per Species (Frequency >= 15)')
    ax.set_xticks(x)
    ax.set_xticklabels(species_list, rotation=45, ha='right')

    # Add a legend
    ax.legend()

    # Display the plot
    plt.tight_layout()
    plt.show()

# File path
file_path = r"C:\Users\wku\Advancedprogram\project\GC_calc-complexity-project\mature.fa"
let7_family_presence(file_path)