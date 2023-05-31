# Home-take-exam-project
#### Question 1: How many species are there in the current miRBase release?
# ANSWERS:
## the code for the current miRBase release and it` path_file is:
import re
from collections import Counter
import matplotlib.pyplot as plt

file_path = r"C:\Users\wku\Advancedprogram\project\GC_calc-complexity-project\mature.fa"

species_codes = []

with open(file_path, 'r') as file:
    for line in file:
        match = re.match(r'^>(\w+)', line)
        if match:
            species_codes.append(match.group(1))

species_counts = Counter(species_codes)
total_species = len(species_counts)
print("Total number of species:", total_species)

# Filter species with frequencies greater than or equal to 220
min_frequency = 220
filtered_species_counts = {species: count for species, count in species_counts.items() if count >= min_frequency}

# Sort filtered species based on microRNA counts in ascending order
sorted_species_counts = sorted(filtered_species_counts.items(), key=lambda x: x[1])

# Extract species and their counts from sorted list
species = [item[0] for item in sorted_species_counts]
counts = [item[1] for item in sorted_species_counts]

# Plotting the bar chart
plt.figure(figsize=(10, 6))  # Set the figure size
plt.bar(species, counts)  # Create the bar plot
plt.xlabel('Species')  # Set the x-axis label
plt.ylabel('Count')  # Set the y-axis label
plt.title('Number of miRNA per Species (Lowest to Highest)')  # Set the plot title
plt.xticks(rotation=90)  # Rotate x-axis labels for better visibility
plt.tight_layout()  # Adjust the layout for better spacing
plt.show()  # Display the plot

# the possible out put number of species in the current miRBase release are:
### Total number of species: 271
#### Task: generate an ordered plot (i.e. from lowest to highest) of number of miRNAs / species
#![Figure_1](https://github.com/aberakenea/Home-take-exam-project/assets/130226484/62702380-cb84-4908-adf6-bcfd95016ee4)

#### Question 2: how many **let-7** miRNAs are there in the current release of miRBase
# ANSWERs:

##  it`s out put is Total number of let-7 miRNAs across all species: 740

#### Question 3: what is the current version of miRBase?
### ANSWERS:
## code for the current version of miRBase and it`s path_file is:
'''
Created on May 20, 2023

@author: ABERA KENEA
'''
import requests

import re



def get_mirbase_version():

    try:

        url = "https://www.mirbase.org/ftp/CURRENT/README"

        response = requests.get(url)

        response.raise_for_status()  # Raise an exception for HTTP errors

        readme_text = response.text

        # Search for the version information in the README file

        version_match = re.search(r"Release (\d+\.\d+)", readme_text)

        if version_match:

            miRBase_version = version_match.group(1)

            return miRBase_version

    except requests.RequestException as e:

        print("An error occurred while making the request:", e)

    except re.error:

        print("An error occurred while searching for the version information.")

    

    return None
# Call the function to retrieve the miRBase version

mirbase_version = get_mirbase_version()



if mirbase_version:

    print("Current version of miRBase:", mirbase_version)

else:

    print("Unable to retrieve the current version of miRBase.")

## the possible out put is:
# Current version of miRBase: 22.1
#### Task: generate a plot to show which let miRNAs are present in each species.
## ANSWER:
## the code to generate a plot to show which let miRNAs are present in each species is the following: 
'''
Created on May 28, 2023

@author: ABERA KENEA
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

## ![Figure_2](https://github.com/aberakenea/Home-take-exam-project/assets/130226484/17063864-7d72-453a-a2b3-a52f8a994cbf)
## Levenshtein Distance
What is the average levenshtein distance for the let-7 miRNAs for each species?
e.g., for human, we have the following let-7 miRNAs

```
>hsa-let-7a-5p
>hsa-let-7a-3p
>hsa-let-7a-2-
>hsa-let-7b-5p
>hsa-let-7b-3p
>hsa-let-7c-5p
>hsa-let-7c-3p
>hsa-let-7d-5p
>hsa-let-7d-3p
>hsa-let-7e-5p
>hsa-let-7e-3p
>hsa-let-7f-5p
```
##ANSWER:
Average Levenshtein distance for let-7 miRNA in species 'cel-': 15.00
Average Levenshtein distance for let-7 miRNA in species 'hsa-': 9.78
Average Levenshtein distance for let-7 miRNA in species 'mmu-': 9.90
Average Levenshtein distance for let-7 miRNA in species 'dme-': 14.00
Average Levenshtein distance for let-7 miRNA in species 'rno-': 9.66
Average Levenshtein distance for let-7 miRNA in species 'gga-': 9.43
Average Levenshtein distance for let-7 miRNA in species 'dre-': 6.35
Average Levenshtein distance for let-7 miRNA in species 'ssc-': 8.27
Average Levenshtein distance for let-7 miRNA in species 'fru-': 3.50
Average Levenshtein distance for let-7 miRNA in species 'tni-': 3.50
Average Levenshtein distance for let-7 miRNA in species 'bta-': 5.44
Average Levenshtein distance for let-7 miRNA in species 'xtr-': 3.14
Average Levenshtein distance for let-7 miRNA in species 'bmo-': 14.00
Average Levenshtein distance for let-7 miRNA in species 'sme-': 12.14
Average Levenshtein distance for let-7 miRNA in species 'mdo-': 9.73
Average Levenshtein distance for let-7 miRNA in species 'oan-': 9.54
Average Levenshtein distance for let-7 miRNA in species 'odi-': 6.50
Average Levenshtein distance for let-7 miRNA in species 'cin-': 10.13
Average Levenshtein distance for let-7 miRNA in species 'csa-': 3.67
Average Levenshtein distance for let-7 miRNA in species 'mml-': 9.80
Average Levenshtein distance for let-7 miRNA in species 'cfa-': 6.29
Average Levenshtein distance for let-7 miRNA in species 'ptr-': 2.82
Average Levenshtein distance for let-7 miRNA in species 'tca-': 15.00
Average Levenshtein distance for let-7 miRNA in species 'bfl-': 11.67
Average Levenshtein distance for let-7 miRNA in species 'sma-': 17.00
Average Levenshtein distance for let-7 miRNA in species 'eca-': 2.27
Average Levenshtein distance for let-7 miRNA in species 'cqu-': 13.00
Average Levenshtein distance for let-7 miRNA in species 'tgu-': 9.50
Average Levenshtein distance for let-7 miRNA in species 'oar-': 3.67
Average Levenshtein distance for let-7 miRNA in species 'ppy-': 2.82
Average Levenshtein distance for let-7 miRNA in species 'pma-': 7.60
Average Levenshtein distance for let-7 miRNA in species 'egr-': 14.00
Average Levenshtein distance for let-7 miRNA in species 'emu-': 14.00
Average Levenshtein distance for let-7 miRNA in species 'asu-': 15.00
Average Levenshtein distance for let-7 miRNA in species 'aca-': 9.47
Average Levenshtein distance for let-7 miRNA in species 'ola-': 6.19
Average Levenshtein distance for let-7 miRNA in species 'sha-': 2.67
Average Levenshtein distance for let-7 miRNA in species 'cgr-': 10.02
Average Levenshtein distance for let-7 miRNA in species 'ggo-': 2.57
Average Levenshtein distance for let-7 miRNA in species 'pol-': 10.13
Average Levenshtein distance for let-7 miRNA in species 'ccr-': 3.30
Average Levenshtein distance for let-7 miRNA in species 'ipu-': 3.13
Average Levenshtein distance for let-7 miRNA in species 'hhi-': 4.00
Average Levenshtein distance for let-7 miRNA in species 'prd-': 16.00
Average Levenshtein distance for let-7 miRNA in species 'bbe-': 12.00
Average Levenshtein distance for let-7 miRNA in species 'ssa-': 9.41
Average Levenshtein distance for let-7 miRNA in species 'efu-': 5.33
Average Levenshtein distance for let-7 miRNA in species 'gsa-': 17.00
Average Levenshtein distance for let-7 miRNA in species 'cpi-': 9.50
Average Levenshtein distance for let-7 miRNA in species 'ami-': 9.58
Average Levenshtein distance for let-7 miRNA in species 'cli-': 9.56
Average Levenshtein distance for let-7 miRNA in species 'pbv-': 9.39
Average Levenshtein distance for let-7 miRNA in species 'chi-': 9.73
Average Levenshtein distance for let-7 miRNA in species 'tch-': 6.80
Average Levenshtein distance for let-7 miRNA in species 'oha-': 9.96
Average Levenshtein distance for let-7 miRNA in species 'cja-': 2.82
Average Levenshtein distance for let-7 miRNA in species 'pal-': 10.01
Average Levenshtein distance for let-7 miRNA in species 'hpo-': 13.00
Average Levenshtein distance for let-7 miRNA in species 'tcf-': 14.00
Average Levenshtein distance for let-7 miRNA in species 'abu-': 2.83
Average Levenshtein distance for let-7 miRNA in species 'mze-': 2.89
Average Levenshtein distance for let-7 miRNA in species 'nbr-': 2.89
Average Levenshtein distance for let-7 miRNA in species 'oni-': 2.78
Average Levenshtein distance for let-7 miRNA in species 'pny-': 2.83
Average Levenshtein distance for let-7 miRNA in species 'gmo-': 9.30
Average Levenshtein distance for let-7 miRNA in species 'dqu-': 15.00
Average Levenshtein distance for let-7 miRNA in species 'pca-': 15.00
Average Levenshtein distance for let-7 miRNA in species 'pte-': 13.00
Average Levenshtein distance for let-7 miRNA in species 'xla-': 9.93
Average Levenshtein distance for let-7 miRNA in species 'cpo-': 9.78
Average Levenshtein distance for let-7 miRNA in species 'dno-': 9.68
Average Levenshtein distance for let-7 miRNA in species 'ocu-': 9.72
Average Levenshtein distance for let-7 miRNA in species 'mle-': 13.00
Average Levenshtein distance for let-7 miRNA in species 'ppa-': 2.57
Average Levenshtein distance for let-7 miRNA in species 'mmr-': 2.10
Average Levenshtein distance for let-7 miRNA in species 'dma-': 2.10
Average Levenshtein distance for let-7 miRNA in species 'nle-': 2.53
Average Levenshtein distance for let-7 miRNA in species 'sbo-': 2.00
Average Levenshtein distance for let-7 miRNA in species 'pha-': 2.57
Average Levenshtein distance for let-7 miRNA in species 'oga-': 2.57
#### Question 4: what is the average Levenshtein distance among all pairs for human?
#ANSWER:
let-7 miRNA codes:
hsa-let-7a-5p MIMAT0000062 Homo sapiens let-7a-5p
hsa-let-7a-3p MIMAT0004481 Homo sapiens let-7a-3p
hsa-let-7a-2-3p MIMAT0010195 Homo sapiens let-7a-2-3p
hsa-let-7b-5p MIMAT0000063 Homo sapiens let-7b-5p
hsa-let-7b-3p MIMAT0004482 Homo sapiens let-7b-3p
hsa-let-7c-5p MIMAT0000064 Homo sapiens let-7c-5p
hsa-let-7c-3p MIMAT0026472 Homo sapiens let-7c-3p
hsa-let-7d-5p MIMAT0000065 Homo sapiens let-7d-5p
hsa-let-7d-3p MIMAT0004484 Homo sapiens let-7d-3p
hsa-let-7e-5p MIMAT0000066 Homo sapiens let-7e-5p
hsa-let-7e-3p MIMAT0004485 Homo sapiens let-7e-3p
hsa-let-7f-5p MIMAT0000067 Homo sapiens let-7f-5p
hsa-let-7f-1-3p MIMAT0004486 Homo sapiens let-7f-1-3p
hsa-let-7f-2-3p MIMAT0004487 Homo sapiens let-7f-2-3p
hsa-let-7g-5p MIMAT0000414 Homo sapiens let-7g-5p
hsa-let-7g-3p MIMAT0004584 Homo sapiens let-7g-3p
hsa-let-7i-5p MIMAT0000415 Homo sapiens let-7i-5p
hsa-let-7i-3p MIMAT0004585 Homo sapiens let-7i-3p
Pairwise Levenshtein distances for hsa-let-7a-5p MIMAT0000062 Homo sapiens let-7a-5p
Levenshtein distance between hsa-let-7a-5p MIMAT0000062 Homo sapiens let-7a-5p and hsa-let-7a-3p MIMAT0004481 Homo sapiens let-7a-3p: 14
Levenshtein distance between hsa-let-7a-5p MIMAT0000062 Homo sapiens let-7a-5p and hsa-let-7a-2-3p MIMAT0010195 Homo sapiens let-7a-2-3p: 15
Levenshtein distance between hsa-let-7a-5p MIMAT0000062 Homo sapiens let-7a-5p and hsa-let-7b-5p MIMAT0000063 Homo sapiens let-7b-5p: 2
Levenshtein distance between hsa-let-7a-5p MIMAT0000062 Homo sapiens let-7a-5p and hsa-let-7b-3p MIMAT0004482 Homo sapiens let-7b-3p: 16
Levenshtein distance between hsa-let-7a-5p MIMAT0000062 Homo sapiens let-7a-5p and hsa-let-7c-5p MIMAT0000064 Homo sapiens let-7c-5p: 1
Levenshtein distance between hsa-let-7a-5p MIMAT0000062 Homo sapiens let-7a-5p and hsa-let-7c-3p MIMAT0026472 Homo sapiens let-7c-3p: 14
Levenshtein distance between hsa-let-7a-5p MIMAT0000062 Homo sapiens let-7a-5p and hsa-let-7d-5p MIMAT0000065 Homo sapiens let-7d-5p: 2
Levenshtein distance between hsa-let-7a-5p MIMAT0000062 Homo sapiens let-7a-5p and hsa-let-7d-3p MIMAT0004484 Homo sapiens let-7d-3p: 15
Levenshtein distance between hsa-let-7a-5p MIMAT0000062 Homo sapiens let-7a-5p and hsa-let-7e-5p MIMAT0000066 Homo sapiens let-7e-5p: 1
Levenshtein distance between hsa-let-7a-5p MIMAT0000062 Homo sapiens let-7a-5p and hsa-let-7e-3p MIMAT0004485 Homo sapiens let-7e-3p: 16
Levenshtein distance between hsa-let-7a-5p MIMAT0000062 Homo sapiens let-7a-5p and hsa-let-7f-5p MIMAT0000067 Homo sapiens let-7f-5p: 1
Levenshtein distance between hsa-let-7a-5p MIMAT0000062 Homo sapiens let-7a-5p and hsa-let-7f-1-3p MIMAT0004486 Homo sapiens let-7f-1-3p: 15
Levenshtein distance between hsa-let-7a-5p MIMAT0000062 Homo sapiens let-7a-5p and hsa-let-7f-2-3p MIMAT0004487 Homo sapiens let-7f-2-3p: 14
Levenshtein distance between hsa-let-7a-5p MIMAT0000062 Homo sapiens let-7a-5p and hsa-let-7g-5p MIMAT0000414 Homo sapiens let-7g-5p: 2
Levenshtein distance between hsa-let-7a-5p MIMAT0000062 Homo sapiens let-7a-5p and hsa-let-7g-3p MIMAT0004584 Homo sapiens let-7g-3p: 15
Levenshtein distance between hsa-let-7a-5p MIMAT0000062 Homo sapiens let-7a-5p and hsa-let-7i-5p MIMAT0000415 Homo sapiens let-7i-5p: 4
Levenshtein distance between hsa-let-7a-5p MIMAT0000062 Homo sapiens let-7a-5p and hsa-let-7i-3p MIMAT0004585 Homo sapiens let-7i-3p: 12
Pairwise Levenshtein distances for hsa-let-7a-3p MIMAT0004481 Homo sapiens let-7a-3p
Levenshtein distance between hsa-let-7a-3p MIMAT0004481 Homo sapiens let-7a-3p and hsa-let-7a-5p MIMAT0000062 Homo sapiens let-7a-5p: 14
Levenshtein distance between hsa-let-7a-3p MIMAT0004481 Homo sapiens let-7a-3p and hsa-let-7a-2-3p MIMAT0010195 Homo sapiens let-7a-2-3p: 7
Levenshtein distance between hsa-let-7a-3p MIMAT0004481 Homo sapiens let-7a-3p and hsa-let-7b-5p MIMAT0000063 Homo sapiens let-7b-5p: 13
Levenshtein distance between hsa-let-7a-3p MIMAT0004481 Homo sapiens let-7a-3p and hsa-let-7b-3p MIMAT0004482 Homo sapiens let-7b-3p: 4
Levenshtein distance between hsa-let-7a-3p MIMAT0004481 Homo sapiens let-7a-3p and hsa-let-7c-5p MIMAT0000064 Homo sapiens let-7c-5p: 14
Levenshtein distance between hsa-let-7a-3p MIMAT0004481 Homo sapiens let-7a-3p and hsa-let-7c-3p MIMAT0026472 Homo sapiens let-7c-3p: 6
Levenshtein distance between hsa-let-7a-3p MIMAT0004481 Homo sapiens let-7a-3p and hsa-let-7d-5p MIMAT0000065 Homo sapiens let-7d-5p: 15
Levenshtein distance between hsa-let-7a-3p MIMAT0004481 Homo sapiens let-7a-3p and hsa-let-7d-3p MIMAT0004484 Homo sapiens let-7d-3p: 5
Levenshtein distance between hsa-let-7a-3p MIMAT0004481 Homo sapiens let-7a-3p and hsa-let-7e-5p MIMAT0000066 Homo sapiens let-7e-5p: 14
Levenshtein distance between hsa-let-7a-3p MIMAT0004481 Homo sapiens let-7a-3p and hsa-let-7e-3p MIMAT0004485 Homo sapiens let-7e-3p: 7
Levenshtein distance between hsa-let-7a-3p MIMAT0004481 Homo sapiens let-7a-3p and hsa-let-7f-5p MIMAT0000067 Homo sapiens let-7f-5p: 13
Levenshtein distance between hsa-let-7a-3p MIMAT0004481 Homo sapiens let-7a-3p and hsa-let-7f-1-3p MIMAT0004486 Homo sapiens let-7f-1-3p: 4
Levenshtein distance between hsa-let-7a-3p MIMAT0004481 Homo sapiens let-7a-3p and hsa-let-7f-2-3p MIMAT0004487 Homo sapiens let-7f-2-3p: 2
Levenshtein distance between hsa-let-7a-3p MIMAT0004481 Homo sapiens let-7a-3p and hsa-let-7g-5p MIMAT0000414 Homo sapiens let-7g-5p: 14
Levenshtein distance between hsa-let-7a-3p MIMAT0004481 Homo sapiens let-7a-3p and hsa-let-7g-3p MIMAT0004584 Homo sapiens let-7g-3p: 6
Levenshtein distance between hsa-let-7a-3p MIMAT0004481 Homo sapiens let-7a-3p and hsa-let-7i-5p MIMAT0000415 Homo sapiens let-7i-5p: 13
Levenshtein distance between hsa-let-7a-3p MIMAT0004481 Homo sapiens let-7a-3p and hsa-let-7i-3p MIMAT0004585 Homo sapiens let-7i-3p: 7
Pairwise Levenshtein distances for hsa-let-7a-2-3p MIMAT0010195 Homo sapiens let-7a-2-3p
Levenshtein distance between hsa-let-7a-2-3p MIMAT0010195 Homo sapiens let-7a-2-3p and hsa-let-7a-5p MIMAT0000062 Homo sapiens let-7a-5p: 15
Levenshtein distance between hsa-let-7a-2-3p MIMAT0010195 Homo sapiens let-7a-2-3p and hsa-let-7a-3p MIMAT0004481 Homo sapiens let-7a-3p: 7
Levenshtein distance between hsa-let-7a-2-3p MIMAT0010195 Homo sapiens let-7a-2-3p and hsa-let-7b-5p MIMAT0000063 Homo sapiens let-7b-5p: 16
Levenshtein distance between hsa-let-7a-2-3p MIMAT0010195 Homo sapiens let-7a-2-3p and hsa-let-7b-3p MIMAT0004482 Homo sapiens let-7b-3p: 6
Levenshtein distance between hsa-let-7a-2-3p MIMAT0010195 Homo sapiens let-7a-2-3p and hsa-let-7c-5p MIMAT0000064 Homo sapiens let-7c-5p: 15
Levenshtein distance between hsa-let-7a-2-3p MIMAT0010195 Homo sapiens let-7a-2-3p and hsa-let-7c-3p MIMAT0026472 Homo sapiens let-7c-3p: 2
Levenshtein distance between hsa-let-7a-2-3p MIMAT0010195 Homo sapiens let-7a-2-3p and hsa-let-7d-5p MIMAT0000065 Homo sapiens let-7d-5p: 14
Levenshtein distance between hsa-let-7a-2-3p MIMAT0010195 Homo sapiens let-7a-2-3p and hsa-let-7d-3p MIMAT0004484 Homo sapiens let-7d-3p: 7
Levenshtein distance between hsa-let-7a-2-3p MIMAT0010195 Homo sapiens let-7a-2-3p and hsa-let-7e-5p MIMAT0000066 Homo sapiens let-7e-5p: 15
Levenshtein distance between hsa-let-7a-2-3p MIMAT0010195 Homo sapiens let-7a-2-3p and hsa-let-7e-3p MIMAT0004485 Homo sapiens let-7e-3p: 2
Levenshtein distance between hsa-let-7a-2-3p MIMAT0010195 Homo sapiens let-7a-2-3p and hsa-let-7f-5p MIMAT0000067 Homo sapiens let-7f-5p: 15
Levenshtein distance between hsa-let-7a-2-3p MIMAT0010195 Homo sapiens let-7a-2-3p and hsa-let-7f-1-3p MIMAT0004486 Homo sapiens let-7f-1-3p: 8
Levenshtein distance between hsa-let-7a-2-3p MIMAT0010195 Homo sapiens let-7a-2-3p and hsa-let-7f-2-3p MIMAT0004487 Homo sapiens let-7f-2-3p: 5
Levenshtein distance between hsa-let-7a-2-3p MIMAT0010195 Homo sapiens let-7a-2-3p and hsa-let-7g-5p MIMAT0000414 Homo sapiens let-7g-5p: 15
Levenshtein distance between hsa-let-7a-2-3p MIMAT0010195 Homo sapiens let-7a-2-3p and hsa-let-7g-3p MIMAT0004584 Homo sapiens let-7g-3p: 6
Levenshtein distance between hsa-let-7a-2-3p MIMAT0010195 Homo sapiens let-7a-2-3p and hsa-let-7i-5p MIMAT0000415 Homo sapiens let-7i-5p: 14
Levenshtein distance between hsa-let-7a-2-3p MIMAT0010195 Homo sapiens let-7a-2-3p and hsa-let-7i-3p MIMAT0004585 Homo sapiens let-7i-3p: 9
Pairwise Levenshtein distances for hsa-let-7b-5p MIMAT0000063 Homo sapiens let-7b-5p
Levenshtein distance between hsa-let-7b-5p MIMAT0000063 Homo sapiens let-7b-5p and hsa-let-7a-5p MIMAT0000062 Homo sapiens let-7a-5p: 2
Levenshtein distance between hsa-let-7b-5p MIMAT0000063 Homo sapiens let-7b-5p and hsa-let-7a-3p MIMAT0004481 Homo sapiens let-7a-3p: 13
Levenshtein distance between hsa-let-7b-5p MIMAT0000063 Homo sapiens let-7b-5p and hsa-let-7a-2-3p MIMAT0010195 Homo sapiens let-7a-2-3p: 16
Levenshtein distance between hsa-let-7b-5p MIMAT0000063 Homo sapiens let-7b-5p and hsa-let-7b-3p MIMAT0004482 Homo sapiens let-7b-3p: 16
Levenshtein distance between hsa-let-7b-5p MIMAT0000063 Homo sapiens let-7b-5p and hsa-let-7c-5p MIMAT0000064 Homo sapiens let-7c-5p: 1
Levenshtein distance between hsa-let-7b-5p MIMAT0000063 Homo sapiens let-7b-5p and hsa-let-7c-3p MIMAT0026472 Homo sapiens let-7c-3p: 15
Levenshtein distance between hsa-let-7b-5p MIMAT0000063 Homo sapiens let-7b-5p and hsa-let-7d-5p MIMAT0000065 Homo sapiens let-7d-5p: 4
Levenshtein distance between hsa-let-7b-5p MIMAT0000063 Homo sapiens let-7b-5p and hsa-let-7d-3p MIMAT0004484 Homo sapiens let-7d-3p: 15
Levenshtein distance between hsa-let-7b-5p MIMAT0000063 Homo sapiens let-7b-5p and hsa-let-7e-5p MIMAT0000066 Homo sapiens let-7e-5p: 3
Levenshtein distance between hsa-let-7b-5p MIMAT0000063 Homo sapiens let-7b-5p and hsa-let-7e-3p MIMAT0004485 Homo sapiens let-7e-3p: 16
Levenshtein distance between hsa-let-7b-5p MIMAT0000063 Homo sapiens let-7b-5p and hsa-let-7f-5p MIMAT0000067 Homo sapiens let-7f-5p: 3
Levenshtein distance between hsa-let-7b-5p MIMAT0000063 Homo sapiens let-7b-5p and hsa-let-7f-1-3p MIMAT0004486 Homo sapiens let-7f-1-3p: 15
Levenshtein distance between hsa-let-7b-5p MIMAT0000063 Homo sapiens let-7b-5p and hsa-let-7f-2-3p MIMAT0004487 Homo sapiens let-7f-2-3p: 14
Levenshtein distance between hsa-let-7b-5p MIMAT0000063 Homo sapiens let-7b-5p and hsa-let-7g-5p MIMAT0000414 Homo sapiens let-7g-5p: 4
Levenshtein distance between hsa-let-7b-5p MIMAT0000063 Homo sapiens let-7b-5p and hsa-let-7g-3p MIMAT0004584 Homo sapiens let-7g-3p: 15
Levenshtein distance between hsa-let-7b-5p MIMAT0000063 Homo sapiens let-7b-5p and hsa-let-7i-5p MIMAT0000415 Homo sapiens let-7i-5p: 3
Levenshtein distance between hsa-let-7b-5p MIMAT0000063 Homo sapiens let-7b-5p and hsa-let-7i-3p MIMAT0004585 Homo sapiens let-7i-3p: 12
Pairwise Levenshtein distances for hsa-let-7b-3p MIMAT0004482 Homo sapiens let-7b-3p
Levenshtein distance between hsa-let-7b-3p MIMAT0004482 Homo sapiens let-7b-3p and hsa-let-7a-5p MIMAT0000062 Homo sapiens let-7a-5p: 16
Levenshtein distance between hsa-let-7b-3p MIMAT0004482 Homo sapiens let-7b-3p and hsa-let-7a-3p MIMAT0004481 Homo sapiens let-7a-3p: 4
Levenshtein distance between hsa-let-7b-3p MIMAT0004482 Homo sapiens let-7b-3p and hsa-let-7a-2-3p MIMAT0010195 Homo sapiens let-7a-2-3p: 6
Levenshtein distance between hsa-let-7b-3p MIMAT0004482 Homo sapiens let-7b-3p and hsa-let-7b-5p MIMAT0000063 Homo sapiens let-7b-5p: 16
Levenshtein distance between hsa-let-7b-3p MIMAT0004482 Homo sapiens let-7b-3p and hsa-let-7c-5p MIMAT0000064 Homo sapiens let-7c-5p: 16
Levenshtein distance between hsa-let-7b-3p MIMAT0004482 Homo sapiens let-7b-3p and hsa-let-7c-3p MIMAT0026472 Homo sapiens let-7c-3p: 5
Levenshtein distance between hsa-let-7b-3p MIMAT0004482 Homo sapiens let-7b-3p and hsa-let-7d-5p MIMAT0000065 Homo sapiens let-7d-5p: 15
Levenshtein distance between hsa-let-7b-3p MIMAT0004482 Homo sapiens let-7b-3p and hsa-let-7d-3p MIMAT0004484 Homo sapiens let-7d-3p: 4
Levenshtein distance between hsa-let-7b-3p MIMAT0004482 Homo sapiens let-7b-3p and hsa-let-7e-5p MIMAT0000066 Homo sapiens let-7e-5p: 16
Levenshtein distance between hsa-let-7b-3p MIMAT0004482 Homo sapiens let-7b-3p and hsa-let-7e-3p MIMAT0004485 Homo sapiens let-7e-3p: 6
Levenshtein distance between hsa-let-7b-3p MIMAT0004482 Homo sapiens let-7b-3p and hsa-let-7f-5p MIMAT0000067 Homo sapiens let-7f-5p: 16
Levenshtein distance between hsa-let-7b-3p MIMAT0004482 Homo sapiens let-7b-3p and hsa-let-7f-1-3p MIMAT0004486 Homo sapiens let-7f-1-3p: 2
Levenshtein distance between hsa-let-7b-3p MIMAT0004482 Homo sapiens let-7b-3p and hsa-let-7f-2-3p MIMAT0004487 Homo sapiens let-7f-2-3p: 4
Levenshtein distance between hsa-let-7b-3p MIMAT0004482 Homo sapiens let-7b-3p and hsa-let-7g-5p MIMAT0000414 Homo sapiens let-7g-5p: 16
Levenshtein distance between hsa-let-7b-3p MIMAT0004482 Homo sapiens let-7b-3p and hsa-let-7g-3p MIMAT0004584 Homo sapiens let-7g-3p: 6
Levenshtein distance between hsa-let-7b-3p MIMAT0004482 Homo sapiens let-7b-3p and hsa-let-7i-5p MIMAT0000415 Homo sapiens let-7i-5p: 15
Levenshtein distance between hsa-let-7b-3p MIMAT0004482 Homo sapiens let-7b-3p and hsa-let-7i-3p MIMAT0004585 Homo sapiens let-7i-3p: 6
Pairwise Levenshtein distances for hsa-let-7c-5p MIMAT0000064 Homo sapiens let-7c-5p
Levenshtein distance between hsa-let-7c-5p MIMAT0000064 Homo sapiens let-7c-5p and hsa-let-7a-5p MIMAT0000062 Homo sapiens let-7a-5p: 1
Levenshtein distance between hsa-let-7c-5p MIMAT0000064 Homo sapiens let-7c-5p and hsa-let-7a-3p MIMAT0004481 Homo sapiens let-7a-3p: 14
Levenshtein distance between hsa-let-7c-5p MIMAT0000064 Homo sapiens let-7c-5p and hsa-let-7a-2-3p MIMAT0010195 Homo sapiens let-7a-2-3p: 15
Levenshtein distance between hsa-let-7c-5p MIMAT0000064 Homo sapiens let-7c-5p and hsa-let-7b-5p MIMAT0000063 Homo sapiens let-7b-5p: 1
Levenshtein distance between hsa-let-7c-5p MIMAT0000064 Homo sapiens let-7c-5p and hsa-let-7b-3p MIMAT0004482 Homo sapiens let-7b-3p: 16
Levenshtein distance between hsa-let-7c-5p MIMAT0000064 Homo sapiens let-7c-5p and hsa-let-7c-3p MIMAT0026472 Homo sapiens let-7c-3p: 14
Levenshtein distance between hsa-let-7c-5p MIMAT0000064 Homo sapiens let-7c-5p and hsa-let-7d-5p MIMAT0000065 Homo sapiens let-7d-5p: 3
Levenshtein distance between hsa-let-7c-5p MIMAT0000064 Homo sapiens let-7c-5p and hsa-let-7d-3p MIMAT0004484 Homo sapiens let-7d-3p: 15
Levenshtein distance between hsa-let-7c-5p MIMAT0000064 Homo sapiens let-7c-5p and hsa-let-7e-5p MIMAT0000066 Homo sapiens let-7e-5p: 2
Levenshtein distance between hsa-let-7c-5p MIMAT0000064 Homo sapiens let-7c-5p and hsa-let-7e-3p MIMAT0004485 Homo sapiens let-7e-3p: 16
Levenshtein distance between hsa-let-7c-5p MIMAT0000064 Homo sapiens let-7c-5p and hsa-let-7f-5p MIMAT0000067 Homo sapiens let-7f-5p: 2
Levenshtein distance between hsa-let-7c-5p MIMAT0000064 Homo sapiens let-7c-5p and hsa-let-7f-1-3p MIMAT0004486 Homo sapiens let-7f-1-3p: 15
Levenshtein distance between hsa-let-7c-5p MIMAT0000064 Homo sapiens let-7c-5p and hsa-let-7f-2-3p MIMAT0004487 Homo sapiens let-7f-2-3p: 14
Levenshtein distance between hsa-let-7c-5p MIMAT0000064 Homo sapiens let-7c-5p and hsa-let-7g-5p MIMAT0000414 Homo sapiens let-7g-5p: 3
Levenshtein distance between hsa-let-7c-5p MIMAT0000064 Homo sapiens let-7c-5p and hsa-let-7g-3p MIMAT0004584 Homo sapiens let-7g-3p: 15
Levenshtein distance between hsa-let-7c-5p MIMAT0000064 Homo sapiens let-7c-5p and hsa-let-7i-5p MIMAT0000415 Homo sapiens let-7i-5p: 4
Levenshtein distance between hsa-let-7c-5p MIMAT0000064 Homo sapiens let-7c-5p and hsa-let-7i-3p MIMAT0004585 Homo sapiens let-7i-3p: 12
Pairwise Levenshtein distances for hsa-let-7c-3p MIMAT0026472 Homo sapiens let-7c-3p
Levenshtein distance between hsa-let-7c-3p MIMAT0026472 Homo sapiens let-7c-3p and hsa-let-7a-5p MIMAT0000062 Homo sapiens let-7a-5p: 14
Levenshtein distance between hsa-let-7c-3p MIMAT0026472 Homo sapiens let-7c-3p and hsa-let-7a-3p MIMAT0004481 Homo sapiens let-7a-3p: 6
Levenshtein distance between hsa-let-7c-3p MIMAT0026472 Homo sapiens let-7c-3p and hsa-let-7a-2-3p MIMAT0010195 Homo sapiens let-7a-2-3p: 2
Levenshtein distance between hsa-let-7c-3p MIMAT0026472 Homo sapiens let-7c-3p and hsa-let-7b-5p MIMAT0000063 Homo sapiens let-7b-5p: 15
Levenshtein distance between hsa-let-7c-3p MIMAT0026472 Homo sapiens let-7c-3p and hsa-let-7b-3p MIMAT0004482 Homo sapiens let-7b-3p: 5
Levenshtein distance between hsa-let-7c-3p MIMAT0026472 Homo sapiens let-7c-3p and hsa-let-7c-5p MIMAT0000064 Homo sapiens let-7c-5p: 14
Levenshtein distance between hsa-let-7c-3p MIMAT0026472 Homo sapiens let-7c-3p and hsa-let-7d-5p MIMAT0000065 Homo sapiens let-7d-5p: 14
Levenshtein distance between hsa-let-7c-3p MIMAT0026472 Homo sapiens let-7c-3p and hsa-let-7d-3p MIMAT0004484 Homo sapiens let-7d-3p: 6
Levenshtein distance between hsa-let-7c-3p MIMAT0026472 Homo sapiens let-7c-3p and hsa-let-7e-5p MIMAT0000066 Homo sapiens let-7e-5p: 14
Levenshtein distance between hsa-let-7c-3p MIMAT0026472 Homo sapiens let-7c-3p and hsa-let-7e-3p MIMAT0004485 Homo sapiens let-7e-3p: 4
Levenshtein distance between hsa-let-7c-3p MIMAT0026472 Homo sapiens let-7c-3p and hsa-let-7f-5p MIMAT0000067 Homo sapiens let-7f-5p: 14
Levenshtein distance between hsa-let-7c-3p MIMAT0026472 Homo sapiens let-7c-3p and hsa-let-7f-1-3p MIMAT0004486 Homo sapiens let-7f-1-3p: 7
Levenshtein distance between hsa-let-7c-3p MIMAT0026472 Homo sapiens let-7c-3p and hsa-let-7f-2-3p MIMAT0004487 Homo sapiens let-7f-2-3p: 6
Levenshtein distance between hsa-let-7c-3p MIMAT0026472 Homo sapiens let-7c-3p and hsa-let-7g-5p MIMAT0000414 Homo sapiens let-7g-5p: 14
Levenshtein distance between hsa-let-7c-3p MIMAT0026472 Homo sapiens let-7c-3p and hsa-let-7g-3p MIMAT0004584 Homo sapiens let-7g-3p: 7
Levenshtein distance between hsa-let-7c-3p MIMAT0026472 Homo sapiens let-7c-3p and hsa-let-7i-5p MIMAT0000415 Homo sapiens let-7i-5p: 13
Levenshtein distance between hsa-let-7c-3p MIMAT0026472 Homo sapiens let-7c-3p and hsa-let-7i-3p MIMAT0004585 Homo sapiens let-7i-3p: 8
Pairwise Levenshtein distances for hsa-let-7d-5p MIMAT0000065 Homo sapiens let-7d-5p
Levenshtein distance between hsa-let-7d-5p MIMAT0000065 Homo sapiens let-7d-5p and hsa-let-7a-5p MIMAT0000062 Homo sapiens let-7a-5p: 2
Levenshtein distance between hsa-let-7d-5p MIMAT0000065 Homo sapiens let-7d-5p and hsa-let-7a-3p MIMAT0004481 Homo sapiens let-7a-3p: 15
Levenshtein distance between hsa-let-7d-5p MIMAT0000065 Homo sapiens let-7d-5p and hsa-let-7a-2-3p MIMAT0010195 Homo sapiens let-7a-2-3p: 14
Levenshtein distance between hsa-let-7d-5p MIMAT0000065 Homo sapiens let-7d-5p and hsa-let-7b-5p MIMAT0000063 Homo sapiens let-7b-5p: 4
Levenshtein distance between hsa-let-7d-5p MIMAT0000065 Homo sapiens let-7d-5p and hsa-let-7b-3p MIMAT0004482 Homo sapiens let-7b-3p: 15
Levenshtein distance between hsa-let-7d-5p MIMAT0000065 Homo sapiens let-7d-5p and hsa-let-7c-5p MIMAT0000064 Homo sapiens let-7c-5p: 3
Levenshtein distance between hsa-let-7d-5p MIMAT0000065 Homo sapiens let-7d-5p and hsa-let-7c-3p MIMAT0026472 Homo sapiens let-7c-3p: 14
Levenshtein distance between hsa-let-7d-5p MIMAT0000065 Homo sapiens let-7d-5p and hsa-let-7d-3p MIMAT0004484 Homo sapiens let-7d-3p: 14
Levenshtein distance between hsa-let-7d-5p MIMAT0000065 Homo sapiens let-7d-5p and hsa-let-7e-5p MIMAT0000066 Homo sapiens let-7e-5p: 3
Levenshtein distance between hsa-let-7d-5p MIMAT0000065 Homo sapiens let-7d-5p and hsa-let-7e-3p MIMAT0004485 Homo sapiens let-7e-3p: 15
Levenshtein distance between hsa-let-7d-5p MIMAT0000065 Homo sapiens let-7d-5p and hsa-let-7f-5p MIMAT0000067 Homo sapiens let-7f-5p: 3
Levenshtein distance between hsa-let-7d-5p MIMAT0000065 Homo sapiens let-7d-5p and hsa-let-7f-1-3p MIMAT0004486 Homo sapiens let-7f-1-3p: 14
Levenshtein distance between hsa-let-7d-5p MIMAT0000065 Homo sapiens let-7d-5p and hsa-let-7f-2-3p MIMAT0004487 Homo sapiens let-7f-2-3p: 15
Levenshtein distance between hsa-let-7d-5p MIMAT0000065 Homo sapiens let-7d-5p and hsa-let-7g-5p MIMAT0000414 Homo sapiens let-7g-5p: 4
Levenshtein distance between hsa-let-7d-5p MIMAT0000065 Homo sapiens let-7d-5p and hsa-let-7g-3p MIMAT0004584 Homo sapiens let-7g-3p: 14
Levenshtein distance between hsa-let-7d-5p MIMAT0000065 Homo sapiens let-7d-5p and hsa-let-7i-5p MIMAT0000415 Homo sapiens let-7i-5p: 6
Levenshtein distance between hsa-let-7d-5p MIMAT0000065 Homo sapiens let-7d-5p and hsa-let-7i-3p MIMAT0004585 Homo sapiens let-7i-3p: 12
Pairwise Levenshtein distances for hsa-let-7d-3p MIMAT0004484 Homo sapiens let-7d-3p
Levenshtein distance between hsa-let-7d-3p MIMAT0004484 Homo sapiens let-7d-3p and hsa-let-7a-5p MIMAT0000062 Homo sapiens let-7a-5p: 15
Levenshtein distance between hsa-let-7d-3p MIMAT0004484 Homo sapiens let-7d-3p and hsa-let-7a-3p MIMAT0004481 Homo sapiens let-7a-3p: 5
Levenshtein distance between hsa-let-7d-3p MIMAT0004484 Homo sapiens let-7d-3p and hsa-let-7a-2-3p MIMAT0010195 Homo sapiens let-7a-2-3p: 7
Levenshtein distance between hsa-let-7d-3p MIMAT0004484 Homo sapiens let-7d-3p and hsa-let-7b-5p MIMAT0000063 Homo sapiens let-7b-5p: 15
Levenshtein distance between hsa-let-7d-3p MIMAT0004484 Homo sapiens let-7d-3p and hsa-let-7b-3p MIMAT0004482 Homo sapiens let-7b-3p: 4
Levenshtein distance between hsa-let-7d-3p MIMAT0004484 Homo sapiens let-7d-3p and hsa-let-7c-5p MIMAT0000064 Homo sapiens let-7c-5p: 15
Levenshtein distance between hsa-let-7d-3p MIMAT0004484 Homo sapiens let-7d-3p and hsa-let-7c-3p MIMAT0026472 Homo sapiens let-7c-3p: 6
Levenshtein distance between hsa-let-7d-3p MIMAT0004484 Homo sapiens let-7d-3p and hsa-let-7d-5p MIMAT0000065 Homo sapiens let-7d-5p: 14
Levenshtein distance between hsa-let-7d-3p MIMAT0004484 Homo sapiens let-7d-3p and hsa-let-7e-5p MIMAT0000066 Homo sapiens let-7e-5p: 14
Levenshtein distance between hsa-let-7d-3p MIMAT0004484 Homo sapiens let-7d-3p and hsa-let-7e-3p MIMAT0004485 Homo sapiens let-7e-3p: 5
Levenshtein distance between hsa-let-7d-3p MIMAT0004484 Homo sapiens let-7d-3p and hsa-let-7f-5p MIMAT0000067 Homo sapiens let-7f-5p: 15
Levenshtein distance between hsa-let-7d-3p MIMAT0004484 Homo sapiens let-7d-3p and hsa-let-7f-1-3p MIMAT0004486 Homo sapiens let-7f-1-3p: 6
Levenshtein distance between hsa-let-7d-3p MIMAT0004484 Homo sapiens let-7d-3p and hsa-let-7f-2-3p MIMAT0004487 Homo sapiens let-7f-2-3p: 6
Levenshtein distance between hsa-let-7d-3p MIMAT0004484 Homo sapiens let-7d-3p and hsa-let-7g-5p MIMAT0000414 Homo sapiens let-7g-5p: 15
Levenshtein distance between hsa-let-7d-3p MIMAT0004484 Homo sapiens let-7d-3p and hsa-let-7g-3p MIMAT0004584 Homo sapiens let-7g-3p: 7
Levenshtein distance between hsa-let-7d-3p MIMAT0004484 Homo sapiens let-7d-3p and hsa-let-7i-5p MIMAT0000415 Homo sapiens let-7i-5p: 13
Levenshtein distance between hsa-let-7d-3p MIMAT0004484 Homo sapiens let-7d-3p and hsa-let-7i-3p MIMAT0004585 Homo sapiens let-7i-3p: 7
Pairwise Levenshtein distances for hsa-let-7e-5p MIMAT0000066 Homo sapiens let-7e-5p
Levenshtein distance between hsa-let-7e-5p MIMAT0000066 Homo sapiens let-7e-5p and hsa-let-7a-5p MIMAT0000062 Homo sapiens let-7a-5p: 1
Levenshtein distance between hsa-let-7e-5p MIMAT0000066 Homo sapiens let-7e-5p and hsa-let-7a-3p MIMAT0004481 Homo sapiens let-7a-3p: 14
Levenshtein distance between hsa-let-7e-5p MIMAT0000066 Homo sapiens let-7e-5p and hsa-let-7a-2-3p MIMAT0010195 Homo sapiens let-7a-2-3p: 15
Levenshtein distance between hsa-let-7e-5p MIMAT0000066 Homo sapiens let-7e-5p and hsa-let-7b-5p MIMAT0000063 Homo sapiens let-7b-5p: 3
Levenshtein distance between hsa-let-7e-5p MIMAT0000066 Homo sapiens let-7e-5p and hsa-let-7b-3p MIMAT0004482 Homo sapiens let-7b-3p: 16
Levenshtein distance between hsa-let-7e-5p MIMAT0000066 Homo sapiens let-7e-5p and hsa-let-7c-5p MIMAT0000064 Homo sapiens let-7c-5p: 2
Levenshtein distance between hsa-let-7e-5p MIMAT0000066 Homo sapiens let-7e-5p and hsa-let-7c-3p MIMAT0026472 Homo sapiens let-7c-3p: 14
Levenshtein distance between hsa-let-7e-5p MIMAT0000066 Homo sapiens let-7e-5p and hsa-let-7d-5p MIMAT0000065 Homo sapiens let-7d-5p: 3
Levenshtein distance between hsa-let-7e-5p MIMAT0000066 Homo sapiens let-7e-5p and hsa-let-7d-3p MIMAT0004484 Homo sapiens let-7d-3p: 14
Levenshtein distance between hsa-let-7e-5p MIMAT0000066 Homo sapiens let-7e-5p and hsa-let-7e-3p MIMAT0004485 Homo sapiens let-7e-3p: 15
Levenshtein distance between hsa-let-7e-5p MIMAT0000066 Homo sapiens let-7e-5p and hsa-let-7f-5p MIMAT0000067 Homo sapiens let-7f-5p: 2
Levenshtein distance between hsa-let-7e-5p MIMAT0000066 Homo sapiens let-7e-5p and hsa-let-7f-1-3p MIMAT0004486 Homo sapiens let-7f-1-3p: 16
Levenshtein distance between hsa-let-7e-5p MIMAT0000066 Homo sapiens let-7e-5p and hsa-let-7f-2-3p MIMAT0004487 Homo sapiens let-7f-2-3p: 15
Levenshtein distance between hsa-let-7e-5p MIMAT0000066 Homo sapiens let-7e-5p and hsa-let-7g-5p MIMAT0000414 Homo sapiens let-7g-5p: 3
Levenshtein distance between hsa-let-7e-5p MIMAT0000066 Homo sapiens let-7e-5p and hsa-let-7g-3p MIMAT0004584 Homo sapiens let-7g-3p: 14
Levenshtein distance between hsa-let-7e-5p MIMAT0000066 Homo sapiens let-7e-5p and hsa-let-7i-5p MIMAT0000415 Homo sapiens let-7i-5p: 5
Levenshtein distance between hsa-let-7e-5p MIMAT0000066 Homo sapiens let-7e-5p and hsa-let-7i-3p MIMAT0004585 Homo sapiens let-7i-3p: 13
Pairwise Levenshtein distances for hsa-let-7e-3p MIMAT0004485 Homo sapiens let-7e-3p
Levenshtein distance between hsa-let-7e-3p MIMAT0004485 Homo sapiens let-7e-3p and hsa-let-7a-5p MIMAT0000062 Homo sapiens let-7a-5p: 16
Levenshtein distance between hsa-let-7e-3p MIMAT0004485 Homo sapiens let-7e-3p and hsa-let-7a-3p MIMAT0004481 Homo sapiens let-7a-3p: 7
Levenshtein distance between hsa-let-7e-3p MIMAT0004485 Homo sapiens let-7e-3p and hsa-let-7a-2-3p MIMAT0010195 Homo sapiens let-7a-2-3p: 2
Levenshtein distance between hsa-let-7e-3p MIMAT0004485 Homo sapiens let-7e-3p and hsa-let-7b-5p MIMAT0000063 Homo sapiens let-7b-5p: 16
Levenshtein distance between hsa-let-7e-3p MIMAT0004485 Homo sapiens let-7e-3p and hsa-let-7b-3p MIMAT0004482 Homo sapiens let-7b-3p: 6
Levenshtein distance between hsa-let-7e-3p MIMAT0004485 Homo sapiens let-7e-3p and hsa-let-7c-5p MIMAT0000064 Homo sapiens let-7c-5p: 16
Levenshtein distance between hsa-let-7e-3p MIMAT0004485 Homo sapiens let-7e-3p and hsa-let-7c-3p MIMAT0026472 Homo sapiens let-7c-3p: 4
Levenshtein distance between hsa-let-7e-3p MIMAT0004485 Homo sapiens let-7e-3p and hsa-let-7d-5p MIMAT0000065 Homo sapiens let-7d-5p: 15
Levenshtein distance between hsa-let-7e-3p MIMAT0004485 Homo sapiens let-7e-3p and hsa-let-7d-3p MIMAT0004484 Homo sapiens let-7d-3p: 5
Levenshtein distance between hsa-let-7e-3p MIMAT0004485 Homo sapiens let-7e-3p and hsa-let-7e-5p MIMAT0000066 Homo sapiens let-7e-5p: 15
Levenshtein distance between hsa-let-7e-3p MIMAT0004485 Homo sapiens let-7e-3p and hsa-let-7f-5p MIMAT0000067 Homo sapiens let-7f-5p: 16
Levenshtein distance between hsa-let-7e-3p MIMAT0004485 Homo sapiens let-7e-3p and hsa-let-7f-1-3p MIMAT0004486 Homo sapiens let-7f-1-3p: 8
Levenshtein distance between hsa-let-7e-3p MIMAT0004485 Homo sapiens let-7e-3p and hsa-let-7f-2-3p MIMAT0004487 Homo sapiens let-7f-2-3p: 5
Levenshtein distance between hsa-let-7e-3p MIMAT0004485 Homo sapiens let-7e-3p and hsa-let-7g-5p MIMAT0000414 Homo sapiens let-7g-5p: 15
Levenshtein distance between hsa-let-7e-3p MIMAT0004485 Homo sapiens let-7e-3p and hsa-let-7g-3p MIMAT0004584 Homo sapiens let-7g-3p: 7
Levenshtein distance between hsa-let-7e-3p MIMAT0004485 Homo sapiens let-7e-3p and hsa-let-7i-5p MIMAT0000415 Homo sapiens let-7i-5p: 15
Levenshtein distance between hsa-let-7e-3p MIMAT0004485 Homo sapiens let-7e-3p and hsa-let-7i-3p MIMAT0004585 Homo sapiens let-7i-3p: 11
Pairwise Levenshtein distances for hsa-let-7f-5p MIMAT0000067 Homo sapiens let-7f-5p
Levenshtein distance between hsa-let-7f-5p MIMAT0000067 Homo sapiens let-7f-5p and hsa-let-7a-5p MIMAT0000062 Homo sapiens let-7a-5p: 1
Levenshtein distance between hsa-let-7f-5p MIMAT0000067 Homo sapiens let-7f-5p and hsa-let-7a-3p MIMAT0004481 Homo sapiens let-7a-3p: 13
Levenshtein distance between hsa-let-7f-5p MIMAT0000067 Homo sapiens let-7f-5p and hsa-let-7a-2-3p MIMAT0010195 Homo sapiens let-7a-2-3p: 15
Levenshtein distance between hsa-let-7f-5p MIMAT0000067 Homo sapiens let-7f-5p and hsa-let-7b-5p MIMAT0000063 Homo sapiens let-7b-5p: 3
Levenshtein distance between hsa-let-7f-5p MIMAT0000067 Homo sapiens let-7f-5p and hsa-let-7b-3p MIMAT0004482 Homo sapiens let-7b-3p: 16
Levenshtein distance between hsa-let-7f-5p MIMAT0000067 Homo sapiens let-7f-5p and hsa-let-7c-5p MIMAT0000064 Homo sapiens let-7c-5p: 2
Levenshtein distance between hsa-let-7f-5p MIMAT0000067 Homo sapiens let-7f-5p and hsa-let-7c-3p MIMAT0026472 Homo sapiens let-7c-3p: 14
Levenshtein distance between hsa-let-7f-5p MIMAT0000067 Homo sapiens let-7f-5p and hsa-let-7d-5p MIMAT0000065 Homo sapiens let-7d-5p: 3
Levenshtein distance between hsa-let-7f-5p MIMAT0000067 Homo sapiens let-7f-5p and hsa-let-7d-3p MIMAT0004484 Homo sapiens let-7d-3p: 15
Levenshtein distance between hsa-let-7f-5p MIMAT0000067 Homo sapiens let-7f-5p and hsa-let-7e-5p MIMAT0000066 Homo sapiens let-7e-5p: 2
Levenshtein distance between hsa-let-7f-5p MIMAT0000067 Homo sapiens let-7f-5p and hsa-let-7e-3p MIMAT0004485 Homo sapiens let-7e-3p: 16
Levenshtein distance between hsa-let-7f-5p MIMAT0000067 Homo sapiens let-7f-5p and hsa-let-7f-1-3p MIMAT0004486 Homo sapiens let-7f-1-3p: 14
Levenshtein distance between hsa-let-7f-5p MIMAT0000067 Homo sapiens let-7f-5p and hsa-let-7f-2-3p MIMAT0004487 Homo sapiens let-7f-2-3p: 13
Levenshtein distance between hsa-let-7f-5p MIMAT0000067 Homo sapiens let-7f-5p and hsa-let-7g-5p MIMAT0000414 Homo sapiens let-7g-5p: 2
Levenshtein distance between hsa-let-7f-5p MIMAT0000067 Homo sapiens let-7f-5p and hsa-let-7g-3p MIMAT0004584 Homo sapiens let-7g-3p: 14
Levenshtein distance between hsa-let-7f-5p MIMAT0000067 Homo sapiens let-7f-5p and hsa-let-7i-5p MIMAT0000415 Homo sapiens let-7i-5p: 4
Levenshtein distance between hsa-let-7f-5p MIMAT0000067 Homo sapiens let-7f-5p and hsa-let-7i-3p MIMAT0004585 Homo sapiens let-7i-3p: 12
Pairwise Levenshtein distances for hsa-let-7f-1-3p MIMAT0004486 Homo sapiens let-7f-1-3p
Levenshtein distance between hsa-let-7f-1-3p MIMAT0004486 Homo sapiens let-7f-1-3p and hsa-let-7a-5p MIMAT0000062 Homo sapiens let-7a-5p: 15
Levenshtein distance between hsa-let-7f-1-3p MIMAT0004486 Homo sapiens let-7f-1-3p and hsa-let-7a-3p MIMAT0004481 Homo sapiens let-7a-3p: 4
Levenshtein distance between hsa-let-7f-1-3p MIMAT0004486 Homo sapiens let-7f-1-3p and hsa-let-7a-2-3p MIMAT0010195 Homo sapiens let-7a-2-3p: 8
Levenshtein distance between hsa-let-7f-1-3p MIMAT0004486 Homo sapiens let-7f-1-3p and hsa-let-7b-5p MIMAT0000063 Homo sapiens let-7b-5p: 15
Levenshtein distance between hsa-let-7f-1-3p MIMAT0004486 Homo sapiens let-7f-1-3p and hsa-let-7b-3p MIMAT0004482 Homo sapiens let-7b-3p: 2
Levenshtein distance between hsa-let-7f-1-3p MIMAT0004486 Homo sapiens let-7f-1-3p and hsa-let-7c-5p MIMAT0000064 Homo sapiens let-7c-5p: 15
Levenshtein distance between hsa-let-7f-1-3p MIMAT0004486 Homo sapiens let-7f-1-3p and hsa-let-7c-3p MIMAT0026472 Homo sapiens let-7c-3p: 7
Levenshtein distance between hsa-let-7f-1-3p MIMAT0004486 Homo sapiens let-7f-1-3p and hsa-let-7d-5p MIMAT0000065 Homo sapiens let-7d-5p: 14
Levenshtein distance between hsa-let-7f-1-3p MIMAT0004486 Homo sapiens let-7f-1-3p and hsa-let-7d-3p MIMAT0004484 Homo sapiens let-7d-3p: 6
Levenshtein distance between hsa-let-7f-1-3p MIMAT0004486 Homo sapiens let-7f-1-3p and hsa-let-7e-5p MIMAT0000066 Homo sapiens let-7e-5p: 16
Levenshtein distance between hsa-let-7f-1-3p MIMAT0004486 Homo sapiens let-7f-1-3p and hsa-let-7e-3p MIMAT0004485 Homo sapiens let-7e-3p: 8
Levenshtein distance between hsa-let-7f-1-3p MIMAT0004486 Homo sapiens let-7f-1-3p and hsa-let-7f-5p MIMAT0000067 Homo sapiens let-7f-5p: 14
Levenshtein distance between hsa-let-7f-1-3p MIMAT0004486 Homo sapiens let-7f-1-3p and hsa-let-7f-2-3p MIMAT0004487 Homo sapiens let-7f-2-3p: 4
Levenshtein distance between hsa-let-7f-1-3p MIMAT0004486 Homo sapiens let-7f-1-3p and hsa-let-7g-5p MIMAT0000414 Homo sapiens let-7g-5p: 16
Levenshtein distance between hsa-let-7f-1-3p MIMAT0004486 Homo sapiens let-7f-1-3p and hsa-let-7g-3p MIMAT0004584 Homo sapiens let-7g-3p: 7
Levenshtein distance between hsa-let-7f-1-3p MIMAT0004486 Homo sapiens let-7f-1-3p and hsa-let-7i-5p MIMAT0000415 Homo sapiens let-7i-5p: 15
Levenshtein distance between hsa-let-7f-1-3p MIMAT0004486 Homo sapiens let-7f-1-3p and hsa-let-7i-3p MIMAT0004585 Homo sapiens let-7i-3p: 7
Pairwise Levenshtein distances for hsa-let-7f-2-3p MIMAT0004487 Homo sapiens let-7f-2-3p
Levenshtein distance between hsa-let-7f-2-3p MIMAT0004487 Homo sapiens let-7f-2-3p and hsa-let-7a-5p MIMAT0000062 Homo sapiens let-7a-5p: 14
Levenshtein distance between hsa-let-7f-2-3p MIMAT0004487 Homo sapiens let-7f-2-3p and hsa-let-7a-3p MIMAT0004481 Homo sapiens let-7a-3p: 2
Levenshtein distance between hsa-let-7f-2-3p MIMAT0004487 Homo sapiens let-7f-2-3p and hsa-let-7a-2-3p MIMAT0010195 Homo sapiens let-7a-2-3p: 5
Levenshtein distance between hsa-let-7f-2-3p MIMAT0004487 Homo sapiens let-7f-2-3p and hsa-let-7b-5p MIMAT0000063 Homo sapiens let-7b-5p: 14
Levenshtein distance between hsa-let-7f-2-3p MIMAT0004487 Homo sapiens let-7f-2-3p and hsa-let-7b-3p MIMAT0004482 Homo sapiens let-7b-3p: 4
Levenshtein distance between hsa-let-7f-2-3p MIMAT0004487 Homo sapiens let-7f-2-3p and hsa-let-7c-5p MIMAT0000064 Homo sapiens let-7c-5p: 14
Levenshtein distance between hsa-let-7f-2-3p MIMAT0004487 Homo sapiens let-7f-2-3p and hsa-let-7c-3p MIMAT0026472 Homo sapiens let-7c-3p: 6
Levenshtein distance between hsa-let-7f-2-3p MIMAT0004487 Homo sapiens let-7f-2-3p and hsa-let-7d-5p MIMAT0000065 Homo sapiens let-7d-5p: 15
Levenshtein distance between hsa-let-7f-2-3p MIMAT0004487 Homo sapiens let-7f-2-3p and hsa-let-7d-3p MIMAT0004484 Homo sapiens let-7d-3p: 6
Levenshtein distance between hsa-let-7f-2-3p MIMAT0004487 Homo sapiens let-7f-2-3p and hsa-let-7e-5p MIMAT0000066 Homo sapiens let-7e-5p: 15
Levenshtein distance between hsa-let-7f-2-3p MIMAT0004487 Homo sapiens let-7f-2-3p and hsa-let-7e-3p MIMAT0004485 Homo sapiens let-7e-3p: 5
Levenshtein distance between hsa-let-7f-2-3p MIMAT0004487 Homo sapiens let-7f-2-3p and hsa-let-7f-5p MIMAT0000067 Homo sapiens let-7f-5p: 13
Levenshtein distance between hsa-let-7f-2-3p MIMAT0004487 Homo sapiens let-7f-2-3p and hsa-let-7f-1-3p MIMAT0004486 Homo sapiens let-7f-1-3p: 4
Levenshtein distance between hsa-let-7f-2-3p MIMAT0004487 Homo sapiens let-7f-2-3p and hsa-let-7g-5p MIMAT0000414 Homo sapiens let-7g-5p: 14
Levenshtein distance between hsa-let-7f-2-3p MIMAT0004487 Homo sapiens let-7f-2-3p and hsa-let-7g-3p MIMAT0004584 Homo sapiens let-7g-3p: 6
Levenshtein distance between hsa-let-7f-2-3p MIMAT0004487 Homo sapiens let-7f-2-3p and hsa-let-7i-5p MIMAT0000415 Homo sapiens let-7i-5p: 13
Levenshtein distance between hsa-let-7f-2-3p MIMAT0004487 Homo sapiens let-7f-2-3p and hsa-let-7i-3p MIMAT0004585 Homo sapiens let-7i-3p: 8
Pairwise Levenshtein distances for hsa-let-7g-5p MIMAT0000414 Homo sapiens let-7g-5p
Levenshtein distance between hsa-let-7g-5p MIMAT0000414 Homo sapiens let-7g-5p and hsa-let-7a-5p MIMAT0000062 Homo sapiens let-7a-5p: 2
Levenshtein distance between hsa-let-7g-5p MIMAT0000414 Homo sapiens let-7g-5p and hsa-let-7a-3p MIMAT0004481 Homo sapiens let-7a-3p: 14
Levenshtein distance between hsa-let-7g-5p MIMAT0000414 Homo sapiens let-7g-5p and hsa-let-7a-2-3p MIMAT0010195 Homo sapiens let-7a-2-3p: 15
Levenshtein distance between hsa-let-7g-5p MIMAT0000414 Homo sapiens let-7g-5p and hsa-let-7b-5p MIMAT0000063 Homo sapiens let-7b-5p: 4
Levenshtein distance between hsa-let-7g-5p MIMAT0000414 Homo sapiens let-7g-5p and hsa-let-7b-3p MIMAT0004482 Homo sapiens let-7b-3p: 16
Levenshtein distance between hsa-let-7g-5p MIMAT0000414 Homo sapiens let-7g-5p and hsa-let-7c-5p MIMAT0000064 Homo sapiens let-7c-5p: 3
Levenshtein distance between hsa-let-7g-5p MIMAT0000414 Homo sapiens let-7g-5p and hsa-let-7c-3p MIMAT0026472 Homo sapiens let-7c-3p: 14
Levenshtein distance between hsa-let-7g-5p MIMAT0000414 Homo sapiens let-7g-5p and hsa-let-7d-5p MIMAT0000065 Homo sapiens let-7d-5p: 4
Levenshtein distance between hsa-let-7g-5p MIMAT0000414 Homo sapiens let-7g-5p and hsa-let-7d-3p MIMAT0004484 Homo sapiens let-7d-3p: 15
Levenshtein distance between hsa-let-7g-5p MIMAT0000414 Homo sapiens let-7g-5p and hsa-let-7e-5p MIMAT0000066 Homo sapiens let-7e-5p: 3
Levenshtein distance between hsa-let-7g-5p MIMAT0000414 Homo sapiens let-7g-5p and hsa-let-7e-3p MIMAT0004485 Homo sapiens let-7e-3p: 15
Levenshtein distance between hsa-let-7g-5p MIMAT0000414 Homo sapiens let-7g-5p and hsa-let-7f-5p MIMAT0000067 Homo sapiens let-7f-5p: 2
Levenshtein distance between hsa-let-7g-5p MIMAT0000414 Homo sapiens let-7g-5p and hsa-let-7f-1-3p MIMAT0004486 Homo sapiens let-7f-1-3p: 16
Levenshtein distance between hsa-let-7g-5p MIMAT0000414 Homo sapiens let-7g-5p and hsa-let-7f-2-3p MIMAT0004487 Homo sapiens let-7f-2-3p: 14
Levenshtein distance between hsa-let-7g-5p MIMAT0000414 Homo sapiens let-7g-5p and hsa-let-7g-3p MIMAT0004584 Homo sapiens let-7g-3p: 16
Levenshtein distance between hsa-let-7g-5p MIMAT0000414 Homo sapiens let-7g-5p and hsa-let-7i-5p MIMAT0000415 Homo sapiens let-7i-5p: 2
Levenshtein distance between hsa-let-7g-5p MIMAT0000414 Homo sapiens let-7g-5p and hsa-let-7i-3p MIMAT0004585 Homo sapiens let-7i-3p: 13
Pairwise Levenshtein distances for hsa-let-7g-3p MIMAT0004584 Homo sapiens let-7g-3p
Levenshtein distance between hsa-let-7g-3p MIMAT0004584 Homo sapiens let-7g-3p and hsa-let-7a-5p MIMAT0000062 Homo sapiens let-7a-5p: 15
Levenshtein distance between hsa-let-7g-3p MIMAT0004584 Homo sapiens let-7g-3p and hsa-let-7a-3p MIMAT0004481 Homo sapiens let-7a-3p: 6
Levenshtein distance between hsa-let-7g-3p MIMAT0004584 Homo sapiens let-7g-3p and hsa-let-7a-2-3p MIMAT0010195 Homo sapiens let-7a-2-3p: 6
Levenshtein distance between hsa-let-7g-3p MIMAT0004584 Homo sapiens let-7g-3p and hsa-let-7b-5p MIMAT0000063 Homo sapiens let-7b-5p: 15
Levenshtein distance between hsa-let-7g-3p MIMAT0004584 Homo sapiens let-7g-3p and hsa-let-7b-3p MIMAT0004482 Homo sapiens let-7b-3p: 6
Levenshtein distance between hsa-let-7g-3p MIMAT0004584 Homo sapiens let-7g-3p and hsa-let-7c-5p MIMAT0000064 Homo sapiens let-7c-5p: 15
Levenshtein distance between hsa-let-7g-3p MIMAT0004584 Homo sapiens let-7g-3p and hsa-let-7c-3p MIMAT0026472 Homo sapiens let-7c-3p: 7
Levenshtein distance between hsa-let-7g-3p MIMAT0004584 Homo sapiens let-7g-3p and hsa-let-7d-5p MIMAT0000065 Homo sapiens let-7d-5p: 14
Levenshtein distance between hsa-let-7g-3p MIMAT0004584 Homo sapiens let-7g-3p and hsa-let-7d-3p MIMAT0004484 Homo sapiens let-7d-3p: 7
Levenshtein distance between hsa-let-7g-3p MIMAT0004584 Homo sapiens let-7g-3p and hsa-let-7e-5p MIMAT0000066 Homo sapiens let-7e-5p: 14
Levenshtein distance between hsa-let-7g-3p MIMAT0004584 Homo sapiens let-7g-3p and hsa-let-7e-3p MIMAT0004485 Homo sapiens let-7e-3p: 7
Levenshtein distance between hsa-let-7g-3p MIMAT0004584 Homo sapiens let-7g-3p and hsa-let-7f-5p MIMAT0000067 Homo sapiens let-7f-5p: 14
Levenshtein distance between hsa-let-7g-3p MIMAT0004584 Homo sapiens let-7g-3p and hsa-let-7f-1-3p MIMAT0004486 Homo sapiens let-7f-1-3p: 7
Levenshtein distance between hsa-let-7g-3p MIMAT0004584 Homo sapiens let-7g-3p and hsa-let-7f-2-3p MIMAT0004487 Homo sapiens let-7f-2-3p: 6
Levenshtein distance between hsa-let-7g-3p MIMAT0004584 Homo sapiens let-7g-3p and hsa-let-7g-5p MIMAT0000414 Homo sapiens let-7g-5p: 16
Levenshtein distance between hsa-let-7g-3p MIMAT0004584 Homo sapiens let-7g-3p and hsa-let-7i-5p MIMAT0000415 Homo sapiens let-7i-5p: 15
Levenshtein distance between hsa-let-7g-3p MIMAT0004584 Homo sapiens let-7g-3p and hsa-let-7i-3p MIMAT0004585 Homo sapiens let-7i-3p: 5
Pairwise Levenshtein distances for hsa-let-7i-5p MIMAT0000415 Homo sapiens let-7i-5p
Levenshtein distance between hsa-let-7i-5p MIMAT0000415 Homo sapiens let-7i-5p and hsa-let-7a-5p MIMAT0000062 Homo sapiens let-7a-5p: 4
Levenshtein distance between hsa-let-7i-5p MIMAT0000415 Homo sapiens let-7i-5p and hsa-let-7a-3p MIMAT0004481 Homo sapiens let-7a-3p: 13
Levenshtein distance between hsa-let-7i-5p MIMAT0000415 Homo sapiens let-7i-5p and hsa-let-7a-2-3p MIMAT0010195 Homo sapiens let-7a-2-3p: 14
Levenshtein distance between hsa-let-7i-5p MIMAT0000415 Homo sapiens let-7i-5p and hsa-let-7b-5p MIMAT0000063 Homo sapiens let-7b-5p: 3
Levenshtein distance between hsa-let-7i-5p MIMAT0000415 Homo sapiens let-7i-5p and hsa-let-7b-3p MIMAT0004482 Homo sapiens let-7b-3p: 15
Levenshtein distance between hsa-let-7i-5p MIMAT0000415 Homo sapiens let-7i-5p and hsa-let-7c-5p MIMAT0000064 Homo sapiens let-7c-5p: 4
Levenshtein distance between hsa-let-7i-5p MIMAT0000415 Homo sapiens let-7i-5p and hsa-let-7c-3p MIMAT0026472 Homo sapiens let-7c-3p: 13
Levenshtein distance between hsa-let-7i-5p MIMAT0000415 Homo sapiens let-7i-5p and hsa-let-7d-5p MIMAT0000065 Homo sapiens let-7d-5p: 6
Levenshtein distance between hsa-let-7i-5p MIMAT0000415 Homo sapiens let-7i-5p and hsa-let-7d-3p MIMAT0004484 Homo sapiens let-7d-3p: 13
Levenshtein distance between hsa-let-7i-5p MIMAT0000415 Homo sapiens let-7i-5p and hsa-let-7e-5p MIMAT0000066 Homo sapiens let-7e-5p: 5
Levenshtein distance between hsa-let-7i-5p MIMAT0000415 Homo sapiens let-7i-5p and hsa-let-7e-3p MIMAT0004485 Homo sapiens let-7e-3p: 15
Levenshtein distance between hsa-let-7i-5p MIMAT0000415 Homo sapiens let-7i-5p and hsa-let-7f-5p MIMAT0000067 Homo sapiens let-7f-5p: 4
Levenshtein distance between hsa-let-7i-5p MIMAT0000415 Homo sapiens let-7i-5p and hsa-let-7f-1-3p MIMAT0004486 Homo sapiens let-7f-1-3p: 15
Levenshtein distance between hsa-let-7i-5p MIMAT0000415 Homo sapiens let-7i-5p and hsa-let-7f-2-3p MIMAT0004487 Homo sapiens let-7f-2-3p: 13
Levenshtein distance between hsa-let-7i-5p MIMAT0000415 Homo sapiens let-7i-5p and hsa-let-7g-5p MIMAT0000414 Homo sapiens let-7g-5p: 2
Levenshtein distance between hsa-let-7i-5p MIMAT0000415 Homo sapiens let-7i-5p and hsa-let-7g-3p MIMAT0004584 Homo sapiens let-7g-3p: 15
Levenshtein distance between hsa-let-7i-5p MIMAT0000415 Homo sapiens let-7i-5p and hsa-let-7i-3p MIMAT0004585 Homo sapiens let-7i-3p: 12
Pairwise Levenshtein distances for hsa-let-7i-3p MIMAT0004585 Homo sapiens let-7i-3p
Levenshtein distance between hsa-let-7i-3p MIMAT0004585 Homo sapiens let-7i-3p and hsa-let-7a-5p MIMAT0000062 Homo sapiens let-7a-5p: 12
Levenshtein distance between hsa-let-7i-3p MIMAT0004585 Homo sapiens let-7i-3p and hsa-let-7a-3p MIMAT0004481 Homo sapiens let-7a-3p: 7
Levenshtein distance between hsa-let-7i-3p MIMAT0004585 Homo sapiens let-7i-3p and hsa-let-7a-2-3p MIMAT0010195 Homo sapiens let-7a-2-3p: 9
Levenshtein distance between hsa-let-7i-3p MIMAT0004585 Homo sapiens let-7i-3p and hsa-let-7b-5p MIMAT0000063 Homo sapiens let-7b-5p: 12
Levenshtein distance between hsa-let-7i-3p MIMAT0004585 Homo sapiens let-7i-3p and hsa-let-7b-3p MIMAT0004482 Homo sapiens let-7b-3p: 6
Levenshtein distance between hsa-let-7i-3p MIMAT0004585 Homo sapiens let-7i-3p and hsa-let-7c-5p MIMAT0000064 Homo sapiens let-7c-5p: 12
Levenshtein distance between hsa-let-7i-3p MIMAT0004585 Homo sapiens let-7i-3p and hsa-let-7c-3p MIMAT0026472 Homo sapiens let-7c-3p: 8
Levenshtein distance between hsa-let-7i-3p MIMAT0004585 Homo sapiens let-7i-3p and hsa-let-7d-5p MIMAT0000065 Homo sapiens let-7d-5p: 12
Levenshtein distance between hsa-let-7i-3p MIMAT0004585 Homo sapiens let-7i-3p and hsa-let-7d-3p MIMAT0004484 Homo sapiens let-7d-3p: 7
Levenshtein distance between hsa-let-7i-3p MIMAT0004585 Homo sapiens let-7i-3p and hsa-let-7e-5p MIMAT0000066 Homo sapiens let-7e-5p: 13
Levenshtein distance between hsa-let-7i-3p MIMAT0004585 Homo sapiens let-7i-3p and hsa-let-7e-3p MIMAT0004485 Homo sapiens let-7e-3p: 11
Levenshtein distance between hsa-let-7i-3p MIMAT0004585 Homo sapiens let-7i-3p and hsa-let-7f-5p MIMAT0000067 Homo sapiens let-7f-5p: 12
Levenshtein distance between hsa-let-7i-3p MIMAT0004585 Homo sapiens let-7i-3p and hsa-let-7f-1-3p MIMAT0004486 Homo sapiens let-7f-1-3p: 7
Levenshtein distance between hsa-let-7i-3p MIMAT0004585 Homo sapiens let-7i-3p and hsa-let-7f-2-3p MIMAT0004487 Homo sapiens let-7f-2-3p: 8
Levenshtein distance between hsa-let-7i-3p MIMAT0004585 Homo sapiens let-7i-3p and hsa-let-7g-5p MIMAT0000414 Homo sapiens let-7g-5p: 13
Levenshtein distance between hsa-let-7i-3p MIMAT0004585 Homo sapiens let-7i-3p and hsa-let-7g-3p MIMAT0004584 Homo sapiens let-7g-3p: 5
Levenshtein distance between hsa-let-7i-3p MIMAT0004585 Homo sapiens let-7i-3p and hsa-let-7i-5p MIMAT0000415 Homo sapiens let-7i-5p: 12

# AND IT`S PLOT IS AS FOLLOWING:
![Figure_3](https://github.com/aberakenea/Home-take-exam-project/assets/130226484/dfb34d05-3163-479a-8955-e852794884f7)

#### Task: repeat this for all species and plot
#ANSWER:
##NOTE: it is possible to do for all species by editing the above code of "hsa" to get their Levenshtein distance between let_7_miRNA and their plots`, but the reason why I did not done is because my pc`s RAM can not run the code for all species and I just done only for hsa, then it is possible to do for all species if the the computer do have enuogh RAM to run the code for all species instead of "hsa". finally we will get their Levenshtein distance between let_7_miRNA and  plot for all species. 

#### Question 5: What is the levenshtein distance for each let-7 miRNA across all species?
## for example, let-7a is present in 121 species, what is the average Levenshtein distance among all pairs?
## ANSWER: 
## the code for let-7e miRNA across all species is:
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
## I have done for code of levenshtein distance of let-7e miRNA across all species and get the following results:
Species: cin-let-7e - pbv-let-7e-5p | Levenshtein Distance: 5
Species: cin-let-7e - pbv-let-7e-3p | Levenshtein Distance: 15
Species: cin-let-7e - chi-let-7e-5p | Levenshtein Distance: 4
Species: cin-let-7e - chi-let-7e-3p | Levenshtein Distance: 16
Species: cin-let-7e - tch-let-7e-5p | Levenshtein Distance: 4
Species: cin-let-7e - oha-let-7e-5p | Levenshtein Distance: 5
Species: cin-let-7e - oha-let-7e-3p | Levenshtein Distance: 15
Species: cin-let-7e - pal-let-7e-5p | Levenshtein Distance: 4
Species: cin-let-7e - pal-let-7e-3p | Levenshtein Distance: 15
Species: cin-let-7e - cgr-let-7e | Levenshtein Distance: 16
Species: cin-let-7e - abu-let-7e | Levenshtein Distance: 5
Species: cin-let-7e - mze-let-7e | Levenshtein Distance: 5
Species: cin-let-7e - nbr-let-7e | Levenshtein Distance: 5
Species: cin-let-7e - oni-let-7e | Levenshtein Distance: 5
Species: cin-let-7e - pny-let-7e | Levenshtein Distance: 5
Species: cin-let-7e - gmo-let-7e-5p | Levenshtein Distance: 5
Species: cin-let-7e - gmo-let-7e-3p | Levenshtein Distance: 15
Species: cin-let-7e - xla-let-7e-5p | Levenshtein Distance: 4
Species: cin-let-7e - xla-let-7e-3p | Levenshtein Distance: 15
Species: cin-let-7e - cpo-let-7e-5p | Levenshtein Distance: 4
Species: cin-let-7e - cpo-let-7e-3p | Levenshtein Distance: 16
Species: cin-let-7e - ppa-let-7e | Levenshtein Distance: 4
Species: cin-let-7e - cja-let-7e | Levenshtein Distance: 4
Species: cin-let-7e - nle-let-7e | Levenshtein Distance: 4
Species: cin-let-7e - sbo-let-7e | Levenshtein Distance: 4
Species: cin-let-7e - pha-let-7e | Levenshtein Distance: 4
Species: cin-let-7e - oga-let-7e | Levenshtein Distance: 4
Species: mml-let-7e-5p - mml-let-7e-3p | Levenshtein Distance: 15
Species: mml-let-7e-5p - cfa-let-7e | Levenshtein Distance: 0
Species: mml-let-7e-5p - ptr-let-7e | Levenshtein Distance: 0
Species: mml-let-7e-5p - eca-let-7e | Levenshtein Distance: 0
Species: mml-let-7e-5p - ssc-let-7e | Levenshtein Distance: 0
Species: mml-let-7e-5p - tgu-let-7e-5p | Levenshtein Distance: 3
Species: mml-let-7e-5p - tgu-let-7e-3p | Levenshtein Distance: 15
Species: mml-let-7e-5p - ppy-let-7e | Levenshtein Distance: 0
Species: mml-let-7e-5p - aca-let-7e-5p | Levenshtein Distance: 3
Species: mml-let-7e-5p - aca-let-7e-3p | Levenshtein Distance: 14
Species: mml-let-7e-5p - ola-let-7e | Levenshtein Distance: 6
Species: mml-let-7e-5p - ggo-let-7e | Levenshtein Distance: 1
Species: mml-let-7e-5p - ipu-let-7e | Levenshtein Distance: 3
Species: mml-let-7e-5p - ssa-let-7e-5p | Levenshtein Distance: 3
Species: mml-let-7e-5p - ssa-let-7e-3p | Levenshtein Distance: 15
Species: mml-let-7e-5p - efu-let-7e | Levenshtein Distance: 2
Species: mml-let-7e-5p - cpi-let-7e-5p | Levenshtein Distance: 3
Species: mml-let-7e-5p - cpi-let-7e-3p | Levenshtein Distance: 13
Species: mml-let-7e-5p - ami-let-7e-5p | Levenshtein Distance: 3
Species: mml-let-7e-5p - ami-let-7e-3p | Levenshtein Distance: 13
Species: mml-let-7e-5p - cli-let-7e-5p | Levenshtein Distance: 3
Species: mml-let-7e-5p - cli-let-7e-3p | Levenshtein Distance: 15
Species: mml-let-7e-5p - pbv-let-7e-5p | Levenshtein Distance: 3
Species: mml-let-7e-5p - pbv-let-7e-3p | Levenshtein Distance: 15
Species: mml-let-7e-5p - chi-let-7e-5p | Levenshtein Distance: 0
Species: mml-let-7e-5p - chi-let-7e-3p | Levenshtein Distance: 15
Species: mml-let-7e-5p - tch-let-7e-5p | Levenshtein Distance: 0
Species: mml-let-7e-5p - oha-let-7e-5p | Levenshtein Distance: 3
Species: mml-let-7e-5p - oha-let-7e-3p | Levenshtein Distance: 14
Species: mml-let-7e-5p - pal-let-7e-5p | Levenshtein Distance: 0
Species: mml-let-7e-5p - pal-let-7e-3p | Levenshtein Distance: 14
Species: mml-let-7e-5p - cgr-let-7e | Levenshtein Distance: 15
Species: mml-let-7e-5p - abu-let-7e | Levenshtein Distance: 3
Species: mml-let-7e-5p - mze-let-7e | Levenshtein Distance: 3
Species: mml-let-7e-5p - nbr-let-7e | Levenshtein Distance: 3
Species: mml-let-7e-5p - oni-let-7e | Levenshtein Distance: 3
Species: mml-let-7e-5p - pny-let-7e | Levenshtein Distance: 3
Species: mml-let-7e-5p - gmo-let-7e-5p | Levenshtein Distance: 3
Species: mml-let-7e-5p - gmo-let-7e-3p | Levenshtein Distance: 15
Species: mml-let-7e-5p - xla-let-7e-5p | Levenshtein Distance: 2
Species: mml-let-7e-5p - xla-let-7e-3p | Levenshtein Distance: 15
Species: mml-let-7e-5p - cpo-let-7e-5p | Levenshtein Distance: 0
Species: mml-let-7e-5p - cpo-let-7e-3p | Levenshtein Distance: 15
Species: mml-let-7e-5p - ppa-let-7e | Levenshtein Distance: 0
Species: mml-let-7e-5p - cja-let-7e | Levenshtein Distance: 0
Species: mml-let-7e-5p - nle-let-7e | Levenshtein Distance: 0
Species: mml-let-7e-5p - sbo-let-7e | Levenshtein Distance: 0
Species: mml-let-7e-5p - pha-let-7e | Levenshtein Distance: 0
Species: mml-let-7e-5p - oga-let-7e | Levenshtein Distance: 0
Species: mml-let-7e-3p - cfa-let-7e | Levenshtein Distance: 15
Species: mml-let-7e-3p - ptr-let-7e | Levenshtein Distance: 15
Species: mml-let-7e-3p - eca-let-7e | Levenshtein Distance: 15
Species: mml-let-7e-3p - ssc-let-7e | Levenshtein Distance: 15
Species: mml-let-7e-3p - tgu-let-7e-5p | Levenshtein Distance: 16
Species: mml-let-7e-3p - tgu-let-7e-3p | Levenshtein Distance: 6
Species: mml-let-7e-3p - ppy-let-7e | Levenshtein Distance: 15
Species: mml-let-7e-3p - aca-let-7e-5p | Levenshtein Distance: 16
Species: mml-let-7e-3p - aca-let-7e-3p | Levenshtein Distance: 10
Species: mml-let-7e-3p - ola-let-7e | Levenshtein Distance: 13
Species: mml-let-7e-3p - ggo-let-7e | Levenshtein Distance: 16
Species: mml-let-7e-3p - ipu-let-7e | Levenshtein Distance: 16
Species: mml-let-7e-3p - ssa-let-7e-5p | Levenshtein Distance: 16
Species: mml-let-7e-3p - ssa-let-7e-3p | Levenshtein Distance: 6
Species: mml-let-7e-3p - efu-let-7e | Levenshtein Distance: 15
Species: mml-let-7e-3p - cpi-let-7e-5p | Levenshtein Distance: 16
Species: mml-let-7e-3p - cpi-let-7e-3p | Levenshtein Distance: 8
Species: mml-let-7e-3p - ami-let-7e-5p | Levenshtein Distance: 16
Species: mml-let-7e-3p - ami-let-7e-3p | Levenshtein Distance: 8
Species: mml-let-7e-3p - cli-let-7e-5p | Levenshtein Distance: 16
Species: mml-let-7e-3p - cli-let-7e-3p | Levenshtein Distance: 6
Species: mml-let-7e-3p - pbv-let-7e-5p | Levenshtein Distance: 16
Species: mml-let-7e-3p - pbv-let-7e-3p | Levenshtein Distance: 6
Species: mml-let-7e-3p - chi-let-7e-5p | Levenshtein Distance: 15
Species: mml-let-7e-3p - chi-let-7e-3p | Levenshtein Distance: 0
Species: mml-let-7e-3p - tch-let-7e-5p | Levenshtein Distance: 15
Species: mml-let-7e-3p - oha-let-7e-5p | Levenshtein Distance: 16
Species: mml-let-7e-3p - oha-let-7e-3p | Levenshtein Distance: 10
Species: mml-let-7e-3p - pal-let-7e-5p | Levenshtein Distance: 15
Species: mml-let-7e-3p - pal-let-7e-3p | Levenshtein Distance: 1
Species: mml-let-7e-3p - cgr-let-7e | Levenshtein Distance: 0
Species: mml-let-7e-3p - abu-let-7e | Levenshtein Distance: 16
Species: mml-let-7e-3p - mze-let-7e | Levenshtein Distance: 16
Species: mml-let-7e-3p - nbr-let-7e | Levenshtein Distance: 16
Species: mml-let-7e-3p - oni-let-7e | Levenshtein Distance: 16
Species: mml-let-7e-3p - pny-let-7e | Levenshtein Distance: 16
Species: mml-let-7e-3p - gmo-let-7e-5p | Levenshtein Distance: 16
Species: mml-let-7e-3p - gmo-let-7e-3p | Levenshtein Distance: 6
Species: mml-let-7e-3p - xla-let-7e-5p | Levenshtein Distance: 16
Species: mml-let-7e-3p - xla-let-7e-3p | Levenshtein Distance: 6
Species: mml-let-7e-3p - cpo-let-7e-5p | Levenshtein Distance: 15
Species: mml-let-7e-3p - cpo-let-7e-3p | Levenshtein Distance: 0
Species: mml-let-7e-3p - ppa-let-7e | Levenshtein Distance: 15
Species: mml-let-7e-3p - cja-let-7e | Levenshtein Distance: 15
Species: mml-let-7e-3p - nle-let-7e | Levenshtein Distance: 15
Species: mml-let-7e-3p - sbo-let-7e | Levenshtein Distance: 15
Species: mml-let-7e-3p - pha-let-7e | Levenshtein Distance: 15
Species: mml-let-7e-3p - oga-let-7e | Levenshtein Distance: 15
Species: cfa-let-7e - ptr-let-7e | Levenshtein Distance: 0
Species: cfa-let-7e - eca-let-7e | Levenshtein Distance: 0
Species: cfa-let-7e - ssc-let-7e | Levenshtein Distance: 0
Species: cfa-let-7e - tgu-let-7e-5p | Levenshtein Distance: 3
Species: cfa-let-7e - tgu-let-7e-3p | Levenshtein Distance: 15
Species: cfa-let-7e - ppy-let-7e | Levenshtein Distance: 0
Species: cfa-let-7e - aca-let-7e-5p | Levenshtein Distance: 3
Species: cfa-let-7e - aca-let-7e-3p | Levenshtein Distance: 14
Species: cfa-let-7e - ola-let-7e | Levenshtein Distance: 6
Species: cfa-let-7e - ggo-let-7e | Levenshtein Distance: 1
Species: cfa-let-7e - ipu-let-7e | Levenshtein Distance: 3
Species: cfa-let-7e - ssa-let-7e-5p | Levenshtein Distance: 3
Species: cfa-let-7e - ssa-let-7e-3p | Levenshtein Distance: 15
Species: cfa-let-7e - efu-let-7e | Levenshtein Distance: 2
Species: cfa-let-7e - cpi-let-7e-5p | Levenshtein Distance: 3
Species: cfa-let-7e - cpi-let-7e-3p | Levenshtein Distance: 13
Species: cfa-let-7e - ami-let-7e-5p | Levenshtein Distance: 3
Species: cfa-let-7e - ami-let-7e-3p | Levenshtein Distance: 13
Species: cfa-let-7e - cli-let-7e-5p | Levenshtein Distance: 3
Species: cfa-let-7e - cli-let-7e-3p | Levenshtein Distance: 15
Species: cfa-let-7e - pbv-let-7e-5p | Levenshtein Distance: 3
Species: cfa-let-7e - pbv-let-7e-3p | Levenshtein Distance: 15
Species: cfa-let-7e - chi-let-7e-5p | Levenshtein Distance: 0
Species: cfa-let-7e - chi-let-7e-3p | Levenshtein Distance: 15
Species: cfa-let-7e - tch-let-7e-5p | Levenshtein Distance: 0
Species: cfa-let-7e - oha-let-7e-5p | Levenshtein Distance: 3
Species: cfa-let-7e - oha-let-7e-3p | Levenshtein Distance: 14
Species: cfa-let-7e - pal-let-7e-5p | Levenshtein Distance: 0
Species: cfa-let-7e - pal-let-7e-3p | Levenshtein Distance: 14
Species: cfa-let-7e - cgr-let-7e | Levenshtein Distance: 15
Species: cfa-let-7e - abu-let-7e | Levenshtein Distance: 3
Species: cfa-let-7e - mze-let-7e | Levenshtein Distance: 3
Species: cfa-let-7e - nbr-let-7e | Levenshtein Distance: 3
Species: cfa-let-7e - oni-let-7e | Levenshtein Distance: 3
Species: cfa-let-7e - pny-let-7e | Levenshtein Distance: 3
Species: cfa-let-7e - gmo-let-7e-5p | Levenshtein Distance: 3
Species: cfa-let-7e - gmo-let-7e-3p | Levenshtein Distance: 15
Species: cfa-let-7e - xla-let-7e-5p | Levenshtein Distance: 2
Species: cfa-let-7e - xla-let-7e-3p | Levenshtein Distance: 15
Species: cfa-let-7e - cpo-let-7e-5p | Levenshtein Distance: 0
Species: cfa-let-7e - cpo-let-7e-3p | Levenshtein Distance: 15
Species: cfa-let-7e - ppa-let-7e | Levenshtein Distance: 0
Species: cfa-let-7e - cja-let-7e | Levenshtein Distance: 0
Species: cfa-let-7e - nle-let-7e | Levenshtein Distance: 0
Species: cfa-let-7e - sbo-let-7e | Levenshtein Distance: 0
Species: cfa-let-7e - pha-let-7e | Levenshtein Distance: 0
Species: cfa-let-7e - oga-let-7e | Levenshtein Distance: 0
Species: ptr-let-7e - eca-let-7e | Levenshtein Distance: 0
Species: ptr-let-7e - ssc-let-7e | Levenshtein Distance: 0
Species: ptr-let-7e - tgu-let-7e-5p | Levenshtein Distance: 3
Species: ptr-let-7e - tgu-let-7e-3p | Levenshtein Distance: 15
Species: ptr-let-7e - ppy-let-7e | Levenshtein Distance: 0
Species: ptr-let-7e - aca-let-7e-5p | Levenshtein Distance: 3
Species: ptr-let-7e - aca-let-7e-3p | Levenshtein Distance: 14
Species: ptr-let-7e - ola-let-7e | Levenshtein Distance: 6
Species: ptr-let-7e - ggo-let-7e | Levenshtein Distance: 1
Species: ptr-let-7e - ipu-let-7e | Levenshtein Distance: 3
Species: ptr-let-7e - ssa-let-7e-5p | Levenshtein Distance: 3
Species: ptr-let-7e - ssa-let-7e-3p | Levenshtein Distance: 15
Species: ptr-let-7e - efu-let-7e | Levenshtein Distance: 2
Species: ptr-let-7e - cpi-let-7e-5p | Levenshtein Distance: 3
Species: ptr-let-7e - cpi-let-7e-3p | Levenshtein Distance: 13
Species: ptr-let-7e - ami-let-7e-5p | Levenshtein Distance: 3
Species: ptr-let-7e - ami-let-7e-3p | Levenshtein Distance: 13
Species: ptr-let-7e - cli-let-7e-5p | Levenshtein Distance: 3
Species: ptr-let-7e - cli-let-7e-3p | Levenshtein Distance: 15
Species: ptr-let-7e - pbv-let-7e-5p | Levenshtein Distance: 3
Species: ptr-let-7e - pbv-let-7e-3p | Levenshtein Distance: 15
Species: ptr-let-7e - chi-let-7e-5p | Levenshtein Distance: 0
Species: ptr-let-7e - chi-let-7e-3p | Levenshtein Distance: 15
Species: ptr-let-7e - tch-let-7e-5p | Levenshtein Distance: 0
Species: ptr-let-7e - oha-let-7e-5p | Levenshtein Distance: 3
Species: ptr-let-7e - oha-let-7e-3p | Levenshtein Distance: 14
Species: ptr-let-7e - pal-let-7e-5p | Levenshtein Distance: 0
Species: ptr-let-7e - pal-let-7e-3p | Levenshtein Distance: 14
Species: ptr-let-7e - cgr-let-7e | Levenshtein Distance: 15
Species: ptr-let-7e - abu-let-7e | Levenshtein Distance: 3
Species: ptr-let-7e - mze-let-7e | Levenshtein Distance: 3
Species: ptr-let-7e - nbr-let-7e | Levenshtein Distance: 3
Species: ptr-let-7e - oni-let-7e | Levenshtein Distance: 3
Species: ptr-let-7e - pny-let-7e | Levenshtein Distance: 3
Species: ptr-let-7e - gmo-let-7e-5p | Levenshtein Distance: 3
Species: ptr-let-7e - gmo-let-7e-3p | Levenshtein Distance: 15
Species: ptr-let-7e - xla-let-7e-5p | Levenshtein Distance: 2
Species: ptr-let-7e - xla-let-7e-3p | Levenshtein Distance: 15
Species: ptr-let-7e - cpo-let-7e-5p | Levenshtein Distance: 0
Species: ptr-let-7e - cpo-let-7e-3p | Levenshtein Distance: 15
Species: ptr-let-7e - ppa-let-7e | Levenshtein Distance: 0
Species: ptr-let-7e - cja-let-7e | Levenshtein Distance: 0
Species: ptr-let-7e - nle-let-7e | Levenshtein Distance: 0
Species: ptr-let-7e - sbo-let-7e | Levenshtein Distance: 0
Species: ptr-let-7e - pha-let-7e | Levenshtein Distance: 0
Species: ptr-let-7e - oga-let-7e | Levenshtein Distance: 0
Species: eca-let-7e - ssc-let-7e | Levenshtein Distance: 0
Species: eca-let-7e - tgu-let-7e-5p | Levenshtein Distance: 3
Species: eca-let-7e - tgu-let-7e-3p | Levenshtein Distance: 15
Species: eca-let-7e - ppy-let-7e | Levenshtein Distance: 0
Species: eca-let-7e - aca-let-7e-5p | Levenshtein Distance: 3
Species: eca-let-7e - aca-let-7e-3p | Levenshtein Distance: 14
Species: eca-let-7e - ola-let-7e | Levenshtein Distance: 6
Species: eca-let-7e - ggo-let-7e | Levenshtein Distance: 1
Species: eca-let-7e - ipu-let-7e | Levenshtein Distance: 3
Species: eca-let-7e - ssa-let-7e-5p | Levenshtein Distance: 3
Species: eca-let-7e - ssa-let-7e-3p | Levenshtein Distance: 15
Species: eca-let-7e - efu-let-7e | Levenshtein Distance: 2
Species: eca-let-7e - cpi-let-7e-5p | Levenshtein Distance: 3
Species: eca-let-7e - cpi-let-7e-3p | Levenshtein Distance: 13
Species: eca-let-7e - ami-let-7e-5p | Levenshtein Distance: 3
Species: eca-let-7e - ami-let-7e-3p | Levenshtein Distance: 13
Species: eca-let-7e - cli-let-7e-5p | Levenshtein Distance: 3
Species: eca-let-7e - cli-let-7e-3p | Levenshtein Distance: 15
Species: eca-let-7e - pbv-let-7e-5p | Levenshtein Distance: 3
Species: eca-let-7e - pbv-let-7e-3p | Levenshtein Distance: 15
Species: eca-let-7e - chi-let-7e-5p | Levenshtein Distance: 0
Species: eca-let-7e - chi-let-7e-3p | Levenshtein Distance: 15
Species: eca-let-7e - tch-let-7e-5p | Levenshtein Distance: 0
Species: eca-let-7e - oha-let-7e-5p | Levenshtein Distance: 3
Species: eca-let-7e - oha-let-7e-3p | Levenshtein Distance: 14
Species: eca-let-7e - pal-let-7e-5p | Levenshtein Distance: 0
Species: eca-let-7e - pal-let-7e-3p | Levenshtein Distance: 14
Species: eca-let-7e - cgr-let-7e | Levenshtein Distance: 15
Species: eca-let-7e - abu-let-7e | Levenshtein Distance: 3
Species: eca-let-7e - mze-let-7e | Levenshtein Distance: 3
Species: eca-let-7e - nbr-let-7e | Levenshtein Distance: 3
Species: eca-let-7e - oni-let-7e | Levenshtein Distance: 3
Species: eca-let-7e - pny-let-7e | Levenshtein Distance: 3
Species: eca-let-7e - gmo-let-7e-5p | Levenshtein Distance: 3
Species: eca-let-7e - gmo-let-7e-3p | Levenshtein Distance: 15
Species: eca-let-7e - xla-let-7e-5p | Levenshtein Distance: 2
Species: eca-let-7e - xla-let-7e-3p | Levenshtein Distance: 15
Species: eca-let-7e - cpo-let-7e-5p | Levenshtein Distance: 0
Species: eca-let-7e - cpo-let-7e-3p | Levenshtein Distance: 15
Species: eca-let-7e - ppa-let-7e | Levenshtein Distance: 0
Species: eca-let-7e - cja-let-7e | Levenshtein Distance: 0
Species: eca-let-7e - nle-let-7e | Levenshtein Distance: 0
Species: eca-let-7e - sbo-let-7e | Levenshtein Distance: 0
Species: eca-let-7e - pha-let-7e | Levenshtein Distance: 0
Species: eca-let-7e - oga-let-7e | Levenshtein Distance: 0
Species: ssc-let-7e - tgu-let-7e-5p | Levenshtein Distance: 3
Species: ssc-let-7e - tgu-let-7e-3p | Levenshtein Distance: 15
Species: ssc-let-7e - ppy-let-7e | Levenshtein Distance: 0
Species: ssc-let-7e - aca-let-7e-5p | Levenshtein Distance: 3
Species: ssc-let-7e - aca-let-7e-3p | Levenshtein Distance: 14
Species: ssc-let-7e - ola-let-7e | Levenshtein Distance: 6
Species: ssc-let-7e - ggo-let-7e | Levenshtein Distance: 1
Species: ssc-let-7e - ipu-let-7e | Levenshtein Distance: 3
Species: ssc-let-7e - ssa-let-7e-5p | Levenshtein Distance: 3
Species: ssc-let-7e - ssa-let-7e-3p | Levenshtein Distance: 15
Species: ssc-let-7e - efu-let-7e | Levenshtein Distance: 2
Species: ssc-let-7e - cpi-let-7e-5p | Levenshtein Distance: 3
Species: ssc-let-7e - cpi-let-7e-3p | Levenshtein Distance: 13
Species: ssc-let-7e - ami-let-7e-5p | Levenshtein Distance: 3
Species: ssc-let-7e - ami-let-7e-3p | Levenshtein Distance: 13
Species: ssc-let-7e - cli-let-7e-5p | Levenshtein Distance: 3
Species: ssc-let-7e - cli-let-7e-3p | Levenshtein Distance: 15
Species: ssc-let-7e - pbv-let-7e-5p | Levenshtein Distance: 3
Species: ssc-let-7e - pbv-let-7e-3p | Levenshtein Distance: 15
Species: ssc-let-7e - chi-let-7e-5p | Levenshtein Distance: 0
Species: ssc-let-7e - chi-let-7e-3p | Levenshtein Distance: 15
Species: ssc-let-7e - tch-let-7e-5p | Levenshtein Distance: 0
Species: ssc-let-7e - oha-let-7e-5p | Levenshtein Distance: 3
Species: ssc-let-7e - oha-let-7e-3p | Levenshtein Distance: 14
Species: ssc-let-7e - pal-let-7e-5p | Levenshtein Distance: 0
Species: ssc-let-7e - pal-let-7e-3p | Levenshtein Distance: 14
Species: ssc-let-7e - cgr-let-7e | Levenshtein Distance: 15
Species: ssc-let-7e - abu-let-7e | Levenshtein Distance: 3
Species: ssc-let-7e - mze-let-7e | Levenshtein Distance: 3
Species: ssc-let-7e - nbr-let-7e | Levenshtein Distance: 3
Species: ssc-let-7e - oni-let-7e | Levenshtein Distance: 3
Species: ssc-let-7e - pny-let-7e | Levenshtein Distance: 3
Species: ssc-let-7e - gmo-let-7e-5p | Levenshtein Distance: 3
Species: ssc-let-7e - gmo-let-7e-3p | Levenshtein Distance: 15
Species: ssc-let-7e - xla-let-7e-5p | Levenshtein Distance: 2
Species: ssc-let-7e - xla-let-7e-3p | Levenshtein Distance: 15
Species: ssc-let-7e - cpo-let-7e-5p | Levenshtein Distance: 0
Species: ssc-let-7e - cpo-let-7e-3p | Levenshtein Distance: 15
Species: ssc-let-7e - ppa-let-7e | Levenshtein Distance: 0
Species: ssc-let-7e - cja-let-7e | Levenshtein Distance: 0
Species: ssc-let-7e - nle-let-7e | Levenshtein Distance: 0
Species: ssc-let-7e - sbo-let-7e | Levenshtein Distance: 0
Species: ssc-let-7e - pha-let-7e | Levenshtein Distance: 0
Species: ssc-let-7e - oga-let-7e | Levenshtein Distance: 0
Species: tgu-let-7e-5p - tgu-let-7e-3p | Levenshtein Distance: 15
Species: tgu-let-7e-5p - ppy-let-7e | Levenshtein Distance: 3
Species: tgu-let-7e-5p - aca-let-7e-5p | Levenshtein Distance: 0
Species: tgu-let-7e-5p - aca-let-7e-3p | Levenshtein Distance: 14
Species: tgu-let-7e-5p - ola-let-7e | Levenshtein Distance: 4
Species: tgu-let-7e-5p - ggo-let-7e | Levenshtein Distance: 4
Species: tgu-let-7e-5p - ipu-let-7e | Levenshtein Distance: 0
Species: tgu-let-7e-5p - ssa-let-7e-5p | Levenshtein Distance: 0
Species: tgu-let-7e-5p - ssa-let-7e-3p | Levenshtein Distance: 15
Species: tgu-let-7e-5p - efu-let-7e | Levenshtein Distance: 5
Species: tgu-let-7e-5p - cpi-let-7e-5p | Levenshtein Distance: 0
Species: tgu-let-7e-5p - cpi-let-7e-3p | Levenshtein Distance: 13
Species: tgu-let-7e-5p - ami-let-7e-5p | Levenshtein Distance: 0
Species: tgu-let-7e-5p - ami-let-7e-3p | Levenshtein Distance: 13
Species: tgu-let-7e-5p - cli-let-7e-5p | Levenshtein Distance: 0
Species: tgu-let-7e-5p - cli-let-7e-3p | Levenshtein Distance: 15
Species: tgu-let-7e-5p - pbv-let-7e-5p | Levenshtein Distance: 0
Species: tgu-let-7e-5p - pbv-let-7e-3p | Levenshtein Distance: 15
Species: tgu-let-7e-5p - chi-let-7e-5p | Levenshtein Distance: 3
Species: tgu-let-7e-5p - chi-let-7e-3p | Levenshtein Distance: 16
Species: tgu-let-7e-5p - tch-let-7e-5p | Levenshtein Distance: 3
Species: tgu-let-7e-5p - oha-let-7e-5p | Levenshtein Distance: 0
Species: tgu-let-7e-5p - oha-let-7e-3p | Levenshtein Distance: 14
Species: tgu-let-7e-5p - pal-let-7e-5p | Levenshtein Distance: 3
Species: tgu-let-7e-5p - pal-let-7e-3p | Levenshtein Distance: 15
Species: tgu-let-7e-5p - cgr-let-7e | Levenshtein Distance: 16
Species: tgu-let-7e-5p - abu-let-7e | Levenshtein Distance: 0
Species: tgu-let-7e-5p - mze-let-7e | Levenshtein Distance: 0
Species: tgu-let-7e-5p - nbr-let-7e | Levenshtein Distance: 0
Species: tgu-let-7e-5p - oni-let-7e | Levenshtein Distance: 0
Species: tgu-let-7e-5p - pny-let-7e | Levenshtein Distance: 0
Species: tgu-let-7e-5p - gmo-let-7e-5p | Levenshtein Distance: 0
Species: tgu-let-7e-5p - gmo-let-7e-3p | Levenshtein Distance: 15
Species: tgu-let-7e-5p - xla-let-7e-5p | Levenshtein Distance: 3
Species: tgu-let-7e-5p - xla-let-7e-3p | Levenshtein Distance: 15
Species: tgu-let-7e-5p - cpo-let-7e-5p | Levenshtein Distance: 3
Species: tgu-let-7e-5p - cpo-let-7e-3p | Levenshtein Distance: 16
Species: tgu-let-7e-5p - ppa-let-7e | Levenshtein Distance: 3
Species: tgu-let-7e-5p - cja-let-7e | Levenshtein Distance: 3
Species: tgu-let-7e-5p - nle-let-7e | Levenshtein Distance: 3
Species: tgu-let-7e-5p - sbo-let-7e | Levenshtein Distance: 3
Species: tgu-let-7e-5p - pha-let-7e | Levenshtein Distance: 3
Species: tgu-let-7e-5p - oga-let-7e | Levenshtein Distance: 3
Species: tgu-let-7e-3p - ppy-let-7e | Levenshtein Distance: 15
Species: tgu-let-7e-3p - aca-let-7e-5p | Levenshtein Distance: 15
Species: tgu-let-7e-3p - aca-let-7e-3p | Levenshtein Distance: 4
Species: tgu-let-7e-3p - ola-let-7e | Levenshtein Distance: 15
Species: tgu-let-7e-3p - ggo-let-7e | Levenshtein Distance: 15
Species: tgu-let-7e-3p - ipu-let-7e | Levenshtein Distance: 15
Species: tgu-let-7e-3p - ssa-let-7e-5p | Levenshtein Distance: 15
Species: tgu-let-7e-3p - ssa-let-7e-3p | Levenshtein Distance: 0
Species: tgu-let-7e-3p - efu-let-7e | Levenshtein Distance: 15
Species: tgu-let-7e-3p - cpi-let-7e-5p | Levenshtein Distance: 15
Species: tgu-let-7e-3p - cpi-let-7e-3p | Levenshtein Distance: 2
Species: tgu-let-7e-3p - ami-let-7e-5p | Levenshtein Distance: 15
Species: tgu-let-7e-3p - ami-let-7e-3p | Levenshtein Distance: 2
Species: tgu-let-7e-3p - cli-let-7e-5p | Levenshtein Distance: 15
Species: tgu-let-7e-3p - cli-let-7e-3p | Levenshtein Distance: 0
Species: tgu-let-7e-3p - pbv-let-7e-5p | Levenshtein Distance: 15
Species: tgu-let-7e-3p - pbv-let-7e-3p | Levenshtein Distance: 0
Species: tgu-let-7e-3p - chi-let-7e-5p | Levenshtein Distance: 15
Species: tgu-let-7e-3p - chi-let-7e-3p | Levenshtein Distance: 6
Species: tgu-let-7e-3p - tch-let-7e-5p | Levenshtein Distance: 15
Species: tgu-let-7e-3p - oha-let-7e-5p | Levenshtein Distance: 15
Species: tgu-let-7e-3p - oha-let-7e-3p | Levenshtein Distance: 4
Species: tgu-let-7e-3p - pal-let-7e-5p | Levenshtein Distance: 15
Species: tgu-let-7e-3p - pal-let-7e-3p | Levenshtein Distance: 7
Species: tgu-let-7e-3p - cgr-let-7e | Levenshtein Distance: 6
Species: tgu-let-7e-3p - abu-let-7e | Levenshtein Distance: 15
Species: tgu-let-7e-3p - mze-let-7e | Levenshtein Distance: 15
Species: tgu-let-7e-3p - nbr-let-7e | Levenshtein Distance: 15
Species: tgu-let-7e-3p - oni-let-7e | Levenshtein Distance: 15
Species: tgu-let-7e-3p - pny-let-7e | Levenshtein Distance: 15
Species: tgu-let-7e-3p - gmo-let-7e-5p | Levenshtein Distance: 15
Species: tgu-let-7e-3p - gmo-let-7e-3p | Levenshtein Distance: 0
Species: tgu-let-7e-3p - xla-let-7e-5p | Levenshtein Distance: 15
Species: tgu-let-7e-3p - xla-let-7e-3p | Levenshtein Distance: 0
Species: tgu-let-7e-3p - cpo-let-7e-5p | Levenshtein Distance: 15
Species: tgu-let-7e-3p - cpo-let-7e-3p | Levenshtein Distance: 6
Species: tgu-let-7e-3p - ppa-let-7e | Levenshtein Distance: 15
Species: tgu-let-7e-3p - cja-let-7e | Levenshtein Distance: 15
Species: tgu-let-7e-3p - nle-let-7e | Levenshtein Distance: 15
Species: tgu-let-7e-3p - sbo-let-7e | Levenshtein Distance: 15
Species: tgu-let-7e-3p - pha-let-7e | Levenshtein Distance: 15
Species: tgu-let-7e-3p - oga-let-7e | Levenshtein Distance: 15
Species: ppy-let-7e - aca-let-7e-5p | Levenshtein Distance: 3
Species: ppy-let-7e - aca-let-7e-3p | Levenshtein Distance: 14
Species: ppy-let-7e - ola-let-7e | Levenshtein Distance: 6
Species: ppy-let-7e - ggo-let-7e | Levenshtein Distance: 1
Species: ppy-let-7e - ipu-let-7e | Levenshtein Distance: 3
Species: ppy-let-7e - ssa-let-7e-5p | Levenshtein Distance: 3
Species: ppy-let-7e - ssa-let-7e-3p | Levenshtein Distance: 15
Species: ppy-let-7e - efu-let-7e | Levenshtein Distance: 2
Species: ppy-let-7e - cpi-let-7e-5p | Levenshtein Distance: 3
Species: ppy-let-7e - cpi-let-7e-3p | Levenshtein Distance: 13
Species: ppy-let-7e - ami-let-7e-5p | Levenshtein Distance: 3
Species: ppy-let-7e - ami-let-7e-3p | Levenshtein Distance: 13
Species: ppy-let-7e - cli-let-7e-5p | Levenshtein Distance: 3
Species: ppy-let-7e - cli-let-7e-3p | Levenshtein Distance: 15
Species: ppy-let-7e - pbv-let-7e-5p | Levenshtein Distance: 3
Species: ppy-let-7e - pbv-let-7e-3p | Levenshtein Distance: 15
Species: ppy-let-7e - chi-let-7e-5p | Levenshtein Distance: 0
Species: ppy-let-7e - chi-let-7e-3p | Levenshtein Distance: 15
Species: ppy-let-7e - tch-let-7e-5p | Levenshtein Distance: 0
Species: ppy-let-7e - oha-let-7e-5p | Levenshtein Distance: 3
Species: ppy-let-7e - oha-let-7e-3p | Levenshtein Distance: 14
Species: ppy-let-7e - pal-let-7e-5p | Levenshtein Distance: 0
Species: ppy-let-7e - pal-let-7e-3p | Levenshtein Distance: 14
Species: ppy-let-7e - cgr-let-7e | Levenshtein Distance: 15
Species: ppy-let-7e - abu-let-7e | Levenshtein Distance: 3
Species: ppy-let-7e - mze-let-7e | Levenshtein Distance: 3
Species: ppy-let-7e - nbr-let-7e | Levenshtein Distance: 3
Species: ppy-let-7e - oni-let-7e | Levenshtein Distance: 3
Species: ppy-let-7e - pny-let-7e | Levenshtein Distance: 3
Species: ppy-let-7e - gmo-let-7e-5p | Levenshtein Distance: 3
Species: ppy-let-7e - gmo-let-7e-3p | Levenshtein Distance: 15
Species: ppy-let-7e - xla-let-7e-5p | Levenshtein Distance: 2
Species: ppy-let-7e - xla-let-7e-3p | Levenshtein Distance: 15
Species: ppy-let-7e - cpo-let-7e-5p | Levenshtein Distance: 0
Species: ppy-let-7e - cpo-let-7e-3p | Levenshtein Distance: 15
Species: ppy-let-7e - ppa-let-7e | Levenshtein Distance: 0
Species: ppy-let-7e - cja-let-7e | Levenshtein Distance: 0
Species: ppy-let-7e - nle-let-7e | Levenshtein Distance: 0
Species: ppy-let-7e - sbo-let-7e | Levenshtein Distance: 0
Species: ppy-let-7e - pha-let-7e | Levenshtein Distance: 0
Species: ppy-let-7e - oga-let-7e | Levenshtein Distance: 0
Species: aca-let-7e-5p - aca-let-7e-3p | Levenshtein Distance: 14
Species: aca-let-7e-5p - ola-let-7e | Levenshtein Distance: 4
Species: aca-let-7e-5p - ggo-let-7e | Levenshtein Distance: 4
Species: aca-let-7e-5p - ipu-let-7e | Levenshtein Distance: 0
Species: aca-let-7e-5p - ssa-let-7e-5p | Levenshtein Distance: 0
Species: aca-let-7e-5p - ssa-let-7e-3p | Levenshtein Distance: 15
Species: aca-let-7e-5p - efu-let-7e | Levenshtein Distance: 5
Species: aca-let-7e-5p - cpi-let-7e-5p | Levenshtein Distance: 0
Species: aca-let-7e-5p - cpi-let-7e-3p | Levenshtein Distance: 13
Species: aca-let-7e-5p - ami-let-7e-5p | Levenshtein Distance: 0
Species: aca-let-7e-5p - ami-let-7e-3p | Levenshtein Distance: 13
Species: aca-let-7e-5p - cli-let-7e-5p | Levenshtein Distance: 0
Species: aca-let-7e-5p - cli-let-7e-3p | Levenshtein Distance: 15
Species: aca-let-7e-5p - pbv-let-7e-5p | Levenshtein Distance: 0
Species: aca-let-7e-5p - pbv-let-7e-3p | Levenshtein Distance: 15
Species: aca-let-7e-5p - chi-let-7e-5p | Levenshtein Distance: 3
Species: aca-let-7e-5p - chi-let-7e-3p | Levenshtein Distance: 16
Species: aca-let-7e-5p - tch-let-7e-5p | Levenshtein Distance: 3
Species: aca-let-7e-5p - oha-let-7e-5p | Levenshtein Distance: 0
Species: aca-let-7e-5p - oha-let-7e-3p | Levenshtein Distance: 14
Species: aca-let-7e-5p - pal-let-7e-5p | Levenshtein Distance: 3
Species: aca-let-7e-5p - pal-let-7e-3p | Levenshtein Distance: 15
Species: aca-let-7e-5p - cgr-let-7e | Levenshtein Distance: 16
Species: aca-let-7e-5p - abu-let-7e | Levenshtein Distance: 0
Species: aca-let-7e-5p - mze-let-7e | Levenshtein Distance: 0
Species: aca-let-7e-5p - nbr-let-7e | Levenshtein Distance: 0
Species: aca-let-7e-5p - oni-let-7e | Levenshtein Distance: 0
Species: aca-let-7e-5p - pny-let-7e | Levenshtein Distance: 0
Species: aca-let-7e-5p - gmo-let-7e-5p | Levenshtein Distance: 0
Species: aca-let-7e-5p - gmo-let-7e-3p | Levenshtein Distance: 15
Species: aca-let-7e-5p - xla-let-7e-5p | Levenshtein Distance: 3
Species: aca-let-7e-5p - xla-let-7e-3p | Levenshtein Distance: 15
Species: aca-let-7e-5p - cpo-let-7e-5p | Levenshtein Distance: 3
Species: aca-let-7e-5p - cpo-let-7e-3p | Levenshtein Distance: 16
Species: aca-let-7e-5p - ppa-let-7e | Levenshtein Distance: 3
Species: aca-let-7e-5p - cja-let-7e | Levenshtein Distance: 3
Species: aca-let-7e-5p - nle-let-7e | Levenshtein Distance: 3
Species: aca-let-7e-5p - sbo-let-7e | Levenshtein Distance: 3
Species: aca-let-7e-5p - pha-let-7e | Levenshtein Distance: 3
Species: aca-let-7e-5p - oga-let-7e | Levenshtein Distance: 3
Species: aca-let-7e-3p - ola-let-7e | Levenshtein Distance: 15
Species: aca-let-7e-3p - ggo-let-7e | Levenshtein Distance: 14
Species: aca-let-7e-3p - ipu-let-7e | Levenshtein Distance: 14
Species: aca-let-7e-3p - ssa-let-7e-5p | Levenshtein Distance: 14
Species: aca-let-7e-3p - ssa-let-7e-3p | Levenshtein Distance: 4
Species: aca-let-7e-3p - efu-let-7e | Levenshtein Distance: 15
Species: aca-let-7e-3p - cpi-let-7e-5p | Levenshtein Distance: 14
Species: aca-let-7e-3p - cpi-let-7e-3p | Levenshtein Distance: 4
Species: aca-let-7e-3p - ami-let-7e-5p | Levenshtein Distance: 14
Species: aca-let-7e-3p - ami-let-7e-3p | Levenshtein Distance: 4
Species: aca-let-7e-3p - cli-let-7e-5p | Levenshtein Distance: 14
Species: aca-let-7e-3p - cli-let-7e-3p | Levenshtein Distance: 4
Species: aca-let-7e-3p - pbv-let-7e-5p | Levenshtein Distance: 14
Species: aca-let-7e-3p - pbv-let-7e-3p | Levenshtein Distance: 4
Species: aca-let-7e-3p - chi-let-7e-5p | Levenshtein Distance: 14
Species: aca-let-7e-3p - chi-let-7e-3p | Levenshtein Distance: 10
Species: aca-let-7e-3p - tch-let-7e-5p | Levenshtein Distance: 14
Species: aca-let-7e-3p - oha-let-7e-5p | Levenshtein Distance: 14
Species: aca-let-7e-3p - oha-let-7e-3p | Levenshtein Distance: 0
Species: aca-let-7e-3p - pal-let-7e-5p | Levenshtein Distance: 14
Species: aca-let-7e-3p - pal-let-7e-3p | Levenshtein Distance: 9
Species: aca-let-7e-3p - cgr-let-7e | Levenshtein Distance: 10
Species: aca-let-7e-3p - abu-let-7e | Levenshtein Distance: 14
Species: aca-let-7e-3p - mze-let-7e | Levenshtein Distance: 14
Species: aca-let-7e-3p - nbr-let-7e | Levenshtein Distance: 14
Species: aca-let-7e-3p - oni-let-7e | Levenshtein Distance: 14
Species: aca-let-7e-3p - pny-let-7e | Levenshtein Distance: 14
Species: aca-let-7e-3p - gmo-let-7e-5p | Levenshtein Distance: 14
Species: aca-let-7e-3p - gmo-let-7e-3p | Levenshtein Distance: 4
Species: aca-let-7e-3p - xla-let-7e-5p | Levenshtein Distance: 14
Species: aca-let-7e-3p - xla-let-7e-3p | Levenshtein Distance: 4
Species: aca-let-7e-3p - cpo-let-7e-5p | Levenshtein Distance: 14
Species: aca-let-7e-3p - cpo-let-7e-3p | Levenshtein Distance: 10
Species: aca-let-7e-3p - ppa-let-7e | Levenshtein Distance: 14
Species: aca-let-7e-3p - cja-let-7e | Levenshtein Distance: 14
Species: aca-let-7e-3p - nle-let-7e | Levenshtein Distance: 14
Species: aca-let-7e-3p - sbo-let-7e | Levenshtein Distance: 14
Species: aca-let-7e-3p - pha-let-7e | Levenshtein Distance: 14
Species: aca-let-7e-3p - oga-let-7e | Levenshtein Distance: 14
Species: ola-let-7e - ggo-let-7e | Levenshtein Distance: 5
Species: ola-let-7e - ipu-let-7e | Levenshtein Distance: 4
Species: ola-let-7e - ssa-let-7e-5p | Levenshtein Distance: 4
Species: ola-let-7e - ssa-let-7e-3p | Levenshtein Distance: 15
Species: ola-let-7e - efu-let-7e | Levenshtein Distance: 8
Species: ola-let-7e - cpi-let-7e-5p | Levenshtein Distance: 4
Species: ola-let-7e - cpi-let-7e-3p | Levenshtein Distance: 13
Species: ola-let-7e - ami-let-7e-5p | Levenshtein Distance: 4
Species: ola-let-7e - ami-let-7e-3p | Levenshtein Distance: 13
Species: ola-let-7e - cli-let-7e-5p | Levenshtein Distance: 4
Species: ola-let-7e - cli-let-7e-3p | Levenshtein Distance: 15
Species: ola-let-7e - pbv-let-7e-5p | Levenshtein Distance: 4
Species: ola-let-7e - pbv-let-7e-3p | Levenshtein Distance: 15
Species: ola-let-7e - chi-let-7e-5p | Levenshtein Distance: 6
Species: ola-let-7e - chi-let-7e-3p | Levenshtein Distance: 13
Species: ola-let-7e - tch-let-7e-5p | Levenshtein Distance: 6
Species: ola-let-7e - oha-let-7e-5p | Levenshtein Distance: 4
Species: ola-let-7e - oha-let-7e-3p | Levenshtein Distance: 15
Species: ola-let-7e - pal-let-7e-5p | Levenshtein Distance: 6
Species: ola-let-7e - pal-let-7e-3p | Levenshtein Distance: 13
Species: ola-let-7e - cgr-let-7e | Levenshtein Distance: 13
Species: ola-let-7e - abu-let-7e | Levenshtein Distance: 4
Species: ola-let-7e - mze-let-7e | Levenshtein Distance: 4
Species: ola-let-7e - nbr-let-7e | Levenshtein Distance: 4
Species: ola-let-7e - oni-let-7e | Levenshtein Distance: 4
Species: ola-let-7e - pny-let-7e | Levenshtein Distance: 4
Species: ola-let-7e - gmo-let-7e-5p | Levenshtein Distance: 4
Species: ola-let-7e - gmo-let-7e-3p | Levenshtein Distance: 15
Species: ola-let-7e - xla-let-7e-5p | Levenshtein Distance: 6
Species: ola-let-7e - xla-let-7e-3p | Levenshtein Distance: 15
Species: ola-let-7e - cpo-let-7e-5p | Levenshtein Distance: 6
Species: ola-let-7e - cpo-let-7e-3p | Levenshtein Distance: 13
Species: ola-let-7e - ppa-let-7e | Levenshtein Distance: 6
Species: ola-let-7e - cja-let-7e | Levenshtein Distance: 6
Species: ola-let-7e - nle-let-7e | Levenshtein Distance: 6
Species: ola-let-7e - sbo-let-7e | Levenshtein Distance: 6
Species: ola-let-7e - pha-let-7e | Levenshtein Distance: 6
Species: ola-let-7e - oga-let-7e | Levenshtein Distance: 6
Species: ggo-let-7e - ipu-let-7e | Levenshtein Distance: 4
Species: ggo-let-7e - ssa-let-7e-5p | Levenshtein Distance: 4
Species: ggo-let-7e - ssa-let-7e-3p | Levenshtein Distance: 15
Species: ggo-let-7e - efu-let-7e | Levenshtein Distance: 3
Species: ggo-let-7e - cpi-let-7e-5p | Levenshtein Distance: 4
Species: ggo-let-7e - cpi-let-7e-3p | Levenshtein Distance: 13
Species: ggo-let-7e - ami-let-7e-5p | Levenshtein Distance: 4
Species: ggo-let-7e - ami-let-7e-3p | Levenshtein Distance: 13
Species: ggo-let-7e - cli-let-7e-5p | Levenshtein Distance: 4
Species: ggo-let-7e - cli-let-7e-3p | Levenshtein Distance: 15
Species: ggo-let-7e - pbv-let-7e-5p | Levenshtein Distance: 4
Species: ggo-let-7e - pbv-let-7e-3p | Levenshtein Distance: 15
Species: ggo-let-7e - chi-let-7e-5p | Levenshtein Distance: 1
Species: ggo-let-7e - chi-let-7e-3p | Levenshtein Distance: 16
Species: ggo-let-7e - tch-let-7e-5p | Levenshtein Distance: 1
Species: ggo-let-7e - oha-let-7e-5p | Levenshtein Distance: 4
Species: ggo-let-7e - oha-let-7e-3p | Levenshtein Distance: 14
Species: ggo-let-7e - pal-let-7e-5p | Levenshtein Distance: 1
Species: ggo-let-7e - pal-let-7e-3p | Levenshtein Distance: 15
Species: ggo-let-7e - cgr-let-7e | Levenshtein Distance: 16
Species: ggo-let-7e - abu-let-7e | Levenshtein Distance: 4
Species: ggo-let-7e - mze-let-7e | Levenshtein Distance: 4
Species: ggo-let-7e - nbr-let-7e | Levenshtein Distance: 4
Species: ggo-let-7e - oni-let-7e | Levenshtein Distance: 4
Species: ggo-let-7e - pny-let-7e | Levenshtein Distance: 4
Species: ggo-let-7e - gmo-let-7e-5p | Levenshtein Distance: 4
Species: ggo-let-7e - gmo-let-7e-3p | Levenshtein Distance: 15
Species: ggo-let-7e - xla-let-7e-5p | Levenshtein Distance: 3
Species: ggo-let-7e - xla-let-7e-3p | Levenshtein Distance: 15
Species: ggo-let-7e - cpo-let-7e-5p | Levenshtein Distance: 1
Species: ggo-let-7e - cpo-let-7e-3p | Levenshtein Distance: 16
Species: ggo-let-7e - ppa-let-7e | Levenshtein Distance: 1
Species: ggo-let-7e - cja-let-7e | Levenshtein Distance: 1
Species: ggo-let-7e - nle-let-7e | Levenshtein Distance: 1
Species: ggo-let-7e - sbo-let-7e | Levenshtein Distance: 1
Species: ggo-let-7e - pha-let-7e | Levenshtein Distance: 1
Species: ggo-let-7e - oga-let-7e | Levenshtein Distance: 1
Species: ipu-let-7e - ssa-let-7e-5p | Levenshtein Distance: 0
Species: ipu-let-7e - ssa-let-7e-3p | Levenshtein Distance: 15
Species: ipu-let-7e - efu-let-7e | Levenshtein Distance: 5
Species: ipu-let-7e - cpi-let-7e-5p | Levenshtein Distance: 0
Species: ipu-let-7e - cpi-let-7e-3p | Levenshtein Distance: 13
Species: ipu-let-7e - ami-let-7e-5p | Levenshtein Distance: 0
Species: ipu-let-7e - ami-let-7e-3p | Levenshtein Distance: 13
Species: ipu-let-7e - cli-let-7e-5p | Levenshtein Distance: 0
Species: ipu-let-7e - cli-let-7e-3p | Levenshtein Distance: 15
Species: ipu-let-7e - pbv-let-7e-5p | Levenshtein Distance: 0
Species: ipu-let-7e - pbv-let-7e-3p | Levenshtein Distance: 15
Species: ipu-let-7e - chi-let-7e-5p | Levenshtein Distance: 3
Species: ipu-let-7e - chi-let-7e-3p | Levenshtein Distance: 16
Species: ipu-let-7e - tch-let-7e-5p | Levenshtein Distance: 3
Species: ipu-let-7e - oha-let-7e-5p | Levenshtein Distance: 0
Species: ipu-let-7e - oha-let-7e-3p | Levenshtein Distance: 14
Species: ipu-let-7e - pal-let-7e-5p | Levenshtein Distance: 3
Species: ipu-let-7e - pal-let-7e-3p | Levenshtein Distance: 15
Species: ipu-let-7e - cgr-let-7e | Levenshtein Distance: 16
Species: ipu-let-7e - abu-let-7e | Levenshtein Distance: 0
Species: ipu-let-7e - mze-let-7e | Levenshtein Distance: 0
Species: ipu-let-7e - nbr-let-7e | Levenshtein Distance: 0
Species: ipu-let-7e - oni-let-7e | Levenshtein Distance: 0
Species: ipu-let-7e - pny-let-7e | Levenshtein Distance: 0
Species: ipu-let-7e - gmo-let-7e-5p | Levenshtein Distance: 0
Species: ipu-let-7e - gmo-let-7e-3p | Levenshtein Distance: 15
Species: ipu-let-7e - xla-let-7e-5p | Levenshtein Distance: 3
Species: ipu-let-7e - xla-let-7e-3p | Levenshtein Distance: 15
Species: ipu-let-7e - cpo-let-7e-5p | Levenshtein Distance: 3
Species: ipu-let-7e - cpo-let-7e-3p | Levenshtein Distance: 16
Species: ipu-let-7e - ppa-let-7e | Levenshtein Distance: 3
Species: ipu-let-7e - cja-let-7e | Levenshtein Distance: 3
Species: ipu-let-7e - nle-let-7e | Levenshtein Distance: 3
Species: ipu-let-7e - sbo-let-7e | Levenshtein Distance: 3
Species: ipu-let-7e - pha-let-7e | Levenshtein Distance: 3
Species: ipu-let-7e - oga-let-7e | Levenshtein Distance: 3
Species: ssa-let-7e-5p - ssa-let-7e-3p | Levenshtein Distance: 15
Species: ssa-let-7e-5p - efu-let-7e | Levenshtein Distance: 5
Species: ssa-let-7e-5p - cpi-let-7e-5p | Levenshtein Distance: 0
Species: ssa-let-7e-5p - cpi-let-7e-3p | Levenshtein Distance: 13
Species: ssa-let-7e-5p - ami-let-7e-5p | Levenshtein Distance: 0
Species: ssa-let-7e-5p - ami-let-7e-3p | Levenshtein Distance: 13
Species: ssa-let-7e-5p - cli-let-7e-5p | Levenshtein Distance: 0
Species: ssa-let-7e-5p - cli-let-7e-3p | Levenshtein Distance: 15
Species: ssa-let-7e-5p - pbv-let-7e-5p | Levenshtein Distance: 0
Species: ssa-let-7e-5p - pbv-let-7e-3p | Levenshtein Distance: 15
Species: ssa-let-7e-5p - chi-let-7e-5p | Levenshtein Distance: 3
Species: ssa-let-7e-5p - chi-let-7e-3p | Levenshtein Distance: 16
Species: ssa-let-7e-5p - tch-let-7e-5p | Levenshtein Distance: 3
Species: ssa-let-7e-5p - oha-let-7e-5p | Levenshtein Distance: 0
Species: ssa-let-7e-5p - oha-let-7e-3p | Levenshtein Distance: 14
Species: ssa-let-7e-5p - pal-let-7e-5p | Levenshtein Distance: 3
Species: ssa-let-7e-5p - pal-let-7e-3p | Levenshtein Distance: 15
Species: ssa-let-7e-5p - cgr-let-7e | Levenshtein Distance: 16
Species: ssa-let-7e-5p - abu-let-7e | Levenshtein Distance: 0
Species: ssa-let-7e-5p - mze-let-7e | Levenshtein Distance: 0
Species: ssa-let-7e-5p - nbr-let-7e | Levenshtein Distance: 0
Species: ssa-let-7e-5p - oni-let-7e | Levenshtein Distance: 0
Species: ssa-let-7e-5p - pny-let-7e | Levenshtein Distance: 0
Species: ssa-let-7e-5p - gmo-let-7e-5p | Levenshtein Distance: 0
Species: ssa-let-7e-5p - gmo-let-7e-3p | Levenshtein Distance: 15
Species: ssa-let-7e-5p - xla-let-7e-5p | Levenshtein Distance: 3
Species: ssa-let-7e-5p - xla-let-7e-3p | Levenshtein Distance: 15
Species: ssa-let-7e-5p - cpo-let-7e-5p | Levenshtein Distance: 3
Species: ssa-let-7e-5p - cpo-let-7e-3p | Levenshtein Distance: 16
Species: ssa-let-7e-5p - ppa-let-7e | Levenshtein Distance: 3
Species: ssa-let-7e-5p - cja-let-7e | Levenshtein Distance: 3
Species: ssa-let-7e-5p - nle-let-7e | Levenshtein Distance: 3
Species: ssa-let-7e-5p - sbo-let-7e | Levenshtein Distance: 3
Species: ssa-let-7e-5p - pha-let-7e | Levenshtein Distance: 3
Species: ssa-let-7e-5p - oga-let-7e | Levenshtein Distance: 3
Species: ssa-let-7e-3p - efu-let-7e | Levenshtein Distance: 15
Species: ssa-let-7e-3p - cpi-let-7e-5p | Levenshtein Distance: 15
Species: ssa-let-7e-3p - cpi-let-7e-3p | Levenshtein Distance: 2
Species: ssa-let-7e-3p - ami-let-7e-5p | Levenshtein Distance: 15
Species: ssa-let-7e-3p - ami-let-7e-3p | Levenshtein Distance: 2
Species: ssa-let-7e-3p - cli-let-7e-5p | Levenshtein Distance: 15
Species: ssa-let-7e-3p - cli-let-7e-3p | Levenshtein Distance: 0
Species: ssa-let-7e-3p - pbv-let-7e-5p | Levenshtein Distance: 15
Species: ssa-let-7e-3p - pbv-let-7e-3p | Levenshtein Distance: 0
Species: ssa-let-7e-3p - chi-let-7e-5p | Levenshtein Distance: 15
Species: ssa-let-7e-3p - chi-let-7e-3p | Levenshtein Distance: 6
Species: ssa-let-7e-3p - tch-let-7e-5p | Levenshtein Distance: 15
Species: ssa-let-7e-3p - oha-let-7e-5p | Levenshtein Distance: 15
Species: ssa-let-7e-3p - oha-let-7e-3p | Levenshtein Distance: 4
Species: ssa-let-7e-3p - pal-let-7e-5p | Levenshtein Distance: 15
Species: ssa-let-7e-3p - pal-let-7e-3p | Levenshtein Distance: 7
Species: ssa-let-7e-3p - cgr-let-7e | Levenshtein Distance: 6
Species: ssa-let-7e-3p - abu-let-7e | Levenshtein Distance: 15
Species: ssa-let-7e-3p - mze-let-7e | Levenshtein Distance: 15
Species: ssa-let-7e-3p - nbr-let-7e | Levenshtein Distance: 15
Species: ssa-let-7e-3p - oni-let-7e | Levenshtein Distance: 15
Species: ssa-let-7e-3p - pny-let-7e | Levenshtein Distance: 15
Species: ssa-let-7e-3p - gmo-let-7e-5p | Levenshtein Distance: 15
Species: ssa-let-7e-3p - gmo-let-7e-3p | Levenshtein Distance: 0
Species: ssa-let-7e-3p - xla-let-7e-5p | Levenshtein Distance: 15
Species: ssa-let-7e-3p - xla-let-7e-3p | Levenshtein Distance: 0
Species: ssa-let-7e-3p - cpo-let-7e-5p | Levenshtein Distance: 15
Species: ssa-let-7e-3p - cpo-let-7e-3p | Levenshtein Distance: 6
Species: ssa-let-7e-3p - ppa-let-7e | Levenshtein Distance: 15
Species: ssa-let-7e-3p - cja-let-7e | Levenshtein Distance: 15
Species: ssa-let-7e-3p - nle-let-7e | Levenshtein Distance: 15
Species: ssa-let-7e-3p - sbo-let-7e | Levenshtein Distance: 15
Species: ssa-let-7e-3p - pha-let-7e | Levenshtein Distance: 15
Species: ssa-let-7e-3p - oga-let-7e | Levenshtein Distance: 15
Species: efu-let-7e - cpi-let-7e-5p | Levenshtein Distance: 5
Species: efu-let-7e - cpi-let-7e-3p | Levenshtein Distance: 15
Species: efu-let-7e - ami-let-7e-5p | Levenshtein Distance: 5
Species: efu-let-7e - ami-let-7e-3p | Levenshtein Distance: 15
Species: efu-let-7e - cli-let-7e-5p | Levenshtein Distance: 5
Species: efu-let-7e - cli-let-7e-3p | Levenshtein Distance: 15
Species: efu-let-7e - pbv-let-7e-5p | Levenshtein Distance: 5
Species: efu-let-7e - pbv-let-7e-3p | Levenshtein Distance: 15
Species: efu-let-7e - chi-let-7e-5p | Levenshtein Distance: 2
Species: efu-let-7e - chi-let-7e-3p | Levenshtein Distance: 15
Species: efu-let-7e - tch-let-7e-5p | Levenshtein Distance: 2
Species: efu-let-7e - oha-let-7e-5p | Levenshtein Distance: 5
Species: efu-let-7e - oha-let-7e-3p | Levenshtein Distance: 15
Species: efu-let-7e - pal-let-7e-5p | Levenshtein Distance: 2
Species: efu-let-7e - pal-let-7e-3p | Levenshtein Distance: 14
Species: efu-let-7e - cgr-let-7e | Levenshtein Distance: 15
Species: efu-let-7e - abu-let-7e | Levenshtein Distance: 5
Species: efu-let-7e - mze-let-7e | Levenshtein Distance: 5
Species: efu-let-7e - nbr-let-7e | Levenshtein Distance: 5
Species: efu-let-7e - oni-let-7e | Levenshtein Distance: 5
Species: efu-let-7e - pny-let-7e | Levenshtein Distance: 5
Species: efu-let-7e - gmo-let-7e-5p | Levenshtein Distance: 5
Species: efu-let-7e - gmo-let-7e-3p | Levenshtein Distance: 15
Species: efu-let-7e - xla-let-7e-5p | Levenshtein Distance: 4
Species: efu-let-7e - xla-let-7e-3p | Levenshtein Distance: 15
Species: efu-let-7e - cpo-let-7e-5p | Levenshtein Distance: 2
Species: efu-let-7e - cpo-let-7e-3p | Levenshtein Distance: 15
Species: efu-let-7e - ppa-let-7e | Levenshtein Distance: 2
Species: efu-let-7e - cja-let-7e | Levenshtein Distance: 2
Species: efu-let-7e - nle-let-7e | Levenshtein Distance: 2
Species: efu-let-7e - sbo-let-7e | Levenshtein Distance: 2
Species: efu-let-7e - pha-let-7e | Levenshtein Distance: 2
Species: efu-let-7e - oga-let-7e | Levenshtein Distance: 2
Species: cpi-let-7e-5p - cpi-let-7e-3p | Levenshtein Distance: 13
Species: cpi-let-7e-5p - ami-let-7e-5p | Levenshtein Distance: 0
Species: cpi-let-7e-5p - ami-let-7e-3p | Levenshtein Distance: 13
Species: cpi-let-7e-5p - cli-let-7e-5p | Levenshtein Distance: 0
Species: cpi-let-7e-5p - cli-let-7e-3p | Levenshtein Distance: 15
Species: cpi-let-7e-5p - pbv-let-7e-5p | Levenshtein Distance: 0
Species: cpi-let-7e-5p - pbv-let-7e-3p | Levenshtein Distance: 15
Species: cpi-let-7e-5p - chi-let-7e-5p | Levenshtein Distance: 3
Species: cpi-let-7e-5p - chi-let-7e-3p | Levenshtein Distance: 16
Species: cpi-let-7e-5p - tch-let-7e-5p | Levenshtein Distance: 3
Species: cpi-let-7e-5p - oha-let-7e-5p | Levenshtein Distance: 0
Species: cpi-let-7e-5p - oha-let-7e-3p | Levenshtein Distance: 14
Species: cpi-let-7e-5p - pal-let-7e-5p | Levenshtein Distance: 3
Species: cpi-let-7e-5p - pal-let-7e-3p | Levenshtein Distance: 15
Species: cpi-let-7e-5p - cgr-let-7e | Levenshtein Distance: 16
Species: cpi-let-7e-5p - abu-let-7e | Levenshtein Distance: 0
Species: cpi-let-7e-5p - mze-let-7e | Levenshtein Distance: 0
Species: cpi-let-7e-5p - nbr-let-7e | Levenshtein Distance: 0
Species: cpi-let-7e-5p - oni-let-7e | Levenshtein Distance: 0
Species: cpi-let-7e-5p - pny-let-7e | Levenshtein Distance: 0
Species: cpi-let-7e-5p - gmo-let-7e-5p | Levenshtein Distance: 0
Species: cpi-let-7e-5p - gmo-let-7e-3p | Levenshtein Distance: 15
Species: cpi-let-7e-5p - xla-let-7e-5p | Levenshtein Distance: 3
Species: cpi-let-7e-5p - xla-let-7e-3p | Levenshtein Distance: 15
Species: cpi-let-7e-5p - cpo-let-7e-5p | Levenshtein Distance: 3
Species: cpi-let-7e-5p - cpo-let-7e-3p | Levenshtein Distance: 16
Species: cpi-let-7e-5p - ppa-let-7e | Levenshtein Distance: 3
Species: cpi-let-7e-5p - cja-let-7e | Levenshtein Distance: 3
Species: cpi-let-7e-5p - nle-let-7e | Levenshtein Distance: 3
Species: cpi-let-7e-5p - sbo-let-7e | Levenshtein Distance: 3
Species: cpi-let-7e-5p - pha-let-7e | Levenshtein Distance: 3
Species: cpi-let-7e-5p - oga-let-7e | Levenshtein Distance: 3
Species: cpi-let-7e-3p - ami-let-7e-5p | Levenshtein Distance: 13
Species: cpi-let-7e-3p - ami-let-7e-3p | Levenshtein Distance: 0
Species: cpi-let-7e-3p - cli-let-7e-5p | Levenshtein Distance: 13
Species: cpi-let-7e-3p - cli-let-7e-3p | Levenshtein Distance: 2
Species: cpi-let-7e-3p - pbv-let-7e-5p | Levenshtein Distance: 13
Species: cpi-let-7e-3p - pbv-let-7e-3p | Levenshtein Distance: 2
Species: cpi-let-7e-3p - chi-let-7e-5p | Levenshtein Distance: 13
Species: cpi-let-7e-3p - chi-let-7e-3p | Levenshtein Distance: 8
Species: cpi-let-7e-3p - tch-let-7e-5p | Levenshtein Distance: 13
Species: cpi-let-7e-3p - oha-let-7e-5p | Levenshtein Distance: 13
Species: cpi-let-7e-3p - oha-let-7e-3p | Levenshtein Distance: 4
Species: cpi-let-7e-3p - pal-let-7e-5p | Levenshtein Distance: 13
Species: cpi-let-7e-3p - pal-let-7e-3p | Levenshtein Distance: 7
Species: cpi-let-7e-3p - cgr-let-7e | Levenshtein Distance: 8
Species: cpi-let-7e-3p - abu-let-7e | Levenshtein Distance: 13
Species: cpi-let-7e-3p - mze-let-7e | Levenshtein Distance: 13
Species: cpi-let-7e-3p - nbr-let-7e | Levenshtein Distance: 13
Species: cpi-let-7e-3p - oni-let-7e | Levenshtein Distance: 13
Species: cpi-let-7e-3p - pny-let-7e | Levenshtein Distance: 13
Species: cpi-let-7e-3p - gmo-let-7e-5p | Levenshtein Distance: 13
Species: cpi-let-7e-3p - gmo-let-7e-3p | Levenshtein Distance: 2
Species: cpi-let-7e-3p - xla-let-7e-5p | Levenshtein Distance: 13
Species: cpi-let-7e-3p - xla-let-7e-3p | Levenshtein Distance: 2
Species: cpi-let-7e-3p - cpo-let-7e-5p | Levenshtein Distance: 13
Species: cpi-let-7e-3p - cpo-let-7e-3p | Levenshtein Distance: 8
Species: cpi-let-7e-3p - ppa-let-7e | Levenshtein Distance: 13
Species: cpi-let-7e-3p - cja-let-7e | Levenshtein Distance: 13
Species: cpi-let-7e-3p - nle-let-7e | Levenshtein Distance: 13
Species: cpi-let-7e-3p - sbo-let-7e | Levenshtein Distance: 13
Species: cpi-let-7e-3p - pha-let-7e | Levenshtein Distance: 13
Species: cpi-let-7e-3p - oga-let-7e | Levenshtein Distance: 13
Species: ami-let-7e-5p - ami-let-7e-3p | Levenshtein Distance: 13
Species: ami-let-7e-5p - cli-let-7e-5p | Levenshtein Distance: 0
Species: ami-let-7e-5p - cli-let-7e-3p | Levenshtein Distance: 15
Species: ami-let-7e-5p - pbv-let-7e-5p | Levenshtein Distance: 0
Species: ami-let-7e-5p - pbv-let-7e-3p | Levenshtein Distance: 15
Species: ami-let-7e-5p - chi-let-7e-5p | Levenshtein Distance: 3
Species: ami-let-7e-5p - chi-let-7e-3p | Levenshtein Distance: 16
Species: ami-let-7e-5p - tch-let-7e-5p | Levenshtein Distance: 3
Species: ami-let-7e-5p - oha-let-7e-5p | Levenshtein Distance: 0
Species: ami-let-7e-5p - oha-let-7e-3p | Levenshtein Distance: 14
Species: ami-let-7e-5p - pal-let-7e-5p | Levenshtein Distance: 3
Species: ami-let-7e-5p - pal-let-7e-3p | Levenshtein Distance: 15
Species: ami-let-7e-5p - cgr-let-7e | Levenshtein Distance: 16
Species: ami-let-7e-5p - abu-let-7e | Levenshtein Distance: 0
Species: ami-let-7e-5p - mze-let-7e | Levenshtein Distance: 0
Species: ami-let-7e-5p - nbr-let-7e | Levenshtein Distance: 0
Species: ami-let-7e-5p - oni-let-7e | Levenshtein Distance: 0
Species: ami-let-7e-5p - pny-let-7e | Levenshtein Distance: 0
Species: ami-let-7e-5p - gmo-let-7e-5p | Levenshtein Distance: 0
Species: ami-let-7e-5p - gmo-let-7e-3p | Levenshtein Distance: 15
Species: ami-let-7e-5p - xla-let-7e-5p | Levenshtein Distance: 3
Species: ami-let-7e-5p - xla-let-7e-3p | Levenshtein Distance: 15
Species: ami-let-7e-5p - cpo-let-7e-5p | Levenshtein Distance: 3
Species: ami-let-7e-5p - cpo-let-7e-3p | Levenshtein Distance: 16
Species: ami-let-7e-5p - ppa-let-7e | Levenshtein Distance: 3
Species: ami-let-7e-5p - cja-let-7e | Levenshtein Distance: 3
Species: ami-let-7e-5p - nle-let-7e | Levenshtein Distance: 3
Species: ami-let-7e-5p - sbo-let-7e | Levenshtein Distance: 3
Species: ami-let-7e-5p - pha-let-7e | Levenshtein Distance: 3
Species: ami-let-7e-5p - oga-let-7e | Levenshtein Distance: 3
Species: ami-let-7e-3p - cli-let-7e-5p | Levenshtein Distance: 13
Species: ami-let-7e-3p - cli-let-7e-3p | Levenshtein Distance: 2
Species: ami-let-7e-3p - pbv-let-7e-5p | Levenshtein Distance: 13
Species: ami-let-7e-3p - pbv-let-7e-3p | Levenshtein Distance: 2
Species: ami-let-7e-3p - chi-let-7e-5p | Levenshtein Distance: 13
Species: ami-let-7e-3p - chi-let-7e-3p | Levenshtein Distance: 8
Species: ami-let-7e-3p - tch-let-7e-5p | Levenshtein Distance: 13
Species: ami-let-7e-3p - oha-let-7e-5p | Levenshtein Distance: 13
Species: ami-let-7e-3p - oha-let-7e-3p | Levenshtein Distance: 4
Species: ami-let-7e-3p - pal-let-7e-5p | Levenshtein Distance: 13
Species: ami-let-7e-3p - pal-let-7e-3p | Levenshtein Distance: 7
Species: ami-let-7e-3p - cgr-let-7e | Levenshtein Distance: 8
Species: ami-let-7e-3p - abu-let-7e | Levenshtein Distance: 13
Species: ami-let-7e-3p - mze-let-7e | Levenshtein Distance: 13
Species: ami-let-7e-3p - nbr-let-7e | Levenshtein Distance: 13
Species: ami-let-7e-3p - oni-let-7e | Levenshtein Distance: 13
Species: ami-let-7e-3p - pny-let-7e | Levenshtein Distance: 13
Species: ami-let-7e-3p - gmo-let-7e-5p | Levenshtein Distance: 13
Species: ami-let-7e-3p - gmo-let-7e-3p | Levenshtein Distance: 2
Species: ami-let-7e-3p - xla-let-7e-5p | Levenshtein Distance: 13
Species: ami-let-7e-3p - xla-let-7e-3p | Levenshtein Distance: 2
Species: ami-let-7e-3p - cpo-let-7e-5p | Levenshtein Distance: 13
Species: ami-let-7e-3p - cpo-let-7e-3p | Levenshtein Distance: 8
Species: ami-let-7e-3p - ppa-let-7e | Levenshtein Distance: 13
Species: ami-let-7e-3p - cja-let-7e | Levenshtein Distance: 13
Species: ami-let-7e-3p - nle-let-7e | Levenshtein Distance: 13
Species: ami-let-7e-3p - sbo-let-7e | Levenshtein Distance: 13
Species: ami-let-7e-3p - pha-let-7e | Levenshtein Distance: 13
Species: ami-let-7e-3p - oga-let-7e | Levenshtein Distance: 13
Species: cli-let-7e-5p - cli-let-7e-3p | Levenshtein Distance: 15
Species: cli-let-7e-5p - pbv-let-7e-5p | Levenshtein Distance: 0
Species: cli-let-7e-5p - pbv-let-7e-3p | Levenshtein Distance: 15
Species: cli-let-7e-5p - chi-let-7e-5p | Levenshtein Distance: 3
Species: cli-let-7e-5p - chi-let-7e-3p | Levenshtein Distance: 16
Species: cli-let-7e-5p - tch-let-7e-5p | Levenshtein Distance: 3
Species: cli-let-7e-5p - oha-let-7e-5p | Levenshtein Distance: 0
Species: cli-let-7e-5p - oha-let-7e-3p | Levenshtein Distance: 14
Species: cli-let-7e-5p - pal-let-7e-5p | Levenshtein Distance: 3
Species: cli-let-7e-5p - pal-let-7e-3p | Levenshtein Distance: 15
Species: cli-let-7e-5p - cgr-let-7e | Levenshtein Distance: 16
Species: cli-let-7e-5p - abu-let-7e | Levenshtein Distance: 0
Species: cli-let-7e-5p - mze-let-7e | Levenshtein Distance: 0
Species: cli-let-7e-5p - nbr-let-7e | Levenshtein Distance: 0
Species: cli-let-7e-5p - oni-let-7e | Levenshtein Distance: 0
Species: cli-let-7e-5p - pny-let-7e | Levenshtein Distance: 0
Species: cli-let-7e-5p - gmo-let-7e-5p | Levenshtein Distance: 0
Species: cli-let-7e-5p - gmo-let-7e-3p | Levenshtein Distance: 15
Species: cli-let-7e-5p - xla-let-7e-5p | Levenshtein Distance: 3
Species: cli-let-7e-5p - xla-let-7e-3p | Levenshtein Distance: 15
Species: cli-let-7e-5p - cpo-let-7e-5p | Levenshtein Distance: 3
Species: cli-let-7e-5p - cpo-let-7e-3p | Levenshtein Distance: 16
Species: cli-let-7e-5p - ppa-let-7e | Levenshtein Distance: 3
Species: cli-let-7e-5p - cja-let-7e | Levenshtein Distance: 3
Species: cli-let-7e-5p - nle-let-7e | Levenshtein Distance: 3
Species: cli-let-7e-5p - sbo-let-7e | Levenshtein Distance: 3
Species: cli-let-7e-5p - pha-let-7e | Levenshtein Distance: 3
Species: cli-let-7e-5p - oga-let-7e | Levenshtein Distance: 3
Species: cli-let-7e-3p - pbv-let-7e-5p | Levenshtein Distance: 15
Species: cli-let-7e-3p - pbv-let-7e-3p | Levenshtein Distance: 0
Species: cli-let-7e-3p - chi-let-7e-5p | Levenshtein Distance: 15
Species: cli-let-7e-3p - chi-let-7e-3p | Levenshtein Distance: 6
Species: cli-let-7e-3p - tch-let-7e-5p | Levenshtein Distance: 15
Species: cli-let-7e-3p - oha-let-7e-5p | Levenshtein Distance: 15
Species: cli-let-7e-3p - oha-let-7e-3p | Levenshtein Distance: 4
Species: cli-let-7e-3p - pal-let-7e-5p | Levenshtein Distance: 15
Species: cli-let-7e-3p - pal-let-7e-3p | Levenshtein Distance: 7
Species: cli-let-7e-3p - cgr-let-7e | Levenshtein Distance: 6
Species: cli-let-7e-3p - abu-let-7e | Levenshtein Distance: 15
Species: cli-let-7e-3p - mze-let-7e | Levenshtein Distance: 15
Species: cli-let-7e-3p - nbr-let-7e | Levenshtein Distance: 15
Species: cli-let-7e-3p - oni-let-7e | Levenshtein Distance: 15
Species: cli-let-7e-3p - pny-let-7e | Levenshtein Distance: 15
Species: cli-let-7e-3p - gmo-let-7e-5p | Levenshtein Distance: 15
Species: cli-let-7e-3p - gmo-let-7e-3p | Levenshtein Distance: 0
Species: cli-let-7e-3p - xla-let-7e-5p | Levenshtein Distance: 15
Species: cli-let-7e-3p - xla-let-7e-3p | Levenshtein Distance: 0
Species: cli-let-7e-3p - cpo-let-7e-5p | Levenshtein Distance: 15
Species: cli-let-7e-3p - cpo-let-7e-3p | Levenshtein Distance: 6
Species: cli-let-7e-3p - ppa-let-7e | Levenshtein Distance: 15
Species: cli-let-7e-3p - cja-let-7e | Levenshtein Distance: 15
Species: cli-let-7e-3p - nle-let-7e | Levenshtein Distance: 15
Species: cli-let-7e-3p - sbo-let-7e | Levenshtein Distance: 15
Species: cli-let-7e-3p - pha-let-7e | Levenshtein Distance: 15
Species: cli-let-7e-3p - oga-let-7e | Levenshtein Distance: 15
Species: pbv-let-7e-5p - pbv-let-7e-3p | Levenshtein Distance: 15
Species: pbv-let-7e-5p - chi-let-7e-5p | Levenshtein Distance: 3
Species: pbv-let-7e-5p - chi-let-7e-3p | Levenshtein Distance: 16
Species: pbv-let-7e-5p - tch-let-7e-5p | Levenshtein Distance: 3
Species: pbv-let-7e-5p - oha-let-7e-5p | Levenshtein Distance: 0
Species: pbv-let-7e-5p - oha-let-7e-3p | Levenshtein Distance: 14
Species: pbv-let-7e-5p - pal-let-7e-5p | Levenshtein Distance: 3
Species: pbv-let-7e-5p - pal-let-7e-3p | Levenshtein Distance: 15
Species: pbv-let-7e-5p - cgr-let-7e | Levenshtein Distance: 16
Species: pbv-let-7e-5p - abu-let-7e | Levenshtein Distance: 0
Species: pbv-let-7e-5p - mze-let-7e | Levenshtein Distance: 0
Species: pbv-let-7e-5p - nbr-let-7e | Levenshtein Distance: 0
Species: pbv-let-7e-5p - oni-let-7e | Levenshtein Distance: 0
Species: pbv-let-7e-5p - pny-let-7e | Levenshtein Distance: 0
Species: pbv-let-7e-5p - gmo-let-7e-5p | Levenshtein Distance: 0
Species: pbv-let-7e-5p - gmo-let-7e-3p | Levenshtein Distance: 15
Species: pbv-let-7e-5p - xla-let-7e-5p | Levenshtein Distance: 3
Species: pbv-let-7e-5p - xla-let-7e-3p | Levenshtein Distance: 15
Species: pbv-let-7e-5p - cpo-let-7e-5p | Levenshtein Distance: 3
Species: pbv-let-7e-5p - cpo-let-7e-3p | Levenshtein Distance: 16
Species: pbv-let-7e-5p - ppa-let-7e | Levenshtein Distance: 3
Species: pbv-let-7e-5p - cja-let-7e | Levenshtein Distance: 3
Species: pbv-let-7e-5p - nle-let-7e | Levenshtein Distance: 3
Species: pbv-let-7e-5p - sbo-let-7e | Levenshtein Distance: 3
Species: pbv-let-7e-5p - pha-let-7e | Levenshtein Distance: 3
Species: pbv-let-7e-5p - oga-let-7e | Levenshtein Distance: 3
Species: pbv-let-7e-3p - chi-let-7e-5p | Levenshtein Distance: 15
Species: pbv-let-7e-3p - chi-let-7e-3p | Levenshtein Distance: 6
Species: pbv-let-7e-3p - tch-let-7e-5p | Levenshtein Distance: 15
Species: pbv-let-7e-3p - oha-let-7e-5p | Levenshtein Distance: 15
Species: pbv-let-7e-3p - oha-let-7e-3p | Levenshtein Distance: 4
Species: pbv-let-7e-3p - pal-let-7e-5p | Levenshtein Distance: 15
Species: pbv-let-7e-3p - pal-let-7e-3p | Levenshtein Distance: 7
Species: pbv-let-7e-3p - cgr-let-7e | Levenshtein Distance: 6
Species: pbv-let-7e-3p - abu-let-7e | Levenshtein Distance: 15
Species: pbv-let-7e-3p - mze-let-7e | Levenshtein Distance: 15
Species: pbv-let-7e-3p - nbr-let-7e | Levenshtein Distance: 15
Species: pbv-let-7e-3p - oni-let-7e | Levenshtein Distance: 15
Species: pbv-let-7e-3p - pny-let-7e | Levenshtein Distance: 15
Species: pbv-let-7e-3p - gmo-let-7e-5p | Levenshtein Distance: 15
Species: pbv-let-7e-3p - gmo-let-7e-3p | Levenshtein Distance: 0
Species: pbv-let-7e-3p - xla-let-7e-5p | Levenshtein Distance: 15
Species: pbv-let-7e-3p - xla-let-7e-3p | Levenshtein Distance: 0
Species: pbv-let-7e-3p - cpo-let-7e-5p | Levenshtein Distance: 15
Species: pbv-let-7e-3p - cpo-let-7e-3p | Levenshtein Distance: 6
Species: pbv-let-7e-3p - ppa-let-7e | Levenshtein Distance: 15
Species: pbv-let-7e-3p - cja-let-7e | Levenshtein Distance: 15
Species: pbv-let-7e-3p - nle-let-7e | Levenshtein Distance: 15
Species: pbv-let-7e-3p - sbo-let-7e | Levenshtein Distance: 15
Species: pbv-let-7e-3p - pha-let-7e | Levenshtein Distance: 15
Species: pbv-let-7e-3p - oga-let-7e | Levenshtein Distance: 15
Species: chi-let-7e-5p - chi-let-7e-3p | Levenshtein Distance: 15
Species: chi-let-7e-5p - tch-let-7e-5p | Levenshtein Distance: 0
Species: chi-let-7e-5p - oha-let-7e-5p | Levenshtein Distance: 3
Species: chi-let-7e-5p - oha-let-7e-3p | Levenshtein Distance: 14
Species: chi-let-7e-5p - pal-let-7e-5p | Levenshtein Distance: 0
Species: chi-let-7e-5p - pal-let-7e-3p | Levenshtein Distance: 14
Species: chi-let-7e-5p - cgr-let-7e | Levenshtein Distance: 15
Species: chi-let-7e-5p - abu-let-7e | Levenshtein Distance: 3
Species: chi-let-7e-5p - mze-let-7e | Levenshtein Distance: 3
Species: chi-let-7e-5p - nbr-let-7e | Levenshtein Distance: 3
Species: chi-let-7e-5p - oni-let-7e | Levenshtein Distance: 3
Species: chi-let-7e-5p - pny-let-7e | Levenshtein Distance: 3
Species: chi-let-7e-5p - gmo-let-7e-5p | Levenshtein Distance: 3
Species: chi-let-7e-5p - gmo-let-7e-3p | Levenshtein Distance: 15
Species: chi-let-7e-5p - xla-let-7e-5p | Levenshtein Distance: 2
Species: chi-let-7e-5p - xla-let-7e-3p | Levenshtein Distance: 15
Species: chi-let-7e-5p - cpo-let-7e-5p | Levenshtein Distance: 0
Species: chi-let-7e-5p - cpo-let-7e-3p | Levenshtein Distance: 15
Species: chi-let-7e-5p - ppa-let-7e | Levenshtein Distance: 0
Species: chi-let-7e-5p - cja-let-7e | Levenshtein Distance: 0
Species: chi-let-7e-5p - nle-let-7e | Levenshtein Distance: 0
Species: chi-let-7e-5p - sbo-let-7e | Levenshtein Distance: 0
Species: chi-let-7e-5p - pha-let-7e | Levenshtein Distance: 0
Species: chi-let-7e-5p - oga-let-7e | Levenshtein Distance: 0
Species: chi-let-7e-3p - tch-let-7e-5p | Levenshtein Distance: 15
Species: chi-let-7e-3p - oha-let-7e-5p | Levenshtein Distance: 16
Species: chi-let-7e-3p - oha-let-7e-3p | Levenshtein Distance: 10
Species: chi-let-7e-3p - pal-let-7e-5p | Levenshtein Distance: 15
Species: chi-let-7e-3p - pal-let-7e-3p | Levenshtein Distance: 1
Species: chi-let-7e-3p - cgr-let-7e | Levenshtein Distance: 0
Species: chi-let-7e-3p - abu-let-7e | Levenshtein Distance: 16
Species: chi-let-7e-3p - mze-let-7e | Levenshtein Distance: 16
Species: chi-let-7e-3p - nbr-let-7e | Levenshtein Distance: 16
Species: chi-let-7e-3p - oni-let-7e | Levenshtein Distance: 16
Species: chi-let-7e-3p - pny-let-7e | Levenshtein Distance: 16
Species: chi-let-7e-3p - gmo-let-7e-5p | Levenshtein Distance: 16
Species: chi-let-7e-3p - gmo-let-7e-3p | Levenshtein Distance: 6
Species: chi-let-7e-3p - xla-let-7e-5p | Levenshtein Distance: 16
Species: chi-let-7e-3p - xla-let-7e-3p | Levenshtein Distance: 6
Species: chi-let-7e-3p - cpo-let-7e-5p | Levenshtein Distance: 15
Species: chi-let-7e-3p - cpo-let-7e-3p | Levenshtein Distance: 0
Species: chi-let-7e-3p - ppa-let-7e | Levenshtein Distance: 15
Species: chi-let-7e-3p - cja-let-7e | Levenshtein Distance: 15
Species: chi-let-7e-3p - nle-let-7e | Levenshtein Distance: 15
Species: chi-let-7e-3p - sbo-let-7e | Levenshtein Distance: 15
Species: chi-let-7e-3p - pha-let-7e | Levenshtein Distance: 15
Species: chi-let-7e-3p - oga-let-7e | Levenshtein Distance: 15
Species: tch-let-7e-5p - oha-let-7e-5p | Levenshtein Distance: 3
Species: tch-let-7e-5p - oha-let-7e-3p | Levenshtein Distance: 14
Species: tch-let-7e-5p - pal-let-7e-5p | Levenshtein Distance: 0
Species: tch-let-7e-5p - pal-let-7e-3p | Levenshtein Distance: 14
Species: tch-let-7e-5p - cgr-let-7e | Levenshtein Distance: 15
Species: tch-let-7e-5p - abu-let-7e | Levenshtein Distance: 3
Species: tch-let-7e-5p - mze-let-7e | Levenshtein Distance: 3
Species: tch-let-7e-5p - nbr-let-7e | Levenshtein Distance: 3
Species: tch-let-7e-5p - oni-let-7e | Levenshtein Distance: 3
Species: tch-let-7e-5p - pny-let-7e | Levenshtein Distance: 3
Species: tch-let-7e-5p - gmo-let-7e-5p | Levenshtein Distance: 3
Species: tch-let-7e-5p - gmo-let-7e-3p | Levenshtein Distance: 15
Species: tch-let-7e-5p - xla-let-7e-5p | Levenshtein Distance: 2
Species: tch-let-7e-5p - xla-let-7e-3p | Levenshtein Distance: 15
Species: tch-let-7e-5p - cpo-let-7e-5p | Levenshtein Distance: 0
Species: tch-let-7e-5p - cpo-let-7e-3p | Levenshtein Distance: 15
Species: tch-let-7e-5p - ppa-let-7e | Levenshtein Distance: 0
Species: tch-let-7e-5p - cja-let-7e | Levenshtein Distance: 0
Species: tch-let-7e-5p - nle-let-7e | Levenshtein Distance: 0
Species: tch-let-7e-5p - sbo-let-7e | Levenshtein Distance: 0
Species: tch-let-7e-5p - pha-let-7e | Levenshtein Distance: 0
Species: tch-let-7e-5p - oga-let-7e | Levenshtein Distance: 0
Species: oha-let-7e-5p - oha-let-7e-3p | Levenshtein Distance: 14
Species: oha-let-7e-5p - pal-let-7e-5p | Levenshtein Distance: 3
Species: oha-let-7e-5p - pal-let-7e-3p | Levenshtein Distance: 15
Species: oha-let-7e-5p - cgr-let-7e | Levenshtein Distance: 16
Species: oha-let-7e-5p - abu-let-7e | Levenshtein Distance: 0
Species: oha-let-7e-5p - mze-let-7e | Levenshtein Distance: 0
Species: oha-let-7e-5p - nbr-let-7e | Levenshtein Distance: 0
Species: oha-let-7e-5p - oni-let-7e | Levenshtein Distance: 0
Species: oha-let-7e-5p - pny-let-7e | Levenshtein Distance: 0
Species: oha-let-7e-5p - gmo-let-7e-5p | Levenshtein Distance: 0
Species: oha-let-7e-5p - gmo-let-7e-3p | Levenshtein Distance: 15
Species: oha-let-7e-5p - xla-let-7e-5p | Levenshtein Distance: 3
Species: oha-let-7e-5p - xla-let-7e-3p | Levenshtein Distance: 15
Species: oha-let-7e-5p - cpo-let-7e-5p | Levenshtein Distance: 3
Species: oha-let-7e-5p - cpo-let-7e-3p | Levenshtein Distance: 16
Species: oha-let-7e-5p - ppa-let-7e | Levenshtein Distance: 3
Species: oha-let-7e-5p - cja-let-7e | Levenshtein Distance: 3
Species: oha-let-7e-5p - nle-let-7e | Levenshtein Distance: 3
Species: oha-let-7e-5p - sbo-let-7e | Levenshtein Distance: 3
Species: oha-let-7e-5p - pha-let-7e | Levenshtein Distance: 3
Species: oha-let-7e-5p - oga-let-7e | Levenshtein Distance: 3
Species: oha-let-7e-3p - pal-let-7e-5p | Levenshtein Distance: 14
Species: oha-let-7e-3p - pal-let-7e-3p | Levenshtein Distance: 9
Species: oha-let-7e-3p - cgr-let-7e | Levenshtein Distance: 10
Species: oha-let-7e-3p - abu-let-7e | Levenshtein Distance: 14
Species: oha-let-7e-3p - mze-let-7e | Levenshtein Distance: 14
Species: oha-let-7e-3p - nbr-let-7e | Levenshtein Distance: 14
Species: oha-let-7e-3p - oni-let-7e | Levenshtein Distance: 14
Species: oha-let-7e-3p - pny-let-7e | Levenshtein Distance: 14
Species: oha-let-7e-3p - gmo-let-7e-5p | Levenshtein Distance: 14
Species: oha-let-7e-3p - gmo-let-7e-3p | Levenshtein Distance: 4
Species: oha-let-7e-3p - xla-let-7e-5p | Levenshtein Distance: 14
Species: oha-let-7e-3p - xla-let-7e-3p | Levenshtein Distance: 4
Species: oha-let-7e-3p - cpo-let-7e-5p | Levenshtein Distance: 14
Species: oha-let-7e-3p - cpo-let-7e-3p | Levenshtein Distance: 10
Species: oha-let-7e-3p - ppa-let-7e | Levenshtein Distance: 14
Species: oha-let-7e-3p - cja-let-7e | Levenshtein Distance: 14
Species: oha-let-7e-3p - nle-let-7e | Levenshtein Distance: 14
Species: oha-let-7e-3p - sbo-let-7e | Levenshtein Distance: 14
Species: oha-let-7e-3p - pha-let-7e | Levenshtein Distance: 14
Species: oha-let-7e-3p - oga-let-7e | Levenshtein Distance: 14
Species: pal-let-7e-5p - pal-let-7e-3p | Levenshtein Distance: 14
Species: pal-let-7e-5p - cgr-let-7e | Levenshtein Distance: 15
Species: pal-let-7e-5p - abu-let-7e | Levenshtein Distance: 3
Species: pal-let-7e-5p - mze-let-7e | Levenshtein Distance: 3
Species: pal-let-7e-5p - nbr-let-7e | Levenshtein Distance: 3
Species: pal-let-7e-5p - oni-let-7e | Levenshtein Distance: 3
Species: pal-let-7e-5p - pny-let-7e | Levenshtein Distance: 3
Species: pal-let-7e-5p - gmo-let-7e-5p | Levenshtein Distance: 3
Species: pal-let-7e-5p - gmo-let-7e-3p | Levenshtein Distance: 15
Species: pal-let-7e-5p - xla-let-7e-5p | Levenshtein Distance: 2
Species: pal-let-7e-5p - xla-let-7e-3p | Levenshtein Distance: 15
Species: pal-let-7e-5p - cpo-let-7e-5p | Levenshtein Distance: 0
Species: pal-let-7e-5p - cpo-let-7e-3p | Levenshtein Distance: 15
Species: pal-let-7e-5p - ppa-let-7e | Levenshtein Distance: 0
Species: pal-let-7e-5p - cja-let-7e | Levenshtein Distance: 0
Species: pal-let-7e-5p - nle-let-7e | Levenshtein Distance: 0
Species: pal-let-7e-5p - sbo-let-7e | Levenshtein Distance: 0
Species: pal-let-7e-5p - pha-let-7e | Levenshtein Distance: 0
Species: pal-let-7e-5p - oga-let-7e | Levenshtein Distance: 0
Species: pal-let-7e-3p - cgr-let-7e | Levenshtein Distance: 1
Species: pal-let-7e-3p - abu-let-7e | Levenshtein Distance: 15
Species: pal-let-7e-3p - mze-let-7e | Levenshtein Distance: 15
Species: pal-let-7e-3p - nbr-let-7e | Levenshtein Distance: 15
Species: pal-let-7e-3p - oni-let-7e | Levenshtein Distance: 15
Species: pal-let-7e-3p - pny-let-7e | Levenshtein Distance: 15
Species: pal-let-7e-3p - gmo-let-7e-5p | Levenshtein Distance: 15
Species: pal-let-7e-3p - gmo-let-7e-3p | Levenshtein Distance: 7
Species: pal-let-7e-3p - xla-let-7e-5p | Levenshtein Distance: 15
Species: pal-let-7e-3p - xla-let-7e-3p | Levenshtein Distance: 7
Species: pal-let-7e-3p - cpo-let-7e-5p | Levenshtein Distance: 14
Species: pal-let-7e-3p - cpo-let-7e-3p | Levenshtein Distance: 1
Species: pal-let-7e-3p - ppa-let-7e | Levenshtein Distance: 14
Species: pal-let-7e-3p - cja-let-7e | Levenshtein Distance: 14
Species: pal-let-7e-3p - nle-let-7e | Levenshtein Distance: 14
Species: pal-let-7e-3p - sbo-let-7e | Levenshtein Distance: 14
Species: pal-let-7e-3p - pha-let-7e | Levenshtein Distance: 14
Species: pal-let-7e-3p - oga-let-7e | Levenshtein Distance: 14
Species: cgr-let-7e - abu-let-7e | Levenshtein Distance: 16
Species: cgr-let-7e - mze-let-7e | Levenshtein Distance: 16
Species: cgr-let-7e - nbr-let-7e | Levenshtein Distance: 16
Species: cgr-let-7e - oni-let-7e | Levenshtein Distance: 16
Species: cgr-let-7e - pny-let-7e | Levenshtein Distance: 16
Species: cgr-let-7e - gmo-let-7e-5p | Levenshtein Distance: 16
Species: cgr-let-7e - gmo-let-7e-3p | Levenshtein Distance: 6
Species: cgr-let-7e - xla-let-7e-5p | Levenshtein Distance: 16
Species: cgr-let-7e - xla-let-7e-3p | Levenshtein Distance: 6
Species: cgr-let-7e - cpo-let-7e-5p | Levenshtein Distance: 15
Species: cgr-let-7e - cpo-let-7e-3p | Levenshtein Distance: 0
Species: cgr-let-7e - ppa-let-7e | Levenshtein Distance: 15
Species: cgr-let-7e - cja-let-7e | Levenshtein Distance: 15
Species: cgr-let-7e - nle-let-7e | Levenshtein Distance: 15
Species: cgr-let-7e - sbo-let-7e | Levenshtein Distance: 15
Species: cgr-let-7e - pha-let-7e | Levenshtein Distance: 15
Species: cgr-let-7e - oga-let-7e | Levenshtein Distance: 15
Species: abu-let-7e - mze-let-7e | Levenshtein Distance: 0
Species: abu-let-7e - nbr-let-7e | Levenshtein Distance: 0
Species: abu-let-7e - oni-let-7e | Levenshtein Distance: 0
Species: abu-let-7e - pny-let-7e | Levenshtein Distance: 0
Species: abu-let-7e - gmo-let-7e-5p | Levenshtein Distance: 0
Species: abu-let-7e - gmo-let-7e-3p | Levenshtein Distance: 15
Species: abu-let-7e - xla-let-7e-5p | Levenshtein Distance: 3
Species: abu-let-7e - xla-let-7e-3p | Levenshtein Distance: 15
Species: abu-let-7e - cpo-let-7e-5p | Levenshtein Distance: 3
Species: abu-let-7e - cpo-let-7e-3p | Levenshtein Distance: 16
Species: abu-let-7e - ppa-let-7e | Levenshtein Distance: 3
Species: abu-let-7e - cja-let-7e | Levenshtein Distance: 3
Species: abu-let-7e - nle-let-7e | Levenshtein Distance: 3
Species: abu-let-7e - sbo-let-7e | Levenshtein Distance: 3
Species: abu-let-7e - pha-let-7e | Levenshtein Distance: 3
Species: abu-let-7e - oga-let-7e | Levenshtein Distance: 3
Species: mze-let-7e - nbr-let-7e | Levenshtein Distance: 0
Species: mze-let-7e - oni-let-7e | Levenshtein Distance: 0
Species: mze-let-7e - pny-let-7e | Levenshtein Distance: 0
Species: mze-let-7e - gmo-let-7e-5p | Levenshtein Distance: 0
Species: mze-let-7e - gmo-let-7e-3p | Levenshtein Distance: 15
Species: mze-let-7e - xla-let-7e-5p | Levenshtein Distance: 3
Species: mze-let-7e - xla-let-7e-3p | Levenshtein Distance: 15
Species: mze-let-7e - cpo-let-7e-5p | Levenshtein Distance: 3
Species: mze-let-7e - cpo-let-7e-3p | Levenshtein Distance: 16
Species: mze-let-7e - ppa-let-7e | Levenshtein Distance: 3
Species: mze-let-7e - cja-let-7e | Levenshtein Distance: 3
Species: mze-let-7e - nle-let-7e | Levenshtein Distance: 3
Species: mze-let-7e - sbo-let-7e | Levenshtein Distance: 3
Species: mze-let-7e - pha-let-7e | Levenshtein Distance: 3
Species: mze-let-7e - oga-let-7e | Levenshtein Distance: 3
Species: nbr-let-7e - oni-let-7e | Levenshtein Distance: 0
Species: nbr-let-7e - pny-let-7e | Levenshtein Distance: 0
Species: nbr-let-7e - gmo-let-7e-5p | Levenshtein Distance: 0
Species: nbr-let-7e - gmo-let-7e-3p | Levenshtein Distance: 15
Species: nbr-let-7e - xla-let-7e-5p | Levenshtein Distance: 3
Species: nbr-let-7e - xla-let-7e-3p | Levenshtein Distance: 15
Species: nbr-let-7e - cpo-let-7e-5p | Levenshtein Distance: 3
Species: nbr-let-7e - cpo-let-7e-3p | Levenshtein Distance: 16
Species: nbr-let-7e - ppa-let-7e | Levenshtein Distance: 3
Species: nbr-let-7e - cja-let-7e | Levenshtein Distance: 3
Species: nbr-let-7e - nle-let-7e | Levenshtein Distance: 3
Species: nbr-let-7e - sbo-let-7e | Levenshtein Distance: 3
Species: nbr-let-7e - pha-let-7e | Levenshtein Distance: 3
Species: nbr-let-7e - oga-let-7e | Levenshtein Distance: 3
Species: oni-let-7e - pny-let-7e | Levenshtein Distance: 0
Species: oni-let-7e - gmo-let-7e-5p | Levenshtein Distance: 0
Species: oni-let-7e - gmo-let-7e-3p | Levenshtein Distance: 15
Species: oni-let-7e - xla-let-7e-5p | Levenshtein Distance: 3
Species: oni-let-7e - xla-let-7e-3p | Levenshtein Distance: 15
Species: oni-let-7e - cpo-let-7e-5p | Levenshtein Distance: 3
Species: oni-let-7e - cpo-let-7e-3p | Levenshtein Distance: 16
Species: oni-let-7e - ppa-let-7e | Levenshtein Distance: 3
Species: oni-let-7e - cja-let-7e | Levenshtein Distance: 3
Species: oni-let-7e - nle-let-7e | Levenshtein Distance: 3
Species: oni-let-7e - sbo-let-7e | Levenshtein Distance: 3
Species: oni-let-7e - pha-let-7e | Levenshtein Distance: 3
Species: oni-let-7e - oga-let-7e | Levenshtein Distance: 3
Species: pny-let-7e - gmo-let-7e-5p | Levenshtein Distance: 0
Species: pny-let-7e - gmo-let-7e-3p | Levenshtein Distance: 15
Species: pny-let-7e - xla-let-7e-5p | Levenshtein Distance: 3
Species: pny-let-7e - xla-let-7e-3p | Levenshtein Distance: 15
Species: pny-let-7e - cpo-let-7e-5p | Levenshtein Distance: 3
Species: pny-let-7e - cpo-let-7e-3p | Levenshtein Distance: 16
Species: pny-let-7e - ppa-let-7e | Levenshtein Distance: 3
Species: pny-let-7e - cja-let-7e | Levenshtein Distance: 3
Species: pny-let-7e - nle-let-7e | Levenshtein Distance: 3
Species: pny-let-7e - sbo-let-7e | Levenshtein Distance: 3
Species: pny-let-7e - pha-let-7e | Levenshtein Distance: 3
Species: pny-let-7e - oga-let-7e | Levenshtein Distance: 3
Species: gmo-let-7e-5p - gmo-let-7e-3p | Levenshtein Distance: 15
Species: gmo-let-7e-5p - xla-let-7e-5p | Levenshtein Distance: 3
Species: gmo-let-7e-5p - xla-let-7e-3p | Levenshtein Distance: 15
Species: gmo-let-7e-5p - cpo-let-7e-5p | Levenshtein Distance: 3
Species: gmo-let-7e-5p - cpo-let-7e-3p | Levenshtein Distance: 16
Species: gmo-let-7e-5p - ppa-let-7e | Levenshtein Distance: 3
Species: gmo-let-7e-5p - cja-let-7e | Levenshtein Distance: 3
Species: gmo-let-7e-5p - nle-let-7e | Levenshtein Distance: 3
Species: gmo-let-7e-5p - sbo-let-7e | Levenshtein Distance: 3
Species: gmo-let-7e-5p - pha-let-7e | Levenshtein Distance: 3
Species: gmo-let-7e-5p - oga-let-7e | Levenshtein Distance: 3
Species: gmo-let-7e-3p - xla-let-7e-5p | Levenshtein Distance: 15
Species: gmo-let-7e-3p - xla-let-7e-3p | Levenshtein Distance: 0
Species: gmo-let-7e-3p - cpo-let-7e-5p | Levenshtein Distance: 15
Species: gmo-let-7e-3p - cpo-let-7e-3p | Levenshtein Distance: 6
Species: gmo-let-7e-3p - ppa-let-7e | Levenshtein Distance: 15
Species: gmo-let-7e-3p - cja-let-7e | Levenshtein Distance: 15
Species: gmo-let-7e-3p - nle-let-7e | Levenshtein Distance: 15
Species: gmo-let-7e-3p - sbo-let-7e | Levenshtein Distance: 15
Species: gmo-let-7e-3p - pha-let-7e | Levenshtein Distance: 15
Species: gmo-let-7e-3p - oga-let-7e | Levenshtein Distance: 15
Species: xla-let-7e-5p - xla-let-7e-3p | Levenshtein Distance: 15
Species: xla-let-7e-5p - cpo-let-7e-5p | Levenshtein Distance: 2
Species: xla-let-7e-5p - cpo-let-7e-3p | Levenshtein Distance: 16
Species: xla-let-7e-5p - ppa-let-7e | Levenshtein Distance: 2
Species: xla-let-7e-5p - cja-let-7e | Levenshtein Distance: 2
Species: xla-let-7e-5p - nle-let-7e | Levenshtein Distance: 2
Species: xla-let-7e-5p - sbo-let-7e | Levenshtein Distance: 2
Species: xla-let-7e-5p - pha-let-7e | Levenshtein Distance: 2
Species: xla-let-7e-5p - oga-let-7e | Levenshtein Distance: 2
Species: xla-let-7e-3p - cpo-let-7e-5p | Levenshtein Distance: 15
Species: xla-let-7e-3p - cpo-let-7e-3p | Levenshtein Distance: 6
Species: xla-let-7e-3p - ppa-let-7e | Levenshtein Distance: 15
Species: xla-let-7e-3p - cja-let-7e | Levenshtein Distance: 15
Species: xla-let-7e-3p - nle-let-7e | Levenshtein Distance: 15
Species: xla-let-7e-3p - sbo-let-7e | Levenshtein Distance: 15
Species: xla-let-7e-3p - pha-let-7e | Levenshtein Distance: 15
Species: xla-let-7e-3p - oga-let-7e | Levenshtein Distance: 15
Species: cpo-let-7e-5p - cpo-let-7e-3p | Levenshtein Distance: 15
Species: cpo-let-7e-5p - ppa-let-7e | Levenshtein Distance: 0
Species: cpo-let-7e-5p - cja-let-7e | Levenshtein Distance: 0
Species: cpo-let-7e-5p - nle-let-7e | Levenshtein Distance: 0
Species: cpo-let-7e-5p - sbo-let-7e | Levenshtein Distance: 0
Species: cpo-let-7e-5p - pha-let-7e | Levenshtein Distance: 0
Species: cpo-let-7e-5p - oga-let-7e | Levenshtein Distance: 0
Species: cpo-let-7e-3p - ppa-let-7e | Levenshtein Distance: 15
Species: cpo-let-7e-3p - cja-let-7e | Levenshtein Distance: 15
Species: cpo-let-7e-3p - nle-let-7e | Levenshtein Distance: 15
Species: cpo-let-7e-3p - sbo-let-7e | Levenshtein Distance: 15
Species: cpo-let-7e-3p - pha-let-7e | Levenshtein Distance: 15
Species: cpo-let-7e-3p - oga-let-7e | Levenshtein Distance: 15
Species: ppa-let-7e - cja-let-7e | Levenshtein Distance: 0
Species: ppa-let-7e - nle-let-7e | Levenshtein Distance: 0
Species: ppa-let-7e - sbo-let-7e | Levenshtein Distance: 0
Species: ppa-let-7e - pha-let-7e | Levenshtein Distance: 0
Species: ppa-let-7e - oga-let-7e | Levenshtein Distance: 0
Species: cja-let-7e - nle-let-7e | Levenshtein Distance: 0
Species: cja-let-7e - sbo-let-7e | Levenshtein Distance: 0
Species: cja-let-7e - pha-let-7e | Levenshtein Distance: 0
Species: cja-let-7e - oga-let-7e | Levenshtein Distance: 0
Species: nle-let-7e - sbo-let-7e | Levenshtein Distance: 0
Species: nle-let-7e - pha-let-7e | Levenshtein Distance: 0
Species: nle-let-7e - oga-let-7e | Levenshtein Distance: 0
Species: sbo-let-7e - pha-let-7e | Levenshtein Distance: 0
Species: sbo-let-7e - oga-let-7e | Levenshtein Distance: 0
Species: pha-let-7e - oga-let-7e | Levenshtein Distance: 0
Total 'let-7e' miRNAs: 64
Average Levenshtein Distance of 'let-7e' miRNAs: 7.68
## the code for plot of let-7e miRNA across all species is the following:
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
## the matplot for let-7e miRNA across all species is the following:
![Figure_5](https://github.com/aberakenea/Home-take-exam-project/assets/130226484/cc1cd506-0956-4586-bca0-895629df91ca)

#### task: repeat for all let-7 miRNAs and plot
## ANSWER:
## NOTE: we can modify the above code of levensthein distance let-7e miRNAs and plot for all let-7 miRNAs and plot accros all species, but still I can not able to run code for all let-7 miRNAs and plots accros all species due to the shortage of my PCs RAM, Dear Profecer SIMON, you can run the code for all let_7 miRNAs and plots hope your computer do have enuogh storage.



