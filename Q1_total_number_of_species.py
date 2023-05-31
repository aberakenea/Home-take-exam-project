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