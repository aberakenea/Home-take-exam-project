'''
Created on May 28, 2023

@author: ABERA KENEA
'''
file_path = r"C:\Users\wku\Advancedprogram\project\GC_calc-complexity-project\mature.fa"

def count_let7_miRNAs(file_path):
    let7_count = 0
   
    with open(file_path, "r") as file:
        lines = file.readlines()
       
        for line in lines:
            if line.startswith(">"):
                miRNA_name = line.strip()[1:]
               
                if "let-7" in miRNA_name:
                    let7_count += 1
   
    return let7_count

let7_count = count_let7_miRNAs(file_path)
print("Total number of let-7 miRNAs across all species:", let7_count)