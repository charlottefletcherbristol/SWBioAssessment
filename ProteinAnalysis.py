#### PROTEIN ANALYSIS TOOL ####

#This code allows for the computation of various physical and chemical parameters for a given protein. These parameters can be deduced solely from the protein sequence. Parameters:
    #protein length
    #molecular weight
    #theoretical isoelectric point (pI)
    #amino acid composition
    #extinction coefficients
    #GRAVY (Grand Average of Hydropathicity)
#Documentation for parameter calculations can be found in the associated README file
    
#IMPORT PROGRAMS
import numpy as np
import math
import matplotlib.pyplot as plt
import Bio                                          #Biopython
from Bio.SeqUtils.ProtParam import ProteinAnalysis  #ProtParam module based on the Expasy Proteomics Server

#INPUT PROTEIN SEQUENCE
protein = input("Input protein sequence:")     #user is prompted to input protein sequence
    
for residue in protein:                         #error checking:
    accepted_chars = "ARNDBCEQZGHILKMFPSTWYV"
    if residue in accepted_chars:              #only one-letter amino     acid code is accepted
        continue
    else:                                      #else, error message displays and code exits
        print("Error: Protein sequence has non-recognised characters. Only one-letter amino acid code accepted.")
        exit()

        
##CALCULATE PARAMETERS##

#Protein length, molecular weight and pI:

print("Number of amino acids:", len(protein))    #print protein length. Len() counts the number of characters
                                                 #in the protein string.

analysed_protein = ProteinAnalysis(protein)                          #use Biopython's ProteinAnalysis function
print("Molecular weight (Da):", analysed_protein.molecular_weight()) #in-built tool to calculate molecular weight
print("Theoretical pI:", analysed_protein.isoelectric_point())       #in-built tool to calculate pI

#Amino acid composition:

Amino_acid_count = analysed_protein.count_amino_acids()        #in-built tool for counting amino acids
print("Amino acid frequencies:", Amino_acid_count)             #tool generates a dictionary

fig = plt.figure()                                             #create bar graph of amino acid frequencies:
ax = fig.add_axes([0,0,1,1])
x = Amino_acid_count.keys()                                    # x-values are dictionary keys, i.e. amino acids
y = Amino_acid_count.values()                                  # y-values are dictionary values, i.e. count data
plt.bar(x, y, color = "silver")
plt.title("Amino Acid Composition", size=12)                   #label axes and give title
plt.xlabel("Amino Acid", size=10)
plt.ylabel("Frequency", size=10)
plt.show()                                                     #display graph in pop-out window
 
#Number of positively and negatively charged amino acids:
                                                                        #positive = R and K
                                                                        #negative = D and E

count_pos = protein.count('R') + protein.count('K')                      #sum the counts of R and K in protein seq
print("Total number of positively charged residues (R + K):", count_pos) #print

count_neg = protein.count('D') + protein.count('E')                      #sum the counts of D and E in protein seq
print("Total number of negatively charged residues (D + E):", count_neg) #print

#Extinction coefficients

Ext_Y = 1490                 #assign appropriate values to extinction variables
Ext_W = 5500 
Ext_Cystine = 125

#calculate extinction coefficient for protein
E_protein = protein.count('Y')*Ext_Y + protein.count('W')*Ext_W + protein.count('C')*0.5*Ext_Cystine
#calculate absorbance of protein
Absorb_protein = E_protein/analysed_protein.molecular_weight() 

#print values
print("Extinction coefficient in units of  M-1 cm-1, at 280 nm measured in water.")
print("Ext. coefficient: ", E_protein)
print("Abs 0.1% (=1 g/l):", Absorb_protein, ", assuming all pairs of cysteine (C) residues form cystines")

#GRAVY (Grand Average of Hydropathicity)

#make a dictionary of the hydropathicity values associated with each amino acid
residue_hydropathicity = {"A" : 1.8, "R" : -4.5, "N" : -3.5, "D" : -3.5, "C" : 2.5, "Q" : -3.5, "E" : -3.5, 
                  "G" : -0.4, "H" : -3.2, "I" : 4.5, "L" : 3.8, "K" : -3.9, "M" : 1.9, "F" : 2.8, 
                  "P" : -1.6, "S" : -0.8, "T" :  -0.7, "W" : -0.9, "Y" : -1.3, "V" : 4.2}

# make a list of hydropathicity scores for each amino acid in the protein sequence
hydropathicity_scores = []                            #create empty list                           
for residue in protein:
        value = residue_hydropathicity.get(residue)   #iterate through every amino acid in the protein, and 
        hydropathicity_scores.append(value)           #append its hydropathicity score to the list 
        
                                                      #GRAVY = mean of this list
print("Grand average of hydropathicity (GRAVY):", np.mean(hydropathicity_scores))

#GRAPH: 5-point rolling average of hydropathicity across protein

#generate rolling average values:
window=5
cut_off=math.floor(int(window)/2)                              #creating the lower cut off value
upper=math.ceil(int(window)/2)                                 #creating the upper cut off value
window_averages=[]                                             #create empty list to contain averages from windows
for i in range(int(cut_off),len(protein)-int(cut_off),1):      #iterate through every score starting from lower cut off                                                                    #value and ending at the upper cul off value
    windows=hydropathicity_scores[i-int(cut_off):i+int(upper)] #defining a window as the position of i-cut-off to i+cut off                                                                #which creates the window
    average=np.mean(windows)                                   #create an average of the window
    window_averages.append(average)                            #append the average to final list of averages used in                                                                      #plotting
    
#plot rolling average values:
plt.axes()
plt.plot(window_averages, '-', c='midnightblue',linewidth=0.8) # plotting the average scores as a navy line, with width of                                                                #0.8
plt.xlabel("Amino Acid Position",size=10)                      # label x-axis
plt.ylabel("Average Hydropathicity",size=10)                   # label y-axis
plt.tick_params(axis='x', labelsize=8)                         # defining axes parameters (for legibility)
plt.tick_params(axis='y',labelsize=8)
plt.title("5-point Rolling Average of Hydropathicity Across Protein ")            # defining the title of the graph
plt.show()                                                     #display graph in pop out window








