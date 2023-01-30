import re

Annotated_Clades = "h5n1_ha_clades_final.txt"
Sequences_Without_Clades = "h5nX_ha.fa"
Output_File_Path = "h5n1_ha_with_clades.fasta"

"""
In this file, we are going to add our clades onto our metadata from our fasta file. We will
take our initial fasta file, and use that data to write a new fasta file where
the metadata contains the clades on the end
"""

clades_dict = {}
lines_for_file = []

with open(Annotated_Clades, "r") as Clades: #here we're grabbing the txt file with our clades in it
    for line in Clades: #for each line in our clades list
        if "CLADE" not in line and "VIRUS STRAIN" not in line: #making sure we're not looking at the headers
            split_line = line.split(" ") #here we're going to split the lines by the spaces
            for item in split_line: #here we're going to go through each item in our line
                if "/"  in item or "Jude" in item: #here we're just looking for the metadata, which should be either the St. Jude strian or have a /
                    strain = ">" + item.strip() #here we're adding the leading carrot and removing any newline characters
                    
                elif item != "":
                    clades_dict[strain] = item.strip() #here we're going to add all the lines which are not just spaces (So only our clade) to our dict

with open(Sequences_Without_Clades, "r") as Sequences: #opening our fasta without sequences
    for Line in Sequences: #for each line
        Line = Line.strip() #removing newline characters
        if Line and Line[0] == ">": #ensuring we are both at a line and that our line starts with a carrot (for the metadata)
            if Line in clades_dict.keys(): #if our strain is already in our clades dict
                clade = clades_dict[Line] #we're going to grab the clade name
                if clade != "UNRECOGNIZABLE": #we're going to ignore unrecognizable clades
                    metadata = Line + "|" + clade #we're going to add the clade to the metadata here
                    lines_for_file.append(metadata) #then we're going to add the metadata to our lines for file list

                else:
                    lines_for_file.append(Line) #adding the sequences to the lines for file (unaltered)

with open(Output_File_Path, "w") as output_file:
    output_file.write('\n'.join(lines_for_file)) #since we removed the newlines earlier we're going to have to add them back here when we write to the new file