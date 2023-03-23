'''
This code is made to annotate 2.3.4.4 clades onto our tree, it can be used to add any clades onto the tree simply change the input file
and some of the parsing rules.
'''

import baltic as bt
import re

input_file_path = "Data/2.3.4.4-0906.fas"
output_path = "Output/2.3.4.4_Guide_Tree.fasta"
full_tree_path = "Data/FastaGuide_ha.fasta"
compiled_output_path = "Output/h5nxLABEL2-3-4-4Anottated_ha.fasta"

#opening our input file which contains the clades we want to add
with open(input_file_path, "r") as guide:
    source_lines = guide.read().splitlines()

#opening our tree fasta where we want to add these clades
with open(full_tree_path, "r") as full_fasta:
    full_lines = full_fasta.read().splitlines()


'''
Here we are going to format our guide file (which has the clades we want to add)
Because each file is formatted in a unique way, this particular piece of code will likely have
to be changed depending on the file you're using as your guide

It's going to ouput two lists,

lines for file, which is just every line in the guide fasta separated by commas
and strains which contains every strain
'''
def format_guide(guide, output):
    lines_for_file = []
    strains = []
    for line in guide: 
        if line and line[0] == ">": #if our line is metadata
            skip = False #initialize our skip variable
            if "}" not in line and "{" not in line: #because in our guide, only strains with the 2.3.4.4 clades have {} I just picked these characters out
                skip = True #if we don't have the 2.3.4.4 reclassification, we don't care about that strain

            if "like" in line: #we don't want any 2.3.4.4-like clades
                skip = True

            Clade_split = line.split("{") #we're gonna split by { so we should get our strain on the left and the new clade on the right
            Date_split = Clade_split[0].split("/") #here we're going to split by / because we need to grab the date here
            unformatted_date = Date_split[-1] #here we're grabbing the actaul date which we will format later
            strain = Date_split #here our strain is the list of strings which came from our date_split where we split by /
            strain.pop() #we're gonna remove the last instance, which should be our date, which should leave us with just the strain

            if "||" in unformatted_date: #if we have || in our date we're going to remove the ||
                Date_split_final = unformatted_date.split("||")
                Date = Date_split_final[-1]

            #this code is just going to replace dates which have no days or months with -XX-XX
            elif unformatted_date[-1] == "_":
                Date_split_final = unformatted_date.split("_")
                if len(Date_split_final[0]) == 4:
                    Date = Date_split_final[0] + "-XX-XX"

                #we have 1 date that is just 14 so I'm going to replace that with 2014-XX-XX
                elif len(unformatted_date[0]) != 4:
                    Date = Date_split_final[0]
                    if Date == "14":
                        Date = "2014-XX-XX"

            #we're gonna grab the clade from our split earlier
            #this is just a bunch of formatting, it will need to be changed based on what file you're using
            clade_final = Clade_split[-1].split("}")
            clade = clade_final[0]   
            strain = "/".join(strain)
            strain = strain.replace("-", "_")
            metadata = strain + "|" + Date + "|" + clade #here we're grabbing our final metadata

            if strain in strains: #if we already have our strain in our strains list
                skip = True #we skip
            else:
                strains.append(strain) #otherwise we're going to add our strain

            if skip == False: #finally, if skip is false we add it to our lines_for_file list
                lines_for_file.append(metadata)

        elif skip == False: #here we're adding the actual sequence if our strain skip = false
            lines_for_file.append(line)

    return lines_for_file, strains


'''
Here we're going to actually write our new fasta file

'''
def append_full(full_fasta, lines_from_guide, output_file_path, strains):
    for line in full_fasta: #for lines in our full trtee
        if line and line[0] == ">": #if our line is metadata
            skip = False #intialize our skip variable
            metadata = line.split("|") #split our metadata by | so we can grab things easier
            strain = metadata[0] #grab the strain
            if strain in strains: #if we already have our strain from our guide file
                skip = True #skip our strain in our full file

            if skip == False: #add to lines
                lines_from_guide.append(line)

        elif skip == False: #add to lines
            lines_from_guide.append(line)

    #write to new file
    with open(output_file_path, "w") as output:
        output.write("\n".join(lines_from_guide))

    return None
            
#formatting our guide
guide_lines, guide_strains = format_guide(source_lines, output_path)

#appending our fasta with the new clades
append_full(full_lines, guide_lines, compiled_output_path, guide_strains)