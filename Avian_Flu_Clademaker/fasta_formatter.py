"""
This file is going to reformat our fasta file so that it can be used by both LABEL and by nextstrain
it is going to write to a new file all the reformatted lines, leaving the original file alone, in
case you need that original file again
"""

import re

Input_File_Path = "Data/h5nX_ha_unformatted.fa" #file where the unformatted sequences are held
Output_File_Path_Nextstrain = "Output/h5n1_ha.fasta" #new file with properly formatted lines


with open(Input_File_Path, "r") as source_file:
    source_lines = source_file.read().splitlines()


"""
This function is going to act as our first pass through our file, it will reformat the metadata
and assign sequences to the metadata in a nested dictionary

output:

Dict = {strain : {date: date_specificity variable, metadata: string, sequence: [sequence], sequence_length: X}}

metadata which is assigned multiple sequences will have values which are a list, while metadata
that is assigned only a single sequence will have a string as its value
"""
def first_pass(source_lines):
    sq_dict = {}
    for line in source_lines: #here we're going to go through each line
        if line and line[0] == ">": #if our line is a line and is metadata
            skip_line = False #by intializing skip to false here, the line will only be skipped if we find a reason to skip it, otherwise, it will be added
            match = re.search("(\/(\d\d)\s)", line) #here we're searching for dates that are only 2 numbers and nothing else

            if match: #if we find a match
                number = int(match.groups()[1]) #here we're going to grab the number we found using regex

                if number < 24: #if our 2-digit date is less than 24 we're going to assume it's from the 2000's
                    line = line.replace(match.groups()[0], "/20" + str(number))
                else: #otherwise we're going to assume it's from the 1900's
                    line = line.replace(match.groups()[0], "/19" + str(number))

            line = line.replace(" |", "|") #remove space between |
            line = line.replace(" ", "_") #replace spaces with underscores
            line = line.replace("/|", "|") #reformat


            segment_checker = line.split("|") #split the metadata based on |

            if len(segment_checker) < 5: #this won't always work, but for our data it will
                line = line[0] + "|" + line[1:len(segment_checker)] # if our metadata has less than 5 segments we can assume it doesn't have a species, and we can add a blank segment
            
            if "|--" in line:
                skip_line = True #if we find no date, we are going to skip the line because we only want metadata with date information

            elif "--" in line and "|--" not in line: #if we have a year but are missing month and day
                line = line.replace("--", "-XX-XX")
                date_specificity = 1 #only contains year natively
            elif "-|" in line: #if we are missing only a day
                line = line.replace("-|", "-XX|")
                date_specificity = 2 # contains year and month natively
            else:
                date_specificity = 3 #contains year, month, and day natively


            strain = segment_checker[0] #here we're going to grab the strain name, which should be the very first segment of our metadata

            if skip_line == False: #we're only going to append our data to our dictionary if it has not been flagged to be skipped
                if strain in sq_dict.keys(): #if strain is already in the dictionary
                    sq_dict[strain].append( #we're adding the new data to the strain that's already there
                        {"date": date_specificity, #we can parse this out later in the second pass
                        "metadata": line, 
                        "sequence": [],
                        "sequence_length": 0})
                else: #if the strain is not already in the dictionary we're going to add it
                    sq_dict.update({strain: []})
                    sq_dict[strain].append(
                        {"date": date_specificity, 
                        "metadata": line, 
                        "sequence":[],
                        "sequence_length": 0})

        else: #this is going to be the sequences, since anything which does not start with a > will be a sequence
            if skip_line == False:
                sq_dict[strain][-1]["sequence"].append(line)
                sq_dict[strain][-1]["sequence_length"] += len(line.strip())

    return sq_dict

"""
Here we are going to go through all the data we just got in our first pass and remove any duplicates
We are going to choose which strain to keep based on the following variables with those at the top
of the list having more value than those lower on the list:

- Date specificity
- Sequence Length

This will write a new file where the metadata has been reformatted and the duplicates have been removed
also returns the lines which are going to be used for the file as a list
"""
def second_pass(first):
    lines_for_file = []
    for strain in first: #for each strain in our dictionary (we are currently 1 dict deep)
        most_specific = 0 #here we are going to intialize our check variable to 0 that way if this is the first sequence it will be added no matter what
        if len(first[strain]) > 1: #this will only trigger for strains which have a duplicate
            for i in range(len(first[strain])):
                date = first[strain][i]["date"]
                """
                The next bit of code gets a little complex so I'm going to break it down here rather then comment on every line:

                So first, we're going to go through each item in our strain. This will include at least 2 sets of metadata, sequences, etc.
                Every time we do, we're going to grab its date_specificity, which we're setting to "date"
                from there, we're going to compare it to most_specific

                FOR THE FIRST SQUENCE:
                    Because most_specific is initialized to zero to start, the following code is always going to choose the first strain as its most specific
                    to start, but this will change 

                Once we reach the next set of metadata, we're going to have a new date. If the new date specificity is higher
                than the current, we'll ammend our index that we want to keep

                If our date_specificity is the same, we're going to pick based on sequence
                If our current sequence is longer, nothing will change, if our new sequence is longer we'll replace our index to keep with that 
                sequence's index.
                
                This will allow us to select which sequence we want to keep (along with its corresponding metadata) based on both date
                and sequence length
                """
                if date > most_specific:
                    most_specific = date
                    index_to_keep = i
                
                elif date == most_specific:
                    if first[strain][index_to_keep]["sequence_length"] < first[strain][i]["sequence_length"]:
                        index_to_keep = i

            sequence = first[strain][index_to_keep]["sequence"]
            metadata = first[strain][index_to_keep]["metadata"]

            lines_for_file.append(metadata) #here we are adding the metadata to our file
            for index in range(len(sequence)): #followed by the sequence, which will be a list of lines
                lines_for_file.append(sequence[index])

        else: #if we only have 1 sequence for a strain, we can just add it right away without having to check for dupes
            metadata = first[strain][0]["metadata"]
            sequence = first[strain][0]["sequence"]
            lines_for_file.append(metadata)
            for indx in range(len(sequence)):
                lines_for_file.append(sequence[indx])

    with open(Output_File_Path_Nextstrain, "w") as output_file:
        #joining the lines with newline characters because we removed them earlier
        output_file.write('\n'.join(lines_for_file)) 

    return lines_for_file

with_dupes = first_pass(source_lines)
second_pass(with_dupes)