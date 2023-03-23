import csv
LABEL_path = "Data/LABEL_Guide.fasta"
Output_Fasta_LABEL = "Output/FastaGuide_ha.fasta"

'''
This file is just to format some other files we want.

'''
def LABEL_guide(path):
    data = {}
    strains = {}
    lines_for_file = []
    skip = None

    with open(path, 'r') as source_file:
        source_lines = source_file.read().splitlines()

    for line in source_lines:
        if line and line[0] == ">":
            if "|" in line:
                splitter = line.split("|")
                split = splitter[1].split("{")
                strain = ">" + split[0]

            else:
                split = line.split("{")
                strain = split[0]

            strain = strain.strip("_")
            clade = split[-1].strip("}")

            DataD = strain.split("/")
            Date = DataD[-1]

            if len(Date) == 2:
                if int(Date) > 23:
                    Date = "19" + Date

                else:
                    Date = "20" + Date
            
            Date = Date + "-XX-XX"
            
            if strain not in data.keys():
                data[strain] = {
                    "clade": clade,
                    "sequence": []}

            if strain in strains.keys():
                strains[strain].append("LABEL")

            elif strain not in strains.keys():
                strains[strain] = ["LABEL"]

            metadata = strain + "|" + Date + "|" + clade
            Strain = strain + "|" + Date + "|"
            if metadata not in lines_for_file and Strain not in lines_for_file:
                lines_for_file.append(metadata)
                skip = False
                
            else:
                skip = True
                
        elif skip == False:
            lines_for_file.append(line)
            data[strain]["sequence"].append(line)

    return strains, data, lines_for_file
