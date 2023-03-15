import csv
LABEL_path = "Data/LABEL_Guide.fasta"
WILEY_Path = "Data/doi.org-10.1111-j.1750-2659.2011.00298.x_guide_tree.fasta"
Output_Path = "Output/compiled_tree.csv"
Output_Fasta_LABEL = "Output/FastaGuide_ha.fasta"
Output_Fasta_WILEY = "Output/WILEY_Guide_ha.fasta"
Output_Compiled = "Output/CompiledGuide_ha.fasta"


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

def WILEY_Guide(path, strains, data, lines_for_file):
    skip = None
    with open(path, 'r') as source_file:
        source_lines = source_file.read().splitlines()

    for line in source_lines:
        if line and line[0] == ">":
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
                strains[strain].append("WILEY_1")

            elif strain not in strains.keys():
                strains[strain] = ["WILEY_1"]

            metadata = strain + "|" + Date + "|" + clade
            Strain = strain + "|" + Date + "|"
            if metadata not in lines_for_file and Strain not in lines_for_file and metadata not in data.keys():
                lines_for_file.append(metadata)
                skip = False
                
            else:
                skip = True
    
        elif skip == False:
            lines_for_file.append(line)
            data[strain]["sequence"].append(line)

    return strains, data, lines_for_file

LStrains, LData, fasta = LABEL_guide(LABEL_path)
#strains, data, LFasta = WILEY_Guide(WILEY_Path, LStrains, LData, fasta)


with open(Output_Fasta_LABEL, "w") as output_file:
    output_file.write("\n".join(fasta))



"""
Sources:

LABEL - LABEL training data
WILEY_1 - https://onlinelibrary.wiley.com/doi/10.1111/j.1750-2659.2011.00298.x

"""
"""
with open(Output_Path, "w", newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter=';', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        for object in data.keys():
            found_in = ",".join(strains[object])
            writer.writerow([object, data[object]["clade"], found_in, data[object]["sequence"]])
"""