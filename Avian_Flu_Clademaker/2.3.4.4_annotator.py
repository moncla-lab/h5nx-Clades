import baltic as bt
import re

input_file_path = "Data/2.3.4.4-0906.fas"
output_path = "Output/2.3.4.4_Guide_Tree.fasta"
full_tree_path = "Data/FastaGuide_ha.fasta"
compiled_output_path = "Output/h5nxLABEL2-3-4-4Anottated_ha.fasta"

with open(input_file_path, "r") as guide:
    source_lines = guide.read().splitlines()

with open(full_tree_path, "r") as full_fasta:
    full_lines = full_fasta.read().splitlines()

def format_guide(guide, output):
    lines_for_file = []
    strains = []
    for line in guide:
        if line and line[0] == ">":
            skip = False
            if "}" not in line and "{" not in line:
                skip = True

            if "like" in line:
                skip = True

            Clade_split = line.split("{")
            Date_split = Clade_split[0].split("/")
            unformatted_date = Date_split[-1]
            strain = Date_split
            strain.pop()

            if "||" in unformatted_date:
                Date_split_final = unformatted_date.split("||")
                Date = Date_split_final[-1]

            elif unformatted_date[-1] == "_":
                Date_split_final = unformatted_date.split("_")
                if len(Date_split_final[0]) == 4:
                    Date = Date_split_final[0] + "-XX-XX"

                elif len(unformatted_date[0]) != 4:
                    Date = Date_split_final[0]
                    if Date == "14":
                        Date = "2014-XX-XX"


            clade_final = Clade_split[-1].split("}")
            clade = clade_final[0]   
            strain = "/".join(strain)
            strain = strain.replace("-", "_")
            metadata = strain + "|" + Date + "|" + clade

            if strain in strains:
                skip = True
            else:
                strains.append(strain)

            if skip == False:
                lines_for_file.append(metadata)

        elif skip == False:
            lines_for_file.append(line)

    return lines_for_file, strains


def append_full(full_fasta, lines_from_guide, output_file_path, strains):
    for line in full_fasta:
        if line and line[0] == ">":
            skip = False
            metadata = line.split("|")
            strain = metadata[0]
            if strain in strains:
                skip = True

            if skip == False:
                lines_from_guide.append(line)

        elif skip == False:
            lines_from_guide.append(line)

    with open(output_file_path, "w") as output:
        output.write("\n".join(lines_from_guide))

    return None
            

guide_lines, guide_strains = format_guide(source_lines, output_path)
append_full(full_lines, guide_lines, compiled_output_path, guide_strains)