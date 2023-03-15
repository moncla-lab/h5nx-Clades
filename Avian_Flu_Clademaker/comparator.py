import re
import baltic as bt
import decimal

import Mutation_Finder as mf
import Mutation_Characterization as mc

Input_Tree_Path_Full = "Data/flu_avian_h5nxFullNonMafft-noDT382-S4220-0813082-005_ha_clock_enabled.JSON"
Input_Tree_Path_LABEL = "Data/LABELGuideTree/flu_avian_FastaGuide_ha.JSON"
Input_Tree_Path_WILEY = "Data/flu_avian_WILEYGuide_ha.JSON"
mytree_Full, mymeta_Full = bt.loadJSON(Input_Tree_Path_LABEL)
mytree_WILEY, mymeta_WILEY = bt.loadJSON(Input_Tree_Path_WILEY)
mytree_LABEL, mymeta_LABEL = bt.loadJSON(Input_Tree_Path_LABEL)
outputWL = "Output/comparison_WILEY_LABEL_percentage.txt"
OutputWF = "Output/comparison_WILEY_FULL_percentage.txt"
OutputLF = "Output/comparison_LABEL_FULL_percentage.txt"
exclude_list = ["A/falcon/Saudi_Arabia/D1795/2005", 'A/chicken/Giza/CAI37/2009']
parent_dict = ['2.1.3.2', '1', '1.1', '2.1.1', 
'2.2', '2.2.1', '2.2.1.1', '2.3.2.1', '2.3.4', '7', 
'2.2.2', '2.3.2', '1.1.1', '6', '2.1.3.1']


allCladesF, leafCladesF = mf.leaf_clades(mytree_Full, exclude_list)
nodeCladesF, uncleanNodesF = mf.node_clades(mytree_Full, leafCladesF)
nodesF = mc.find_nodes(nodeCladesF, mytree_Full, parent_dict)
HAF = mc.node_HA_muts(mytree_Full, nodesF)
nucF = mc.node_nuc_muts(mytree_Full, nodesF)

allCladesW, leafCladesW = mf.leaf_clades(mytree_WILEY, exclude_list)
nodeCladesW, uncleanNodesW = mf.node_clades(mytree_WILEY, leafCladesW)
nodesW = mc.find_nodes(nodeCladesW, mytree_WILEY, parent_dict)
HAW = mc.node_HA_muts(mytree_WILEY, nodesW)
nucW = mc.node_nuc_muts(mytree_WILEY, nodesW)

allCladesL, leafCladesL = mf.leaf_clades(mytree_LABEL, exclude_list)
nodeCladesL, uncleanNodesL = mf.node_clades(mytree_LABEL, leafCladesL)
nodesL = mc.find_nodes(nodeCladesL, mytree_LABEL, parent_dict)
HAL = mc.node_HA_muts(mytree_LABEL, nodesL)
nucL = mc.node_nuc_muts(mytree_LABEL, nodesL)

def comparison(HaF, NucF, HaG, NucG):
    compare = {}
    shared = {}
    for clade in HaF.keys():
        similar = []
        all = []
        shared[clade] = {"HA": [], "Nuc": []}
        compare[clade] = {"HA": None, "Nuc": None}
        if clade in HaG.keys() and clade in HaF.keys():
            for mutation in HaF[clade]:
                if mutation in HaG[clade]:
                    similar.append(mutation)
                    all.append(mutation)

                elif mutation not in all: 
                    all.append(mutation)

            for mutation in HaG[clade]:
                if mutation in all or mutation in similar:
                    pass

                else:
                    all.append(mutation)

        else:
            compare[clade]["HA"] = 0

        if compare[clade]["HA"] != 0:
            shared[clade]["HA"] = similar
            number = len(similar)/len(all)
            number = round(number, 5)
            compare[clade]["HA"] = number
    
    for clade in NucF.keys():
        similar = []
        all = []
        if clade not in compare.keys():
            compare[clade] = {"HA": None, "Nuc": None}
            shared[clade] = {"HA": [], "Nuc": []}

        if clade in NucG.keys() and clade in NucF.keys():
            for mutation in NucF[clade]:
                if mutation in NucG[clade]:
                    similar.append(mutation)
                    all.append(mutation)

                else:
                    all.append(mutation)

            for mutation in NucG[clade]:
                if mutation in all or mutation in similar:
                    pass
                else:
                    all.append(mutation)

        else:
            compare[clade]["Nuc"] = 0

        if compare[clade]["Nuc"] != 0:
            shared[clade]["Nuc"] = similar
            number = len(similar)/len(all)
            number = round(number, 5)
            compare[clade]["Nuc"] = number

    return compare, shared

def compile(shared, shared_2):
    compiled = {}
    for clade in shared.keys():
        compiled[clade] = {"HA": [], "Nuc": []}
        if len(shared[clade]["HA"]) != 0:
            for HAmutation in shared[clade]["HA"]:
                if HAmutation in shared_2[clade]["HA"]:
                    compiled[clade]["HA"].append(HAmutation)
        if len(shared[clade]["Nuc"]) != 0:
            for Nmutation in shared[clade]["Nuc"]:
                if Nmutation in shared_2[clade]["Nuc"]:
                    compiled[clade]["Nuc"].append(Nmutation)

    return compiled

percentages, shared = comparison(HAF, nucF, HAL, nucL)
Percs, shared_2 = comparison(HAF, nucF, HAL, nucL)

compiled = compile(shared, shared_2)


with open(OutputLF, "w") as output_file:
    output_file.write("Clade \t\t HA \t\t Nucleotide\n")
    for clade in Percs.keys():
        HA = Percs[clade]["HA"]
        Nuc = Percs[clade]["Nuc"]
        if HA != None:
            HA = round(100*HA, 2)
        if Nuc != None:
            Nuc = round(100*Nuc, 2)
        output_file.write(f"{clade}\t\t%{HA}\t\t%{Nuc}\n")


"""
def write_to_file(x, y, z):
    with open(x, "w") as file:
        file.write("Clade\tGene\tSite\tAlt")
    return None
"""


"""
with open(outputAT, "w") as output_file:
    output_file.write("Clade \t\t HA \t\t Nucleotide\n")
    for clade in shared.keys():
        HA = "\n\t\t".join(compiled[clade]["HA"])
        Nuc = "\n\t\t\t\t".join(compiled[clade]["Nuc"])
        output_file.write(f"{clade}\t\t{HA}\t\t{Nuc}\n\n")
"""
