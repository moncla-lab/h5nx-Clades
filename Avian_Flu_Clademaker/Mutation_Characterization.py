import re
import baltic as bt

import Mutation_Finder as mf

Input_Tree_Path = "Data/flu_avian_h5nxLABEL2-3-4-4Anottated_ha.json"
#Input_Tree_Path = "Data/LABELGuideTree/flu_avian_FastaGuide_ha.JSON"
mytree, mymeta = bt.loadJSON(Input_Tree_Path)
output_file_path = "Output/defining_clades_LABEL.txt"
#output_file_path = "Output/unique_clades_guide.txt"
'''
exclude_list = ["A/falcon/Saudi_Arabia/D1795/2005", 'A/chicken/Giza/CAI37/2009', 'A/duck/China/E319-2/03', 
'A/whooper_swan/Mongolia/21/2010', 'A/Mallard/WI/944/82', 'A/turkey/TX/1-40-82/1982', 'A/chicken/Pekanbaru/BPPVRII-91-3/2011',
 'A/Fujian/1/2007', 'A/Hong_Kong/156/97', 'A/goose/Guangdong/08/2005', 'A/goose/Jilin/hb/2003', 'A/Fujian/1/2007']
'''
exclude_list = []
parent_dict = ['2.1.3.2', '1', '1.1', '2.1.1', 
'2.2', '2.2.1', '2.2.1.1', '2.3.2.1', '2.3.4', '7', 
'2.2.2', '2.3.2', '6', '2.1.3.1', '2.1.2']

 
def find_nodes(NC, tree, parents):
    output = {}
    node_length = {}
    for item in tree.Objects:
        if item.branchType == "node":
            length = len(item.leaves)
            node_length[item.traits["name"]] = length

    for node in NC:
        if len(NC[node]) == 1 and NC[node][0] not in parents:
            clade = NC[node][0]
            if "Am_nonGsGD" not in clade and "EA_nonGsGD" not in clade and "like" not in clade:
                if clade in output.keys():
                    if node_length[output[clade]] < node_length[node]:
                        output[clade] = node
                
                else:
                    output[clade] = node

        if len(NC[node]) > 1:
            skip = False
            shortest = ['clade', 2000000]
            for clade in NC[node]:
                if len(clade) < shortest[-1]:
                    shortest = [clade, len(clade)]

                elif len(clade) == shortest[-1]:
                    if clade == '2.1.1':
                        shortest = [clade, len(clade)]

                    #elif clade == '1.1.1':
                        #shortest = [clade, len(clade)]

                    elif clade == '6':
                        shortest = [clade, len(clade)]

                    elif clade == '2.1.3.1':
                        shortest = [clade, len(clade)]

                    elif clade == '3':
                        shortest = [clade, len(clade)]

                    elif clade == '2.1.2':
                        shortest = [clade, len(clade)]

            if shortest[0] in parents:
                for clade in NC[node]:
                    if "like" in clade and shortest != "3":
                        skip = True

                    elif shortest[0] == '2.1.1':
                        if '2.1' not in clade[0:3]:
                            skip = True
                    
                    #elif shortest[0] == '1.1.1':
                        #if '1.1' not in clade[0:4]:
                           # skip = True
                    
                    elif shortest[0] == '2.2':
                        if '2.2' not in clade[0:3] and '2.3' not in clade:
                            skip = True
                    
                    elif shortest[0] == '6':
                        if '6' not in clade[0:3] and '7' not in clade[0]:
                            skip = True

                    elif shortest[0] == '2.1.3.1':
                        if '2.1.3' not in clade[0:6]:
                            skip = True

                    elif shortest[0] == '3':
                        if '4' not in clade[0] and '3' not in clade[0]:
                            skip = True

                    elif shortest[0] not in clade[0:shortest[-1]]:
                        skip = True

                    elif shortest[0] not in parents:
                        skip = True

                    elif shortest[0] == '2.1.2':
                        if '2.1.2' not in clade and '2.1.3' not in clade:
                            skip = True

                if skip == False:
                    if shortest[0] in output.keys():
                        if shortest[-1] > node_length[output[shortest[0]]]:
                            output[shortest[0]] = node

                    else:
                        output[shortest[0]] = node
    return output

def node_HA_muts(Tree, nf):
    node_HAmuts = {}
    
    for object in Tree.Objects:
        if object.branchType == "node":
            mutations = object.traits["branch_attrs"]["mutations"]
            name = object.traits["name"]
            for item in nf.keys():
                if nf[item] == name:
                    if "HA" in mutations:
                        node_HAmuts[item] = object.traits["branch_attrs"]["mutations"]["HA"]

    return node_HAmuts

def node_nuc_muts(Tree, nf):
    node_nucmuts = {}
    
    for object in Tree.Objects:
        if object.branchType == "node":
            mutations = object.traits["branch_attrs"]["mutations"]
            name = object.traits["name"]
            for item in nf.keys():
                if nf[item] == name:
                    if "nuc" in mutations:
                        node_nucmuts[item] = object.traits["branch_attrs"]["mutations"]["nuc"]
    return node_nucmuts


allClades, leafClades = mf.leaf_clades(mytree, exclude_list)
nodeClades, unclean = mf.node_clades(mytree, leafClades)
nodes = find_nodes(nodeClades, mytree, parent_dict)
print(nodes)
#HA = node_HA_muts(mytree, nodes)
#nuc = node_nuc_muts(mytree, nodes)






