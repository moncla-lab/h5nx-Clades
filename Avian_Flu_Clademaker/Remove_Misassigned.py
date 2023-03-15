#File to remove single missassigned leaves

import re
import baltic as bt
import Mutation_Finder as mf
import Mutation_Characterization as mc

tree_path = "Data/flu_avian_h5FullNonMafftGsGDGroupedClock4noDT382-S4220-0813082-005-2.3.4.4_ha.JSON"
exclude_list = []
mytree, mymeta = bt.loadJSON(tree_path)

def remove_misassigned(cleanNodes, rawNodes, tree):
    misassigned = []
    for node in cleanNodes:
        if len(cleanNodes[node]) == 2:
            for item in tree.Objects:
                if item.traits["name"] == node:
                    print(item.leaves)
    return misassigned


allClades, Leaf_Clades  = mf.leaf_clades(mytree, exclude_list)
cleanNodes, rawNodes = mf.node_clades(mytree, Leaf_Clades)
remove_misassigned(cleanNodes, rawNodes, mytree)
