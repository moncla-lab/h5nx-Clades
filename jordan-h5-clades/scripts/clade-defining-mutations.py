# import packages
import json
import importlib.util
import pandas as pd
import argparse


parser = argparse.ArgumentParser()

parser.add_argument('--baltic', type=str, help='path to baltic.py (https://github.com/evogytis/baltic)', required=True)
parser.add_argument('--tree', type=str, help='path to tree json file', required=True)
parser.add_argument('--clade_mutations', type=str, help='path to output clade defining mutations tsv file', required=True)
parser.add_argument('--excluded_tips', type=str, help='optional path to txt file containing any tips to exclude when determining clade defining mutations', required=False)
parser.add_argument('--outgroup_clades', type=str, help='optional path to txt file containing any outgroup clades to ignore when determining clade defining mutations', required=False)
parser.add_argument('--init_mutations', type=str, help='optional path to tsv file containing any mutations to manually include, preformatted as required for augur clades with headers', required=False)

args = parser.parse_args()

baltic = args.baltic
tree = args.tree
clade_mutations = args.clade_mutations
excluded_tips = args.excluded_tips
outgroup_clades = args.outgroup_clades
init_mutations = args.init_mutations

# if --excluded_tips arg provided, convert to list
if excluded_tips:
    with open(excluded_tips, 'r') as f:
        excluded_tips_list = f.read().split('\n')
else:
    excluded_tips_list = []

# if --outgroup_clades arg provided, convert to list
if outgroup_clades:
    with open(outgroup_clades, 'r') as f:
        outgroup_clades_list = f.read().split('\n')
else:
    outgroup_clades_list = []

# if --init_mutations arg provided, open as df
if init_mutations:
    init_df = pd.read_csv(init_mutations, sep='\t', header=None)
else:
    init_df = None





# define functions

def load_module(name, path):
    module_spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(module_spec)
    module_spec.loader.exec_module(module)
    return module

def bt_read_in_tree_json(input_tree):
    '''read in a tree in json format'''
    
    with open(input_tree) as json_file:
        json_tree = json.load(json_file)
    # json_translation = {'absoluteTime':lambda k: k.traits['node_attrs']['num_date']['value'],'name':'name'} ## allows baltic to find correct attributes in JSON, height and name are required at a minimum
    json_translation = {'length':lambda k: k.traits['node_attrs']['div'],'name':'name'} ## allows baltic to find correct attributes in JSON, height and name are required at a minimum
    bt_tree, meta = bt.loadJSON(json_tree, json_translation)
    return json_tree, bt_tree

def get_clade_inheritance(tree, clade_relationships = None, parent = None):
    '''return a dictionary of clade parent-child relationships, stored as {child: parent}'''
    
    if clade_relationships is None:
        clade_relationships = {}
        
    if tree.branchType == 'node':
        if tree in bt_lcas.values():
            clade = [cl for cl in bt_lcas.keys() if bt_lcas[cl] == tree][0]
            clade_relationships[clade] = parent
            parent = clade
        children = tree.children
        for child in children:
            get_clade_inheritance(child, clade_relationships, parent)
    return clade_relationships


# get parent-child clade relationships
def make_tree_from_relationships(clade_relationships_dict):
    '''return a tree (as nested dicts) of clade relationships'''
    
    # convert dictionary into a list of (child, parent) tuples
    relationships_list = [(k, v) for k, v in clade_relationships_dict.items()]
    
    # convert list of tuples into a nested dict tree
    roots = set()
    mapping = {}
    for child, parent in relationships_list:
        childitem = mapping.get(child, None)
        if childitem is None:
            childitem = {}
            mapping[child] = childitem
        else:
            roots.discard(child)
        parentitem = mapping.get(parent, None)
        if parentitem is None:
            mapping[parent] = {child:childitem}
            roots.add(parent)
        else:
            parentitem[child] = childitem
    
    # return the nested dict tree, excluding the basal None 'parent'
    return {i:mapping[i] for i in sorted(roots)}[None]

            
def get_clade_heights(tree, clade_heights = None, i = None):
    '''return a dictionary of clade "heights", where the most basal clade has a
    height of 0, its child clades a height of 1, grandchild clades 2, etc.'''
    
    clade_heights = clade_heights or {}
    i = i or 0
    for root in tree:
        if i in clade_heights:
            clade_heights[i].append(root)
        else:
            clade_heights[i] = [root]
        if tree[root] != {}:
            get_clade_heights(tree[root], clade_heights, i+1)
    return clade_heights


def return_nt_muts_on_branch(branch):
    '''return list of nucleotide mutations on a given branch'''
    
    nt_muts = []
    if 'branch_attrs' in branch.traits:
        if 'mutations' in branch.traits['branch_attrs']:
            if 'nuc' in branch.traits['branch_attrs']['mutations']:
                nt_muts = branch.traits['branch_attrs']['mutations']['nuc']
    return nt_muts


def return_ha_muts_on_branch(branch):
    '''return list of ha amino acid mutations on a given branch'''
    
    ha_muts = []
    if 'branch_attrs' in branch.traits:
        if 'mutations' in branch.traits['branch_attrs']:
            if 'HA' in branch.traits['branch_attrs']['mutations']:
                ha_muts = branch.traits['branch_attrs']['mutations']['HA'] 
    return ha_muts


def return_all_muts_between_lcas(starting_node, ending_node, muts = None, i = None):
    '''return dictionary containing lists of all nt and ha mutations that occur between two nodes'''
    
    muts = muts or {'nuc': [], 'ha': []}
    i = i or 0
    
    # all leaves that descend from the ending node
    # will use these to determine if the ending node descends from the node currently being analyzed
    end_leaves = ending_node.leaves
    
    # set an empty list of mutations and enumerate the children of the starting node; children can be tips or nodes
    children = starting_node.children
    
    for child in children:        
        """if the child is a leaf:
            move on, too far down the tree (we want a node, not terminal tip)"""
        """if the child is an internal node:
            first, test whether that child node contains the target tips in its children.
            child.leaves will output a list of the names of all tips descending from that node. if not, pass. 
            if the node does contain the target end tip in its leaves, keep traversing down that node recursively, 
            collecting mutations as you go"""

        if child.branchType == "leaf":
            # if child is a leaf, we went too far
            pass
        
        elif child.branchType == "node":
            # if child is a node, check if it's the ending node
            ## if it is, add the branch muts and return the list of muts
            ## if it's not the ending node, check if the ending node is in its children
            ### if so, add muts and keep going
            ### if not, we are not on the path to the ending node. pass
            if child == ending_node:
                # found ending node
                nt_muts = return_nt_muts_on_branch(child)
                muts['nuc'].append(nt_muts)
                ha_muts = return_ha_muts_on_branch(child)
                muts['ha'].append(ha_muts)
                
            elif all(leaf in child.leaves for leaf in end_leaves):
                # found node on path to ending node
                nt_muts = return_nt_muts_on_branch(child)
                muts['nuc'].append(nt_muts)
                ha_muts = return_ha_muts_on_branch(child)
                muts['ha'].append(ha_muts)
                
                # continue iterating down the path
                return_all_muts_between_lcas(child, ending_node, muts, i)
                
            else:
                # node not on correct path
                pass
    
    # flatten the list so that you don't have nested lists
    nt_flat_list = [item for sublist in muts['nuc'] for item in sublist]
    ha_flat_list = [item for sublist in muts['ha'] for item in sublist]
    flat_list = [nt_flat_list, ha_flat_list]
    return flat_list


# load baltic
bt = load_module('bt', baltic)





## get clade defining mutations for all clades, excluding -like clades

## exclude any tips which appear to be mislabled (e.g., the tip clusters with a different clade)
## doing so will ensure that the called LCA for the clade will be the LCA of all sequences which are clustering well

## need to add in a better filter for polybasic cleavage site mutations
    ## currently, just filtering based on GsGd polybasic cleavage site range
    ## however, 2.3.2.1 and 2.3.4.4 trees use a different reference sequence
    ## could pass in reference sequence, find the cleavage site, and assign filter to the appropriate aa/nt range

json_tree, bt_tree = bt_read_in_tree_json(tree)

# make dictionary with each clade as a key and a list of all leaves (as baltic objects) of that clade as values
clade_bt_leaves = {}
for leaf in bt_tree.getExternal():
    name = leaf.traits['name']
    if name in excluded_tips_list:
        print('excluding', name)
    elif 'clade' in leaf.traits:
        clade = leaf.traits['clade']
        if clade not in outgroup_clades_list:
            if not clade in clade_bt_leaves:
                clade_bt_leaves[clade]  = [leaf]
            else:
                clade_bt_leaves[clade].append(leaf)
        # previously ignored any '-like' clades and assigned manually, but think it's best to include automatically
#         if '-like' not in clade:
#             if not clade in clade_bt_leaves:
#                 clade_bt_leaves[clade]  = [leaf]
#             else:
#                 clade_bt_leaves[clade].append(leaf)

# find the branch that is the last common ancestor (lca) to all leaves in each clade
# store it in a clade:lca dict
bt_lcas = {}  
for clade, leaf_list in clade_bt_leaves.items():
    try:
        ancestor = bt_tree.commonAncestor(leaf_list)
        # sometimes baltic commonAncestor returns an 'empty' node when the true lca should be the root of the tree
        # check if the ancestor is empty (traits is empty dict {}) — if so, assign the tree root as lca instead
        if ancestor.traits != {}:
            bt_lcas[clade] = ancestor
        else:
            bt_lcas[clade] = bt_tree.root
    except AssertionError:
        # baltic commonAncestor throws error if there are < 2 tips in list
        # if so, can't find an lca — pull mutations from the branch instead
        print('Clade', clade, 'has only one member and thus no LCA can be found')
        bt_lcas[clade] = leaf_list[0]

# pull the mutations from each lca branch
# store them in a clade:mutations dict
bt_lcas_mutations = {}
for clade, lca in bt_lcas.items():
    try:
        bt_lcas_mutations[clade] = lca.traits['branch_attrs']['mutations']
    except KeyError:
        print('LCA of clade', clade, 'has no mutations on the branch and will not be included')
        

# determine parent-child clade relationships
clade_relationships = {}
clade_relationships = get_clade_inheritance(bt_tree.root, clade_relationships)
clade_relationships_tree = make_tree_from_relationships(clade_relationships)
clade_heights = get_clade_heights(clade_relationships_tree)

# dict of mutations that occur between the LCA of a parent clade and the LCA of a clade of interest
# e.g., lca_to_lca_muts['2.3.2.1d'] = [[nucmut1, nucmut2, ...], [hamut1, hamut2, ...]]
## where nucmut1, nucmut2, ... are all the nt mutations that occurred between 2.3.2.1c LCA and 2.3.2.1d LCA
## and hamut1, hamut2, ... are the ha mutations that occurred between the LCAs
lca_to_lca_muts = {}
lca_to_lca_muts_positions = {}

for i in [j for j in clade_heights if j > 0]:
    for clade in clade_heights[i]:
        parent = clade_relationships[clade]
        if clade in bt_lcas and parent in bt_lcas:
            lca1 = bt_lcas[parent]
            lca2 = bt_lcas[clade]
            lca_to_lca_muts[clade] = return_all_muts_between_lcas(lca1, lca2)
            nt_muts = lca_to_lca_muts[clade][0]
            ha_muts = lca_to_lca_muts[clade][1]
            nt_muts_positions = [mut[1:-1] for mut in nt_muts]
            ha_muts_positions = [mut[1:-1] for mut in ha_muts]
            
            # if site toggles multiple times between lcas, need to remove all but the final mutation
            # for instance, if parent lca has nt 100A and child lca has 100T, but 100C is an intermediate,
            # then we will see both '100C' and '100T' as lca-to-lca muts but we only want the final 100T
            # thus, need to remove the first n-1 muts of n muts with identical position values
            rep_nt_muts_positions = set(pos for pos in nt_muts_positions if nt_muts_positions.count(pos) > 1)
            rep_ha_muts_positions = set(pos for pos in ha_muts_positions if ha_muts_positions.count(pos) > 1)

            for nt_pos in rep_nt_muts_positions:
                rep_nt_muts = [mut for mut in nt_muts if mut[1:-1] == nt_pos]
                for mut in rep_nt_muts[:-1]:
                    nt_muts.remove(mut)
            for ha_pos in rep_ha_muts_positions:
                rep_ha_muts = [mut for mut in ha_muts if mut[1:-1] == ha_pos]
                for mut in rep_ha_muts[:-1]:
                    ha_muts.remove(mut)
            
            for nt_mut in [nt_mut for nt_mut in nt_muts if nt_mut in bt_lcas_mutations[clade]['nuc']]:
                nt_muts.remove(nt_mut)
                
            if 'HA' in bt_lcas_mutations[clade]:
                for ha_mut in [ha_mut for ha_mut in ha_muts if ha_mut in bt_lcas_mutations[clade]['HA']]:
                    ha_muts.remove(ha_mut)
            
            nt_muts_positions = [mut[1:-1] for mut in nt_muts]
            ha_muts_positions = [mut[1:-1] for mut in ha_muts]
            
            lca_to_lca_muts_positions[clade] = [nt_muts_positions, ha_muts_positions]
            

# convert the clade:mutations dict to a df in the correct format for augur clades
clade_data = ['clade'] # clade that is being defined
muttype_data = ['gene'] # type of mutation (nuc or HA)
mutsite_data = ['site'] # site of mutation (bp or aa position)
mut_data = ['alt'] # alternative 
unique_data = ['unique'] # unique mutation (not propagated from parental clade)


#for i in [j for j in clade_heights if j > 0]:
for i in clade_heights:
    for clade in clade_heights[i]:
        # get mutations from dict and append data to lists as appropriate
        mutations = bt_lcas_mutations[clade]        
        if clade_relationships[clade] != None and clade_relationships[clade] != '0':
            clade_data.append(clade)
            muttype_data.append('clade')
            mutsite_data.append(clade_relationships[clade])
            mut_data.append('')
            unique_data.append(True)
        if 'nuc' in mutations:
            for mutation in [mut for mut in mutations['nuc'] if 'N' not in mut and '-' not in mut]:
                # if mutation site is between 1036 and 1059 (inclusive) it is part of polybasic cleavage site -- ignore
                if int(mutation[1:-1]) not in range(1036, 1060):
                    clade_data.append(clade)
                    muttype_data.append('nuc')
                    mutsite_data.append(mutation[1:-1])
                    mut_data.append(mutation[-1])
                    unique_data.append(True)
        if 'HA' in mutations:
            for mutation in [mut for mut in mutations['HA'] if 'X' not in mut and '-' not in mut]:
                # if mutation site is between 339 and 346 (inclusive) it is part of polybasic cleavage site -- ignore
                if int(mutation[1:-1]) not in range(339, 347):
                    clade_data.append(clade)
                    muttype_data.append('HA')
                    mutsite_data.append(mutation[1:-1])
                    mut_data.append(mutation[-1])
                    unique_data.append(True)

        # explicitly add mutations propagated from parent clade to lists
        parent = clade_relationships[clade]
        while parent != None:            
            parent_mutations = bt_lcas_mutations[parent]

            if 'nuc' in parent_mutations:
                for mutation in [mut for mut in parent_mutations['nuc'] if 'N' not in mut and '-' not in mut]:
                    # if mutation site is between 1036 and 1059 (inclusive) it is part of polybasic cleavage site -- ignore
                    if int(mutation[1:-1]) not in range(1036, 1060):
                        if mutation[1:-1] not in lca_to_lca_muts_positions[clade][0]:
                            clade_data.append(clade)
                            muttype_data.append('nuc')
                            mutsite_data.append(mutation[1:-1])
                            mut_data.append(mutation[-1])
                            unique_data.append(False)
                        else:
                            corr_mutation = [mut for mut in lca_to_lca_muts[clade][0] if mut[1:-1]==mutation[1:-1]][0]
                            clade_data.append(clade)
                            muttype_data.append('nuc')
                            mutsite_data.append(corr_mutation[1:-1])
                            mut_data.append(corr_mutation[-1])
                            unique_data.append(True)
            if 'HA' in parent_mutations:
                for mutation in [mut for mut in parent_mutations['HA'] if 'X' not in mut and '-' not in mut]:
                    # if mutation site is between 339 and 346 (inclusive) it is part of polybasic cleavage site -- ignore
                    if int(mutation[1:-1]) not in range(339, 347):
                        if mutation[1:-1] not in lca_to_lca_muts_positions[clade][1]:
                            clade_data.append(clade)
                            muttype_data.append('HA')
                            mutsite_data.append(mutation[1:-1])
                            mut_data.append(mutation[-1])
                            unique_data.append(False)
                        else:
                            corr_mutation = [mut for mut in lca_to_lca_muts[clade][1] if mut[1:-1]==mutation[1:-1]][0]
                            clade_data.append(clade)
                            muttype_data.append('HA')
                            mutsite_data.append(corr_mutation[1:-1])
                            mut_data.append(corr_mutation[-1])
                            unique_data.append(True)
            parent = clade_relationships[parent]






df = pd.DataFrame(list(zip(clade_data, muttype_data, mutsite_data, mut_data, unique_data)))
df = df[(df[4]==True) | (df[4]=='unique')].drop(columns=4)

if init_df is not None:
    df = pd.concat([init_df, df.iloc[1:]])

df.to_csv(clade_mutations, sep="\t", index=False, header=False)