"""
This file takes a tree JSON and pulls out the unique mutations (from both nodes and leaves) for each clade
and prints them as a readable.

An parent_clades list should be initialized to the clades you want to be your highest level clades -- the ones we want mutations to define

The input_tree_path can be changed to whatever you need and the overarching clades should be changed based on how you need them
"""

import re
import baltic as bt

#Input_Tree_Path = "Data/flu_avian_h5nxFullNonMafft-noDT382-S4220-0813082-005_ha_clock_enabled.JSON"
Input_Tree_Path = "Data/flu_avian_h5nxLABEL2-3-4-4Anottated_ha.json"
mytree, mymeta = bt.loadJSON(Input_Tree_Path)
parent_clades = ["Am_nonGsGD", "EA_nonGsGD", "0", "1", "2.1", "2.2", "2.3.1/2", "2.3.3/4", "3", "4", "5", "6", "7", "8", "9", "1-8-9-like", "2-like"]
#output_file_path = "Output/unique_clades_full.txt"
output_file_path = "Output/defining_clades_LABEL.txt"
'''
exclude_list = ["A/falcon/Saudi_Arabia/D1795/2005", 'A/chicken/Giza/CAI37/2009', 'A/duck/China/E319-2/03', 
'A/whooper_swan/Mongolia/21/2010', 'A/Mallard/WI/944/82', 'A/turkey/TX/1-40-82/1982', 'A/chicken/Pekanbaru/BPPVRII-91-3/2011',
 'A/Fujian/1/2007', 'A/Hong_Kong/156/97', 'A/goose/Guangdong/08/2005', 'A/goose/Jilin/hb/2003', 'A/Fujian/1/2007']
'''
exclude_list = []
"""
Here we want to get all the clades that are present in our tree and assign the leaves to them
This function takes the tree as an input variable and returns a list of all clades found in the 
tree along with a dictionary that looks like the following:

LeafClades = {'leaf_name': 'clade', 'leaf_name': 'clade', ...}
"""
def leaf_clades(Tree, exclude):
    cladesAll = [] #the list of all clades present in the tree
    LeafClades = {} #the dictionary which will pair leaves to their clades 

    for item in Tree.Objects: #for each item
        if item.branchType == "leaf": #if the item is a leaf
            if item.traits["name"] not in exclude:
                if "h5_label_clade" in item.traits.keys():
                    clade = item.traits["h5_label_clade"] #variable clade is being set to the leaf's clade
                    if clade != 'UNRECOGNIZABLE':
                        LeafClades[item.traits["name"]] = clade #here we are writing the LeafClades dictionary, setting the leaf name to the key and the clade to the value
                                                            #this will make a dictionary where both the key and the value are strings
                    if clade not in cladesAll and clade != "UNRECOGNIZABLE": #making sure the clade does not already appear in our complete list of clades -- we only want each clade to show up once
                        cladesAll.append(clade) #if it doesn't already show up, add it
    

    return cladesAll, LeafClades


"""
Here we define the categories for our overarching clades, we want to assign each leaf's clade to an overarching value for which we will
determine unique mutations later on.

Since not every individual child clade needs to be defined as its own family, we can just define characteristics for overall families.

the input is the list of clades we determined from our leaf_clades() function along with a key of clades which we want to assign as our parent clades.

Returns a dictionary that looks like the following:
hierarchical_clades = {'parent_clade': ['child_clade', 'child_clade', 'child_clade'...], 'parent_clade': ['child_clade', 'child_clade', 'child_clade'...], ...}
"""
def clade_heirarchy(all_clades, parentClades):
    hierarchical_clades = {} #initalizing the dictionary which we will return later
    for parent in parentClades: #for each parent clade, we're going to set a key to which we will later assign a list of all children clades
        hierarchical_clades[parent] = [] 

    for clade in all_clades: #for each clade in our list of all clades
        if len(clade) == 1 or "non" in clade: #if our clade is any of the non_GsGD clades, or is just a single number then that clade is equivalent to its parent clade
            hierarchical_clades[clade].append(clade) #the reason we are going to have a key with a value which is the same as the key is for later, we are going to have to pull keys from values and 
                                                    #having some keys which have no values would make that very difficult


        elif clade[0] in "13456789": #if the first character of our clade is in this list, it will just be appended to the family that its first character denotes
            hierarchical_clades[clade[0]].append(clade) 
        elif "2.3.1" in clade or "2.3.2" in clade: #here we assign all 2.3.1 and 2.3.2 clades to the 2.3.1/2 overarching family
            hierarchical_clades["2.3.1/2"].append(clade)
        elif "2.1" in clade[0:3]:
            hierarchical_clades["2.1"].append(clade) #here we assign all 2.1 clades to the 2.1 clade
        elif "2.2" in clade[0:3]:
            hierarchical_clades["2.2"].append(clade) #here we assign all the 2.2 clades to the 2.2 clade
        elif "2.3.3" in clade or "2.3.4" in clade: #here we define all the 2.3.3 and 2.3.4 clades to the 2.3.3/4 parent clade
            hierarchical_clades["2.3.3/4"].append(clade)

    return hierarchical_clades #here we are going to return our heirarchical list of clades, which contains the parent clade as the key with all children clades as a list of values


"""
Here we want to pull out all the nucleotide mutations for each leaf and assign it to the parent clade of that leaf. 
This is going to return a nested dictionary that looks like the following:

mutations = {'parent_clade': [['mutation', 'mutation'...], ['mutation', 'mutation', 'mutation'...]], ...}
The dictionary values will be lists of lists, where the parent list contains lists of mutations from each leaf for that clade
"""
def leaf_nuc_muts(hierarchical_list, treeJSON):
    mutations = {} #initialize the mutation dictionary so we can call it later
    for k in treeJSON.Objects: #here we're going to go through every leaf and do the following code for each leaf we find
        skip = False
        if k.branchType == "leaf": #find all leaves
            clade = k.traits["h5_label_clade"] #set the clade as the leaf's clade
            if "like" in clade:
                skip = True
            mut_dict = k.traits["branch_attrs"]["mutations"] #grab the mutations here as a shorter variable so we don't have to write the long one over and over
            if 'nuc' in mut_dict.keys(): #checking if there are nucleotide mutations
                nucleotide_muts = mut_dict['nuc']

            for key, value in hierarchical_list.items(): #here we need to grab the parent clade given the child clade so we are going to grab the key and values from our hierarchical dict.
                if clade in value and skip == False:
                    mutation_clade = key #setting our clade as the parent clade

            if mutation_clade not in mutations.keys() and skip == False: #here we need to check if the clade is already in the dict, if isn't we can just set the mutations to the clade
                mutations[mutation_clade] = [nucleotide_muts]
            
            elif skip == False: #if the clade is already present, we don't want to overwrite the mutations we already have, so we're going to append it
                mutations[mutation_clade].append(nucleotide_muts)

    return mutations #here we return the dictionary of mutations


"""
This function is going to give each node a list of clades derived from its children leaves. This means every node will have the clades of its children.
Later we can use this to determine which node is the defining node for each clade.

Here we are assigning the nodes the children clades, rather than the parents clades. The reason we do this is because we can easily assign parent clades later when we need them,
but it would be much more difficult to assign child clades given parents, and we would essentially have to run this code again.

The output for this function looks like this:
clean_nodes = {'node_name': ['clade', 'clade', 'clade', ...], ...}
It will only return 1 of each clade found in that node's children. No duplicates.
"""
def node_clades(treeJSON, leaf_clades):
    node_clade = {} #initializing the variable, just like usual
    for k in treeJSON.Objects: 
        if k.branchType == "node": #here we are going to run through each node
            name = k.traits["name"] #grab the node name as a string
            node_clade[name] = [] #initialize the dictionary
            for leaf in k.leaves: #here we are going to go through each leaf that is a child of our node
                if leaf in leaf_clades.keys():
                    node_clade[name].append(leaf_clades[leaf]) #here we are going to grab the clade that is assigned to that leaf and assign it to the node *** this will assign duplicates ***


            clean_nodes = {} #here we are going to remove duplicate clades because it will make things prettier
            for node in node_clade: #for each node
                clean_nodes[node] = [] #initialize the new dictionary with the node
                for i in range(len(node_clade[node])): #here we are going to parse through each clade (we are using a length index here because not all clades have the same length)
                    clade = node_clade[node][i] #setting the clade to the clade we just found
                    if clade not in clean_nodes[node]: #here we are only going to add the clade to our dictionary if we don't already have it there
                        clean_nodes[node].append(clade)

    return clean_nodes, node_clade


"""
Here we are going to output all of the nucleotide mutations found in our nodes, 
given our hierarchical list of clades, our tree, and our list of clades that we assigned to our nodes earlier

this function will output a dictionary that looks like this:
node_nMuts = {'parent_clade': [[mutation, mutation, ...], ['mutation', 'mutation', ...]]}
"""
def node_nuc_muts(hList, Tree, node_clades):
    node_nMuts = {} #initialize our dictionary
    for object in Tree.Objects:
        if object.branchType == "node": #for each node in our tree
            mutations_full = object.traits["branch_attrs"]["mutations"] #set a variable to this directory so we don't have to type it out again
            name = object.traits["name"] #grab the node name
            clade = node_clades[name] #set the clade given our name -- this is the child clade and we will have to get the parent clade later on
            if "nuc" in mutations_full: #if we have nucleotide mutations
                for item in clade: #for clade in our list of clades

                    #because the following is nested in our 'for' loop above, it will all repeat for each clade in our list of clades
                    for key, value in hList.items(): #set the clade to the parent clade
                        if item in value:
                            fClade = key

                    if fClade not in node_nMuts.keys(): #if we don't already have that clade, add it and the mutations
                        node_nMuts[fClade] = object.traits["branch_attrs"]["mutations"]["nuc"]
                
                    else: #if we do already have that clade, simply append the mutation list we already have
                        node_nMuts[fClade].append(object.traits["branch_attrs"]["mutations"]["nuc"])


    return node_nMuts


"""
Here we are going to take the two lists we prevously generation, the leaf nucleotide mutation dict and the node nucleotide mutation dict

The dictionary we are going to output looks like this:
all_nMuts = {'clade': ['mutation', 'mutation', ...], 'clade': ...}

Here we will have no duplicates the the values for the dictionary will be a single list, not a list of lists like we had previously. 
The clade keys will be the parent clades
"""
def all_nuc_muts(nodes, leaves):
    all_nMuts = {} #initializing our dictionary -- like always

    #we're going to go through leaves first and then will go through nodes
    for key in leaves:
        if key not in all_nMuts.keys(): #we're going to add the clades to the dictionary if they are no already there
            all_nMuts[key] = []

        for item in range(len(leaves[key])): #here we are going to go through our clade mutations for the clade we pulled out in our first for loop
            if isinstance(leaves[key][item], list): #if our first key is a list, we are going to have to go through it again
                for string in range(len(leaves[key][item])):
                    if leaves[key][item][string] not in all_nMuts[key]: #if our mutation is already present, we don't want to add it again
                        all_nMuts[key].append(leaves[key][item][string])
            else: #if our index is not a list, we can just add things to our mutation list
                if leaves[key][item] not in all_nMuts[key]: # if it's already in our clade's mutation list, we don't want it and will skip it, otherwise, we'll append it to our list
                    all_nMuts[key].append(leaves[key][item])
        
#now we have a dictionary where the keys are the parent clades and the values are a single list of mutations
#we're now going to go through this same code for the nodes

    #nodes pass
    for key_2 in nodes:
        if key_2 not in all_nMuts.keys():
            all_nMuts[key_2] = [] #adding all node clades to the dictionary
        for index in range(len(nodes[key_2])): 
            if isinstance(nodes[key_2][index], list): #checking if we're in a list or not
                for string in range(len(nodes[key_2][index])): #if we're a list, we're going to go through each item in the list seperately
                    if nodes[key_2][index][string] not in all_nMuts[key_2]:
                        all_nMuts[key_2].append(nodes[key_2][index][string])
            else: #if we're not a list we can just add all the mutations
                if nodes[key_2][index] not in all_nMuts[key_2]:
                    all_nMuts[key_2].append(nodes[key_2][index])

    return all_nMuts


"""
Here we want to pull out all the HA mutations for each leaf and assign it to the parent clade of that leaf. 
This is going to return a nested dictionary that looks like the following:

mutations = {'parent_clade': [['mutation', 'mutation'...], ['mutation', 'mutation', 'mutation'...]], ...}
The dictionary values will be lists of lists, where the parent list contains lists of mutations from each leaf for that clade
"""
def leaf_HA_muts(hierarchical_list, treeJSON):
    mutations = {}

    for k in treeJSON.Objects:
        if k.branchType == "leaf": #running through each leaf
            clade = k.traits["h5_label_clade"] #grabbing its clade (child, not parent)
            mut_dict = k.traits["branch_attrs"]["mutations"] #grabbing its mutation trait
            if "HA" in mut_dict.keys(): # if it has HA mutations, we're gonna set our path to that so we don't have to type it over and over again
                HA_muts = mut_dict["HA"]

                mutation_clade = "" #initializing our mutation_clade variable -- this isn't super necessary but it makes me feel better
                for key, value in hierarchical_list.items(): #get the parent clade given the child clade we pulled out earlier
                    if clade in value:
                        mutation_clade = key
                
                if mutation_clade not in mutations.keys(): #here we're going to add our HA mutations to a dictionary where the parent clade is the key
                    mutations[mutation_clade] = [HA_muts] #if there isn't already a key for our clade we add it
                else:
                    mutations[mutation_clade].append(HA_muts) #if there is already a key we just append our new mutation list

    return mutations


"""
Here we are going to output all of the HA mutations found in our nodes, 
given our hierarchical list of clades, our tree, and our list of clades that we assigned to our nodes earlier

this function will output a dictionary that looks like this:
node_nMuts = {'parent_clade': [[mutation, mutation, ...], ['mutation', 'mutation', ...]]}
"""
def node_HA_muts(hList, Tree, node_clades):
    node_HAMuts = {} #initialize our dictionary
    for object in Tree.Objects:
        if object.branchType == "node": #run through each node
            mutations_full = object.traits["branch_attrs"]["mutations"]
            name = object.traits["name"] #grab the name of the node
            clade = node_clades[name] #grab the child clade
            if "HA" in mutations_full: #if we have HA mutations
                for item in clade:
                    for key, value in hList.items():
                        if item in value:
                            fClade = key #grabbing the parent clade and replacing our child clade with that

                    if fClade not in node_HAMuts.keys(): #adding the mutations if we don't have the key
                        node_HAMuts[fClade] = object.traits["branch_attrs"]["mutations"]["HA"]
                
                    else: #appending the mutations if we do have that key already
                        node_HAMuts[fClade].append(object.traits["branch_attrs"]["mutations"]["HA"])

    return node_HAMuts


"""
Here we are going to take the two lists we prevously generation, the leaf HA mutation dict and the node HA mutation dict

The dictionary we are going to output looks like this:
all_HAmuts = {'clade': ['mutation', 'mutation', ...], 'clade': ...}

Here we will have no duplicates the the values for the dictionary will be a single list, not a list of lists like we had previously. 
The clade keys will be the parent clades
"""
def all_HA_muts(nodes, leaves):
    all_HAMuts = {} #initialize our dictionary to be returned later

    for key in leaves: #finding the clade in our leaves list
        if key not in all_HAMuts.keys():
            all_HAMuts[key] = [] #if we don't already have that clade, add it
        for item in range(len(leaves[key])): #running through each of the values for our clade
            if isinstance(leaves[key][item], list): #if our value is a list
                for string in range(len(leaves[key][item])): #run through each item in the list
                    if leaves[key][item][string] not in all_HAMuts[key]: #if our mutation isn't already there, we add it
                        all_HAMuts[key].append(leaves[key][item][string])
            else: #if our value is not a list
                if leaves[key][item] not in all_HAMuts[key]: #we add the mutation if it's not already there
                    all_HAMuts[key].append(leaves[key][item])
        
    #the following code is the exact same code as for the leaves
    for key_2 in nodes:
        if key_2 not in all_HAMuts.keys():
            all_HAMuts[key_2] = []
        for index in range(len(nodes[key_2])):
            if isinstance(nodes[key_2][index], list):
                for string in range(len(nodes[key_2][index])):
                    if nodes[key_2][index][string] not in all_HAMuts[key_2]:
                        all_HAMuts[key_2].append(nodes[key_2][index][string])

            else:
                if nodes[key_2][index] not in all_HAMuts[key_2]:
                    all_HAMuts[key_2].append(nodes[key_2][index])

    return all_HAMuts #return our dictionary of all HA mutations


"""
Here we are going to define all the mutations which are unique to a clade:
    in other words, only mutations which show up in that clade, and that clade only,
    will be assigned to the dictionary we return. If a mutation shows up in any other clade, we will remove it

Takes two dictionaries, the dictionary of all HA mutations and the dict of all nuc mutations

Returns two dictionaries that look like the following:

unique_HA = {'clade': ['mutation', 'mutation',...],...}
unique_nuc = {'clade': ['mutation', 'mutation',...],...}
"""
def unique_muts(HA, nuc):
    #initializing our dictionaries to be returned
    unique_HA = {}
    unique_nuc = {}

    for key in HA.keys(): #because we should only have 1 instance of each clade
        unique_HA[key] = [] #we can simply initalize our list at the start

        for mutation in HA[key]: #for each mutation in HA mutations
            skip = False #initialize our skip variable for later
            for clade in HA.keys(): #for clade in our keys
                if clade != key: # as long as our clade is not the key we are appending to
                        skip = True #if it is present we're going to skip it
            if "-" in mutation or "x" in mutation or "X" in mutation or "*" in mutation:
                skip = True
            if skip == False: #then, if skip is False, we're going to append our mutation to the list
                unique_HA[key].append(mutation)
    #the way the skip variable works is for each mutation it is by default
    #assigned to false, then if it is found in any other clade it is assigned to True
    #because we assign to false before we start going in to the check portion of the code
    #it will only be assigned false by default when we first find the mutation, and no other time
    
    #the following is the exact same code as above but for nuc rather than HA
    for key in nuc.keys():
        unique_nuc[key] = []

        for mutation in nuc[key]:
            skip = False
            for clade in nuc.keys():
                if clade != key:
                    if mutation in nuc[clade]:
                        skip = True

            if skip == False:
                unique_nuc[key].append(mutation)

    return unique_HA, unique_nuc

"""
Here we find all the nodes which define our clades:
here is how we are defining them:
We want the first node in the tree (furthest to the left) which has leaves
that are all the same clade. This means any leaf descending from that node should be 
from the same clade. 

To define this, we are going to look at the list of clades our nodes were defined earlier
and find all the nodes which only have a single clade defined. Once we have nodes with only
single clade leaves, we can take the node that is furthest to the left (which should have the lowest number)

This function will return a dictionary where the clade is the key and the single defining node is the value
"""
"""
def node_finder(nodeClades, hierarchical_list):
    defining_nodes = {} #initialize our variable
    for node in nodeClades.keys(): #for each node
        if len(nodeClades[node]) == 1: #if we have a single clade
            clade = nodeClades[node][0] #we are going to set our clade variable to the clade we just found
            for key, value in hierarchical_list.items(): #then, we are going to get the parent clade because we currently have the child clade
                if clade in value:
                    new_clade = key

            if new_clade in defining_nodes.keys(): #then, if we already have our clade in the dictionary
                checker = node.split("_") #we're gonna split the node so we can get the number on the end
                checker_int = re.search("^0+(\d+)", checker[1]).group(1) #we're gonna grab the number without the leading zeros
                checker_int = int(checker_int)

                current = defining_nodes[new_clade] #grab the number for the clade without the leading zeros
                current = current.split("_")
                current_int = re.search("^0+(\d+)", checker[1]).group(1)
                current_int = int(current_int)

                if checker_int < current_int: #if our new node has a lower number than our current node, we replace the current node with our new node
                    defining_nodes[new_clade] = node

            else:
                defining_nodes[new_clade] = node #if we don't already have our clade in our dict, we're going to add it with our node as the value.

    return defining_nodes
"""
"""
#get leaf clades:
allClades, leafClades = leaf_clades(mytree)

#make hierarchical list
hierarchy = clade_heirarchy(allClades, parent_clades)

#get the clades for each node
nodeClades = node_clades(mytree, leafClades)

#get nucleotide mutations for each leaf
lnm = leaf_nuc_muts(hierarchy, mytree)

#get nucleotide mutations for each node
nnm = node_nuc_muts(hierarchy, mytree, nodeClades)

#combine all nucleotide mutations into a single dictionary with the clade as the key
allNuc = all_nuc_muts(lnm, nnm)

#get all the leaf HA mutations
lHAm = leaf_HA_muts(hierarchy, mytree)

#get all node HA mutations
nHAm = node_HA_muts(hierarchy, mytree, nodeClades)

#combine leaf and node HA mutations into a single dictionary with the clade as the key
allHA = all_HA_muts(nHAm, lHAm)

#get just the unique mutations
HA, Nuc = unique_muts(allHA, allNuc)


#find the nodes which define each clade
Nodes = node_finder(nodeClades, hierarchy)

#print the results we got from unique_muts in a pretty way
pretty_print(HA, Nuc, Nodes, output_file_path)

"""

#allClades, leafClades = leaf_clades(mytree, exclude_list)