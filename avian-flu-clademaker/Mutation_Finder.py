"""
This file takes a tree JSON and pulls out the unique mutations (from both nodes and leaves) for each clade
and prints them as a readable.

An parent_clades list should be initialized to the clades you want to be your highest level clades -- the ones we want mutations to define

The input_tree_path can be changed to whatever you need and the overarching clades should be changed based on how you need them
"""

import re
import baltic as bt

Input_Tree_Path = "Data/flu_avian_h5nxLABEL2-3-4-4Anottated_ha.json"
mytree, mymeta = bt.loadJSON(Input_Tree_Path)
parent_clades = ["Am_nonGsGD", "EA_nonGsGD", "0", "1", "2.1", "2.2", "2.3.1/2", "2.3.3/4", "3", "4", "5", "6", "7", "8", "9", "1-8-9-like", "2-like"]
output_file_path = "Output/defining_clades_LABEL.txt"
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
                    if mutation in HA[clade]:
                        skip = True #if it is present we're going to skip it
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