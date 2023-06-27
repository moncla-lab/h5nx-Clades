"""
This file is intended to run on its own, it can be run from the terminal or used here simply by running it. The current settings should be
good to run with a JSON file assuming you have the file named correctly and in the right location

Data should be kept in the Data file and all files will be output to the Output file
"""


import re
import baltic as bt
import argparse
parser = argparse.ArgumentParser()

import Mutation_Finder as mf

parser.add_argument('--input_tree', type=str, help='path to input tree in JSON format')
parser.add_argument('--output_file', type=str, help='name and path to output file')
parser.add_argument('--exclude_file', type=str, help='file including strains to exclude')
parser.add_argument('--metadata_clade_column_name', type=str, help='the title of the metadata column containing the clade assignment')

args = parser.parse_args()
Input_Tree_Path = args.input_tree
output_file_path = args.output_file
exclude_file = args.exclude_file
clade_column_name = args.metadata_clade_column_name


mytree, mymeta = bt.loadJSON(Input_Tree_Path)
exclude_list = ['A/breeder_duck/Korea/H158/2014'] #this strain is behaving weirdly in our tree, so this code is just going to have us ignore it
                                                #if any other strains are acting weird, just ad them here in quotations and comma separated
parent_list = ['2.1.3.2', '1', '1.1',
'2.2', '2.2.1', '2.2.1.1', '2.3.2.1', '2.3.4', '7',
'2.2.2', '2.3.2', '6', '2.1.3.1', '2.1.2', '2.3.4.4']

"""
This dictionary is important because when we are outputting our tsv file we need our clades to be able to inherit from other clades
so this dictionary allows us to quickly output which clades descend from what without having to manually write it every time. If you find
that your clades are related differently, simply change this code. The format is child:parent, where the parent is the clade which our child
inherits from
"""
relationships = {
    '0': None,
    '1': None,
    '1.1': '1',
    '1.1.1': '1.1',
    '1.1.2': '1.1',
    '2.4': None,
    '2.5': None,
    '2.1.2': '2.1.1',
    '2.1.3': '2.1.2',
    '2.1.1': None,
    '2.1.3.1': '2.1.2',
    '2.1.3.2': '2.1.2',
    '2.1.3.2a': '2.1.3.2',
    '2.1.3.2b': '2.1.3.2',
    '2.1.3.3': '2.1.2',
    '2.2': None,
    '2.2.2': '2.2',
    '2.2.1': '2.2',
    '2.2.1.1': '2.2.1',
    '2.2.1.1a': '2.2.1.1',
    '2.2.1.2': '2.2.1',
    '2.2.2.1': '2.2',
    '2.3.1': None,
    '2.3.2': None,
    '2.3.2.1': '2.3.2',
    '2.3.2.1a': '2.3.2.1',
    '2.3.2.1b': '2.3.2.1',
    '2.3.2.1c': '2.3.2.1',
    '2.3.2.1d': '2.3.2.1',
    '2.3.2.1e': '2.3.2.1',
    '2.3.2.1f': '2.3.2.1',
    '2.3.2.1g': '2.3.2.1',
    '2.3.3': None,
    '2.3.4': None,
    '2.3.4.1': '2.3.4',
    '2.3.4.2': '2.3.4',
    '2.3.4.3': '2.3.4',
    '2.3.4.4': '2.3.4',
    '2.3.4.4a': '2.3.4.4',
    '2.3.4.4b': '2.3.4.4',
    '2.3.4.4c': '2.3.4.4',
    '2.3.4.4d': '2.3.4.4',
    '2.3.4.4e': '2.3.4.4',
    '2.3.4.4f': '2.3.4.4',
    '2.3.4.4g': '2.3.4.4',
    '2.3.4.4h': '2.3.4.4',
    '3': None,
    '4': None,
    '7': None,
    '6': None,
    '7.1': '7',
    '7.2': '7',
    '5': None,
    '9': None,
    '8': None
}


"""
This function is going to run through our tree and get nodes that define our clades based on the following criteria:
1. Nodes must have leaves that are only of a single clade, or of a single clade and its children
2. Nodes must be the most ancestral node that meets the critera in point 1 (so it must have the longest list of leaves)
3. If there is no node that defines a clade, the most ancestral leaf defines it (such as with 2.1.3.2b)

The parent list (which is the list above) decides which of our clades have children

it outputs a dictionary which looks like this:
    output = {'clade': 'node', 'clade': 'node', ...}

it also outputs a list called outliers, which do not have a node to define them and instead only have a single leaf
    outliers = {'strain': 'clade', ...}

"""
def find_nodes(NC, tree, parents):
    output = {} #initialize our variables here
    node_length = {} #node length is going to be important later as this is how we're going to find the most ancestral node
    for item in tree.Objects: #for each item in our tree
        if item.branchType == "node": #if it's a node
            length = len(item.leaves) #we're going to grab the list of leaves that node has and then get its length
            node_length[item.traits["name"]] = length #here we're going to make a dictionary where we can grab the length of our node at a later date
                                                        #since the node will be the key, and the length its value
    ####### THIS IS FOR CLADES WITH NO CHILDREN #######
    for node in NC: #for node in our list of nodes (from node_clades)
        if len(NC[node]) == 1 and NC[node][0] not in parents: #if our node only has a single clade and is not in our list of clades with children
            clade = NC[node][0] #our clade is the clade in our list (since it should only be 1 here)
            if "Am_nonGsGD" not in clade and "EA_nonGsGD" not in clade and "like" not in clade: #if our clade is not any of the ones listed and does not have "like in it"
                if clade in output.keys(): #if clade is already in our output
                    if node_length[output[clade]] < node_length[node]: #if our node is longer than the node we already have
                        output[clade] = node #we're going to replace our node in our output dictionary

                else:
                    output[clade] = node #otherwise, if we don't already have our clade in our output we're going to add it with our current node


    ####### THIS IS FOR CLADES WITH CHILDREN #######
        if len(NC[node]) > 1: #this is if our node has more than a single clade
            skip = False #initialize our skip variable -- I know skip isn't the best, but I really like this method of flagging things
            shortest = ['clade', 2000000] #we're going to intialize our shortest variable which will allow us to select by length later
            for clade in NC[node]: #for clade in our node's list of clades
                if len(clade) < shortest[-1]: #if the length of our clade is shorter than the shortest (which for the first one it always will be)
                    shortest = [clade, len(clade)] #we're going to replace the shortest variable with our clade, and it's length

                elif len(clade) == shortest[-1]: #if we have two clades taht are of equal length we're going to run through the next several elif statements
                    if clade == '6': #if our clade of equal length is 6 we grab 6
                        shortest = [clade, len(clade)]

                    elif clade == '2.1.3.1': #if our clade of equal length is 2.1.3.1 we grab 2.1.3.1
                        shortest = [clade, len(clade)]
                    #etc etc. The reason we do this, is because some of our clades are the same length (character-wise) as their children,
                    #so I wrote this code to grab the parent of a group of clades if multiple clades have the same length
                    elif clade == '3':
                        shortest = [clade, len(clade)]

                    elif clade == '2.1.2':
                        shortest = [clade, len(clade)]

            if shortest[0] in parents: #next, if our shortest is in our parents list (which it should pretty much always be if we're grabbing the right clade)
                for clade in NC[node]: #we go through each clade in our node
                    #here we're going to go through the other clades in our node, and if they are not the children or the clade itself,
                    #we're going to skip that node. I'll go through the first two as an example:


                    if "like" in clade and shortest != "3": #if like is in our clade and our shortest isn't 3, we're going to skip it
                        skip = True #this one is unique because 3 has a 3-like that we want to catch just because of the way our tree is, if you use a different tree
                                    #then this can be changed

                    elif shortest[0] == '2.2': #if our shortest is 2.2
                        if '2.2' not in clade[0:3] and '2.3' not in clade: #if 2.2 or 2.3 are not in our clade
                            skip = True #we skip. This is because 2.2 gives rise to all 2.2's and some 2.3's so we want to catch all of those
                            #the rest of this code is just the same as the logic above but unique clades for each one

                    elif shortest[0] == '6':
                        if '6' not in clade[0:3] and '7' not in clade[0]:
                            skip = True

                    elif shortest[0] == '2.1.2':
                        if '2.1.2' not in clade and '2.1.3' not in clade:
                            skip = True


                    elif shortest[0] == '2.1.3.1':
                        if '2.1.3.1' not in clade and '2.1.3.3' not in clade:
                            skip = True

                    elif shortest[0] == '3':
                        if '4' not in clade[0] and '3' not in clade[0]:
                            skip = True

                    elif shortest[0] not in clade[0:shortest[-1]]:
                        skip = True

                    elif shortest[0] not in parents:
                        skip = True

                #finally, if through all that our skip variable is still false,
                if skip == False:
                    if shortest[0] in output.keys(): #we check to see if our shortest clade is already in our keys
                        #the reason we do this is because our first chunk of code which only grabs nodes with a single clade is going to trigger
                        #for all clades, whether we want it to or not, so we're going to overwrite those incorrect assignments here
                        if shortest[-1] > node_length[output[shortest[0]]]:
                            output[shortest[0]] = node

                    else: #if our clade is not already in our dict, we're going to add it here
                        output[shortest[0]] = node
    #so here we have a clade that only has a single leaf on our tree, so we're just going to grab that strain name and set that as our clade-defining leaf instead of a node
    outliers = {"A/chicken/Indonesia/D10014/2010": '2.1.3.2b'}
    ####for different trees, this could be changed
    return output, outliers

"""
Here we are going to find the HA mutations for each clade, it is going to output a dictionary that looks like this:

node_HAmuts = {'node': ['mutation', 'mutation', 'mutation', ...], ...}
Will also append outliers with 'strain' = ['mutation', 'mutation', ...]

"""
def node_HA_muts(Tree, nf, outliers):
    node_HAmuts = {} #initialize our variable

    for object in Tree.Objects:
        if object.branchType == "node":
            mutations = object.traits["branch_attrs"]["mutations"] #grabbing our mutations so we don't have to write this out later
            name = object.traits["name"] #grabbing the name of our node
            for item in nf.keys(): #here we're going to find if our node is in our nodes dict
                if nf[item] == name: #if our node is in our nodes dict
                    if "HA" in mutations: #grab the HA mutations
                        node_HAmuts[item] = mutations["HA"] #define our node in our new dictionary where the value is our HA mutations


    for strain in outliers.keys(): #this is essentially the same code as above but we're looking only for a single strain and grabbing that
        for item in Tree.Objects:
            if item.traits["name"] == strain:
                mutations = object.traits["branch_attrs"]["mutations"]
                if "HA" in mutations:
                    node_HAmuts[outliers[strain]] = mutations["HA"]

    return node_HAmuts


'''
Here we are going to find the nucleotide mutations for each clade, it is going to output a dictionary that looks like this:

node_HAmuts = {'node': ['mutation', 'mutation', 'mutation', ...], ...}
Will also append outliers with 'strain' = ['mutation', 'mutation', ...]

'''
def node_nuc_muts(Tree, nf, outliers):
    node_nucmuts = {}
    #this is the same code as node_HA_muts so I'm no going to go through it as thoroughly

    for object in Tree.Objects:
        if object.branchType == "node":
            mutations = object.traits["branch_attrs"]["mutations"] #set varaible for later
            name = object.traits["name"] #grab node name
            for item in nf.keys():
                if nf[item] == name: #check if our node name is in our keys
                    if "nuc" in mutations: #if we have nucleotide mutations
                        node_nucmuts[item] =  mutations["nuc"] #set new dict. value


    for strain in outliers.keys():
        for item in Tree.Objects:
            if item.traits["name"] == strain:
                mutations = object.traits["branch_attrs"]["mutations"]
                if "nuc" in mutations:
                    node_nucmuts[outliers[strain]] = mutations["nuc"]

    return node_nucmuts

'''
This function is going to output our tsv file, it is not going to return anything and instead outputs a txt file which
can be used in augur clades
'''
def file_writer(HA, nuc, output_file_path, relationships):
    with open(output_file_path, "w") as tsv:
        tsv.write("clade\tgene\tsite\talt\n")
        for clade in nuc.keys():
            if relationships[clade] != None:
                tsv.write(f"{clade}\tclade\t{relationships[clade]}\n")
            for mutation in nuc[clade]:
                if "-" in mutation or "x" in mutation or "X" in mutation or "n" in mutation or "N" in mutation:
                    pass

                else:
                    tsv.write(f"{clade}\tnuc\t{mutation[1:-1]}\t{mutation[-1]}\n")

            if clade in HA.keys():
                for mutation in HA[clade]:
                    if "-" in mutation or "x" in mutation or "X" in mutation or "n" in mutation or "N" in mutation:
                        pass

                    else:
                        tsv.write(f"{clade}\tHA\t{mutation[1:-1]}\t{mutation[-1]}\n")

            tsv.write("\n")

        #this used to output some text that told us which clades we found no mutations for, but it is no readable by augur clades, so I removed
        #it for now but if you need to just visually see which clades aren't working, these 3 lines will output that info
        #for clade in relationships.keys():
            #if clade not in nuc.keys() and clade not in HA.keys():
                #tsv.write(f"No unique mutations for clade {clade} found")
    #since this function just writes our new file, we don't need to return anything
    return None


#get our clades and leafclades
allClades, leafClades = mf.leaf_clades(mytree, exclude_list, clade_column_name)

#getting our list of clades for each node, also getting our unclean nodes which is not essential right now but may be helpful when
#defining a function to ignore misassigned tips at some point
nodeClades, unclean = mf.node_clades(mytree, leafClades)

#getting our defining nodes, plus any outliers
nodes, leaves = find_nodes(nodeClades, mytree, parent_list)

#getting our HA mutations
HA = node_HA_muts(mytree, nodes, leaves)

#getting our nucleotide mutations
nuc = node_nuc_muts(mytree, nodes, leaves)

#filtering our mutations so only unique ones are returned (unique being they're found in only one clade)
unique_HA, unique_Nuc = mf.unique_muts(HA, nuc)

#this will write our tsv
file_writer(unique_HA, unique_Nuc, output_file_path, relationships)
