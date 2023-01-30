This is the documentation which is a work-in-progress for the Avian-flu Clademaker project

This is a pipeline where we are taking a fast file from a source such as GenBank and using it to define nodes that branch off into a single clade.

The way the pipeline works is the following:

To start you will need a fasta file from GenBank where the formatting is the following:
The search set should be the followingL

Host: Any
Protein: HA
Subtype H: 5
Subtype N: Any

Fasta defline:
>{strain}|{host}|{accession}|{year}-{month}-{day}|{country}

be sure there is a carrot ">" at the front

Once you have that downloaded, it should be added to the directory and named h5nX_ha.fa (or any other name as long as you change the name of the input file in fasta_formatter.py, clade_annotater.py and mutation_finder.py)

Now run fasta_formatter.py It should output a file named h5n1_ha.fasta
***It will not make a new file if there is already a file with this name***

This will give you a file that is usable with nextstrain and with the clade_annotater.py
Next you should run LABEL on your file, and grab the list of clades it outputs as a file, naming it: h5n1_ha_clades_final.txt **or anything else as long as you rename it in clade_annotater.py


Run the clade_annotater.py with your newly generated file. This will give you an output file that is a fasta file similar to the one you initially downloaded, but which contains the clades for each strain added to the end of its metadata lines

From here you can run the fasta with Nextstrain to generate a tree which should have the option to color by clade. ***THIS IS IMPORTANT***
Without the ability to sort by clade the mutation_finder.py file will not work

Once you have this tree, grab the JSON file and drag into the same directory as mutation_finder.py
Run Mutation_finder.py

This will print out a list of your unique mutations for each clade as well as the defining node for each clade. STILL WORKING ON AN OUTPUT FILE FOR THIS


