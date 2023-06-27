Avian flu clademaker README
######################################################################
The order in which the files should be run is the following:
------- BEFORE NEXTSTRAIN -------
unformatted fasta file --> fasta_formatter.py --> formatted fasta file
Run through LABEL (if still using LABEL) otherwise once you have some other txt of clades
formatted fasta file + clade.txt --> clade_annotator.py --> fasta file with clades
fasta file with clades --> 2.3.4.4_annotator.py --> fasta file with 2.3.4.4 clades
------- AFTER NEXTSTRAIN -------
JSON --> mutation_characterization.py --> mutations.txt
######################################################################

NOTE: Mutation_Finder.py is important but is NOT meant to run on its own. It is pretty much only meant to run from Mutation_Characterization.py

The first thing that needs to be done is getting a fasta file with some data. The fasta_formatter.py is designed for data from GenBank, but, there is some code already written in reference_tree_maker.py that should reformatt the LABEL guide tree

##### The reference_tree_maker.py file is not commented at this moment, but largely is older code that is not used, it is just an example of string parsing for the LABEL guide tree #####

Once you have a fasta, run fasta_formatter.py, which should output a file (named whatever you like) in 
/Output

Once you have that output file, go ahead and run it through LABEL or augur clades,
 whatever way you want to get a text file with clade assignments. 
 Next, run it through clade_annotator.py. IMPORTANT NOTE -- clade_annotator.py is made specifically 
 for LABEL assignments and will need to be either changed or completely re-written if another assignment tool is used.

clade_annotator.py will output a fasta that has your strains and sequence data with clades annotated at the end in the format metadata|clade

From here, move your output file to /Data and then run 2.3.4.4 annotator if you have additional clades you want to add. 
There is already a 2.3.4.4 guide file in /Data as the baseline package for this code, 
BUT, if there are other clades you want to add it is pretty simple to change that file, 
just re-write the string parsing to more closely align with whatever file you have rather than the guide file that is already in there.

At this point you should take the fasta file you've made through nextrstrain to produce a JSON file, 
which is required for the remainder of the code.

The most important thing about this code is that the mutation_finder code is written very explicitly for
each JSON file, and will need changes based on what your JSON file looks like. The biggest problem
currently is that there is no automatic way to remove unassigned tips from the mutation finder code,
instead, strains must be automatically added to the exclude list that is initialized at the start of the file.

### This would be a great thing to work on in the future, as the larger the fasta file the more likely it is nextstrain is going to missassign tips ###

The mutation finder, written as is, only needs to be run (with the input and output paths changed for your specific files) and will give you an output tsv of your clades with their unique mutations.

All of this code can be changed based on your input fasta, the output you want, your JSON file, etc. 
it is completely changeable and will require some personalization based on your input data.

Importantly, .py files take inputs from /Data and output files into /Output
