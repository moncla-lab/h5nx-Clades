# README 

This repository contains code for developing a method for defining and assigning clades for H5 viruses. The goal of this project is to develop a methodology similar to that employed for SARS-CoV-2 and seasonal influenza viruses, where clades are defined by identifying monophyletic clusters of viruses on the tree that share a set of clade-defining mutations. The goal is to generate a set of mutations that can be used to define clades, use those to assign clades, and develop a dataset that can be incorporated into NextClade. 

## Requirements: 

1. [baltic](https://github.com/evogytis/baltic) for tree parsing. Install with `pip install baltic`

2. [Nextstrain](https://docs.nextstrain.org/en/latest/install.html)

## Examples and resources:  
1. [Seasonal flu](https://github.com/nextstrain/seasonal-flu/blob/master/config/nextstrain_clades_h3n2_ha.tsv) clade-defining mutations format

2. [`augur_clades.py`](https://github.com/nextstrain/augur/blob/master/augur/clades.py)

## Input Data: 

There are a bunch of different input data files for this. 

* `H5 reference alignment.fas`: this is an alignment of ~400 sequences including the full breadth of the H5 clades. This was emailed to us from Todd Davis. It also includes the recent 2.3.4.4 and 2.3.2.1 splits. This could be the basis of a full tree clade assignment set. 

* `2.3.4.4_update_v1.fas`: this is an alignment of the new 2.3.4.4 virus splits from Tommy Lam. This is from an email to my Penn email on March 16, 2023, and indicated that he and Todd were going to publish a paper on this sometime soon making these splits public. So this is a good reference dataset for the 2.3.4.4bs. 

* `2.3.2.1c_update_v1.fas`: this is an alignment of the new 2.3.2.1c virus splits from Tommy Lam. This is from an email to my Penn email on March 16, 2023, and indicated that he and Todd were going to publish a paper on this sometime soon making these splits public. So this is a good reference dataset for the 2.3.4.4bs. 

## Notes:

Currently, the pipeline is as follows: 

1. Using the reference tree with clade-annotated sequences, run a Nextstrain build to generate a tree in JSON format. It is critical that this is done in a way that retains the reference sequence in the alignment, otherwise the mutations will all be meaningless. Currently, this updated set of sequences from Todd Davis is in `Nextstrain-build-with-updated-clades/data/h5nx-all-clades-reference-2023-04-13.fasta`. The results are in `Nextstrain-build-with-updated-clades/auspice`

2. Using that JSON output, run Sonja's script for calling mutations for each clade: `python avian-flu-clademaker/Mutation_Characterization.py --input_tree Nextstrain-Build-with-updated-clades/auspice/flu_avian_h5nx_ha.json --output_file Nextstrain-Build-with-updated-clades/clade-defining-mutations.txt --metadata_clade_column_name clade`. This will generate `clade-defining-mutations.txt`

3. Rerun with augur clades. The syntax is written out in `Nextstrain-build-with-updated-clades/Snakefile`

