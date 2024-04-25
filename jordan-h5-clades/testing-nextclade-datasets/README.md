#### There are three folders, corresponding to each of the H5Nx Nextclade datasets to be generated:
- all H5Nx clades (/h5nx/)
- 2.3.2.1 and descendant clades (/2321/)
- 2.3.4.4 and descendant clades (/2344/)
<br></br>
#### Within each folder, there are two folders containing:
- files needed to build the Nextstrain tree (/_clade_/build/)
  - the Snakemake file used for the Nextstrain build (/_clade_/builds/Snakefile)
  - the config folder used for the build, including additional config specific to the clade-defining-mutations script (/_clade_/builds/config/)
  - the output clade-defining mutations used to assign clades to internal nodes (/_clade_/builds/clade_defining_mutations_h5nx_ha.tsv)
- files that make up the Nextclade dataset (/_clade_/files/)
  - See the [Nextclade documentation](https://docs.nextstrain.org/projects/nextclade/en/stable/user/datasets.html) for more information on each file
<br></br>
#### Nextclade supports direct links to datasets publicly available on Github
To use a custom dataset, use [clades.nextstrain.org/?dataset-url=_???_](https://clades.nextstrain.org/?dataset-url=???) and replace _???_ with a link to the Github folder

[H5Nx (all clades) Nextclade](https://clades.nextstrain.org/?dataset-url=https://github.com/moncla-lab/h5nx-Clades/tree/main/jordan-h5-clades/testing-nextclade-datasets/h5nx/files)

[H5Nx clade 2.3.2.1](https://clades.nextstrain.org/?dataset-url=https://github.com/moncla-lab/h5nx-Clades/tree/main/jordan-h5-clades/testing-nextclade-datasets/2321/files)

[H5Nx clade 2.3.4.4](https://clades.nextstrain.org/?dataset-url=https://github.com/moncla-lab/h5nx-Clades/tree/main/jordan-h5-clades/testing-nextclade-datasets/2344/files)
