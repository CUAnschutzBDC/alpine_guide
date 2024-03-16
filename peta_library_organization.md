# Peta library Orgaization

I have been working to make the peta library space easy to understand and find everything you need.

The base directory is `/pl/active/Anschutz_BDC`

## resources
I will try to put anything I think may be helpful to you here. Everything in this directory is backed up to the peta archive allocation and aws glacial storage.

Resources include:

### ref
This will be any genome information. Feel free to add your own genome to this folder if you can't find it already so everyone else can have access.

`annotation` - includes all GTF and BED files
`cellranger` - includes cellranger formatted references
`genome` - includes fasta files for genomes
`indicies` - includes indicies necessary for aligners like `star` and `minimap2`

Each of the above directories have the exact same structure. First is the organism, second is the genome build. Please keep this structure so people can easily see the genomes available for their organism.

### single_cell_references
These are single cell references that I've found helpful. Many have seurat objects so you can immediately explore the data. Some have also been made into `clustifyr` references for cell type naming.

### singularity
I will put my singularity images here. These are helpful for reproducability,

### snakemake
These are all of my snakemake pipelines. I'm slowly updating them for use on Alpine

### tutorials
Different tutorials I've taught, including this one!

## analysis

This is where all of your analysis can be stored. To keep it organized, I recommed you making a directory for your PI and then making a direcotory for yourself within that.

For example, Maria, who is in Lori's lab would have this directory `/pl/active/Anschutz_BDC/analysis/sussel/maria`. She can then do all of her analysis there without accidentally working on someone else's data.

## data

This is where I put all the raw fastq files. Each PI has their own folder. I then use the name of the directory given by the sequencing core. Everything in this directory is backed up to the archive allocation and also to aws glacial storage.

Because of this backup system, please tar all fastq and pod5 files produced by nanopore runs as these small files will not work in the archive storage.

## bin

This is where I will keep all of the executable files related to packages I think may be important. You can add this to your path in your `.bashrc` if you find these packages helpful.

## packages
This is where I put packages I have installed from places outside of conda. For example, you can find Cellranger here
