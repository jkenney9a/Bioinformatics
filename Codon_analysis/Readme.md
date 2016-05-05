#Code for analyzing codon frequency usage

This script uses Bipython to take in a list of gene names and a given codon map, and returns the following for each gene:
* Average decoding time
* Normalized average decoding time (normalized to most frequent/fastest codon for a given amino acid)
* The ratio of difference between minimum and maximim protein lengths (relative to max protein length) downloaded from NCBI for a given gene
* The range of the protein lengths downloaded from NCBI for a given gene
* The median protein length of all protein lengths downloaded from NCBI for a given gene

The mRNA sequences are downloaded from the NCBI nucleotide database using Biopython using the following search terms:
  * *gene* + "[Gene Name] AND " + *organism* + "[Organism] AND mRNA[Filter] AND RefSeq[Filter] NOT PREDICTED"

###Inputs:
######Items in bold are required every time the script is run. Other items are only required if FASTA files for genes of interest need to be downloaded.

* **--gene_list** or **--list**
  * A list of gene names to be analyzed or can be in the format of "***GN=<gene name>***" that is the common output
of proteomic experimental data.

* **--output**
  * The output file name

* **--codon_map** or **--map**
  * A CSV file mapping codons to some value such as decoding time or relative frequency
  * At the moment the file must be in the order of: *codon*, *AA abbreviation*, *value*
  
* **--organism**
  * The organism that corresponds to the gene list. 
  * This should be put in quotes (e.g, "mus musculus")
  
* --download
  * True or False. Whether or not to downlaod data from NCBI. This is necessary to obtain the FASTA files for each gene of interest for analysis.
  * If FASTA files have already been downloaded, the script will check the FASTA directory and if the files alread exist, it will not download them again.
  * If FASTA files have not been downloaded, the script will create a directory and download the FASTA files into that directory.

* --email
  * If you are going to download data from NCBI (e.g, the FASTA gene files) then you need to provide your email address
  
