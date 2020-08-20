# Bin_metabolome_analysis
Some scripts to analyze the metabolome of the bin genome.

# Usage
User have to download the all metabolite dataset from Human Metabolome Database (HMDB) and put it in the references folder (https://hmdb.ca/downloads). 
The input data should be tsv format (See demo in the inputs folder). The first column shoud be the orf ID of the bin genome, and the second column represents the KO number of the orf.  
Run the python script to collect the metabolome data from KEGG and HMDB, which will generate a csv format data.  
After the info collection, use the R script to summarize the results.  

# Examples
Python Bin_ko_map_Compounds.py -I ../inputs/6462_u_ko.txt -O ../outputs/6462.csv

# References
Wishart, D. S. et al. HMDB 4.0: the human metabolome database for 2018. Nucleic Acids Res. 46, D608-D617 (2017).  
Ogata, H., . et al. KEGG: Kyoto Encyclopedia of Genes and Genomes. Nucleic Acids Res. 27, 29-34 (2000).
