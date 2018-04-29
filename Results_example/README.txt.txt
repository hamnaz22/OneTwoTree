
This Directory can be downloaded at the end of each run of OTT from the site:
Under Results you can find the 'Download all output' button that will download 
a sipped summary_dir of your OTT job.

The main files in this directory:
-----------------------------------------------------------------------------
1. 1524727033-concat-aligned.fasta / msa_nexus.txt : MSA file (fasta and nexus format).
2. Result_Tree_1524727033.tre / Result_Tree_1524727033.pdf: Phylogeny (Newick format / PDF file).
3. Result_Tree_1524727033_NotUltra.tre: in case the use chose the ultrametric option then we also provide the original version of the phylogeny, before calibration.
4. All_Clusters_data.csv: a table with all clusters chosen for the final alignment, length, number of taxa included and type.
5. AccessionsMatrix.csv: Accessions matrix according to cluster distribution.
6. FinalSpeciesList.txt / Final_TaxId_vs_Name.txt: Taxa list included in the final tree.
7. Added_LargeTaxIds.txt: list of taxIds (accessions) with more than 1000 sequences that were added to the analysis.
8. FinalStatus.txt: the final status of the job submitted on the server.
9. RemovedTaxa.txt: (still under development) list of taxa that was not included in the final alignment and the reason for not being included.
10. Result_Tree_1524727033.tre.r8s (not in use).
11. node_date_config.txt
12. summary_file.txt: summary data including initial/final species number and outgroup species seleced by OTT.
13. userInput.txt: list of taxa as entered by the user in OTT input feild.
14. web_partition.txt: partition file for the phylogeny presented on the web.
