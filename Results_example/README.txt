This Directory can be downloaded at the end of each run of OTT from the site.
Under the "Results" section in the results page (e.g., http://onetwotree.tau.ac.il/results.html?jobId=1524727033&jobTitle=Narcissus) you will find the 'Download all output' button that will download a zipped directory of your OTT job.
The main files in this directory:

-----------------------------------------------------------------------------

1. 1524727033-concat-aligned.fasta / msa_nexus.txt : MSA file (fasta / nexus format).
2. Result_Tree_1524727033.tre / Result_Tree_1524727033.pdf: Phylogeny (Newick format / PDF file).
3. Result_Tree_1524727033_NotUltra.tre: in case the user chooses to reconstruct an ultrametric tree (set "Divergence time estimation" to "yes") OTT provides the original phylogeny, before calibration, as well.
4. All_Clusters_data.csv: a table with all clusters chosen for the final alignment, length, number of taxa included and type of each  marker.
5. AccessionsMatrix.csv: Cluster vs. taxa matrix (the accessions in each cluster per taxa).
6. FinalSpeciesList.txt / Final_TaxId_vs_Name.txt: Taxa list included in the final tree.
7. Added_LargeTaxIds.txt: list of taxIds (accessions) with more than 1000 sequences that were added to the analysis (these accessions are not part of the default pipeline, but rather are added to it after the clustering stage).
8. FinalStatus.txt: the final status of the job submitted on the server.
9. RemovedTaxa.txt: (under development) partial list of taxa that were not included in the final alignment and the reason for not being included (will be fully incorporated in future releases).
10. Result_Tree_1524727033.tre.r8s (not in use).
11. node_date_config.txt : if time calibration was requested the configuration file of the node dating process is shown.
12. summary_file.txt: summary data including initial/final species number and outgroup species seleced by OTT.
13. userInput.txt: list of taxa as entered by the user in OTT input feild.
14. web_partition.txt: partition file or the markers used to reconstruct the phylogeny.
