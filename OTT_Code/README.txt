The code behind OneTwoTree web is mainly written in python. OneTwoTree.py creates the main Job###.sh file 
that will be submitted to our server's queue. Once the job is completed the user will be notified by email 
that the job was completed and the results are available online and ready to be downloaded (zipped dir similar to the Results_example: https://github.com/MayroseLab/OneTwoTree/tree/master/Results_example).

The code is arranged in 3 main directories:
clustering: perl files for performing the main clustering stage of OTT.
ott_objects_defs: python files holding main definitions of classes and functions used by OTT.
ott_scripts: python/R/perl scripts used in our OTT flow.
webServer_files: perl files handling web server parameters.
