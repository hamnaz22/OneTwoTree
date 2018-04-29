The code behind OneTwoTree web is mainly written in python. OneTwoTree.py creates the main Job###.sh file 
that will be submitted to our server's queue. Once the job is completed the user will be notified by email 
that his job was completed and the results are available online and ready to be downloaded (zipped dir similar 
to the Results_example: https://github.com/MayroseLab/OneTwoTree/tree/master/Results_example).

The code is splited into 3 main directories:
clustering: perl files for perfoeming the main clustering stage of OTT.
ott_objects_defs: python files holding main definitions of classes and functions used by OTT.
ott_scripts: python/R/perl scripts used in out OTT flow.
webServer_files: perl files handeleing web server parameters and email notice to the user.
