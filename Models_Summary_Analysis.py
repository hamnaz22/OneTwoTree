__author__ = 'ItayM3'

# Import the os module, for the os.walk function
import os
import csv
import sys
import time
import mmap
import shutil
#from Tree_Globals import *


# python src_MD/Models_Summary_Analysis.py ../GeneBank/Genera_Groups_CCDBLists/Angiosperms_genera.txt


SCRIPTS_DIR = "/groups/itay_mayrose/michaldrori/scripts/"
PRE_TREE_DIR = 'PreTree_Dir'
PHY_EXT = '_phy.phy'
JMT_EXT = '_phy.phy.jmt'
FASTA_FORMAT = "fasta"
PHYLIP_FORMAT = "phylip-relaxed"

MSA_OUTPUT_PATH = '/groups/itay_mayrose/michaldrori/MSA_JmodelTree/output_MD/'
TABLE_HEADERS = "Model\s+-lnL\s+K\s+TEST\s+delta\s+weight\s+cumWeight(\s+uDelta)?"
TREE_LINE = "Tree for the best TEST model = "

DEBUG_PRINT = 1

WORKING_DIR = '/groups/itay_mayrose/michaldrori/MSA_JmodelTree/output_MD/'


def calc_precentage(genusName,NameToRes_file):
    Total_Names = -1
    Exact_Name = 0
    if os.path.exists(NameToRes_file):
        with open(NameToRes_file,'r') as f:
            for line in f:
                if genusName in line:
                    Exact_Name+=1
                    Total_Names+=1
                else:
                    Total_Names+=1
    else:
        return 99999
    return float(Exact_Name)/float(Total_Names)*100



def returnSpeciesNumber(genusName):
    Num = 0
    file_name = WORKING_DIR + genusName + '/concat/' + genusName + '-concat-aligned.fasta'         # Erica-concat-aligned.fasta
    if os.path.exists(file_name):
        with open(file_name,'r') as f:
            for line in f:
                if '>' in line:
                    Num+=1
    else:
        return 12321
    return Num


def returnSpeciesNumberSeqs(genusName):
    #Count the number of species in the input allseq file and compare to the number in the alignment:
    #Input file:
    species_withNoGenusName=[]
    species_names=[]
    organs_allseq_file = WORKING_DIR + genusName + '/' + genusName + '-allseq.fasta'
    if os.path.exists(organs_allseq_file):
        with open(organs_allseq_file,'r') as f:
            for line in f:
                if '>gi' in line:
                    parsed_line = line.split('|')
                    species_names.append(parsed_line[5])
                    if not genusName in parsed_line[5]:
                        species_withNoGenusName.append(parsed_line[5])
    NumSpeciesInitial = (len(set(species_names)))
    species_withNoGenusName=set(species_withNoGenusName)

    species_names_final=[]
    align_file = WORKING_DIR + genusName + '/concat/' + genusName + '-concat-aligned.fasta'
    if os.path.exists(align_file):
        with open(align_file,'r') as f:
            for line in f:
                if '>' in line:
                    parsed_line = line.split('>')
                    species_names_final.append(parsed_line[1])
                    #print(species_names_final)
    species_names_final = (len(set(species_names_final)))
    #print(NumSpeciesInitial)
    #print(species_names_final)
    #print("___________________________")
    return (NumSpeciesInitial, species_names_final, species_withNoGenusName)





    return NumSpeciesInitial

    #file_name = WORKING_DIR + genusName + '/concat/' + genusName + '-concat-aligned.fasta'
    #NumOfSpecies = 0
    #Num = 0
    #if os.path.exists(file_name):
    #    with open(file_name,'r') as f:
    #        for line in f:
    #            if '>' in line:
    #                NumOfSpecies+=1
    #return NumOfSpecies


def copyDirectory(src, dest):
    try:
        shutil.copytree(src, dest)
    # Directories are the same
    except shutil.Error as e:
        print('Directory not copied. Error: %s' % e)
        return 1
    # Any error saying that the directory doesn't exist
    except OSError as e:
        print('Directory not copied. Error: %s' % e)
        return 1
    if os.path.exists(src):
        #rm_command = "rm " + "-rf " + src
        os.remove(src)
    return 0




def Read_CSV (file_name, f_out):

    DebugFlag = 0
    PloidyTot = 0
    row_num = 0
    TaxonName = []
    PolyPloidy = []
    Robust_Score = []
    Robust_Score_minus = []
    Simulation_Score = []
    Simulation_Score_minus = []
    ChromNumList = []
    PloidyCount = 0
    Mean_val = 0
    PolyPloidyNum = 0
    start_idx = 0
    NoChromCount = 0
    NAPloidyInf = 0
    PolyLower_Flag=0
    TaxonName_Main = "MAIN"
    LowPolyPloid=0
    MinPolyploidy=0

    if(DebugFlag == 1):
        f_out.write("The CSV file working on is: %s\n" % file_name)
    f = open(file_name)
    for row in csv.reader(f):
        if row[0]=="Taxon":
            continue
        else:
            if row != 0:
                TaxonName.append(row[start_idx])
                if row[start_idx+1]=='x':
                    ChromNum = 999
                    NoChromCount+=1
                else:
                    ChromNum = int(row[start_idx+1])
                    Robust_Score_minus.append(float(row[start_idx+3]))
                    Simulation_Score_minus.append(float(row[start_idx+4]))
                    #print(mean(row))
                if row[start_idx+2]=='NA':
                    Ploidy = -1
                    NAPloidyInf+=1
                else:
                    Ploidy = float(row[start_idx+2])
                if (Ploidy == 0) & (ChromNum != 999):
                    PloidyTot += ChromNum   #Sum of all Ploidy taxons for average calculation
                    PloidyCount +=1
                    ChromNumList.append(int(row[start_idx+1]))
                if (Ploidy == 1) & (ChromNum != 999):
                    PolyPloidy.append(int(row[start_idx+1]))
                    PolyPloidyNum+=1
                Robust_Score.append(float(row[start_idx+3]))
                Simulation_Score.append(float(row[start_idx+4]))
    row_num=row_num+1
    # EndOf file read

    TaxonName_Main=TaxonName[0].split("_")[0]

    f_out.write("%s, " % TaxonName_Main)
    f_out.write("%d, " % len(TaxonName))
    if(PloidyCount==0):
        f_out.write("NO DEPLOIDS!!!,")
    else:
        Mean_val=PloidyTot/PloidyCount
        #print ("Mean value of Ploidy is:", Mean_val)
        f_out.write("%.02f," % Mean_val)
    f_out.write("%d -> %.2f, " % (NoChromCount, (100*(NoChromCount/len(TaxonName)))))
    f_out.write("%d -> %.2f, " % (NAPloidyInf, (100*(NAPloidyInf/len(TaxonName)))))
    f_out.write("%.2f, " % (sum(Robust_Score)/len(TaxonName)))
    f_out.write("%.2f, " % (sum(Robust_Score_minus)/(len(TaxonName)-NoChromCount)))
    f_out.write("%.2f, " % (sum(Simulation_Score)/len(TaxonName)))
    f_out.write("%.2f, " % (sum(Simulation_Score_minus)/(len(TaxonName)-NoChromCount)))

    if(all_print_flag==1):
        f_out.write("Total Taxon Number in file: %d, " % len(TaxonName))
        f_out.write("Taxon with X Chromosome count: %d -> %.2f, " % (NoChromCount, (100*(NoChromCount/len(TaxonName)))))
        f_out.write("Taxon with NA Ploidy inference: %d -> %.2f, " % (NAPloidyInf, (100*(NAPloidyInf/len(TaxonName)))))
        #f_out.write("Phylogeny robustness score AND Simulation reliability score\n")
        f_out.write("PhylogenyRobustScore avrg: %.2f, " % (sum(Robust_Score)/len(TaxonName)))
        f_out.write("PhylogenyRobustScore avrg - No Chrom data: %.2f, " % (sum(Robust_Score_minus)/(len(TaxonName)-NoChromCount)))
        f_out.write("SimulReliabilityScore avrg: %.2f, " % (sum(Simulation_Score)/len(TaxonName)))
        f_out.write("SimulReliabilityScore avrg - No Chrom data: %.2f, " % (sum(Simulation_Score_minus)/(len(TaxonName)-NoChromCount)))


    # We might have a problem

    if not PolyPloidy:
        f_out.write("No Polyploids,")
    else:
        if not ChromNumList:
            f_out.write("No Deploids,")
        else:
            MinPolyploidy=min(PolyPloidy)
            Max_Ploidy=max(ChromNumList)
            if(MinPolyploidy<Max_Ploidy):
                f_out.write("YES (Ploidy %02f PolyPloidy %02f)" % (Max_Ploidy,MinPolyploidy))
            else:
                f_out.write("NO,")

    f.close()
    if(f.closed == False):
        print ("!!!!!!!!!  ERROR in Closing file   !!!!!!!!", f)
    return

def readModelsSummary (file_name, f_out):
# BASE_NUM
# BASE_NUM_DUPL
# CONST_RATE
# CONST_RATE_DEMI
# CONST_RATE_DEMI_SET
# CONST_RATE_NO_DUPL

    DebugPrint = 0

    BASE_NUM_count=0
    BASE_NUM_DUPL_count=0
    CONST_RATE_count=0
    CONST_RATE_DEMI_count=0
    CONST_RATE_DEMI_SET_count=0
    CONST_RATE_NO_DUPL_count=0

    NumOfTrees=0

    if(DebugPrint == 1):
        f_out.write("ModelsSummary is working on: %s" % (file_name))
    with open(file_name) as f:
        content = f.readlines()
        for i in range(len(content)):
            #print (content[i])
            if "tree" in content[i] and "BASE_NUM_DUPL" in content[i]:
                BASE_NUM_DUPL_count+=1
                NumOfTrees+=1
                continue
            if "tree" in content[i] and "BASE_NUM" in content[i]:
                BASE_NUM_count+=1
                NumOfTrees+=1
                continue
            if "tree" in content[i] and "CONST_RATE_DEMI_SET" in content[i]:
                CONST_RATE_DEMI_SET_count+=1
                NumOfTrees+=1
                continue
            if "tree" in content[i] and "CONST_RATE_DEMI" in content[i]:
                CONST_RATE_DEMI_count+=1
                NumOfTrees+=1
                continue
            if "tree" in content[i] and "CONST_RATE_NO_DUPL" in content[i]:
                CONST_RATE_NO_DUPL_count+=1
                NumOfTrees+=1
                continue
            if "tree" in content[i] and "CONST_RATE" in content[i]:
                CONST_RATE_count+=1
                NumOfTrees+=1
                continue
    #f_out.write("NUM of repetitions for each Model: (Trees num = '%d'\n" % NumOfTrees)
    #f_out.write("+++++++++++++++++++++++++++++++++++++++++++++++\n")

    #f_out.write("BASE_NUM,BASE_NUM_DUPL,CONST_RATE,CONST_RATE_DEMI,CONST_RATE_DEMI_SET,CONST_RATE_NO_DUPL\n")

    f_out.write("%.2f ," % (BASE_NUM_count))
    f_out.write("%.2f ," % (BASE_NUM_DUPL_count))
    f_out.write("%.2f ," % (CONST_RATE_count))
    f_out.write("%.2f ," % (CONST_RATE_DEMI_count))
    f_out.write("%.2f ," % (CONST_RATE_DEMI_SET_count))
    f_out.write("%.2f \n" % (CONST_RATE_NO_DUPL_count))

#   if BASE_NUM_count != 0:
#       f_out.write("%.2f ," % (BASE_NUM_count))
#   if BASE_NUM_DUPL_count != 0:
#       f_out.write("%.2f ," % (BASE_NUM_DUPL_count))
#   if CONST_RATE_count != 0:
#       f_out.write("%.2f ," % (CONST_RATE_count))
#   if CONST_RATE_DEMI_count != 0:
#       f_out.write("%.2f ," % (CONST_RATE_DEMI_count))
#   if CONST_RATE_DEMI_SET_count != 0:
#       f_out.write("%.2f ," % (CONST_RATE_DEMI_SET_count))
#   if CONST_RATE_NO_DUPL_count != 0:
#       f_out.write("%.2f \n" % (CONST_RATE_NO_DUPL_count))

  #fwrite.write( "INSERT INTO `wp_users` (`user_login`, `user_name`) " )
  #fwrite.write( "VALUES ('%s','%s');\n" % (safe_name,safe_login) )

    f.close()
    if(f.closed == False):
        print ("!!!!!!!!!  ERROR in Closing file   !!!!!!!!", f)
    return


def WorkOnPathList (file_name, f_out):


    try:
        f = open(file_name, 'r')
    except IOError:
        f_out.write('cannot open', f)

    with open(file_name) as f:
        PathList = f.readlines()
    for pathName in PathList:
        Read_CSV(pathName, f_out)
        #f_out.write("\n")

        return

############################################################################

def check_tree_status(genusName):
    if os.path.exists(MSA_OUTPUT_PATH+"/" + genusName+ "/concat/mb.out.con.tre"):   # MrBayes Done!!!
        if not os.path.exists(MSA_OUTPUT_PATH+"/" + genusName+ "/concat/parsemb_trees.tre"):       # R failed!!!
            return 'Failed_R: no parsemb_trees.tre file '
        return 'Ready For ChromEvol'                                                                                # Ready For ChromEvol
    else:
        return 'NO Tree'                                                                                 # No Tree!!!  -  Need to run MB or check if failed

############################################################################

def check_chromevol_status(genusName):
#/groups/itay_mayrose/michaldrori/MSA_JmodelTree/output/Croton/Croton_Chromevol/chromevol_out
    if os.path.exists(MSA_OUTPUT_PATH+"/" + genusName + "/" + genusName + "_Chromevol/chromevol_out/ploidy.csv"):
        return 'yes'
    else:
        if os.path.exists(MSA_OUTPUT_PATH+"/" + genusName + "/" + genusName + "_Chromevol/chromevol_out/NoPloidyInference_moreThan50Trees"):
            return 'noPloidy'
        if not os.path.exists(MSA_OUTPUT_PATH+"/" + genusName + "/" + genusName + "_Chromevol/chromevol_out/infer/"):
            if os.path.exists(MSA_OUTPUT_PATH+"/" + genusName + "/" + genusName + "_Chromevol/chromevol_out/No_Infer_directory_MissingCounts.txt"):
                return 'No_Counts'
        else:
            return 'none'




############################################################################
#
#           M A I N
#
############################################################################
file_genera_list=sys.argv[1]   # Genera names for analysis


#rootDir = "/groups/itay_mayrose/share/ploidb/chromevolout-11-02-15/"
#rootDir = "."
# Set the directory you want to start from
#rootDir = '../../itaymay/share/ploidb/chromevolout-11-01-15/'
#rootDir = input("Please enter root dir:")

#rootDir = "/groups/itay_mayrose/share/ploidb/chromevolout-11-02-15/Ardisia/"
# at the beginning:
start_time = time.time()

# Chromevol Analysis
f_out=open('/groups/itay_mayrose/michaldrori/MSA_JmodelTree/Analysis_Files/Analysis_LOG_md.txt', 'w')
f_passed_genera=open('/groups/itay_mayrose/michaldrori/MSA_JmodelTree/Analysis_Files/Passed_AfterAlign_md.txt', 'w')
f_out_chromevol = open('/groups/itay_mayrose/michaldrori/MSA_JmodelTree/Analysis_Files/ChromEvol_Analysis_md.csv', 'w')
f_out_chromevol.write("Taxon_Name, Taxon Number , Mean value - Ploidy , Unknown ChromosomeNum , NA Ploidy inference , PhylogenyRobustScore avrg , PhylogenyRobustScore avrg - No Chrom data , SimulReliabilityScore avrg , SimulReliabilityScore avrg - No Chrom data , PolyPloidy LOWER than Deploidy, BASE_NUM,BASE_NUM_DUPL , CONST_RATE,CONST_RATE_DEMI , CONST_RATE_DEMI_SET , CONST_RATE_NO_DUPL\n") #HEADLINE

all_print_flag = 0



#Alignment/Tree Analysis:
Fail_Flag=0
f = open(file_genera_list)
list_passed_Align = []
#Counters for statistics:
Less_then5Species_Cnt=0
NoOutgroup_Cnt=0
NoSpecies_Cnt=0
NoData_inClusters_Cnt=0
OrthoMcl_issue_Cnt=0
StillRunning_Cnt=0
Total_run_Cnt=0
Passed_Alignment = 0

for row in csv.reader(f):
    genusName=row[0]
    f_out.write(genusName)
    # Verify The right Taxon name was identified:
    NameToRes_file = MSA_OUTPUT_PATH + "/" + genusName + "/names_to_resolve.csv"
    # Check Number of Species - after Alignment Vs. All species that have seqs
    Align_Debug_file = MSA_OUTPUT_PATH + "/" + genusName + "/" + genusName + "-orgnames-allseq.fasta"
#    list_numbers=[]
#    Num = 0
#    if os.path.exists(Align_Debug_file):
#        with open(Align_Debug_file,'r') as f:
#            for line in f:
#                if '>gi|' in line:
#                    parse_line = line.split('|')
#                    if parse_line[3] not in list_numbers:
#                        list_numbers.append(parse_line[3])
#            print('\n')
#            print(list_numbers)
#    SpeciesNum_Align = returnSpeciesNumber(genusName)
    #f_out.write('% s : Number Of Species in Alignment %d \n' % (genusName,SpeciesNum_Align))
    # Check number in debug file - > length of species_list
    Align_Debug_file = (MSA_OUTPUT_PATH + "/" + genusName + "/" + genusName + "-debug.log")
    SpeciesNum_Seqs =  11 #returnSpeciesNumberSeqs(Align_Debug_file)
    #f_out.write('% s : Number Of Species with seqs %d \n' % (genusName,SpeciesNum_Seqs))
    concat_dir = (MSA_OUTPUT_PATH + "/" + genusName + "/concat/")
    if not os.path.exists(MSA_OUTPUT_PATH + "/" + genusName):
        f_out.write("," + "No Alignment data - Missing directory for this genus")
    else:
        Total_run_Cnt+=1
        if not os.path.exists(Align_Debug_file):
            f_out.write("," + "FAILED ALIGN: No debug file")
    ##    else:
    ##        align_stop_file = MSA_OUTPUT_PATH + "/" + genusName + "/pre_tree_stopped"
    ##        if os.path.exists(align_stop_file):
    ##            f_out.write("," + "pre_tree_stopped: ")
    ##            with open(align_stop_file, "r") as stop_f:
    ##                for line in stop_f:
    ##                    f_out.write(line)
        else:
            if os.path.exists(Align_Debug_file):
                f_align = open(Align_Debug_file, 'rb', 0)
                s = mmap.mmap(f_align.fileno(), 0, access=mmap.ACCESS_READ)
                if s.find(b'Finished execution') != -1:
                     if s.find(b'Outgroup list is EMPTY !!!') != -1:
                        f_out.write("," + "FAILED ALIGN:  No Outgroup")
                        NoOutgroup_Cnt+=1
                        Fail_Flag=1
                     else:
                        f_out.write("," + "PASSED Alignment,")
                        Passed_Alignment+=1
                        list_passed_Align.append(genusName)
                        foreign_Species=[]
                        #TaxonNameOutOfAll = calc_precentage(genusName,NameToRes_file)
                        #f_out.write(', # %f ,' % TaxonNameOutOfAll)
                        (SpeciesNum_Init, SpeciesNum_final,foreign_Species)  =  returnSpeciesNumberSeqs(genusName) # Exclude OUT-GROUP
                        if not foreign_Species:
                            f_out.write(' Species# before %d,after %d, ' % (SpeciesNum_Init,(SpeciesNum_final-1)))  # -1 to remove Outgroup !!!
                        else:
                            join_foreignList=(', '.join(foreign_Species))
                            f_out.write(' Species# before %d,after %d, foreign[%s]' % (SpeciesNum_Init,(SpeciesNum_final-1),join_foreignList))  # -1 to remove Outgroup !!!
                        if check_tree_status(genusName)=='Ready For ChromEvol':
                            f_out.write(" After ChromEvol, ")
                            if check_chromevol_status(genusName)=='yes':
                                f_out.write(" Found Ploidy File, ")
                            else:
                                if check_chromevol_status(genusName)=='noPloidy':
                                    f_out.write(" Over 50% NO_DUPL models, ")
                                else:
                                    if check_chromevol_status(genusName)=='No_Counts':
                                        f_out.write(" Not enough Counts, ")
                                    else:
                                        if os.path.exists(MSA_OUTPUT_PATH+"/" + genusName + "/" + genusName + "_Chromevol/chromevol_out"):
                                            f_out.write(" Failed ChromEvol!!!, ")
                        else:
                            f_out.write(check_tree_status(genusName))
                # Check failure status options:
                else:
                     if s.find(b'No species found - STOPPING') != -1:
                        f_out.write("," + "FAILED ALIGN:  No species found")
                        NoSpecies_Cnt+=1
                        Fail_Flag=1
                     elif s.find(b'Less than 5 species found - STOPPING') != -1:
                        f_out.write("," + "FAILED ALIGN:  Less than 5 species found")
                        Less_then5Species_Cnt+=1
                        Fail_Flag=1
                     else:
                        if s.find(b'Not enough data in clusters - STOPPING') != -1:
                            f_out.write("," + "FAILED ALIGN: Not enough data in clusters")
                            NoData_inClusters_Cnt+=1
                            Fail_Flag=1
                        else:
                            #if s.find(b'Found 0 results for db=taxonomy term='+genusName) != -1:
                            if s.find(b'QQQ___') != -1:
                                 f_out.write("," + "FAILED ALIGN: Found 0 results for db=taxonomy term="+genusName)
                                 Fail_Flag=1
                            else:
                                if s.find(b'Found 0 candidates for outgroup') != -1:
                                    f_out.write("," + "FAILED ALIGN: Found 0 candidates for outgroup")
                                    Fail_Flag=1
                                else:
                                    clustering_seqs_path = (MSA_OUTPUT_PATH + "/" + genusName + "/clustering_seqs")
                                    if os.path.exists(clustering_seqs_path):
                                        f_out.write(", Still running......")
                                        StillRunning_Cnt+=1
                                    else:
                                        f_out.write("," + "FAILED ALIGN: Check for OrthoMCL Failure")
                                        OrthoMcl_issue_Cnt+=1
                                    #for file in os.listdir(MSA_OUTPUT_PATH+"/" + genusName + "/"):
                                    #    if file.endswith(".ER"):
                                    #        f_out.write(file)
                                    #        with open(MSA_OUTPUT_PATH+"/" + genusName + "/" + file, 'r') as f_ER:
                                    #            for line in f_ER:
                                    #                if 'Errno' in line:
                                    #                    f_out.write(genusName + "," + "FAILED ALIGN: %s" % line)
                                    #                    Fail_Flag=1
                                        #else:
                                         #   f_out.write(genusName + "," + "Unknown Status, maybe still running")
    src = MSA_OUTPUT_PATH + "/" + genusName
    dst = MSA_OUTPUT_PATH + "/FAILED_ALIGN/"
#    if Fail_Flag == 1:
#            if not os.path.exists(dst):
#                shutil.move(src, dst)
#            os.system("rm -rf " + src)
#            print("!!! Removed folde: " , src)
#            if os.path.exists(MSA_OUTPUT_PATH + "/FAILED_ALIGN/" + genusName):
#                f_out.write(genusName + "," + "FAILED ALIGNMENT")
#            else:
#                f_out.write(genusName + "," + "UNKNOWN STATUS")

#        if os.path.exists(src):
#            if os.listdir(src) == []:
#                print("Empty Folder: " + src)
#                os.system("rm -rf " + src)
#        if not os.path.exists(src):
#            if not os.path.exists(MSA_OUTPUT_PATH + "/FAILED_ALIGN/" + genusName):
#                f_out.write(genusName + "," + "No DIR for this genus!!!")

    # Number os Species with seqs Vs. Number of species included in alignment:
    #f_out.write(", SpeciesNum_Seqs = %d" % SpeciesNum_Seqs)
    #f_out.write(", SpeciesNum_Align = %d" % SpeciesNum_Align)
    f_out.write("\n")
    Fail_Flag=0
    # Chromevol analysis
    #/groups/itay_mayrose/michaldrori/MSA_JmodelTree/output/Hedychium/Hedychium_Chromevol/chromevol_out

    pathName_ploidy=(MSA_OUTPUT_PATH + "/" + genusName + "/" + genusName + "_Chromevol/chromevol_out/ploidy.csv")
    if os.path.exists(pathName_ploidy):
        Read_CSV(pathName_ploidy, f_out_chromevol)
        f_out_chromevol.write("\n")


#print("\n".join(list_passed_Align))
f_passed_genera.write("\n".join(list_passed_Align))
print("====================================\n")
print("No Species Found: %d" % NoSpecies_Cnt)
print("Less_then5Species_Cnt: %d" % Less_then5Species_Cnt)
print("NoOutgroup_Cnt: %d " % NoOutgroup_Cnt)
print("NoData_inClusters_Cnt: %d" % NoData_inClusters_Cnt)
print("OrthoMcl_issue_Cnt: %d" % OrthoMcl_issue_Cnt)
print("Passed_Alignment: %d" % Passed_Alignment)
print("StillRunning_Cnt: %d" % StillRunning_Cnt)
print("====================================\n")
print("Total_run_Cnt: %d" % Total_run_Cnt)

sys.exit(0)


with open("./PathList.txt") as f:
    PathList = f.readlines()
for pathName in PathList:
    pathName=pathName.rstrip('\n')
    Read_CSV(pathName, f_out_chromevol)
    f_out_chromevol.write("\n")


try:
    WorkOnPathList ("./PathList.txt", f_out_chromevol)
except:
    f_out_chromevol.write("No Such file: PathList.txt, in current Directory!!!")


print("%f seconds" % (time.time() - start_time))


##list_names = ['michal','Drori','michal2','Drori2','michal3','Drori3','michal4','Drori4','michal5','Drori5','michal6','Drori6']
##chank_species_lists=[list_names[x:x+2] for x in range(0, len(list_names), 2)]
##print(list_names)
##for name_chank in chank_species_lists:
##    print(name_chank)
##
##sys.exit(0)
##
##i=0
##