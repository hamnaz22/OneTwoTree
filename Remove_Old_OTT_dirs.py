from glob import glob
import os
import subprocess

month_dict={'Jan':1, 'Feb':2, 'Mar':3, 'Apr':4, 'May':5, 'Jun':6, 'Jul':7, 'Aug':8, 'Sep':9, 'Oct':10, 'Nov':11, 'Dec':12}
Exclude_list=['1490596072', '1489627082', '1490596049', '1489087958', '1489627039', '1490595938', '1489914703', '1489314188', '1489407001']
CURRENT_MONTH = 6

for dir_OTT in glob("/bioseq/data/results/oneTwoTree/*/"):
    flag_skip = 'No'
    for exclude_dir in Exclude_list:
        if exclude_dir in dir_OTT:
            print("SKIP %s" % exclude_dir)
            flag_skip = 'Yes'
    if flag_skip != 'Yes':
        #data_split=((os.system("ls -ldc %s" %dir_OTT)).split(' '))
        data_split=(subprocess.check_output("ls -ldc %s" %dir_OTT, shell=True)).decode("utf-8")
        data_split_=data_split.split(' ')
        #check if older than one month:
        if month_dict[data_split_[5]] + 1 < CURRENT_MONTH:
            print("Remove dir: %s" % dir_OTT)
            subprocess.check_output("rm -rf %s" %dir_OTT, shell=True)
        #print(data_split)
       #
       #print(data_split_[4])
       #print(data_split_[5])
       #print(data_split_[6])

