import sqlite3
import pandas as pd
import os
from ploidbCommon import *


DB_Dir = ploidb_config['general']['DB_DIR']
seqs_db_path = DB_Dir+ '/GenbankDB_grp.db'
fasta_file_for_blast = DB_Dir + '/BlastDB/Blast_fst.fst'
grp_names = ['mam','inv','pln'] #sys.argv[3]


Accession_Id_df = dict()
SeqsData_df = dict()

print('** collect data from db')# connect to DB
conn_db = sqlite3.connect(seqs_db_path)
curs_db = conn_db.cursor()
count = 1
query=''
grp_num = len(grp_names)
for grp in grp_names:
    print('Collect data for group %s frm genabnk database:' %grp)
    query = 'select TaxonId,Accesion,organism,sequenceData from Genbank_seqs_%s' %grp

    print(query)
    curs_db.execute(query)
    rows_db = curs_db.fetchall()
    print('** Done select db')

    # Insert data into csv file:
    headers_names_list = [ "TaxonId","Accesion","organism","sequenceData"]
    df = pd.DataFrame(rows_db, columns=headers_names_list)

    TaxonId_df = df['TaxonId']
    Accession_Id_df = df['Accesion']
    organism_df = df['organism']
    SeqsData_df = df['sequenceData']

    df_count = 0
    df_length = len(Accession_Id_df)

    with open(fasta_file_for_blast, 'a') as f_fasta_db:
        while df_count < df_length:
            if df_count % 100000 == 0:
                print(df_count)
            organism_underscore=str(organism_df[df_count]).replace(' ','_')
            description_str = '>' + str(TaxonId_df[df_count])+'|'+ str(Accession_Id_df[df_count])+'|'+ organism_underscore
            f_fasta_db.write(description_str+'\n')
            f_fasta_db.write(str(SeqsData_df[df_count]+'\n'))
            df_count+=1

conn_db.close()

print('Done createing fasta file, now create the blast database...')
os.system('cd %s' % DB_Dir)
os.system('makeblastdb -in Blast_fst.fst -parse_seqids -dbtype nucl')
print('Done !!!')
