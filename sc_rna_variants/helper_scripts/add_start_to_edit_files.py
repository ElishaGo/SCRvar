import pandas as pd

# path to editing DB file
rep_db_path = '/home/eligol/Documents/01_WIS/scrarevar/data/DB_edit_snp/human/TABLE1_hg38.txt'

df = pd.read_csv(rep_db_path, sep = '\t')
df = df.iloc[:,:-3]

# create an end column
df['start'] = df['Position'] - 1
#df['end'] = df['Position'] + 1  #change end

# add the new column and reorder the table
columns = df.columns.tolist()
reorder = [0,-1,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21]
#reorder = [0,1,-1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21]  #change end
columns = [columns[i] for i in reorder]

# apply the new column order
df = df[columns]
#df.rename(columns = {'Region':'#chrom', '#Region':'#chrom', 'Position': 'start'}, inplace=True) #change end
df.rename(columns = {'Region':'#chrom', '#Region':'#chrom', 'Position': 'end'}, inplace=True)

# separate the DB into repetitive and non_repetitive tables
df_rep = df[df['type'] == 'REP']
df_nonrep = df[df['type'] == 'NONREP']

df_rep.to_csv('/home/eligol/Documents/01_WIS/scrarevar/data/DB_edit_snp/human_end_increase_1/edit_rep_end_add_1.bed', sep = '\t',index = False)
df_nonrep.to_csv('/home/eligol/Documents/01_WIS/scrarevar/data/DB_edit_snp/human_end_increase_1/edit_nonrep_end_add_1.bed', sep = '\t',index = False)
