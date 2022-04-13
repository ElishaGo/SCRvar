# TODO: Document, put into function, take parameters out, create simple argparser
# Important to document coordinate system 

import pandas as pd

bed_file_path = '/home/eligol/Documents/01_WIS/scrarevar/data/DB_edit_snp/human/edit_rep.bed'
coordinate_to_filter_path = '/home/eligol/Documents/01_WIS/scrarevar/data/DB_edit_snp/human/coordinate_to_filter_edit_rep_bed.txt'

# load bed file
bed_file = pd.read_csv(bed_file_path, sep='\t')
# load coordination to filter
with open(coordinate_to_filter_path, 'r') as f:
    coor_to_filt = f.readlines()

# process the coordinates to be in one lise, and filter unnecessary strings
coor_to_filt = [' '.join(x).replace('\n', '').replace('>','').replace('(-)','').replace('(+)','') for x in zip(coor_to_filt[0::2], coor_to_filt[1::2])]
# keep only not 'A' references in filter list
coor_to_filt = [c for c in coor_to_filt if c[-1] not in ['a', 'A']]

# convert coordinates to dataframe
df_coor_to_filt = pd.DataFrame(columns=bed_file.columns[:4])
df_coor_to_filt[['#chrom', 'start']] = pd.Series(coor_to_filt).str.split(':', expand=True)
df_coor_to_filt[['start', 'end']] = df_coor_to_filt['start'].str.split('-', expand=True)
df_coor_to_filt[['end', 'Ref']] = df_coor_to_filt['end'].str.split(' ', expand=True)
df_coor_to_filt[['start', 'end']] = df_coor_to_filt[['start', 'end']].astype('int')

# find indices to drop
idx_to_drop = pd.merge(bed_file.reset_index(), df_coor_to_filt.drop('Ref', axis=1), how='inner').set_index('index').index
# drop indices
bed_file_filtered = bed_file.drop(idx_to_drop).set_index('#chrom')

bed_file_filtered.to_csv('/home/eligol/Documents/01_WIS/scrarevar/data/DB_edit_snp/human/edit_rep_filtered.bed', sep='\t')


