#%%
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
from pandas_plink import read_rel

#%% Importing functions 

# Dictionnary for changing well id
WellName_to_WE00No_dict = {'A001': 'WE00001', 'A002': 'WE00002', 'A003': 'WE00003', 'A004': 'WE00004', 'A005': 'WE00005', 'A006': 'WE00006', 'A007': 'WE00007', 'A008': 'WE00008', 'A009': 'WE00009', 'A010': 'WE00010', 'A011': 'WE00011', 'A012': 'WE00012', 'B012': 'WE00013', 'B011': 'WE00014', 'B010': 'WE00015', 'B009': 'WE00016', 'B008': 'WE00017', 'B007': 'WE00018', 'B006': 'WE00019', 'B005': 'WE00020', 'B004': 'WE00021', 'B003': 'WE00022', 'B002': 'WE00023', 'B001': 'WE00024', 'C001': 'WE00025', 'C002': 'WE00026', 'C003': 'WE00027', 'C004': 'WE00028', 'C005': 'WE00029', 'C006': 'WE00030', 'C007': 'WE00031', 'C008': 'WE00032', 'C009': 'WE00033', 'C010': 'WE00034', 'C011': 'WE00035', 'C012': 'WE00036', 'D012': 'WE00037', 'D011': 'WE00038', 'D010': 'WE00039', 'D009': 'WE00040', 'D008': 'WE00041', 'D007': 'WE00042', 'D006': 'WE00043', 'D005': 'WE00044', 'D004': 'WE00045', 'D003': 'WE00046', 'D002': 'WE00047', 'D001': 'WE00048', 'E001': 'WE00049', 'E002': 'WE00050', 'E003': 'WE00051', 'E004': 'WE00052', 'E005': 'WE00053', 'E006': 'WE00054', 'E007': 'WE00055', 'E008': 'WE00056', 'E009': 'WE00057', 'E010': 'WE00058', 'E011': 'WE00059', 'E012': 'WE00060', 'F012': 'WE00061', 'F011': 'WE00062', 'F010': 'WE00063', 'F009': 'WE00064', 'F008': 'WE00065', 'F007': 'WE00066', 'F006': 'WE00067', 'F005': 'WE00068', 'F004': 'WE00069', 'F003': 'WE00070', 'F002': 'WE00071', 'F001': 'WE00072', 'G001': 'WE00073', 'G002': 'WE00074', 'G003': 'WE00075', 'G004': 'WE00076', 'G005': 'WE00077', 'G006': 'WE00078', 'G007': 'WE00079', 'G008': 'WE00080', 'G009': 'WE00081', 'G010': 'WE00082', 'G011': 'WE00083', 'G012': 'WE00084', 'H012': 'WE00085', 'H011': 'WE00086', 'H010': 'WE00087', 'H009': 'WE00088', 'H008': 'WE00089', 'H007': 'WE00090', 'H006': 'WE00091', 'H005': 'WE00092', 'H004': 'WE00093', 'H003': 'WE00094', 'H002': 'WE00095', 'H001': 'WE00096'}
WE00No_to_WellName_dict = {'WE00001': 'A001', 'WE00002': 'A002', 'WE00003': 'A003', 'WE00004': 'A004', 'WE00005': 'A005', 'WE00006': 'A006', 'WE00007': 'A007', 'WE00008': 'A008', 'WE00009': 'A009', 'WE00010': 'A010', 'WE00011': 'A011', 'WE00012': 'A012', 'WE00013': 'B012', 'WE00014': 'B011', 'WE00015': 'B010', 'WE00016': 'B009', 'WE00017': 'B008', 'WE00018': 'B007', 'WE00019': 'B006', 'WE00020': 'B005', 'WE00021': 'B004', 'WE00022': 'B003', 'WE00023': 'B002', 'WE00024': 'B001', 'WE00025': 'C001', 'WE00026': 'C002', 'WE00027': 'C003', 'WE00028': 'C004', 'WE00029': 'C005', 'WE00030': 'C006', 'WE00031': 'C007', 'WE00032': 'C008', 'WE00033': 'C009', 'WE00034': 'C010', 'WE00035': 'C011', 'WE00036': 'C012', 'WE00037': 'D012', 'WE00038': 'D011', 'WE00039': 'D010', 'WE00040': 'D009', 'WE00041': 'D008', 'WE00042': 'D007', 'WE00043': 'D006', 'WE00044': 'D005', 'WE00045': 'D004', 'WE00046': 'D003', 'WE00047': 'D002', 'WE00048': 'D001', 'WE00049': 'E001', 'WE00050': 'E002', 'WE00051': 'E003', 'WE00052': 'E004', 'WE00053': 'E005', 'WE00054': 'E006', 'WE00055': 'E007', 'WE00056': 'E008', 'WE00057': 'E009', 'WE00058': 'E010', 'WE00059': 'E011', 'WE00060': 'E012', 'WE00061': 'F012', 'WE00062': 'F011', 'WE00063': 'F010', 'WE00064': 'F009', 'WE00065': 'F008', 'WE00066': 'F007', 'WE00067': 'F006', 'WE00068': 'F005', 'WE00069': 'F004', 'WE00070': 'F003', 'WE00071': 'F002', 'WE00072': 'F001', 'WE00073': 'G001', 'WE00074': 'G002', 'WE00075': 'G003', 'WE00076': 'G004', 'WE00077': 'G005', 'WE00078': 'G006', 'WE00079': 'G007', 'WE00080': 'G008', 'WE00081': 'G009', 'WE00082': 'G010', 'WE00083': 'G011', 'WE00084': 'G012', 'WE00085': 'H012', 'WE00086': 'H011', 'WE00087': 'H010', 'WE00088': 'H009', 'WE00089': 'H008', 'WE00090': 'H007', 'WE00091': 'H006', 'WE00092': 'H005', 'WE00093': 'H004', 'WE00094': 'H003', 'WE00095': 'H002', 'WE00096': 'H001'}

def relationship_matrix(df, legend='', serie_legend=None, legend_palette='husl', name='clustermap.png'):
    if len(serie_legend) != 0:
        lut = dict(zip(serie_legend.astype(str).unique(), sns.color_palette(legend_palette, serie_legend.astype(str).nunique())))
        mapping_dict = serie_legend.to_dict()
        row_colors = df.index.map(mapping_dict).astype(str).map(lut)
        col_colors = df.columns.map(mapping_dict).astype(str).map(lut)
    else:
        row_colors=None
        col_colors=None

    clust= sns.clustermap(df, method='ward', row_colors=row_colors, col_colors=col_colors, cmap='viridis')
    if legend:
        for label in serie_legend.astype(str).unique():
            clust.ax_col_dendrogram.bar(0, 0,
                                        color=lut[label],
                                        label=label,
                                        linewidth=0)
            clust.ax_col_dendrogram.legend(title=legend,
                                            loc="right",
                                            ncol=1,
                                            bbox_to_anchor=(1.6, -1.6))
    clust.savefig(name)


def main(matrix, matrix_ids=None, info=None):
    # Loading relationship matrix and its ids
    if not matrix_ids: 
        matrix_ids = matrix.replace('.bin', '.id')
    print(matrix_ids)
    arr = np.fromfile(matrix)
    # arr = read_rel(matrix, matrix_ids)
    print(arr)
    ids = pd.read_table(matrix_ids)
    # df_rel = pd.DataFrame(arr, columns=ids.iloc[:,0], index=ids.iloc[:,0])
    df = pd.DataFrame(arr)

    
    # Table to link IDs
    # ids['EMBL_file_ID'] = 'AA' + ids['#IID'].str.split('AA', expand=True)[2]
    # ids['id'] = ids['EMBL_file_ID'].map(id_dict)
    
    print(df.head())
    
    # meta = pd.read_csv('/home/fanny/Documents/Work/HeartRate_Tox/data/normalized_format/EMBL_sample_linkage.csv')
    # id_dict = meta.set_index('EMBL_file_ID')['EMBL_ID'].to_dict()
    # # %%
    
    # # Plotting by strain ID
    # relationship_matrix(df_rel, legend='strain_ID', serie_legend=meta.set_index('EMBL_ID')['strain_ID'], legend_palette='husl', name='strain_ID_clustermap.png')
    
    # # %%
    
    # # Plotting by strain ID
    # relationship_matrix(df_rel, legend='EMBL_seq_batch', serie_legend=meta.set_index('EMBL_ID')['EMBL_seq_batch'], legend_palette='husl', name='EMBL_seq_batch_clustermap.png')

    # # %% Loading coverage
    
    # cov = pd.read_csv('/home/fanny/Documents/Work/HeartRate_Tox/data/normalized_format/sequencing_data.csv')
    # cov_dict = cov.set_index('file_ID')['cov'].to_dict()
    # cov['cov'] = 1/ cov['cov']
    
    # cov['bin'] = pd.cut(cov['cov'], bins=[0,0.01,0.05,0.1,0.5,1,2,5,10,20,50,100,500,1000]).astype(str)
    
    # relationship_matrix(df_rel, legend='bin', serie_legend=cov.sort_values('cov').set_index('ID')['bin'], legend_palette='hls', name='coverage_clustermap.png')
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Script to plot the relationship matrix')
    parser.add_argument('matrix', type=str, help='File containing the relationship matrix issued from plink')
    parser.add_argument('-i', '--matrix_ids', type=str, help='File containing the ids')
    args = parser.parse_args()
    main(**vars(args))
