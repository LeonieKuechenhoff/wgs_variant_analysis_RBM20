
# Author: Leonie Küchenhoff
# Date: Sep 19th, 2022
# Script to save filtered variants after quality filtering steps


import numpy as np
import pandas as pd
import os


# define basedir and paths to files
basedir = '/g/steinmetz/project/leonie_crispr/03_data/01_heartproject/snakemake_vcf/'
os.chdir(basedir)
known_var_dir = 'known_variants_vcfs/overlap/'
paths = ['merged/txt_files/HLT279.merged.remdup.txt',
        'merged/txt_files/HLT282.merged.remdup.txt',
        'merged/txt_files/HLT450.merged.remdup.txt']
names = ['279', '282', '450']

# read in files and save as dict with sample nr. as name
file_dict = {}
for i in zip(paths, names):
    df = pd.read_csv(i[0], delimiter = ',')
    file_dict[i[1]] = df

# for each file, add info wether variant was annotated before
for name in names:
    file = file_dict[name]
    #read in files with known variants
    dbsnp = pd.read_csv(f'{known_var_dir}/{name}_dbsnp_overlap.txt', delimiter = '_', names = ['chr', 'pos', 'ref', 'alt'])
    mgp = pd.read_csv(f'{known_var_dir}/{name}_mgp_overlap.txt', delimiter = '_', names = ['chr', 'pos', 'ref', 'alt'])
    dbsnp['known_dbsnp'] = True
    mgp['known_mgp'] = True

    # add info from known variants to variant lists
    file = pd.merge(file, dbsnp, how = 'left', on = ['chr', 'pos', 'ref', 'alt'])
    print(file.shape)
    file = pd.merge(file, mgp, how = 'left', on = ['chr', 'pos', 'ref', 'alt'])

    # drop duplicates that are caused by the merging step
    file = file.drop_duplicates()

    # alter dataframes
    # 1. pick one af from mutect or haplotype caller

    # first, pick ad value from mutect2
    file.loc[:,'ad_l'] = file.loc[:,'AD|l_mu']
    file.loc[:,'ad_h'] = file.loc[:,'AD|h_mu']
    file.loc[:,'ad_t'] = file.loc[:,'AD|t_mu']
    # replace missing ad values from mutect with those from hc
    file.loc[(file['ad_l'] == '.')|(file['ad_l'] == './.'),'ad_l'] = file.loc[(file['ad_l'] == '.')|(file['ad_l'] == './.'),'AD|l_hc']
    file.loc[(file['ad_h'] == '.')|(file['ad_h'] == './.'),'ad_h'] = file.loc[(file['ad_h'] == '.')|(file['ad_h'] == './.'),'AD|h_hc']
    file.loc[(file['ad_t'] == '.')|(file['ad_t'] == './.'),'ad_t'] = file.loc[(file['ad_t'] == '.')|(file['ad_t'] == './.'),'AD|t_hc']

    # 2. correct gt of mutect with 0 alternative reads to 0,0
    # First, prepare data to be able alter it
    mutect_ad = file[['AD|h_mu','AD|l_mu','AD|t_mu']]
    mutect_ad[mutect_ad == '.'] = '0,0'
    mutect_ad[mutect_ad == './.'] = '0,0'
    allel1_ad = mutect_ad.applymap(lambda x: int(x.split(',')[0])).to_numpy()
    allel2_ad = mutect_ad.applymap(lambda x: int(x.split(',')[1])).to_numpy()
    # gt in variants with no alternative allele will be changed to 0/0 
    # (mutect2 gives all samples a gt 0f 0/1 of at least one sample contains alt allele)
    zero_alt = allel2_ad == 0

    mouse_gt = file[['GT|h_hc','GT|l_hc','GT|t_hc','GT|h_mu','GT|l_mu','GT|t_mu', 'GT|h_ls','GT|l_ls','GT|t_ls']]
    allel1 = mouse_gt.applymap(lambda x: str(x[0])).to_numpy()
    allel2 = mouse_gt.applymap(lambda x: str(x[2])).to_numpy()

    # correct gt of mutect 2
    allel1[:,3:6][zero_alt] = '0'
    allel2[:,3:6][zero_alt] = '0'
    
    newgt = np.array([x1 +'/'+ x2 for x1,x2 in zip(allel1,allel2)])
    
    file.loc[:,['GT|h_hc','GT|l_hc','GT|t_hc','GT|h_mu','GT|l_mu','GT|t_mu', 'GT|h_ls','GT|l_ls','GT|t_ls']] = newgt
    # save file
    file.to_csv(f'merged/txt_files/ad_hc_mu/HLT{name}.correctedgt.dbinfo.txt', sep = '\t', index = False)
    file_dict[name] = file


def get_sets(df_gt):
    '''
    function to get info where gt is  0/1 or 1/1 in the shape of a binary mask
    returns binary mask and indices of locations where variant is present
    '''
    # separate info from GT into two arrays (before: '0,1' - after: allel1 = [0], allel2 = [1])
    allel1 = df_gt.applymap(lambda x: x[0]).to_numpy()
    allel2 = df_gt.applymap(lambda x: x[2]).to_numpy()

    # replace '.' with -1 to ba able to use integers
    allel1[allel1 == '.'] = -1
    allel2[allel2 == '.'] = -1

    # transform arrays of strings to integer arrays
    allel1 = allel1.astype(int)
    allel2 = allel2.astype(int)

    # find all gt that are eiher 1/0 or 0/1 or 1/1 and create binary mask
    gt = allel1 + allel2
    binary = gt > 0

    # create ascending list of length of table that will be later used for venn diagrams
    setlist = np.arange(0,len(gt))
    return binary, setlist

def get_af(df, col1 = 'ad_h', col2 = 'ad_l', col3 = 'ad_t'):
    mouse_ad = df[[col1, col2, col3]]
    allel1 = mouse_ad.applymap(lambda x: int(x.split(',')[0])).to_numpy()
    allel2 = mouse_ad.applymap(lambda x: int(x.split(',')[1])).to_numpy()
    af = allel2 / (allel1 + allel2)
    af = np.nan_to_num(af, 0)
    return af

venn_dict = {}
for i in names:
    # get GT information from complete list
    mouse_gt = file_dict[i][['GT|h_hc','GT|l_hc','GT|t_hc','GT|h_mu','GT|l_mu','GT|t_mu', 'GT|h_ls','GT|l_ls','GT|t_ls']]
    mouse_ad = file_dict[i][['AD|h_hc','AD|l_hc','AD|t_hc','AD|h_mu','AD|l_mu','AD|t_mu']]
    # get GT & AD information from complete list that exlcuded known variants
    mouse_unknown = file_dict[i][
        (file_dict[i]['known_dbsnp'].isna() & file_dict[i]['known_mgp'].isna())
        ]
    mouse_gt_unknown = mouse_unknown[['GT|h_hc','GT|l_hc','GT|t_hc','GT|h_mu','GT|l_mu','GT|t_mu', 'GT|h_ls','GT|l_ls','GT|t_ls']]
    # get info where gt is  0/1 or 1/1 in the shape of a binary mask      
    binary_all, setlist_all = get_sets(mouse_gt)
    binary_uk, setlist_uk = get_sets(mouse_gt_unknown)

    # summarize information in dictionary for easier access
    venn_dict[i] = [binary_all, setlist_all, binary_uk, setlist_uk]

# save file of variant caller overlap - 2 way overlap
filtered_read_count_alta = {}
for i in names:
    df = file_dict[i]  
    # reminder: venn dict = [binary_all, setlist_all, binary_uk, setlist_uk] 
    # separate binary arrays where heart gt is either 1/0 or 0/1 or 1/1 into lists
    setlist_all = venn_dict[i][1]
    heart_hc = venn_dict[i][0][:,0]
    heart_mu = venn_dict[i][0][:,3]
    heart_ls = venn_dict[i][0][:,6]

    liver_hc = venn_dict[i][0][:,1]
    liver_mu = venn_dict[i][0][:,4]
    liver_ls = venn_dict[i][0][:,7]

    tail_hc = venn_dict[i][0][:,2]
    tail_mu = venn_dict[i][0][:,5]
    tail_ls = venn_dict[i][0][:,8]
    
    df['heart_gt_sum'] = np.array(heart_mu, dtype = int) + np.array(heart_hc, dtype = int) + np.array(heart_ls, dtype = int)
    df['liver_gt_sum'] = np.array(liver_mu, dtype = int) + np.array(liver_hc, dtype = int) + np.array(liver_ls, dtype = int)
    df['tail_gt_sum'] = np.array(tail_mu, dtype = int) + np.array(tail_hc, dtype = int) + np.array(tail_ls, dtype = int)

    # get index of overlaps
    heart = setlist_all[(heart_mu & heart_hc) | (heart_mu & heart_ls) | (heart_ls & heart_hc)]
    liver = setlist_all[(liver_mu & liver_hc) | (liver_mu & liver_ls) | (liver_ls & liver_hc)]
    tail = setlist_all[(tail_mu & tail_hc) | (tail_mu & tail_ls) | (tail_ls & tail_hc)]

    # create df that only contain heart / liver / tail overlaps
    heart_df = df.iloc[heart,:]
    heart_df['heart'] = True
    
    liver_df = df.iloc[liver,:]
    liver_df['liver'] = True

    tail_df = df.iloc[tail,:]
    tail_df['tail'] = True

    # perfrom an outer join on dataframes (in two steps as function only takes two df)
    complete_df = pd.merge(heart_df, liver_df, on = df.columns.to_list(), how = 'outer')
    complete_df = pd.merge(complete_df, tail_df, on = df.columns.to_list(), how = 'outer').fillna(False)

    #Have at least 5 reads total per variant per tissue
    mouse_ad = complete_df[['ad_h','ad_l','ad_t']]
    allel1 = mouse_ad.applymap(lambda x: int(x.split(',')[0])).to_numpy()
    allel2 = mouse_ad.applymap(lambda x: int(x.split(',')[1])).to_numpy()
    summed = (allel1 + allel2)
    bool = summed >=5
    all_above = np.all(bool, axis = 1)
    filtered_df = complete_df[all_above]

    #Have at least 2 reads of alt alleles across tissues
    mouse_ad = filtered_df[['ad_h','ad_l','ad_t']]
    allel2 = mouse_ad.applymap(lambda x: int(x.split(',')[1])).to_numpy()

    keep_variants = allel2.sum(axis = 1) > 1
    filtered_df = filtered_df[keep_variants]
    filtered_read_count_alta[i] = filtered_df
    filtered_df.to_csv(f'merged/txt_files/ad_hc_mu/HLT{i}.overlap_three_readabove4_above2_1ed.txt', sep = '\t', index = False)


# for each filter, divide data into somatic vs. germline variants
known_dict = {}
for i in names:
    # tissue specific variants that are unknown
    uk = filtered_read_count_alta[i][(filtered_read_count_alta[i]['known_dbsnp'] == False) & (filtered_read_count_alta[i]['known_mgp'] == False)]

    heart_only = uk[(uk['heart'] == True)& (uk['liver'] != True)& (uk['tail'] != True)]
    af = get_af(heart_only) 
    heart_spec = heart_only.iloc[np.where((af[:,1]==0) &(af[:,2]==0))]

    liver_only = uk[(uk['heart'] != True)& (uk['liver'] == True)& (uk['tail'] != True)]
    af = get_af(liver_only) 
    liver_spec = liver_only.iloc[np.where((af[:,0]==0) &(af[:,2]==0))]

    tail_only = uk[(uk['heart'] != True)& (uk['liver'] != True)& (uk['tail'] == True)]
    af = get_af(tail_only) 
    tail_spec = tail_only.iloc[np.where((af[:,0]==0) &(af[:,1]==0))]
    
    specific = pd.concat([heart_spec, liver_spec, tail_spec])

    # common variants from all three tissues
    hlt = filtered_read_count_alta[i][(filtered_read_count_alta[i]['heart'] == True)& (filtered_read_count_alta[i]['liver'] == True)& (filtered_read_count_alta[i]['tail'] == True)]

    specific.to_csv(f'merged/txt_files/ad_hc_mu/HLT{i}.specific.txt', sep = '\t', index = False)
    hlt.to_csv(f'merged/txt_files/ad_hc_mu/HLT{i}.tisoverlap.txt', sep = '\t', index = False)