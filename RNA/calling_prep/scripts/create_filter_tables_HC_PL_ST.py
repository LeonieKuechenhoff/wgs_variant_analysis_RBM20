# Author: Leonie Küchenhoff
# Date: Sep 19th, 2022
# Script to save filtered variants after quality filtering steps

import numpy as np
import pandas as pd
import argparse

def get_sets(df_gt):
    '''
    function to get info where gt is  0/1 or 1/1 in the shape of a binary mask
    returns binary mask and indices of locations where variant is present
    '''
    df_gt = df_gt.replace('.','./.')
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
    setlist = np.arange(0,gt.shape[0])
    return binary, setlist

def get_af(df, col1 = 'ad_h', col2 = 'ad_l'):
    mouse_ad = df[[col1, col2]]
    allel1 = mouse_ad.applymap(lambda x: int(x.split(',')[0])).to_numpy()
    allel2 = mouse_ad.applymap(lambda x: int(x.split(',')[1])).to_numpy()
    af = allel2 / (allel1 + allel2)
    af = np.nan_to_num(af, 0)
    return af

def duplicate(df):
    '''
    There are very few duplicate variants present in the files which have different annotations in the two duplicating lines. 
    This function will edit these rows in a way that they will be changed to the genotype containing the alt allel.
    For example, if one row has GT = 0/0 and one row has 1/0, the genotype 1/0 will be used for both rows. 
    '''
    # first duplicated row
    first = df.duplicated(subset = ['chr', 'pos', 'ref', 'alt'], keep = 'first')
    # second duplicated row
    second = df.duplicated(subset = ['chr', 'pos', 'ref', 'alt'], keep = 'last')

    subset = df.iloc[np.where(first == True)].iloc[:, -6:]

    allel1 = subset.applymap(lambda x: x[0], na_action = 'ignore').to_numpy()
    allel2 = subset.applymap(lambda x: x[2], na_action = 'ignore').to_numpy()
    # replace '.' with 0 to ba able to use integers
    allel1[allel1 == '.'] = 0
    allel2[allel2 == '.'] = 0
    # transform arrays of strings to integer arrays
    allel1 = allel1.astype(int)
    allel2 = allel2.astype(int)

    subset = df.iloc[np.where(second == True)].iloc[:, -6:]
    subset[subset == '.'] = '0,0'
    allel3 = subset.applymap(lambda x: x[0]).to_numpy()
    allel4 = subset.applymap(lambda x: x[2]).to_numpy()
    # replace '.' with 0 to ba able to use integers
    allel3[allel3 == '.'] = 0
    allel4[allel4 == '.'] = 0
    # transform arrays of strings to integer arrays
    allel3 = allel3.astype(int)
    allel4 = allel4.astype(int)

    # add alleles of duplicate rows 
    # in the upcoming analysis, it makes no difference wether gt is 1/1, 0/1, 1/0, or 2/0,.. therefore they are simply added in this step
    new_allel1 = allel1 + allel3
    new_allel2 = allel2 + allel4
    # go back to gt annotation
    gt = np.char.add(new_allel1.astype('U1'), '/')
    gt = np.char.add(gt, new_allel2.astype('U1'))
    # replace gt in original df
    df.iloc[np.where(first == True)[0],-6:] = gt
    df.iloc[np.where(second == True)[0],-6:] = gt
    # drop duplicates
    df = df.drop_duplicates(subset = ['chr', 'pos', 'ref', 'alt'])

    return df

def read_three_varcaller(strelka_path, plat_path, hc_path, name):
    '''
    Function to read in data from three different varcallers.
    As all varcallers have slightly different outputs, there are read separatley
    '''

    # strelka output requires exclusion of falsely formatted GT entires
    strelka = pd.read_csv(strelka_path, delimiter = '\t', names = ['chr','pos','ref','alt','AD|h_st','AD|l_st','GT|h_st','GT|l_st'], header = 0)
    strelka = strelka[strelka['chr'].isin([
        'chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10',
        'chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19',
        'chrX', 'chrY'])]
    strelka = strelka[(strelka['GT|h_st']!= '1') & (strelka['GT|h_st']!= '0')]
    strelka = strelka[(strelka['GT|l_st']!= '1') & (strelka['GT|l_st']!= '0')]
    
    strelka = strelka[['chr','pos','ref','alt','GT|h_st','GT|l_st','AD|h_st','AD|l_st']]

    # platypus output requires calculation of AD
    plat = pd.read_csv(plat_path, delimiter = '\t', names = ['chr','pos','ref','alt','NV_h','NV_l','NR_h', 'NR_l','GT|h_pl','GT|l_pl'], header = 0)
    plat = plat[plat['chr'].isin([
        'chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10',
        'chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19',
        'chrX', 'chrY'])]
    # only pick monoallelic variants
    plat = plat.astype({'NR_h':str, 'NR_l':str})
    plat['multi']=plat['NR_h'].str.count(r'(,)').to_numpy() + plat['NR_l'].str.count(r'(,)').to_numpy()
    plat = plat[plat['multi'] == 0]
    plat = plat.astype({'NR_h': np.int32, 'NV_h': np.int32, 'NR_l': np.int32, 'NV_l': np.int32})
    ref_h = plat['NR_h']-plat['NV_h']
    ref_l = plat['NR_l']-plat['NV_l']
    alt_h = plat['NV_h']
    alt_l = plat['NV_l']
    plat['AD|h_pl'] =list(map("{0[0]},{0[1]}".format, zip(ref_h, alt_h)))
    plat['AD|l_pl'] =list(map("{0[0]},{0[1]}".format, zip(ref_l, alt_l)))
    
    plat = plat[['chr','pos','ref','alt','GT|h_pl','GT|l_pl','AD|h_pl','AD|l_pl']]
    

    # output haplotype caller
    hc = pd.read_csv(hc_path, delimiter = '\t', names = ['chr','pos','ref','alt','AD|h_hc','AD|l_hc','GT|h_hc','GT|l_hc'], header = 0)
    hc = hc[['chr','pos','ref','alt','GT|h_hc','GT|l_hc','AD|h_hc','AD|l_hc']]
    hc = hc[hc['chr'].isin([
        'chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10',
        'chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19',
        'chrX', 'chrY'])]
    merged = hc.merge(plat,  on = ['chr', 'pos', 'ref', 'alt'], how = 'outer')
    merged = merged.merge(strelka,  on = ['chr', 'pos', 'ref', 'alt'], how = 'outer').fillna('0,0')

    print(merged.head())
    return merged


def main(names, knownvar, input, output, dbsnp, mgp):
    
    dbsnp_df = pd.read_csv(dbsnp, delimiter='\t', header = None, names = ['chr', 'pos', 'ref', 'alt'])
    mgp_df = pd.read_csv(mgp, delimiter='\t', header = None, names = ['chr', 'pos', 'ref', 'alt'])
    dbsnp_df['known_dbsnp'] = True
    mgp_df['known_mgp'] = True

    for name in names:
        print(name)
        strelka = f'{input}/{name}.st.txt'
        plat = f'{input}/{name}.plat.txt'
        hc =f'{input}/{name}.hc.txt'
        merged = read_three_varcaller(strelka, plat, hc, name)
        sample_df_dup = duplicate(merged)

        # get intersectionss of database entries
        # add info from known variants to variant lists
        file = pd.merge(sample_df_dup, dbsnp_df, how = 'left', on = ['chr', 'pos', 'ref', 'alt'])
        print(file.head())
        file = pd.merge(file, mgp_df, how = 'left', on = ['chr', 'pos', 'ref', 'alt'])
        print(file.head())

        file = file.drop_duplicates()
        
        # alter dataframes
        # pick one af from platypus or haplotype caller

        # first, pick ad value from platypus
        file.loc[:,'ad_l'] = file.loc[:,'AD|l_pl']
        file.loc[:,'ad_h'] = file.loc[:,'AD|h_pl']
        # replace missing ad values from platypus with those from hc
        file.loc[(file['ad_l'] == '0/0'),'ad_l'] = file.loc[(file['ad_l'] == '0/0'),'AD|l_hc']
        file.loc[(file['ad_l'] == '0/0'),'ad_h'] = file.loc[(file['ad_h'] == '0/0'),'AD|h_hc']

        file.to_csv(f'{output}/{name}.correctedgt.dbinfo.txt', sep = '\t', index = False)
        

        # get GT information from complete list
        mouse_gt = file[['GT|h_hc','GT|l_hc','GT|h_pl','GT|l_pl', 'GT|h_st','GT|l_st']]

        # get GT information from complete list that exlcuded known variants
        mouse_unknown = file[
            (file['known_dbsnp'].isna() & file['known_mgp'].isna())
            ]
        mouse_gt_unknown = mouse_unknown[['GT|h_hc','GT|l_hc','GT|h_pl','GT|l_pl', 'GT|h_st','GT|l_st']]

        # get info where gt is  0/1 or 1/1 in the shape of a binary mask      
        binary_all, setlist_all = get_sets(mouse_gt)
        binary_uk, setlist_uk = get_sets(mouse_gt_unknown)

        # summarize information in dictionary for easier access
        venn = [binary_all, setlist_all, binary_uk, setlist_uk]

        # save file of variant caller overlap - 2 way overlap

        # reminder: venn dict = [binary_all, setlist_all, binary_uk, setlist_uk] 
        # separate binary arrays where heart gt is either 1/0 or 0/1 or 1/1 into lists
        setlist_all = venn[1]
        heart_hc = venn[0][:,0]
        heart_pl = venn[0][:,2]
        heart_st =venn[0][:,4]

        liver_hc = venn[0][:,1]
        liver_pl = venn[0][:,3]
        liver_st = venn[0][:,5]

        file['heart_gt_sum'] = np.array(heart_pl, dtype = int) + np.array(heart_hc, dtype = int) + np.array(heart_st, dtype = int)
        file['liver_gt_sum'] = np.array(liver_pl, dtype = int) + np.array(liver_hc, dtype = int) + np.array(liver_st, dtype = int)

        # get index of overlaps
        heart = setlist_all[(heart_pl & heart_hc) | (heart_pl & heart_st) | (heart_st & heart_hc)]
        liver = setlist_all[(liver_pl & liver_hc) | (liver_pl & liver_st) | (liver_st & liver_hc)]

        # create df that only contain heart / liver overlaps
        heart_df = file.iloc[heart,:]
        heart_df['heart'] = True
        
        liver_df = file.iloc[liver,:]
        liver_df['liver'] = True


        # perfrom an outer join on dataframes
        complete_df = pd.merge(heart_df, liver_df, on = file.columns.to_list(), how = 'outer').fillna(False)
        #complete_df.to_csv(f'{output}/{name}.overlap_two.txt', sep = '\t', index = False)



        #Have at least 5 reads total per variant 
        mouse_ad = complete_df[['ad_h','ad_l']]
        allel1 = mouse_ad.applymap(lambda x: int(x.split(',')[0])).to_numpy()
        allel2 = mouse_ad.applymap(lambda x: int(x.split(',')[1])).to_numpy()
        summed = (allel1 + allel2)
        bool = summed >=5

        bool2 = complete_df[['heart','liver']].to_numpy()
        #combine info from GT and number of reads
        present = [bool*bool2]
        # only mark variants as present if there were at least 5 reads measured on that position
        complete_df['heart']= present[0][:,0]
        complete_df['liver']= present[0][:,1]

        #Have at least 2 reads of alt alleles across tissues
        mouse_ad = complete_df[['ad_h','ad_l']]
        allel2 = mouse_ad.applymap(lambda x: int(x.split(',')[1])).to_numpy()

        keep_variants = allel2.sum(axis = 1) > 1
        filtered_df = complete_df[keep_variants]
        filtered_df.to_csv(f'{output}/{name}.overlap_two_ed.txt', sep = '\t', index = False)


        # tissue specific variants that are unknown
        uk = filtered_df[(filtered_df['known_dbsnp'] == False) & (filtered_df['known_mgp'] == False)]

        heart_only = uk[(uk['heart'] == True)& (uk['liver'] != True)]
        af = get_af(heart_only) 
        heart_spec = heart_only.iloc[np.where((af[:,1]==0))]

        liver_only = uk[(uk['heart'] != True)& (uk['liver'] == True)]
        af = get_af(liver_only) 
        liver_spec = liver_only.iloc[np.where((af[:,0]==0))]

        
        specific = pd.concat([heart_spec, liver_spec])

        # common variants from all three tissues
        hlt = filtered_df[(filtered_df['heart'] == True)& (filtered_df['liver'] == True)]

        specific.to_csv(f'{output}/{name}.specific.txt', sep = '\t', index = False)
        hlt.to_csv(f'{output}/{name}.tisoverlap.txt', sep = '\t', index = False)
    open(f'{output}/myfile.txt', "x")



if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='')
    parser.add_argument('-n', '--names', help = 'Specify sample names', type = str, nargs = "+")
    parser.add_argument('-k', '--knownvar', help='Specify paths to known variant file from dbsnp', type = str)
    parser.add_argument('-i', '--input', help='Specify paths to input directory. Files must be in format SAMPLE_NAME.VARCALLER.txt. Varcallers: hc/plat/st', type = str)
    parser.add_argument('-o', '--output', help='Specify dir to output', type = str,
                        default = 'merged/txt_files/ad_hc_mu')
    parser.add_argument('-d', '--dbsnp', help='Specify path to known variant file from dbsnp', type = str,
                        default = '/g/steinmetz/project/leonie_crispr/03_data/01_heartproject/known_variants_vcfs/normalized/sorted/dbsnp.txt')
    parser.add_argument('-m', '--mgp', help='Specify path to known variant file from mgp', type = str,
                        default = '/g/steinmetz/project/leonie_crispr/03_data/01_heartproject/known_variants_vcfs/normalized/sorted/mgp.txt')
    args = parser.parse_args()
    #print(len(args.known))
    main(**vars(args))

#python create_filter_tables_HC_PL_ST.py --names HL029_pbs_R HL033_nrch_R HL011_pbs HL012_nrch HL013_nrch HL014_nrch HL279_spry HL321_pbs HL333_pbs HL450_spry HL283_spry HL028_pbs_R HL030_nrch_R HL032_pbs_R HL036_nrch_R --input /g/steinmetz/project/leonie_crispr/03_data/02_rnaseq/snakemake/variant_caller_out/vartest/txt_files --output /g/steinmetz/project/leonie_crispr/03_data/02_rnaseq/snakemake/filtered_tables



#python create_filter_tables_HC_PL_ST.py --names HL030_nrch_R HL032_pbs_R --input /g/steinmetz/project/leonie_crispr/03_data/02_rnaseq/snakemake/variant_caller_out/vartest/txt_files --output /g/steinmetz/project/leonie_crispr/03_data/02_rnaseq/snakemake/filtered_tables


#python create_filter_tables_HC_PL_ST.py --names HL030_nrch_R --input /g/steinmetz/project/leonie_crispr/03_data/02_rnaseq/snakemake/variant_caller_out/vartest/txt_files --output /g/steinmetz/project/leonie_crispr/temp