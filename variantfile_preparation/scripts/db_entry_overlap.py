# Author: Leonie Küchenhoff
# Date: Sep 19th, 2022
# Script to extract overlapping sites with sites annotated in databases. 
# Trims down database entries to those that overlap to avoid reading in all database entries in the following analysis.
# Also corrects genotype of duplicate alleles.

import pandas as pd
import argparse
import numpy as np



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
    subset = df.iloc[np.where(first == True)].iloc[:, -9:]
    allel1 = subset.applymap(lambda x: x[0], na_action = 'ignore').to_numpy()
    allel2 = subset.applymap(lambda x: x[2], na_action = 'ignore').to_numpy()

    # replace '.' with 0 to be able to use integers
    allel1[allel1 == '.'] = 0
    allel2[allel2 == '.'] = 0

    # transform arrays of strings to integer arrays
    allel1 = allel1.astype(int)
    allel2 = allel2.astype(int)
    subset = df.iloc[np.where(second == True)].iloc[:, -9:]
    allel3 = subset.applymap(lambda x: x[0]).to_numpy()
    allel4 = subset.applymap(lambda x: x[2]).to_numpy()

    # replace '.' with 0 to ba able to use integers
    allel3[allel3 == '.'] = 0
    allel4[allel4 == '.'] = 0

    # transform arrays of strings to integer arrays
    allel3 = allel3.astype(int)
    allel4 = allel4.astype(int)

    # add alleles of duplicate rows 
    # in the upcoming analysis, it makes no difference wether gt is 1/1, 0/1, 1/0, or 2/0,
    # therefore they are simply added in this step
    new_allel1 = allel1 + allel3
    new_allel2 = allel2 + allel4

    # go back to gt annotation
    gt = np.char.add(new_allel1.astype('U1'), '/')
    gt = np.char.add(gt, new_allel2.astype('U1'))

    # replace gt in original df
    df.iloc[np.where(first == True)[0],-9:] = gt
    df.iloc[np.where(second == True)[0],-9:] = gt

    # drop duplicates
    df = df.drop_duplicates(subset = ['chr', 'pos', 'ref', 'alt'])
    return df


def main(sample, dbsnp,mgp, output):
    # read in sample and database variants
    sample_df = pd.read_csv(sample, delimiter='\t')
    dbsnp_df = pd.read_csv(dbsnp, delimiter='\t', header = None, names = ['chr', 'pos', 'ref', 'alt'])
    mgp_df = pd.read_csv(mgp, delimiter='\t', header = None, names = ['chr', 'pos', 'ref', 'alt'])

    # correct annotation for duplicates
    sample_df_dup = duplicate(sample_df)

    # get intersectionss of database entries
    isec_snp = pd.merge(dbsnp_df, sample_df_dup , how = 'inner', on = ['chr', 'pos', 'ref', 'alt'])
    isec_mgp = pd.merge(mgp_df, sample_df_dup , how = 'inner', on = ['chr', 'pos', 'ref', 'alt'])

    # save variants that were corrected
    sample_df_dup.to_csv(f'{sample.split(".")[0]}.merged.remdup.txt')

    # save overlapping database entries
    isec_snp[['chr', 'pos', 'ref', 'alt']].to_csv(f'{output}_dbsnp_overlap.txt', index = False, sep = '_', header = False)
    isec_mgp[['chr', 'pos', 'ref', 'alt']].to_csv(f'{output}_mgp_overlap.txt', index = False, sep = '_', header = False)



if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='')
    parser.add_argument('-s', '--sample', help='Specify path to vcf sample file', type = str)
    parser.add_argument('-d', '--dbsnp', help='Specify path to known variant file from dbsnp', type = str)
    parser.add_argument('-m', '--mgp', help='Specify path to known variant file from mgp', type = str)
    parser.add_argument('-o', '--output', help='Specify path to output prefix', type = str)
    args = parser.parse_args()
    main(**vars(args))