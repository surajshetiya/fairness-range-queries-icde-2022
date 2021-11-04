import numpy as np
import matplotlib.pyplot as plt
from pylab import *
import pandas as pd
import collections
import csv
import math
import json
import seaborn as sns
import re


def distQ_Qnew(data):
    algo_dataset_dist = pd.DataFrame()
    for i,row in data.iterrows():
        val_newQ_dataset = []
        if row['card_true_tot_Q'] != "CC IS ALREADY SATISFIED" and row['card_true_tot_Q'] != "REWRITING IS NOT POSSIBLE":
            print(row['newQ'])
            m = re.match(r'.*FROM\s(.*)\sWHERE(.*)', row['newQ'])
            for cond in re.split(' AND | OR ', m.group(2)):
                c = re.match(r'([a-zA-Z\d\_]+)\s*(\<|\>|\<\=|\>\=)\s*([\-]?[\d]+[\.]?[\d]*)', cond.strip())
                val_newQ_dataset += [float(c.group(3).strip())]
            # print(val_newQ_dataset)

        # NORMALIZATION
            val_newQ_dataset_norm = normalized_query(row, val_newQ_dataset)
            # print('valori normalizzati:',val_newQ_dataset_norm)

        # COMPUTATION OF DISTANCE BETWEEN POINTS
            dist = 0.0
            for ii in val_newQ_dataset_norm:
                dist += ii**2
                #print(dist)
            dist = sqrt(dist)
            #print(dist)
            row['dist_Qnew-Q'] = dist
            algo_dataset_dist = algo_dataset_dist.append(row)

    return algo_dataset_dist


def normalized_query(row, val_newQ_sample):
# NORMALIZATION
    val_newQ_sample_norm = []
    for ii, jj, val in zip(json.loads(row['val_min_max_dataset']), json.loads(row['output_prep']), val_newQ_sample):

        if jj['op'] in ['<', '<=']:
            curr_min = float(jj['val_orig'])
            curr_max = float(ii['max'])
            val_newQ_sample_norm +=  [(val - curr_min)/(curr_max-curr_min)]
            print(curr_min, curr_max)
        else:
            curr_min = float(ii['min'])
            curr_max = float(jj['val_orig'])
            val_newQ_sample_norm +=  [(curr_max - val)/(curr_max-curr_min)]
    # print('point coordinates corresponding to normalized Qnew', val_newQ_sample_norm)
    return val_newQ_sample_norm




def main():

    ## FOR COMPUTING THE SAMPLE-BASED MEASURES WE NEED THE RESULTS OBTAINED ON BOTH THE DATASET AND THE SAMPLE ##
    # please use the output of 2_calc_measures_sol.py (for both the sample and the dataset): at the end of that file we compute the distance between the original query and the rewritten one
    # if you do not want to compute the solution measures and you just want to compute the sample-based measures you need to uncomment the following code lines and comment lines

    # algo_dataset0 = pd.read_csv('resultsCRBasePI_dataset.csv')
    # algo_sample0 = pd.read_csv('resultsCRBasePI.csv')
    #
    # # DATASET # computation of the distance between the input query and the rewritten one for the DATASET
    # algo_dataset_new = pd.DataFrame()
    # for i,row in algo_dataset0.iterrows():
    #     if (row['card_true_tot_Q'] != 'REWRITING IS NOT POSSIBLE') & (row['card_true_tot_Q'] != "CC IS ALREADY SATISFIED"):
    #         row['card_true_sa_newQ'] = int(row['card_true_sa_newQ'][1:-1])
    #         algo_dataset_new = algo_dataset_new.append(row)
    #
    # algo_dataset = distQ_Qnew(algo_dataset_new)
    #
    # # SAMPLE # computation of the distance between the input query and the rewritten one for the SAMPLE
    # algo_sample_new = pd.DataFrame()
    # for i,row in algo_sample0.iterrows():
    #     if (row['card_true_tot_Q'] != 'REWRITING IS NOT POSSIBLE') & (row['card_true_tot_Q'] != "CC IS ALREADY SATISFIED"):
    #         row['card_true_sa_newQ'] = int(row['card_true_sa_newQ'][1:-1])
    #         algo_sample_new = algo_sample_new.append(row)
    #
    # algo_sample = distQ_Qnew(algo_sample_new)


    #### INPUT #### insert here the name of your files ##
    algo_dataset = pd.read_csv('resultsCRBasePI_dataset_measures.csv')
    algo_sample = pd.read_csv('resultsCRBasePI_measures.csv')
    #######

    # measures will be computed only in the case the rewritten query exists
    algo_sample = algo_sample[algo_sample['card_true_tot_Q'] != "CC IS ALREADY SATISFIED"]
    algo_sample = algo_sample[algo_sample['card_true_tot_Q'] != "REWRITING IS NOT POSSIBLE"]


    algo_sample2 = pd.DataFrame()
    for i,row in algo_sample.iterrows():
        algo_dataset_tmp = algo_dataset[algo_dataset['n_bin'] == row['n_bin']]
        algo_dataset_tmp = algo_dataset_tmp[algo_dataset_tmp['query'] == row['query']]
        algo_dataset_tmp = algo_dataset_tmp[algo_dataset_tmp['preprocessing'] == row['preprocessing']]

        row['card_true_tot_Q'] = algo_dataset_tmp['card_est_tot_Q'].values[0]
        row['card_true_sa_Q'] = algo_dataset_tmp['card_est_sa_Q'].values[0]

        ##############################
        ######### MINIMALITY #########
        ##############################
        # according to further analysis we slightly change the definition of minimality
        row['minimality'] = abs(int(row['card_true_tot_newQ']) - algo_dataset_tmp['card_true_tot_newQ'].values[0])/algo_dataset_tmp['card_est_tot_Q'].values[0]
        row['minimality_new'] = abs(int(row['card_true_tot_newQ']) - algo_dataset_tmp['card_true_tot_newQ'].values[0])/algo_dataset_tmp['card_true_tot_newQ'].values[0]

        # here we check if the rewritten query satifies the coverage constraint also on the instance
        if row['card_true_sa_newQ'] < row['CC']:
            row['cc_satisfied_on_I'] = 0
        else:
            row['cc_satisfied_on_I'] = 1


        ##############################
        ######### PROXIMITY ##########
        ##############################
        # dist_Qnew-Q is the distance between the input query and the rewritten one: the computation has been done in the file sample-based measures)
        row['dist_Qnew-Q_dataset'] = algo_dataset_tmp['dist_Qnew-Q'].values[0]

        row['proximity'] = abs(row['dist_Qnew-Q'] - row['dist_Qnew-Q_dataset'])


        ##############################
        ##### SOLUTION DISTANCE ######
        ##############################
        #### query sample
        val_newQ_sample = []
        m = re.match(r'.*FROM\s(.*)\sWHERE(.*)', row['newQ'])
        for cond in re.split(' AND | OR ', m.group(2)):
            c = re.match(r'([a-zA-Z\d\_]+)\s*(\<|\>|\<\=|\>\=)\s*([\-]?[\d]+[\.]?[\d]*)', cond.strip())
            val_newQ_sample += [float(c.group(3).strip())]

        ### query dataset
        val_newQ_dataset = []
        m2 = re.match(r'.*FROM\s(.*)\sWHERE(.*)', algo_dataset_tmp['newQ'].values[0])
        for cond in re.split(' AND | OR ', m2.group(2)):
            c = re.match(r'([a-zA-Z\d\_]+)\s*(\<|\>|\<\=|\>\=)\s*([\-]?[\d]+[\.]?[\d]*)', cond.strip())
            val_newQ_dataset += [float(c.group(3).strip())]

        # NORMALIZATION
        val_newQ_dataset_norm = normalized_query(row, val_newQ_dataset)
        val_newQ_sample_norm = normalized_query(row, val_newQ_sample)
        #print(val_newQ_sample, val_newQ_dataset, val_newQ_sample_norm, val_newQ_dataset_norm)

        # COMPUTATION OF DISTANCE BETWEEN POINTS
        dist = 0.0
        for ii in range(0,len(val_newQ_dataset_norm)):
            dist += (val_newQ_dataset_norm[ii] - val_newQ_sample_norm[ii])**2
        dist = sqrt(dist)
        print(dist, i)
        row['solution_distance'] = dist

        algo_sample2 = algo_sample2.append(row)

    ## the new file with the sample-based measures
    algo_sample2.to_csv('results_measures.csv', index=False)


if __name__ == '__main__':
    main()
