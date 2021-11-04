import numpy as np
import matplotlib.pyplot as plt
from pylab import *
import pandas as pd
import collections
import csv
import math
import json
import itertools



def calc_min(bins_list):
    min_val = abs(bins_list[1] - bins_list[0])
    for i in range(2, len(bins_list)):
        if abs(bins_list[i] - bins_list[i-1]) < min_val:
            min_val = abs(bins_list[i] - bins_list[i-1])
    return min_val

def calc_max(bins_list):
    max_val = abs(bins_list[1] - bins_list[0])
    for i in range(2, len(bins_list)):
        if abs(bins_list[i] - bins_list[i-1]) > max_val:
            max_val = abs(bins_list[i] - bins_list[i-1])
    return max_val

def calc_mean(bins_list):
    sum = 0

    for i in range(1, len(bins_list)):
        sum +=  abs(bins_list[i] - bins_list[i-1])
    return sum/(len(bins_list)-1)


# function for computing the max/min diagonals
def calc_min_max_mean_tot(json_obj, fun, n_bin, val_dataset):
    celle = []
    sum = 0
    new_bins_json = []
    for i in json_obj:
        celle += [range(0,len(set(i['bins'])))]

        curr_min = min(i['bins'])
        curr_max = max(i['bins'])

        new_bins = [0]
        if curr_min != curr_max:
            for j in range(1, len(i['bins'])):
                new_bins += [abs(((i['bins'][j] - curr_min)/(curr_max-curr_min)) - ((i['bins'][j-1] - curr_min)/(curr_max-curr_min)))]

        new_bins_json += [new_bins]

    to_visit = list(itertools.product(*celle))[1:]

    sols = []
    for i in to_visit:
        count = 0
        for j,e in enumerate(i):
            count += new_bins_json[j][e]**2
        sols += [math.sqrt(count)]
    return fun(sols)/math.sqrt(len(json_obj))


# function for computing the diagonal of the solution
def calc_min_max_mean_tot_sol(json_obj, sol, n_bin, val_dataset):
    sum = 0
    new_sol = re.split(', ',sol[1:-1])
    num = 0
    for i in json_obj:
        curr_min = min(i['bins'])
        curr_max = max(i['bins'])

        new_bins = [0]
        if curr_min != curr_max:
            for j in range(1, len(i['bins'])):
                new_bins += [abs(((i['bins'][j] - curr_min)/(curr_max-curr_min)) - ((i['bins'][j-1] - curr_min)/(curr_max-curr_min)))]
        else:
            for j in range(1, len(i['bins'])):
                new_bins += [0]

        idx_sol= int(new_sol[num])
        # print(idx_sol, new_bins)
        sol_val = new_bins[idx_sol]
        num += 1

        sum += sol_val**2
    return math.sqrt(sum)/math.sqrt(len(json_obj))


def proximity(json_obj, sol, n_bin, val_dataset):
    new_sol = re.split(', ',sol[1:-1])
    norm_vals_sol = []
    num = 0
    for i in json_obj:
        # print('sooool ', i, new_sol)
        curr_min = min(i['bins'])
        curr_max = max(i['bins'])
        # print('min e max ', curr_min, curr_max)

        norm_vals_bins = []
        if curr_min != curr_max:
            for j in range(0, len(i['bins'])):
                if i['op'] in ['<', '<=']:
                    norm_vals_bins += [((i['bins'][j] - curr_min)/(curr_max-curr_min))]
                else:
                    norm_vals_bins += [((curr_max - i['bins'][j])/(curr_max-curr_min))]
        else:
            # norm_vals_bins = [0] #qui ho sostituito questa riga on le due sotto
            for j in range(0, len(i['bins'])):
                norm_vals_bins += [0]
        # print(norm_vals_sol, norm_vals_bins)

        idx_sol= int(new_sol[num])
        norm_vals_sol += [norm_vals_bins[idx_sol]]
        # print(norm_vals_sol)
        num += 1

    # distance computation between Q and Qnew
    sum = 0.0
    for ii in range(0, len(norm_vals_sol)):
        sum += norm_vals_sol[ii]**2
    return math.sqrt(sum)/math.sqrt(len(json_obj))


def distQ_Qnew(data):
    algo_dataset_dist = pd.DataFrame()
    for i,row in data.iterrows():
        val_newQ_dataset = []
        if row['card_true_tot_Q'] != "CC IS ALREADY SATISFIED" and row['card_true_tot_Q'] != "REWRITING IS NOT POSSIBLE":
            # print(row['newQ'])
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
    ## FOR COMPUTING THE SAMPLE-BASED MEASURES WE NEED THE RESULTS OBTAINED FROM 1_coverage_rewriting_sql.py ##
    #### INPUT #### insert here the name of your files ##
    algo = pd.read_csv('resultsCRBasePI.csv')
    # algo = pd.read_csv('resultsCRBasePI_dataset.csv')
    # print(algo.columns)

    algo_new = pd.DataFrame()
    for i,row in algo.iterrows():
        if (row['card_true_tot_Q'] != 'REWRITING IS NOT POSSIBLE') & (row['card_true_tot_Q'] != "CC IS ALREADY SATISFIED"):
            # print(row['card_est_tot_newQ'], row['card_est_as_Q'], row['card_est_as_Q'][1:-2] ,row['CC'])

            ###############################
            ##### GRID-BASED MEASURES #####
            ###############################
            row['diag_min'] = calc_min_max_mean_tot(json.loads(row['output_prep']), min, row['n_bin'], '')
            row['diag_max'] = calc_min_max_mean_tot(json.loads(row['output_prep']), max, row['n_bin'],  '')

            ###############################
            ### SOLUTION-BASED MEASURES ###
            ###############################
            row['relaxation_degree'] =  (row['card_est_tot_newQ'] - row['card_est_tot_Q'])/row['card_est_tot_Q']
            row['diag_sol'] = calc_min_max_mean_tot_sol(json.loads(row['output_prep']), str(row['solution']), row['n_bin'],'')
            row['proximity'] = proximity(json.loads(row['output_prep']), str(row['solution']), row['n_bin'],'')

            # the following update is needed for the computation of the distance between the orginal query and the rewritten query (few lines below)
            row['card_true_sa_newQ'] = int(row['card_true_sa_newQ'][1:-1])
            algo_new = algo_new.append(row)

    # COMPUTATION OF THE DISTANCE BETWEEN THE INPUT QUERY AND THE REWRITTEN QUERY (this measure is needed for the sample-based measures)
    algo_new2 = distQ_Qnew(algo_new)

    ## the new file with some measures
    algo_new2.to_csv('resultsCRBasePI_measures.csv', index=False)





if __name__ == '__main__':
    main()
