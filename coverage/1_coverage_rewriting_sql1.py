import psycopg2
import pdb
import re
import itertools
import numpy as np
import sys
import math
from functools import reduce
import pandas as pd
import pandas.io.sql as sqlio
from scipy.stats import norm
import json
import time
import math

sys.setrecursionlimit(1000000000)
np.set_printoptions(suppress=True)
np.set_printoptions(threshold=sys.maxsize)

def get_similarity(c, q1, q2):
    cond1 = q1.split('WHERE')[1]
    cond2 = q2.split('WHERE')[1]
    cond1_pred = cond1.split('AND')
    cond2_pred = cond2.split('AND')
    with c.cursor() as cursor:
        cursor.execute('SELECT COUNT(*) FROM urban_gb WHERE (' + cond1_pred[0] + ' OR ' + cond2_pred[0] + ') AND ' + '(' + cond1_pred[1] + ' OR ' + cond2_pred[1] + ') AND ' + '(' + cond1_pred[2] + ' OR ' + cond2_pred[2] + ') AND ' + '(' + cond1_pred[3] + ' OR ' + cond2_pred[3] + ') GROUP BY sa')
        union_values = cursor.fetchall()
        cursor.execute('SELECT COUNT(*) FROM urban_gb WHERE (' + cond1_pred[0] + ' AND ' + cond2_pred[0] + ') AND ' + '(' + cond1_pred[1] + ' AND ' + cond2_pred[1] + ') AND ' + '(' + cond1_pred[2] + ' AND ' + cond2_pred[2] + ') AND ' + '(' + cond1_pred[3] + ' AND ' + cond2_pred[3] + ') GROUP BY sa')
        intersection_values = cursor.fetchall()
        return (intersection_values[0][0] + intersection_values[1][0])/(union_values[0][0] + union_values[1][0])

## functions for the sample
def compute_sample_size(MCE, CL):
    return math.ceil(((norm.ppf(1 - (1 - CL)/2)**2)*(0.5**2))/MCE**2)

def count_tuples(c, table):
    with c.cursor() as cursor:
        cursor.execute('SELECT COUNT(*) FROM ' + table)
        return cursor.fetchall()[0][0]

def do_sample(c, table, sample_size): #, attributes
    sql = 'SELECT * FROM ' + table + ' ORDER BY RANDOM() LIMIT ' + str(sample_size)
    # sql = 'SELECT ' + ', '.join(attributes) + ' FROM ' + table + ' ORDER BY RANDOM() LIMIT ' + str(sample_size)
    return sqlio.read_sql_query(sql, c)

# these functions have been created for tests
# for comparing results sometimes we fix the sample
def do_synthetic_sample_(c, n):
    sql = 'SELECT * FROM sample_'+ n
    return sqlio.read_sql_query(sql, c)

def do_synthetic_sample_dataset(c, n):
    sql = 'SELECT * FROM '+ n
    return sqlio.read_sql_query(sql, c)

# for comparing results we need additional info
def min_max_dataset(c, table, attrs):
    with c.cursor() as cursor:
        min_max_attrs = ', '.join(list(map(lambda x: 'MIN(' + x + '), MAX(' + x + ') ', attrs)))
        cursor.execute('SELECT ' + min_max_attrs + ' FROM ' + table)
        query_result = cursor.fetchall()
        result = []
        count = 0
        for a in attrs:
            result += [{'attr': a,'min': float(query_result[0][count]), 'max': float(query_result[0][count + 1])}]
            count += 2
        return result
################

# function that returns both the total cardinality and the cardinality for each value of the sensitive attribute
def cardinality_estimation(sample, table_size, q, k):
    sample_result_g = sample.query(q)
    est_sel_g = float(sample_result_g.shape[0]) /float(sample.shape[0])
    est_card_g = int(est_sel_g * table_size)

    sens_values_count = []
    for cond in k:
        condition = (' & ').join([sens_attr + ' == \'' + cond['value'][x] + '\'' for x, sens_attr in enumerate(cond['AS'])])
        sens_values_count.append(int(table_size * float(sample_result_g.query(condition).shape[0]) / float(sample.shape[0])))
    return (est_card_g, tuple(sens_values_count))

def get_query(relax_attributes, op, index):
    cond = []
    for j, v in enumerate(relax_attributes):
        cond += [v['attr'] + ' ' + v['op'] + ' ' + str(v['bins'][index[j]]) ]
    cond_str = op.join(cond)
    return cond_str

def get_query_result(relax_attributes, op, index):
    cond = []
    for j, v in enumerate(relax_attributes):
        if v['min_bin'] == v['max_bin']:
            cond += [v['attr'] + ' ' + v['op'] + ' ' + str(v['val_orig']) ]
        else:
            cond += [v['attr'] + ' ' + v['op'] + ' ' + str(v['bins'][index[j]]) ]
    cond_str = op.join(cond)
    return cond_str


def compute_distance(t):
    return np.linalg.norm(t)

def unique(x):
    return list(dict.fromkeys(x))

def chunks(lst, n):
    for i in range(0, len(lst), n):
        yield lst[i:i + n]


# the following functions check if a point is locked or not based on the definition of "upper right"and "lower left" proposed in the papers
def compare_tuples(t1, t2, fun):
    return all(map(fun, t1, t2))

def is_locked(t, lock):
    return any(map(lambda x: compare_tuples(t, x, lambda a, b: a >= b), lock))

def can_not_have_solutions(t, dont_have_solutions):
    return any(map(lambda x: compare_tuples(t, x, lambda a, b: a <= b), dont_have_solutions))

def is_locked_axes(t, lock):
    return list(filter(lambda x: compare_tuples(t, x[1], lambda a, b: a >= b), lock))

def can_not_have_solutions_axes(t, dont_have_solutions):
    return list(filter(lambda x: compare_tuples(t, x[1], lambda a, b: a <= b), dont_have_solutions))


# DIAGONAL PRUNING
# this function returns a point that satisfies the coverage constraint, one that does not satisfy and the sum of the computed estimates (just considering the diagonal)
def diagonal_search(k, sample, card_tot, relax_attributes, boolean_op):
    # high_min and high_max are needed in case the attributes are discretized with different number of bins
    high_min = min(len(x['bins']) for x in relax_attributes)
    high_max = max(len(x['bins']) for x in relax_attributes)

    # low and high are the first and the last point in the diagonal
    high = tuple([high_min-1]*len(relax_attributes))
    low = tuple([0]*len(relax_attributes))
    count_stime = 0

    query = get_query(relax_attributes, boolean_op, low)
    res = (cardinality_estimation(sample, card_tot, query, [*k]), low)
    ## res is a tuple containing both the estimations and the index of the considered point ( (est_tot, (est_sa)), idx )
    # --> in res[0] we have the estimations - res[0][0]: total cardinality estimation, res[0][1]: cardinality estimation for the sensitive attribute
    # --> in res[1] we have the index of the point

    count_stime = count_stime + 1
    if all([x >= int(k[i]['num']) for i, x in enumerate(res[0][1])]):
        return ((res[0][1], res[0][0]), low), None, count_stime

    query = get_query(relax_attributes, boolean_op, high)
    res = (cardinality_estimation(sample, card_tot, query, [*k]), high)
    count_stime = count_stime + 1

    # if the attributes are discretized with different number of bins we consider the point as top right as possible
    if any([x < int(k[i]['num']) for i, x in enumerate(res[0][1])]):
        if high_min != high_max:
            new_high = tuple([len(x['bins'])-1 for x in relax_attributes])
            query = get_query(relax_attributes, boolean_op, new_high)
            new_res = (cardinality_estimation(sample, card_tot, query, [*k]), new_high)
            count_stime = count_stime + 1
            # if the new high doesn't satisfy the constraint, we return new_high
            if any([x < int(k[i]['num']) for i, x in enumerate(new_res[0][1])]):
                return None, ((new_res[0][1], new_res[0][0]), tuple(new_res[1])), count_stime
            else:
                return ((new_res[0][1], new_res[0][0]), new_high), ((res[0][1], res[0][0]), tuple(res[1])), count_stime

        return None, ((res[0][1], res[0][0]), tuple(res[1])), count_stime
    last = res

    # if the furthest point on the diagonal satisfies the coverage constraint and the lowest point does not satisfy, we look for the first point on the diagonal that satisfies the constraint through binary search
    while any(x <= y for x, y in zip(low, high)):
        # mid = l + (r - l)/2;
        mid = tuple(a - b for a, b in zip(high, low))
        mid = tuple(int(i / 2) for i in mid)
        mid = tuple(a + b for a, b in zip(low, mid))
        query = get_query(relax_attributes, boolean_op, mid)
        res = (cardinality_estimation(sample, card_tot, query, [*k]), mid)
        count_stime = count_stime + 1

        if all([x >= int(k[i]['num']) for i, x in enumerate(res[0][1])]):
            high = tuple(i - 1 for i in mid)
            last = tuple(res)
        else:
            low = tuple(i + 1 for i in mid)

    if any([x < int(k[i]['num']) for i, x in enumerate(last[0][1])]):
        return None, ((last[0][1], last[0][0]), tuple(last[1])), count_stime
    else:
        last_1 = tuple(x-1 for x in last[1])
        query = get_query(relax_attributes, boolean_op, last_1)
        res = (cardinality_estimation(sample, card_tot, query, [*k]), last_1)
        count_stime = count_stime + 1
        return ((last[0][1], last[0][0]), last[1]), ((res[0][1], res[0][0]), tuple(res[1])), count_stime


# DIMENSIONAL PRUNING
# this function returns a list of points that satisfy the coverage constraint, a list of points that do not satisfy it and the sum of the computed estimates
def axes_search(k, sample, card_tot, relax_attributes, boolean_op):
    locks = []
    no_sols_tot = []
    count_stime = 0
    for i, v in enumerate(relax_attributes):
        no_sols = []
        j = 0
        # we look for the points that satisfy the coverage constraint and those that do not satisfy it for each axis (with binary search)
        while j < min([len(x['bins']) for x in relax_attributes]):
            high = [j]*len(relax_attributes)
            high[i] = len(relax_attributes[i]['bins']) - 1

            loock = is_locked_axes(high, locks)
            if len(loock) > 0:
                high[i] = min(map(lambda x: x[1][i], loock)) - 1

            high = tuple(high)
            low = [j]*len(relax_attributes)
            low[i] = 0

            no_sool = can_not_have_solutions_axes(low, no_sols_tot)
            if len(no_sool) > 0:
                low[i] = max(map(lambda x: x[1][i], no_sool)) + 1

            if low[i] >= len(relax_attributes[i]['bins']):
                break

            low = tuple(low)
            query = get_query(relax_attributes, boolean_op, low)
            res = (cardinality_estimation(sample, card_tot, query, [*k]), low)
            ## res is a tuple containing both the estimations and the index of the considered point ( (est_tot, (est_sa)), idx )
            # --> in res[0] we have the estimations - res[0][0]: total cardinality estimation, res[0][1]: cardinality estimation for the sensitive attribute
            # --> in res[1] we have the index of the point
            count_stime = count_stime + 1

            if all([x >= int(k[i]['num']) for i, x in enumerate(res[0][1])]):
                locks.append(((res[0][1], res[0][0]), low))
                break

            query = get_query(relax_attributes, boolean_op, high)
            res = (cardinality_estimation(sample, card_tot, query, [*k]), high)
            count_stime = count_stime + 1
            if any([x < int(k[i]['num']) for i, x in enumerate(res[0][1])]):
                no_sols = [((res[0][1], res[0][0]), high)]
                j += 1
                continue

            last = res
            while low[i] < high[i]:
                # mid = l + (r - l)/2;
                mid = [j] * len(relax_attributes)
                mid[i] = int((high[i] - low[i])/2)
                mid[i] = low[i] + mid[i]
                mid = tuple(mid)
                query = get_query(relax_attributes, boolean_op, mid)
                res = (cardinality_estimation(sample, card_tot, query, [*k]), mid)
                count_stime = count_stime + 1

                if all([x >= int(k[i]['num']) for i, x in enumerate(res[0][1])]):
                    high = [j]*len(relax_attributes)
                    high[i] = mid[i] - 1
                    high = tuple(high)
                    last = tuple(res)
                else:
                    low = [j]*len(relax_attributes)
                    low[i] = mid[i] + 1
                    low = tuple(low)
            if all([x >= int(k[i]['num']) for i, x in enumerate(last[0][1])]):
                old_res = ((last[0][1], last[0][0]), last[1])
                idx_lock = list(last[1])
                idx_lock[i] -= 1
                idx_lock = tuple(idx_lock)
                query = get_query(relax_attributes, boolean_op, idx_lock)
                res = (cardinality_estimation(sample, card_tot, query, [*k]), idx_lock)
                count_stime += 1
                while idx_lock[i] > 0 and all([x >= int(k[i]['num']) for i, x in enumerate(res[0][1])]):
                    old_res = ((res[0][1], res[0][0]), res[1])
                    idx_lock = list(idx_lock)
                    idx_lock[i] -= 1
                    idx_lock = tuple(idx_lock)
                    query = get_query(relax_attributes, boolean_op, idx_lock)
                    res = (cardinality_estimation(sample, card_tot, query, [*k]), idx_lock)
                    count_stime += 1
                locks.append(old_res)
                if idx_lock[i] >= 0:
                    no_sols.append(((res[0][1], res[0][0]), res[1]))
            j += 1
        no_sols_tot += no_sols
    return locks, no_sols_tot, count_stime


# COVERAGE-BASED REWRITING (CRBaseI if iter = True, otherwise CRBase)
# this function returns the point that satisfies the coverage constraint by minimazing the total cardinality and the distance from the origin, and the sum of the computed estimates
def minimum_index(k, locks, nosols, sample, card_tot, relax_attributes, boolean_op, n_bin, iter):
    # in LOCK and DONT_HAVE_SOLUTIONS we save the indexes of the points (no cardinality estimation)
    lock = [x[1] for x in locks]
    lock = unique(lock)

    min_lock = locks[0]
    for a in locks[1:]:
        if a[0][1] < min_lock[0][1]:
            min_lock = a

    # MINIMUMS = ((CARDINALITY SA, CARDINALITY TOT), DISTANCE, INDEX)
    minimums = [(min_lock[0], compute_distance(min_lock[1]), min_lock[1])]

    dont_have_solutions = [x[1] for x in nosols]
    dont_have_solutions = unique(dont_have_solutions)

    min_i = [-1]*len(relax_attributes)
    max_i = [len(x['bins']) for x in relax_attributes]

    # we use CURRENT_FRACTION for manage iteration
    if iter:
        current_fraction = n_bin / 2
    else:
        current_fraction = 1

    count = 0
    while current_fraction > 0:
        # We check that the points in the locking set are not in the "upper right" of another point in the same set.
        # Analogous for points in the no_solution set: we check that the points in this set are not in the "lower left" of another point in the same set.
        new_lock = []
        for lk in lock:
            if not is_locked(lk, new_lock):
                new_lock = list(filter(lambda x: not is_locked(x, [lk]), new_lock))
                new_lock.append(lk)
        lock = new_lock

        new_nosols = []
        for lk in dont_have_solutions:
            if not can_not_have_solutions(lk, new_nosols):
                new_nosols = list(filter(lambda x: not can_not_have_solutions(x, [lk]), new_nosols))
                new_nosols.append(lk)
        dont_have_solutions = new_nosols

        for x in lock:
            if x.count(0) == len(x) - 1:
                new_max = [(i, v) for i, v in enumerate(x) if v > 0][0]
                max_i[new_max[0]] = new_max[1]

        for x in dont_have_solutions:
            if [a==len(b['bins'])-1 for a,b in zip(x, relax_attributes)].count(True) == len(x) - 1:
                new_min = [(i, v) for i, v in enumerate(x) if v < len(relax_attributes[i]['bins']) - 1][0]
                min_i[new_min[0]] = new_min[1]

        # creation of the search space
        min_max_idx = [unique([0]+ [x - 1 for x in list(filter(lambda a: a % current_fraction == 0, range(1, len(i['bins']) + 1)))]) for i in relax_attributes]
        # print(min_max_idx)
        min_max_idx = [list(filter(lambda y: y > min_i[i] and y < max_i[i], x)) for i, x in enumerate(min_max_idx)]
        min_max_idx = pd.MultiIndex.from_product(min_max_idx, names=['i_' + str(x) for x in range(0, len(relax_attributes))])
        to_visit = pd.DataFrame(index=min_max_idx)

        for chunk in chunks(lock, 10):
            q_lock = ') & ('.join([' | '.join(['i_' + str(j) + ' < ' + str(y) for j, y in enumerate(x)]) for x in chunk])
            if q_lock != '':
                q_lock = '(' + q_lock + ')'
            to_visit = to_visit.query(q_lock)

        for chunk in chunks(dont_have_solutions, 10):
            q_dont_have_solutions = ') & ('.join([' | '.join(['i_' + str(j) + ' > ' + str(y) for j, y in enumerate(x)]) for x in chunk])
            if q_dont_have_solutions != '':
                q_dont_have_solutions = '(' + q_dont_have_solutions + ')'
            to_visit = to_visit.query(q_dont_have_solutions)

        # to_visit is the list of points to visit, already sorted by the distance
        if len(to_visit) > 0:
            for idx_p,p in enumerate(relax_attributes):
                to_visit['v_'+str(idx_p)] = to_visit.index.map(lambda x: abs(p['op_dist'] - ((p['bins'][x[idx_p]] - p['min_bin'])/ (p['max_bin'] - p['min_bin']))) if (p['max_bin'] - p['min_bin']) != 0 else 0)
            to_visit['dist'] = to_visit.eval('sqrt(' + ' + '.join(['v_' + str(x) + '*' + 'v_' + str(x) for x in range(0, len(relax_attributes))]) + ')')
            to_visit.sort_values(by='dist', inplace=True)

        current_fraction = int(current_fraction/2)

        # As long as to_visit contains elements, for each of these we check if it can be our solution as explained in the paper
        while len(to_visit) > 0:
            u = to_visit.iloc[0]
            u_i = to_visit.index[0]
            u = (u['dist'], u_i)
            to_visit = to_visit.iloc[1:]

            count = count + 1
            query = get_query(relax_attributes, boolean_op, u[1])
            est = cardinality_estimation(sample, card_tot, query, [*k])
            est = (est[1], est[0])

            if all([x >= int(k[i]['num']) for i, x in enumerate(est[0])]):
                if est[1] < minimums[0][0][1] or (est[1] == minimums[0][0][1] and u[0] < minimums[0][1]):
                    minimums = [(est, u[0], u[1])]
                elif est[1] == minimums[0][0][1] and u[0] == minimums[0][1]:
                    minimums.append((est, u[0], u[1]))
                lock = [x for x in lock if not is_locked(x, [u[1]])]
                lock.append(u[1])

                q_lock = ' | '.join(['i_' + str(j) + ' < ' + str(y) for j, y in enumerate(u[1])])
                to_visit = to_visit.query(q_lock)

            elif est[1] >= minimums[0][0][1]:
                lock = [x for x in lock if not is_locked(x, [u[1]])]
                lock.append(u[1])
                q_lock = ' | '.join(['i_' + str(j) + ' < ' + str(y) for j, y in enumerate(u[1])])
                to_visit = to_visit.query(q_lock)
                dont_have_solutions = [x for x in dont_have_solutions if not can_not_have_solutions(x, [u[1]])]
                dont_have_solutions.append(u[1])
            else:
                dont_have_solutions = [x for x in dont_have_solutions if not can_not_have_solutions(x, [u[1]])]
                dont_have_solutions.append(u[1])

    return minimums, count


# function for pre-processing
def min_max(sample, attrs):
    result = []
    for a in attrs:
        result += [{'attr': a,'min': sample[a].min(), 'max': sample[a].max()}]
    return result

def compute_bins_qcut(relax_attributes, n_bin, sample):
    for j, v in enumerate(relax_attributes):
        if v['op'] in ['<', '<=']:
            relax_attributes[j]['op_dist'] = 0
            if sample[v['attr']][sample[v['attr']] >= float(v['val'])].unique().size > 1:
                _,dim_bin = pd.qcut(sample[v['attr']][sample[v['attr']] >= float(v['val'])], q = n_bin-1, duplicates='drop', retbins=True)
                dim_bin = dim_bin.tolist()
            else:
                #dim_bin = [float(v['val']), v['max']]
                dim_bin = [float(v['val']), max(v['max'], float(v['val']))]
        else:
            relax_attributes[j]['op_dist'] = 1
            if sample[v['attr']][sample[v['attr']] <= float(v['val'])].unique().size > 1:
                _,dim_bin = pd.qcut(sample[v['attr']][sample[v['attr']] <= float(v['val'])], q = n_bin-1, duplicates='drop', retbins=True)
                dim_bin = dim_bin.tolist()
            else:
                #dim_bin = [v['min'], float(v['val'])]
                dim_bin = [min(v['min'], float(v['val'])), float(v['val'])]

            dim_bin = dim_bin[::-1]

        dim_bin[0] = float(v['val'])

        if len(set(dim_bin)) == 1:
            relax_attributes[j]['bins'] = dim_bin + dim_bin
        else:
            relax_attributes[j]['bins'] = dim_bin
    return relax_attributes

def compute_bins_distinct(relax_attributes, n_bin, sample):
    for j, v in enumerate(relax_attributes):
        #contare distinct val
        distinct_val = len(sample[v['attr']].unique())
        n_bin_attr = min(n_bin, distinct_val)

        if v['op'] in ['<', '<=']:
            relax_attributes[j]['op_dist'] = 0
            step = (v['max'] - float(v['val'])) / (n_bin_attr - 1)
            dim_bin = [float(v['val']) + step *i for i in range(0, n_bin_attr)]
        else:
            relax_attributes[j]['op_dist'] = 1
            step = (float(v['val']) - v['min']) / (n_bin_attr - 1)
            dim_bin = [float(v['val']) - step *i for i in range(0, n_bin_attr)]

        if len(set(dim_bin)) == 1:
            relax_attributes[j]['bins'] = dim_bin + dim_bin
        else:
            relax_attributes[j]['bins'] = dim_bin
    return relax_attributes

def compute_bins_base(relax_attributes, n_bin, sample):
    for j, v in enumerate(relax_attributes):
        #contare distinct val
        n_bin_attr = min(n_bin, len(sample[v['attr']]))

        if v['op'] in ['<', '<=']:
            relax_attributes[j]['op_dist'] = 0
            step = (float(v['max']) - float(v['val'])) / (n_bin_attr - 1)
            dim_bin = [float(v['val']) + step *i for i in range(0, n_bin_attr)]
        else:
            relax_attributes[j]['op_dist'] = 1
            step = (float(v['val']) - float(v['min'])) / (n_bin_attr - 1)
            dim_bin = [float(v['val']) - step *i for i in range(0, n_bin_attr)]

        if len(set(dim_bin)) == 1:
            relax_attributes[j]['bins'] = dim_bin + dim_bin
        else:
            relax_attributes[j]['bins'] = dim_bin
    return relax_attributes


def cmg_query_ref(sample, relax_attributes, boolean_op, table_size, sample_size, n_bin, sens_attr_list, k, prep, table_name,c, pruning, iter):

    orig_query = boolean_op.join([x['attr'] + ' ' + x['op'] + ' ' + x['val_orig'] for x in relax_attributes])
    #print(orig_query)
    init_card_est = cardinality_estimation(sample, table_size, orig_query,k)
    #print('The cardinality estimation for the input query is : ', init_card_est[0] , '\n the cardinality estimation for the protected group is: ', init_card_est[1] , '\n')

    ##### this part is for tests  #####
    orig_q = 'SELECT * FROM ' + table_name + ' WHERE ' + orig_query.replace('&', 'AND')
    # computation of the true cardinality for the input query
    card_vera_tot_Q = count_tuples(c, table_name + ' WHERE ' + orig_query.replace('&', 'AND'))
    card_vera_sa_Q = []
    for cond in k:
        condition = (' AND ').join([sens_attr + ' = \'' + cond['value'][x] + '\'' for x, sens_attr in enumerate(cond['AS'])])
        card_vera_sa_Q.append(count_tuples(c,  table_name + ' WHERE ('+ orig_query.replace('&', 'AND')  + ') AND (' + condition + ')'))
    # print(card_vera_tot_Q, card_vera_sa_Q)
    # selectivity_Q_est = init_card_est[0]/table_size
    # selectivity_Q_vera = card_vera_tot_Q/table_size
    #####

    # list of cardinality estimation for the given coverage constraints
    sens_attr_ests = list(init_card_est[1])


    ####################################################################################
    #################################### PREPROCESSING #################################
    ####################################################################################

    for i, mm in enumerate(min_max(sample, list(map(lambda x: x['attr'], relax_attributes)))):
        relax_attributes[i]['min'] = float(mm['min'])
        relax_attributes[i]['max'] = float(mm['max'])

        if float(relax_attributes[i]['val_orig']) > relax_attributes[i]['max']:
            relax_attributes[i]['val'] = str(relax_attributes[i]['max'])
        if float(relax_attributes[i]['val_orig']) < relax_attributes[i]['min']:
            relax_attributes[i]['val'] = str(relax_attributes[i]['min'])

    # we consider three types of binning approaches, but we test just two of them:
    # 'base' = 'equi-width' without consider distinct values, 'distinct' = 'equi-width' considering distinct values, 'qcut' = 'equi-depth' considering distinct values
    start_time_prep = time.time()
    if prep == 'base':
        relax_attributes = compute_bins_base(relax_attributes, n_bin, sample)
    if prep == 'distinct':
        relax_attributes = compute_bins_distinct(relax_attributes, n_bin, sample)
    if prep == 'qcut':
        relax_attributes = compute_bins_qcut(relax_attributes, n_bin, sample)
    end_time_prep = time.time()
    time_preprocessing = end_time_prep - start_time_prep

    for i, m in enumerate(relax_attributes):
        relax_attributes[i]['min_bin'] = min(relax_attributes[i]['bins'])
        relax_attributes[i]['max_bin'] = max(relax_attributes[i]['bins'])

    # if the coverage constraint is already satisfied in the original query we do not need the rewriting
    if all([sens_attr_ests[i] >= int(k[i]['num']) for i, x in enumerate([*k])]):
        print('Coverage constraint already satisfied')

        test_result_temp = [len(relax_attributes),
        n_bin,
        sample_size,
        orig_q,
        'CC IS ALREADY SATISFIED','','','','','','','','','','','','','','','','','','','','','','','']
        return test_result_temp


    ####################################################################################
    #################################### PROCESSING ####################################
    ####################################################################################

    locks = [] # in locks we save the points that satisfy the coverage constraint (both the cardinality estimations and the index of the points)
    nosols = [] # in nosols we save the points that do not satisfy the coverage constraint (both the cardinality estimations and the index of the points)
    # similar terminology is used in functions for both pruning and coverage-based rewriting algorithm

    start_time_pruning = time.time()
    if pruning:
        locks_axes, nosols_axes, count_stime_locks_search = axes_search(k, sample, table_size, relax_attributes, boolean_op)
        locks = locks_axes
        nosols = nosols_axes

        l_node, nosols_diag, count_stime_diagonal_search = diagonal_search(k, sample, table_size, relax_attributes, boolean_op)
        if l_node != None:
            locks.append(l_node)
        if nosols_diag != None:
            nosols.append(nosols_diag)
        count_stime_pruning = count_stime_locks_search + count_stime_diagonal_search
    else:
        # If we do not consider pruning rules, we estimate the cardinality of the nearest point and the farthest point.
        # According to the coverage satisfaction, we add them to the corresponding set (nosols or locks)
        count_stime_pruning = 0
        low = tuple([0]*len(relax_attributes))
        high = tuple([len(x['bins'])-1 for x in relax_attributes])
        query_ = get_query(relax_attributes, boolean_op, low)
        res = (cardinality_estimation(sample, table_size, query_, [*k]), low)
        count_stime_pruning = count_stime_pruning + 1
        if all([x >= int(k[i]['num']) for i, x in enumerate(res[0][1])]):
            locks.append(((res[0][1], res[0][0]), low))
        else:
            query_ = get_query(relax_attributes, boolean_op, high)
            new_res = (cardinality_estimation(sample, table_size, query_, [*k]), high)
            count_stime_pruning = count_stime_pruning + 1
            if all([x < int(k[i]['num']) for i, x in enumerate(new_res[0][1])]):
                nosols.append(((new_res[0][1], new_res[0][0]), tuple(new_res[1])))
            else:
                locks.append(((new_res[0][1], new_res[0][0]), tuple(new_res[1])))
                nosols.append(((res[0][1], res[0][0]), low))

    end_time_pruning = time.time()
    time_pruning = end_time_pruning - start_time_pruning

    # If there are no points in the locking set the rewriting is not possible
    if locks == []:
        print('not possible')
        test_result_temp = [len(relax_attributes),
         n_bin,
         sample_size,
         orig_q,
         'REWRITING IS NOT POSSIBLE','','','','','','','','','','','','','','','','','','','','','','','']
        return  test_result_temp

    ## we call the rewriting algorithm
    start_time_algo = time.time()
    min_idx, count_algo = minimum_index(k, locks, nosols, sample, table_size, relax_attributes, boolean_op, n_bin, iter)
    end_time_algo = time.time()
    time_algo = end_time_algo - start_time_algo
    count_stime_tot = count_algo + count_stime_pruning

    new_query = (get_query_result(relax_attributes, boolean_op, min_idx[0][2])).replace(' & ', ' AND ').replace(' | ', ' OR ')
    new_q = 'SELECT * FROM ' + table_name + ' WHERE ' + new_query
    #print('The rewritten query is :\n' + new_q + '\n')

    ##### this part is for the test #####
    # computation of the true cardinality
    card_vera_tot_newQ = count_tuples(c, table_name + ' WHERE ' + new_query)
    card_vera_sa_newQ = []
    for cond in k:
        condition = (' AND ').join([sens_attr + ' = \'' + cond['value'][x] + '\'' for x, sens_attr in enumerate(cond['AS'])])
        card_vera_sa_newQ.append(count_tuples(c,  table_name + ' WHERE ('+ new_query.replace('&', 'AND')  + ') AND (' + condition + ')'))
    #####

    test_result_temp = [len(relax_attributes),
     n_bin,
     sample_size,
     orig_q,
     card_vera_tot_Q,
     card_vera_sa_Q,
     init_card_est[1],
     init_card_est[0],
     k[0]['num'],
     min_idx[0][0][0],
     min_idx[0][0][1],
     json.dumps(relax_attributes), #output preprocessing
     time_preprocessing,
     time_pruning,
     time_algo,
     count_algo,
     count_stime_locks_search,
     count_stime_diagonal_search,
     count_stime_pruning,
     count_stime_tot,
     min_idx[0][2], #solution
     new_q,
     card_vera_tot_newQ,
     card_vera_sa_newQ,
     # selectivity_Q_est, selectivity_Q_vera,
     locks_axes, nosols_axes]

    return test_result_temp


def main():

    #### YOU WILL FIND INPUT DATA IN LINES 654 (for setting the coverage constraint), 701 (for specifying the query) AND 732 (for setting the connection to the database)

    ## SETTING FOR THE PROCESSING
    pruning = True # pruning = False for CRBase, otherwise pruning = True (for CRBaseP and CRBaseIP)
    iter = True # iter = False for CRBase, otherwise iter = True (for CRBaseI and CRBaseIP)

    ## COVERAGE CONSTRAINT DEFINITION
    # the coverage constraint should be in the following form: CC = [{'AS': ['gender'], 'value': ['Female'], 'num': '100'}]
    # for tests we insert the value for 'num' below
    # KEep 11 retries to find a fair range
    CC = [{'AS':['sa']}]
    # example considering more coverage constraints and more than one sensitive attribute:
    # CC = [{'AS': ['gender'], 'value': ['Female'], 'num': '100'}, {'AS': ['gender', 'nationality'], 'value': ['Female', 'Spain'], 'num': '50'}, {'AS': ['gender', 'nationality'], 'value': ['Male', 'Equador'], 'num': '50'}]

    # we extract the information we need from the given coverage constraint and we reorganize them
    sens_attr_cc = []
    for i in range(0, len(CC)):
        for ii in range(0, len(CC[i]['AS'])):
            sens_attr_cc += [CC[i]['AS'][ii]]
    sens_attr_list = list(set(sens_attr_cc))
    #print('Sensitive attribute and coverage constraints: ', sens_attr_list, CC)

    # dataframe for saving the results
    test_result_csv= pd.DataFrame(columns=['n_condizioni',
    'n_bin',
    'sample_size',
    'query',
    'card_true_tot_Q', # this is computed for tests
    'card_true_sa_Q', # this is computed for tests
    'card_est_sa_Q',
    'card_est_tot_Q',
    'CC',
    'card_est_sa_newQ',
    'card_est_tot_newQ',
    'output_prep' ,
    'time_preprocessing',
    'time_pruning',
    'time_algo',
    'count_algo',
    'countPruningDim','countPruningDiag',
    'countPruning',
    'count_stime_tot',
    'solution',
    'newQ',
    'card_true_tot_newQ', # this is computed for tests
    'card_true_sa_newQ', # this is computed for tests
    # 'selectivity_Q_est','selectivity_Q_true', # this is computed for tests
    'locksPruningDim', 'noSolPruningDim',
    'time_sample', 'preprocessing',
    'val_min_max_dataset' # this is computed for tests
    ])
    best_sim = -1.
    best_range = None
    # we consider different types of pre-processing: distinct = equi-width, qcut = equi-depth
    for prep in ['distinct', 'qcut']:
        # for tests we read the queries from a file
        #with open('urbanGB.txt', 'r') as file:
             #for query in file:
                #query = 'SELECT * FROM urban_gb WHERE x > 1500 AND x < 55 AND y < 35 AND y > '
                #query = 'SELECT * FROM urban_gb WHERE x > -0.375799 AND x < -0.191955 AND y > 51.5971 AND y < 51.6302'
                query = 'SELECT * FROM urban_gb WHERE  x > -0.190301 AND x < -0.107706 AND y > 51.4105 AND y < 51.4613'
                query = 'SELECT * FROM urban_gb WHERE  x > -0.178049 AND x < -0.13706 AND y > 51.5199 AND y < 51.6027'


                #print('Compute the disparity in query : {}'.format(query))
                q_parts = query.split('WHERE') # Assume WHERE all caps
                # we set the connection to the database
                c = psycopg2.connect(host= 'localhost', port= 5432, user= 'postgres', password= 'JumpPointer', database= 'db_2d')
                cursor = c.cursor()
                cursor.execute('SELECT sa,COUNT(*) FROM urban_gb WHERE ' + q_parts[1] + ' GROUP BY sa')
                results = cursor.fetchall()
                count_dict = {}
                count_dict[results[0][0]] = results[0][1]
                count_dict[results[1][0]] = results[1][1]
                true_weight = 1
                false_weight = 2
                delta = math.ceil(6514*(results[0][1] + results[1][1])*0.05/10000)
                disp = true_weight * count_dict['true'] - false_weight * count_dict['false']
                if abs(disp) <= delta:
                    # Already fair
                    print("Already fair")
                    return
                if disp > 0:
                     CC[0]['value'] = ['false']
                else:
                     CC[0]['value'] = ['true']
                # we extract the information we need from the query and we reorganize them
                print('Query originale:\n\n' + query + '\n')
                m = re.match(r'.*FROM\s(.*)\sWHERE(.*)', query)
                sens_attr_table = m.group(1)

                relax_attributes = []
                boolean_op = ''

                if (' OR ' in m.group(2) or ' or ' in m.group(2)) and (' AND ' in m.group(2) or ' and ' in m.group(2)):
                    print('Conditions must be all ANDs or all ORs!')
                    exit(1)
                elif ' OR ' in m.group(2):
                    boolean_op = ' | '
                elif ' AND ' in m.group(2):
                    boolean_op = ' & '

                for cond in re.split(' AND | OR ', m.group(2)):
                    c = re.match(r'([a-zA-Z\d\_]+)\s*(\<|\>|\<\=|\>\=)\s*([\-]?[\d]+[\.]?[\d]*)', cond.strip())
                    relax_attributes.append(
                    {
                        'attr': c.group(1).strip(),
                        'op': c.group(2).strip(),
                        'val_orig': c.group(3).strip(),
                        'val': c.group(3).strip()
                    })

                all_attributes = list(set([x['attr'] for x in relax_attributes])) + [sens_attr_list[0]]

                # we set the connection to the database
                c = psycopg2.connect(host= 'localhost', port= 5432, user= 'postgres', password= 'JumpPointer', database= 'db_2d')


                ##### code for test (we save this info for the computation of measures) #####
                val_min_max_dataset = min_max_dataset(c, sens_attr_table, list(map(lambda x: x['attr'], relax_attributes)))

                ##############################################################################
                ############################# SAMPLE GENERATION ##############################
                ##############################################################################
                table_size = count_tuples(c, sens_attr_table)

                ## general approach
                sample_size = compute_sample_size(0.01, 0.99)
                start_time_sample = time.time()
                sample = do_sample(c, sens_attr_table, sample_size)
                # sample = do_synthetic_sample_(c, sample_size)  ## for tests, when needed, we consider a fixed sample
                # sample = do_synthetic_sample_dataset(c, sens_attr_table) ## for tests we also executed directly on the dataset
                end_time_sample = time.time()
                time_sample = end_time_sample - start_time_sample
                # sample_size = len(sample) ## this is used when tests are done on the dataset


                ##### code for test (for setting the coverage constraint) #####
                ## cardinality estimation of the input query for specifying the coverage constraint
                orig_query = boolean_op.join([x['attr'] + ' ' + x['op'] + ' ' + x['val_orig'] for x in relax_attributes])

                # sens_attr_values = distinct_values(c, sens_attr_table, sens_attr_list[0])
                init_card_est = cardinality_estimation(sample, table_size, orig_query, CC)
                card_Q_as = init_card_est[1][0]
                card_Q_tot = init_card_est[0]

                for n_bin in [4,8,16]:
                    print("Using bin size {}".format(n_bin))
                    # here we can define the coverage constraint for tests. In this example we consider qf+20%
                    """
                    for perc in [0.2]: #[0.05,0.1,0.2,0.3]
                        for i, e in enumerate(CC):
                            if init_card_est[1][i] != 0:
                                CC[i]['num'] = int(int(init_card_est[1][i]) + int(init_card_est[1][i])*perc)
                            else:
                                CC[i]['num'] = int(1 + 1*perc)
                        print(CC)

                        test_result = cmg_query_ref(sample, relax_attributes, boolean_op, table_size, sample_size, n_bin, sens_attr_list, CC, prep, sens_attr_table,c, pruning, iter)
                        test_result += [time_sample, prep, json.dumps(val_min_max_dataset)]
                        s=pd.Series(test_result, index=test_result_csv.columns)
                        test_result_csv=test_result_csv.append(s, ignore_index=True)
                    """
                    if disp > delta:
                        # Unfair because its too true
                        # Start with bare minimum and add 10 more in increments of 10
                        start_cc = math.ceil( abs(disp - delta)/false_weight)
                        for i in range(11):
                            CC[0]['num'] = count_dict['false'] + start_cc + 10*i
                            test_result = cmg_query_ref(sample, relax_attributes, boolean_op, table_size, sample_size, n_bin, sens_attr_list, CC, prep, sens_attr_table,c, pruning, iter)
                            if test_result[4] == 'REWRITING IS NOT POSSIBLE':
                                continue
                            if test_result[4] != 'CC IS ALREADY SATISFIED':
                                new_query = test_result[21]
                                cursor = c.cursor()
                                cursor.execute('SELECT sa,COUNT(*) FROM urban_gb WHERE ' + new_query.split('WHERE')[1]  + ' GROUP BY sa')
                                results = cursor.fetchall()
                                count_dict_temp = {}
                                count_dict_temp[results[0][0]] = results[0][1]
                                count_dict_temp[results[1][0]] = results[1][1]
                                disp = true_weight * count_dict_temp['true'] - false_weight * count_dict_temp['false']
                                if disp > delta:
                                    continue
                                print("Fair range found")
                            sim_val = get_similarity(c, new_query, query)
                            if sim_val > best_sim:
                                best_sim = sim_val
                                best_range = new_query
                            if test_result[4] == 'CC IS ALREADY SATISFIED':
                                break

                            test_result += [time_sample, prep, json.dumps(val_min_max_dataset)]
                            s=pd.Series(test_result, index=test_result_csv.columns)
                            test_result_csv=test_result_csv.append(s, ignore_index=True)
                            # Check for fairness
                    else:
                        # Unfair because its too false
                        start_cc = math.ceil( abs(disp - delta)/true_weight)
                        for i in range(11):
                            CC[0]['num'] = count_dict['true'] + start_cc + 10*i
                            # Start with bare minimum and add 10 more in increments of 10
                            test_result = cmg_query_ref(sample, relax_attributes, boolean_op, table_size, sample_size, n_bin, sens_attr_list, CC, prep, sens_attr_table,c, pruning, iter)
                            if test_result[4] == 'REWRITING IS NOT POSSIBLE':
                                continue
                            if test_result[4] != 'CC IS ALREADY SATISFIED': 
                                new_query = test_result[21]
                                cursor = c.cursor()
                                cursor.execute('SELECT sa,COUNT(*) FROM urban_gb WHERE ' + new_query.split('WHERE')[1]  + ' GROUP BY sa')
                                results = cursor.fetchall()
                                count_dict_temp = {}
                                count_dict_temp[results[0][0]] = results[0][1]
                                count_dict_temp[results[1][0]] = results[1][1]
                                disp = true_weight * count_dict_temp['true'] - false_weight * count_dict_temp['false']
                                if disp > delta:
                                    continue
                            print("Fair range found")
                            sim_val = get_similarity(c, new_query, query)
                            if sim_val > best_sim:
                                best_sim = sim_val
                                best_range = new_query
                            if test_result[4] == 'CC IS ALREADY SATISFIED':
                                break
                            test_result += [time_sample, prep, json.dumps(val_min_max_dataset)]
                            s=pd.Series(test_result, index=test_result_csv.columns)
                            test_result_csv=test_result_csv.append(s, ignore_index=True)
                            # Check for fairness
                    """
                    test_result = cmg_query_ref(sample, relax_attributes, boolean_op, table_size, sample_size, n_bin, sens_attr_list, CC, prep, sens_attr_table,c, pruning, iter)
                    test_result += [time_sample, prep, json.dumps(val_min_max_dataset)]
                    s=pd.Series(test_result, index=test_result_csv.columns)
                    test_result_csv=test_result_csv.append(s, ignore_index=True)
                    """
                    # Check for fairness
                    

    test_result_csv.to_csv('resultsCRBasePI.csv', index=False)
    if best_sim > 0.:
        print(best_range)
        print(best_sim)
    # test_result_csv.to_csv('resultsCRBasePI_dataset.csv', index=False)

    # writer=pd.ExcelWriter('resultsCRBasePI.xlsx')
    # test_result_csv.to_excel(writer)
    # writer.save()


if __name__ == '__main__':
    main()
