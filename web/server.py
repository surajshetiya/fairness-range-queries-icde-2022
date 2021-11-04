from flask import Flask
from flask import render_template
from configparser import ConfigParser
import psycopg2
import psycopg2.extras
from flask import request
from werkzeug.utils import secure_filename
import pandas as pd
import time

def get_table_names(db_conn):
    db_cursor = db_conn.cursor()
    db_cursor.execute("SELECT table_schema, table_name FROM information_schema.tables WHERE (table_schema = 'public') ORDER BY table_schema, table_name;")
    list_tables = db_cursor.fetchall()
    db_cursor.close()
    return list(map(lambda x: x[1], list_tables))
    
def drop_table(db_conn, table_name):
    db_cursor = db_conn.cursor()
    db_cursor.execute("DROP TABLE IF EXISTS " + table_name)
    db_conn.commit()
    db_cursor.close()

db_file = 'uniform_1k.csv'

UPLOAD_FOLDER = 'downloads'
# Implementation is for unweighted case
blue_weight = 1
red_weight = -1

config_filename = 'database.ini'
section = 'postgresql'
naive_table = 'naive'
preprocess_table = 'preprocess'

parser = ConfigParser()
parser.read(config_filename)

params = dict(parser.items(section))
"""
params = {}
db_params = parser.items(section)

for db_param in db_params:
    params[db_param[0]] = db_param[1]
"""

print('Connecting to the PostgreSQL database...')
conn = psycopg2.connect(**params)

table_names = get_table_names(conn)
print("Existing Table names are : {}".format(' '.join(table_names)))

print("Droppping table if it exists")
drop_table(conn, naive_table)
drop_table(conn, preprocess_table)


# Create new table
db_cursor = conn.cursor()
db_cursor.execute("CREATE TABLE " + naive_table + " ( value real NOT NULL, sa boolean NOT NULL, PRIMARY KEY (value) );")
conn.commit()
db_cursor.close()

# Create new table
db_cursor = conn.cursor()
db_cursor.execute("CREATE TABLE " + preprocess_table + " ( value real NOT NULL, sa boolean NOT NULL, index Integer NOT NULL, jfb INTEGER, jfr INTEGER, jbb INTEGER, jbr INTEGER, cumltv_sum INTEGER NOT NULL, PRIMARY KEY (index));")
conn.commit()
db_cursor.close()

# List tables
table_names = get_table_names(conn)
print("New Table names are : {}".format(' '.join(table_names)))

# Create new index on table
db_cursor = conn.cursor()
db_cursor.execute("CREATE INDEX IF NOT EXISTS " + preprocess_table + "_value ON " + preprocess_table + "(value);")
conn.commit()
db_cursor.close()

def naive_query(db_conn, start_range, end_range, threshold):
    db_cursor = db_conn.cursor()
    db_cursor.execute("select n1.value , n2.value from naive n1 inner join naive n2 on n1.value < n2.value where (n1.value < {} and n2.value > {});".format(start_range, end_range))
    list_pairs = db_cursor.fetchall()
    db_cursor.close()
    #print("Found {} pairs".format(len(list_pairs)))
    db_cursor = db_conn.cursor()
    range_best = None
    sim_best = 0.
    for pair_i in list_pairs:
        db_cursor.execute("select sa,count(*) from naive where value >= {} and value <= {} group by sa;".format(pair_i[0], pair_i[1]))
        sa_count = db_cursor.fetchall()
        if abs(sa_count[0][1] - sa_count[1][1]) < threshold:
            db_cursor.execute("select count(*) from  {} where value >= {} and value <= {};".format(naive_table, max(pair_i[0], start_range), min(pair_i[1], end_range)))
            inter_pts = db_cursor.fetchall()
            db_cursor.execute("select count(*) from {} where value >= {} and value <= {};".format(naive_table, min(pair_i[0], start_range), max(pair_i[1], end_range)))
            union_pts = db_cursor.fetchall()
            if sim_best < (inter_pts[0][0]/union_pts[0][0]):
                sim_best = inter_pts[0][0]/union_pts[0][0]
                range_best = [pair_i[0], pair_i[1]]
    #print("Sim best is {}".format(sim_best))
    #print("Range best is {}".format(range_best))
    db_cursor.close()

def get_point_with_id(db_conn, id_val):
    db_cursor = db_conn.cursor()
    db_cursor.execute("select * from {} where index={}".format(preprocess_table, id_val))
    point = db_cursor.fetchall()
    db_cursor.close()
    return point

def get_similarity(input_start, input_end, sim_start, sim_end):
    union = max(sim_end[2], input_end[2]) - min(sim_start[2], input_end[2]) + 1
    intersect = min(sim_end[2], input_end[2]) - max(sim_start[2], input_end[2]) + 1
    return intersect/union

def jp_query(db_conn, start_range, end_range, threshold):
    db_cursor = db_conn.cursor()
    db_cursor.execute("select * from {} where value < {} order by value DESC limit 1".format(preprocess_table, start_range))
    point_start = db_cursor.fetchall()[0]
    print(point_start)
    db_cursor.execute("select * from {} where value <= {} order by value DESC limit 1".format(preprocess_table, end_range))
    point_end = db_cursor.fetchall()[0]
    print(point_end)
    sim_best = 0.
    range_best = None
    db_cursor.close()
    if abs(point_start[7]-point_end[7]) < threshold:
        # Input is a fair range
        sim_best = 1.
        range_best = (start_range, end_range)
        #print("Sim best is {}".format(sim_best))
        #print("Range best is {}".format(range_best))
        db_cursor.close()
        return
    #print("Iniitial disparity : {}".format(abs(point_start[7]-point_end[7])))
    stack_fwd_expanse = [point_end]
    stack_fwd_shrink = [point_end]
    point_start_exact = get_point_with_id(db_conn, point_start[2] + 1)[0]
    stack_bwd_expanse = [point_start_exact]
    stack_bwd_shrink = [point_start_exact]
    range_best=[]
    disparity = disp = point_start[7]-point_end[7]
    start_point = point_start
    end_point = point_end
    while abs(disparity) > threshold-1:
        point = None
        loc = 4
        add_val = blue_weight
        if disparity > 0:
            loc = 3
            add_val = red_weight
        if stack_fwd_expanse[-1][loc] == -1:
            break
        point = get_point_with_id(db_conn, stack_fwd_expanse[-1][loc])[0]
        disparity += add_val
        #print("Stack fwd expanse{} disparity {}".format(point, disparity))
        stack_fwd_expanse.append(point)
    disparity = disp
    while abs(disparity) > threshold-1:
        point = None
        loc = 6
        add_val = blue_weight
        if disparity > 0:
            add_val = red_weight
            loc = 5
        if stack_bwd_expanse[-1][loc] == -1:
            break
        point = get_point_with_id(db_conn, stack_bwd_expanse[-1][loc])[0]
        disparity += add_val
        #print("Stack bwd expanse {} disparity {}".format(point, disparity))
        stack_bwd_expanse.append(point)
    disparity = disp
    while abs(disparity) > threshold-1:
        point = None
        loc = 5
        sub_val = red_weight
        if disparity > 0:
            sub_val = blue_weight
            loc = 6
        temp_point = get_point_with_id(db_conn, stack_fwd_shrink[-1][2] + 1)[0]
        if temp_point[loc] == -1:
            break
        point = get_point_with_id(db_conn, temp_point[loc] - 1)[0]
        if point[2] < point_start_exact[2]:
            break
        #print("Stack fwd shrink {} disparity {}".format(point, disparity))
        disparity -= sub_val
        stack_fwd_shrink.append(point)
    disparity = disp
    while abs(disparity) > threshold-1:
        point = None
        loc = 3
        sub_val = red_weight
        if disparity > 0:
            sub_val = blue_weight
            loc = 4
        temp_point = get_point_with_id(db_conn, stack_bwd_shrink[-1][2] - 1)[0]
        if temp_point[loc] == -1:
            break
        point = get_point_with_id(db_conn, temp_point[loc] + 1)[0]
        #print("Stack bwd shrink {} disparity {}".format(point, disparity))
        if point[2] > point_end[2]:
            break
        disparity -= sub_val
        stack_bwd_shrink.append(point)
    # column order - [0 : 'value', 1 : 'sa', 2: 'index', 3: 'jfb', 4 : 'jfr', 5 : 'jbb', 6 : 'jbr', 7 : 'cumltv_sum']
    # Pair stack_bwd_shrink with stack_fwd_shrink and stack_fwd_expanse
    # Pair stack_bwd_expanse with stack_fwd_shrink and stack_fwd_expanse
    for i in range(len(stack_bwd_shrink)):
        if len(stack_fwd_shrink) > abs(disp) - i:
            # Range is stack_bwd_shrink[i] , stack_fwd_shrink[abs(disp)-i]
            sim_val = get_similarity(point_start_exact, point_end, stack_bwd_shrink[i], stack_fwd_shrink[abs(disp)-i])
            if sim_val > sim_best: 
                sim_best = sim_val
                range_best = [stack_bwd_shrink[i], stack_fwd_shrink[abs(disp)-i]]
        if len(stack_fwd_expanse) > abs(disp) - i:
            # Range is stack_bwd_shrink[i] , stack_fwd_expanse[abs(disp)-i]
            sim_val = get_similarity(point_start_exact, point_end, stack_bwd_shrink[i], stack_fwd_expanse[abs(disp)-i])
            if sim_val > sim_best:
                sim_best = sim_val
                range_best = [stack_bwd_shrink[i], stack_fwd_expanse[abs(disp)-i]]
    for i in range(len(stack_bwd_expanse)):
        if len(stack_fwd_shrink) > abs(disp) - i:
            # Range is stack_bwd_expanse[i] , stack_fwd_shrink[abs(disp)-i]
            sim_val = get_similarity(point_start_exact, point_end, stack_bwd_expanse[i], stack_fwd_shrink[abs(disp)-i])
            if sim_val > sim_best:
                sim_best = sim_val
                range_best = [stack_bwd_expanse[i], stack_fwd_shrink[abs(disp)-i]]
        if len(stack_fwd_expanse) > abs(disp) - i:
            # Range is stack_bwd_expanse[i] , stack_fwd_expanse[abs(disp)-i]
            sim_val = get_similarity(point_start_exact, point_end, stack_bwd_expanse[i], stack_fwd_expanse[abs(disp)-i])
            if sim_val > sim_best:
                sim_best = sim_val
                range_best = [stack_bwd_expanse[i], stack_fwd_expanse[abs(disp)-i]]
    #print("Sim best is {}".format(sim_best))
    #print("Range best is {}".format(range_best))
    db_cursor.close()
        

def preprocess(list_val):
    # column order - ['value', 'sa', 'index', 'jfb', 'jfr', 'jbb', 'jbr', 'cumltv_sum']
    sorted_list = sorted(list_val, key=lambda x: x[0])
    cumltv_sum = 0

    # Generate jump pointers
    # Forward jump pointers
    jp_dict = {}
    for i in range(len(sorted_list)):
        if sorted_list[i][1] == 1:
            # Blue
            cumltv_sum += 1
        else:
            # Red
            cumltv_sum -= 1
        sorted_list[i].extend([i, -1, -1, -1, -1, cumltv_sum])
        # 0 for red and 1 for blue
        if cumltv_sum in jp_dict:
            for ind in jp_dict[cumltv_sum][0]:
                sorted_list[ind][4] = i
            for ind in jp_dict[cumltv_sum][1]:
                sorted_list[ind][3] = i
            del jp_dict[cumltv_sum]
        if cumltv_sum + 1 not in jp_dict:
            jp_dict[cumltv_sum + 1] = [[], []]
        jp_dict[cumltv_sum + 1][1].append(i)
        if cumltv_sum - 1 not in jp_dict:
            jp_dict[cumltv_sum - 1] = [[], []]
        jp_dict[cumltv_sum - 1][0].append(i)
    # Backward jump pointers
    cumltv_sum = 0
    jp_dict = {}
    j = len(sorted_list)
    for i in range(len(sorted_list)):
        # Get index in j
        j -= 1
        if sorted_list[j][1] == 1:
            # Blue
            cumltv_sum += 1
        else:
            # Red
            cumltv_sum -= 1
        if cumltv_sum in jp_dict:
            for ind in jp_dict[cumltv_sum][0]:
                sorted_list[ind][6] = j
            for ind in jp_dict[cumltv_sum][1]:
                sorted_list[ind][5] = j
            del jp_dict[cumltv_sum]
        if cumltv_sum + 1 not in jp_dict:
            jp_dict[cumltv_sum + 1] = [[], []]
        jp_dict[cumltv_sum + 1][1].append(j)
        if cumltv_sum - 1 not in jp_dict:
            jp_dict[cumltv_sum - 1] = [[], []]
        jp_dict[cumltv_sum - 1][0].append(j)
    return sorted_list
        

def process_file(filename, db_conn):
    file_path = UPLOAD_FOLDER + '/' + filename
    df = pd.read_csv(file_path, sep=' ', header=None)
    df_columns = [0, 1]
    columns = ['value', 'sa']
    values = "VALUES({})".format(",".join(["%s" for _ in df_columns]))
    insert_stmt = "INSERT INTO {} ({}) {}".format(naive_table, ','.join(columns), values)
    db_cursor = db_conn.cursor()
    list_val = df.values.tolist()
    for i in range(len(list_val)):
        list_val[i][1] = list_val[i][1] > 0.5
    psycopg2.extras.execute_batch(db_cursor, insert_stmt, list_val)
    db_conn.commit()
    db_cursor.close()

    print("Insert completed - naive table")

    # Preprocess procedure
    print("Preprocess start")

    # Measure time taken
    start = time.process_time()
    preprocessed_data = preprocess(list_val)
    time_taken = time.process_time() - start

    pp_columns = ['value', 'sa', 'index', 'jfb', 'jfr', 'jbb', 'jbr', 'cumltv_sum']
    pp_values = "VALUES({})".format(",".join(["%s" for _ in pp_columns]))
    pp_insert_stmt = "INSERT INTO {} ({}) {}".format(preprocess_table, ','.join(pp_columns), pp_values)
    db_cursor = db_conn.cursor()
    psycopg2.extras.execute_batch(db_cursor, pp_insert_stmt, preprocessed_data)
    db_conn.commit()
    db_cursor.close()

    return time_taken


print(process_file(db_file, conn))
naive_start = time.process_time()
naive_query(conn, 0.35, 0.62, 2)
naive_time_taken = time.process_time() - naive_start
print("time taken by naive algorithm is {}".format(naive_time_taken))
jp_start = time.process_time()
jp_query(conn, 0.35, 0.62, 1)
jp_time_taken = time.process_time() - jp_start
print("time taken by jump pointer algorithm is {}".format(jp_time_taken))

"""
TODO:
1. Naive test interface
2. Time measure and report interface
3. Efficient data structure create interface
4. Time measure for efficient jump pointer creation
5. Query processing
6. Batch insert of input into table
"""

app = Flask(__name__)
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER

@app.route('/upload', methods=['POST'])
def upload_file():
    if request.method == 'POST':
        file = request.files['upload_file']
        file.save(f"{secure_filename(file.filename)}")

@app.route("/")
def hello_world():
    return "<p>Hello, World!</p>"
