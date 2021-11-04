import psycopg2
import psycopg2.extras
import pandas as pd

#table_name = 'urban_gb'
table_name = 'uniform'

def drop_table(db_conn, table_name):
    db_cursor = db_conn.cursor()
    db_cursor.execute("DROP TABLE IF EXISTS " + table_name)
    db_conn.commit()
    db_cursor.close()

db_file = 'uniform.csv'
#db_file = '10k_new.csv'

params = {'host':'localhost', 'database':'db_2d', 'user':'postgres', 'password':'JumpPointer'}

conn = psycopg2.connect(**params)

drop_table(conn, table_name)

# Create new table
db_cursor = conn.cursor()
db_cursor.execute("CREATE TABLE " + table_name +" ( x real NOT NULL, y real NOT NULL, sa boolean NOT NULL, PRIMARY KEY (x, y));")
conn.commit()
db_cursor.close()

df = pd.read_csv(db_file, sep=' ', header=None)
df_columns = [0, 1, 2]
columns = ['x', 'y', 'sa']
values = "VALUES({})".format(",".join(["%s" for _ in df_columns]))
insert_stmt = "INSERT INTO {} ({}) {}".format(table_name, ','.join(columns), values)
db_cursor = conn.cursor()
list_val = df.values.tolist()
for i in range(len(list_val)):
    if table_name == 'urban_gb':
        list_val[i][2] = list_val[i][2] > 1.5
    if table_name == 'uniform':
        list_val[i][2] = list_val[i][2] > 0.5
psycopg2.extras.execute_batch(db_cursor, insert_stmt, list_val)
conn.commit()
db_cursor.close()
