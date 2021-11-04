import numpy as np
import pandas as pd
import sys
df = pd.read_csv(sys.argv[1], header=None, delimiter=" ", float_precision="round_trip")
inp = (df[0] >= float(sys.argv[2])) & (df[0] <= float(sys.argv[3])) & (df[1] >= float(sys.argv[4])) & (df[1] <= float(sys.argv[5]))
out = (df[0] >= float(sys.argv[6])) & (df[0] <= float(sys.argv[7])) & (df[1] >= float(sys.argv[8])) & (df[1] <= float(sys.argv[9]))
print("Disparity is {}".format(abs(np.sum(out)-2*np.sum(df[out][2]))))
print("Intersection is {}".format(np.sum(inp & out)))
print("Union is {}".format(np.sum(inp | out)))
print("Jaccard similarity is  is {}".format(np.sum(inp & out)/np.sum(inp | out)))
