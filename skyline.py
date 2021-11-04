"""
Algorithm
===============
1. Data strcuture to create - Store points in r-tree for filtering during query. 
2. Lookup segment tree data structure for finding next point outside range. Or can use sorted per dimension data structure during pre processing for both removal and addition.
3. Preprocessing:
	a. Create R-tree data strcuture.
	b. Create either segment tree or sorted data structure.
4. Query processing:
	a. Create a heap.
	b. Add input range to heap.
	c. Compute neighbours and add each neighbour to heap along with blue and red counts.
		1. Look up dimension wise data structure/ segment tree data structure to find out neighbours.
		2. Get bounding box and lookup r-tree data structure to find points.
		3. Create a skyline out of the points.
		4. Add all neighbouring ranges to heap.
	d. Return from function once a fair range is found.
"""
from rtree import index
from itertools import chain, combinations
from operator import itemgetter, attrgetter
import heapq
import time
from sortedcontainers import SortedList
from fractions import Fraction
import math
import pdb

def readFile(filename):
    in_file = open (filename, "r")
    # Read list of lines
    out = [] # list to save lines
    line = in_file.readline()
    index = 0
    while line:
        point = line.strip().split(" ")
        class_val = int(point[-1])
        point = list(map(float, point[:-1]))
        point.append(1 == class_val)
        point.insert(0, index)
        index += 1
        out.append(point)
        # Read next line
        line = in_file.readline()
    
    # Close file 
    in_file.close()
    
    # out - List of 
    return out

def get_skyline(points, function):
    # Points set with id, bounds, object attributes
    sky = []
    dims = list(range(len(function)))
    if function[0] == 1:
        func_ = lambda x: x["bounds"][0]
    else:
        func_ = lambda x: -x["bounds"][0]
    sorted_points = SortedList(points, key=func_)
    for p in sorted_points:
        dominated = False
        for p_ in sky:
            total = 0
            for index in dims:
                # Function represents the quadrant we are dealing with
                # the opposite value of the function value is our weight
                if function[index]:
                    if p_["bounds"][index] > p["bounds"][index]:
                        total += 1
                        break
                else:
                    if p_["bounds"][index] < p["bounds"][index]:
                        total += 1
                        break
            if total == 0:
                dominated = True
                break
        if not dominated:
            sky.append(p)
    return sky

def get_neighbours(current_range, points_in_input_range, sorted_lists, rtree):
    # Point format : (point, fairness)
    # Range : (left bottom point, right top point)
    dims = int(len(current_range[0])/2)
    no_points = len(sorted_lists[0])
    output = set()

    # ranges - Find the external boundary surrounding the rectangle which consists of one additional point
    ranges = []
    # pdb.set_trace()
    for start_end in range(2): # 0 - start, 1 - end
        if start_end:
            default = math.inf
        else:
            default = -math.inf
        for d in range(dims):
            loc = dims * start_end + d
            value_pos = current_range[0][loc]
            # Create dummy point for bisect
            dummy_point = list(current_range[0][dims:])
            dummy_point[d] = value_pos
            dummy_point.insert(0, -1) # Push in dummy id of -1
            dummy_point.append(False) # Dummy color setting to red
            location = sorted_lists[d].bisect_left(dummy_point)
            if start_end == 1:
                upd_val = 1
            else:
                upd_val = -1
                if location > 0 and sorted_lists[d][location][d+1] > value_pos:
                    location -= 1
            while location >= 0 and location < no_points:
                found = True
                for d_ in range(dims):
                    if d == d_:
                        continue
                    # sorted_lists has an offset of 1
                    if sorted_lists[d][location][d_+1] < current_range[0][d_] or sorted_lists[d][location][d_+1] > current_range[0][d_ + dims]:
                        found = False
                        break
                if found:
                    ranges.append(sorted_lists[d][location][d+1])
                    break
                location += upd_val
            if len(ranges) <= loc:
                ranges.append(default)
    functions = list(chain((dims*(0,),), (l[0] * (0,) + sum(((1,) + (i-j-1) * (0,) for i, j in zip(l[1:], l[:-1])), ()) + (1,) + (dims-l[-1]-1)*(0,) for k in range(1,dims+1) for l in combinations(range(dims), k))))
    for func_ in functions:
        # 0 is start side, 1 is end side
        # TODO: 
        range_val = [None] * (dims*2)
        for dimension, f_ in enumerate(func_):
            if f_:
                range_val[dimension] = current_range[0][dimension + dims]
                range_val[dimension + dims] = ranges[dimension + dims]
            else:
                range_val[dimension] = ranges[dimension]
                range_val[dimension + dims] = current_range[0][dimension]
        # 3 properties - bounds, id, object
        try:
            points = [n for n in rtree.intersection(tuple(range_val), objects=True)]
        except:
            #pdb.set_trace()
            pass
        for ind_points in range(len(points)):
            correct_bounds = []
            for i_x in range(dims):
                correct_bounds.append(points[ind_points].bounds[2*i_x])
            points[ind_points] = {
                "id": points[ind_points].id,
                "object": points[ind_points].object,
                "bounds": correct_bounds
            }
        # pdb.set_trace()
        sky = get_skyline(points, func_)
        for point in sky:
            new_range = [None] * (dims*2)
            for dimension, f_ in enumerate(func_):
                if f_:
                    new_range[dimension] = current_range[0][dimension]
                    new_range[dimension + dims] = point["bounds"][dimension]
                else:
                    new_range[dimension] = point["bounds"][dimension]
                    new_range[dimension + dims] = current_range[0][dimension + dims]
            # Fix new range before pushing and check if the point is not within the input range
            if point["id"] not in points_in_input_range:
                new_fairness = current_range[1]
                if point["object"]:
                    new_fairness += blue_weight
                else:
                    new_fairness -= red_weight
                output.add(tuple([tuple(new_range), new_fairness]))
        # Add points along axis to range
        for ind, val_ in enumerate(ranges):
            if abs(val_) != math.inf:
                new_range = list(current_range[0])
                new_range[ind] = val_
                new_range_ = list(new_range)
                new_range_[(ind+dims)%(2*dims)] = current_range[0][ind]
                new_range = tuple(new_range)
                new_point_ = [n.object for n in rtree.intersection(new_range_, objects=True)]
                if len(new_point_) != 1:
                    #pdb.set_trace()
                    pass
                new_fairness = current_range[1]
                if new_point_[0]:
                    new_fairness += blue_weight
                else:
                    new_fairness -= red_weight
                output.add(tuple([new_range, new_fairness]))
    return output

def shrink_range(current_range, points_in_input_range, sorted_lists):
    # Point format : (point, total_count)
    # Range : (left bottom point, right top point)
    dims = int(len(current_range[0])/2)
    no_points = len(sorted_lists[0])
    output = set()
    for d in range(dims):
        for start_end in range(2): # 0 - start, 1 - end
            loc = dims * start_end + d
            value_pos = current_range[0][loc]
            # Create dummy point for bisect
            dummy_point = list(current_range[0][dims:])
            dummy_point[d] = value_pos
            dummy_point.insert(0, -1) # Push in dummy id of -1
            dummy_point.append(False) # Dummy color setting to red
            location = sorted_lists[d].bisect_left(dummy_point)
            if start_end == 1:
                upd_val = -1
            else:
                upd_val = 1
            while location >= 0 and location < no_points:
                if sorted_lists[d][location][0] in points_in_input_range:
                    # Found the next point location
                    range_ = list(current_range[0])
                    range_[loc] = sorted_lists[d][location][1 + d]
                    output.add(tuple([tuple(range_), current_range[1]]))
                    break
                location += upd_val
    return output

dataset = readFile("data/10k_urban_gb.csv")

# True - blue; False - red

blue_weight = 2
red_weight = 1
diff_val = 5
diff_val = 10

# Last point in each data point stands for the colour
dims = len(dataset[0])-2

start_time = time.perf_counter()
idx = index.Index()
if dims > 2:
    p = index.Property()
    p.dimension = dims
    idx = index.Index('kd_index',properties=p)

# 2. Add all points to r-tree. Use bool as the id
for d in dataset:
    # A point in r tree implementation of library is the same point repeated twice
    point = tuple(chain(*[d[1:-1], d[1:-1]]))
    idx.insert(d[0], point, d[-1])

# 3. Sorted list instead of a segment tree
sorted_lists = []
for d in range(dims):
    sorted_lists.append(SortedList(dataset, key=itemgetter(d+1)))

preprocess_end = time.perf_counter()

# 4. Shortest path using heap 
input_range = tuple([-0.015, 51.45, 0, 51.51])


# Add each point that belongs to range in set
points_set = set()
total_count = 0
for d in dataset:
    add = True
    for dim in range(dims):
        if d[dim+1] < input_range[dim] or d[dim+1] > input_range[dim + dims]:
            add = False
            break
    if add:
        points_set.add(d[0])
        if d[-1]:
            total_count += blue_weight
        else:
            total_count -= red_weight

print("Input range is {} with epsilon value of {} and blue weight as {}, red weight as {}.".format(input_range, diff_val, blue_weight, red_weight))
print("Fairness difference stands at {}".format(total_count))
heap_ = SortedList(key = lambda x: x[0])
heap_.add([Fraction(1,1), len(points_set), len(points_set), set([(input_range, total_count)])])

completed = False
output_range = None
output_similarity = -1
dummy = 1
while len(heap_) > 0 and not completed:
    # Remove last element, largest element
    frac, numer, denom, elems = heap_.pop()
    if dummy % 100 == 0:
        print(frac)
    dummy += 1
    # TODO:
    # Process each element
    frac_add = Fraction(numer, denom+1)
    frac_remove = Fraction(numer-1, denom)
    for rectange in elems:
        # Check if rectange is fair
        if abs(rectange[1]) < diff_val:
            # Closest fair range
            completed = True
            output_range = rectange
            output_similarity = frac
            break

        # rectange - (point-(), fairness)
        ranges_ = get_neighbours(rectange, points_set, sorted_lists, idx)
        # ranges_ - list of (range, total_count)
        # if heap_[heap_.bisect_left(frac_add)][0] != frac_add:
        #     heap_.add([frac_add, numer, denom+1, ranges_])
        # else:
        #     set_ = heap_[heap_.bisect_left(frac_add)][3]
        #     for r in ranges_:
        #         set_.add(r)
        heap_.add([frac_add, numer, denom+1, ranges_])
        # Check neighbourhood ranges and add to range
        ranges_ = shrink_range(rectange, points_set, sorted_lists)
        # ranges_ - list of (range, total_count)     
        # if heap_[heap_.bisect_left(frac_remove)][0] != frac_remove:
        #     heap_.add([frac_remove, numer-1, denom, ranges_])
        # else:
        #     set_ = heap_[heap_.bisect_left(frac_remove)][3]
        #     for r in ranges_:
        #         set_.add(r)
        heap_.add([frac_remove, numer-1, denom, ranges_])

end_time = time.perf_counter()
print("Preprocess time {:.3f} seconds".format(preprocess_end - start_time))
print("Query processing time {:.3f} seconds".format(end_time - preprocess_end))
print("Total time {:.3f} seconds".format(end_time - start_time))
print("Output range is {} and output similarity is {}".format(output_range, output_similarity))
