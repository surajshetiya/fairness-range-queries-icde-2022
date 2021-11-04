#include<iostream>
#include<vector>
#include <algorithm>
#include <random>
#include<unordered_map>
#include <set>
#include <cmath>
#include <chrono> 
#include <sstream>
#include <string>
#include <fstream>

/*
=====
TODO:
=====
1. Sort along each dimension and maintain sorted array - Done
2. Write code to expand along each dimension
3. Update intersection using binary search - Done
4. Find all fair ranges using the fair range method - Done
    a. Add cumulative count for intersection as well - Done
    b. Find the Jaccard similarity using cumulative count - Done
*/

using namespace std;

class Range {
    public:
    vector<double> start, end;

    Range(vector<double> st, vector<double> e) {
        start = st;
        end = e;
    }
};

class Point {
    public:
    
    int id;
    vector<double> X;
    bool color; // color = { true: blue, false: red }
    bool intersect; // intersect with query rectangle
    int sort_dim;

    Point(vector<double> coordinates, bool col, int identifier) {
        X = coordinates;
        color = col;
        id = identifier;
        sort_dim = X.size() - 1;
        intersect = false;
    }

    bool operator < (const Point& other) const
    {
        // Use last dimension for comparision, sorting
        return X[sort_dim] < other.X[sort_dim];
    }
};

void print_range(Range r) {
    cout << "Range is : (";
    for(int i=0; i<r.start.size(); ++i) {
        if(i != 0) {
            cout << " , ";
        }
        cout << r.start.at(i);
    }
    cout << ") to (";
    for(int i=0; i<r.end.size(); ++i) {
        if(i != 0) {
            cout << " , ";
        }
        cout << r.end.at(i);
    }
    cout << ")" << endl;
}

// Number of dimensions
int dims=2, diff_val=5, points_in_db = 10000;
int blue_weight = 2, red_weight = -1;
vector<double> st, e;
Range fair_range=Range(st, e); // Place holder
double best_similarity=0.;
int input_range_size=0; // Will be updated later when we get the query
double similarity_expected = 0.4;
int opt_intersect, opt_union;


vector<Point> generate_database(int count, int dim) {
    vector<Point> db;
    for(int i=0; i< count;++i) {
        vector<double> pnt;
        for(int j=0; j< dim;++j) {
            pnt.push_back(1000.*double(rand())/double(RAND_MAX));
        }
        db.push_back(Point(pnt, rand() % 2 == 0, i));
    }
    return db;
}

// vector<Point> database = generate_database(points_in_db, dims);
vector<Point> database;

void read_urban_gb_file_global_db(string filename) {
    ifstream file_handler(filename);
    if (!file_handler.is_open())  // check file is open, quit if not
    {
        std::cerr << "failed to open file\n";
        return;
    }

    double latitude, longitude;
    int no_vehicle_collusion, id=0;
    string line;
    while (getline(file_handler, line))
    {
        istringstream iss(line);
        if (!(iss >> latitude >> longitude >> no_vehicle_collusion)) {
            cout << "Issues parsing file!" << endl;
            break; 
        } // error
        vector<double> X;
        X.push_back(latitude);
        X.push_back(longitude);
        // there are a total of 3112 1s and 6888 others. So a scale factor of blues=2, reds =1
        Point p(X, no_vehicle_collusion == 1, id++); 
        database.push_back(p);
    }
    points_in_db = database.size();
}

bool compare_range(Range r1, Range r2) {
    if(r1.start.size() != r2.start.size() || r1.end.size() != r2.end.size()) {
        return true;
    }
    for(int i=0; i < r1.start.size(); ++i) {
        if(r1.start.at(i) != r2.start.at(i)) {
            return true;
        }
    }
    for(int i=0; i < r1.end.size(); ++i) {
        if(r1.end.at(i) != r2.end.at(i)) {
            return true;
        }
    }
    return false;
}

int color_to_value(bool color) {
    return color ? blue_weight : red_weight; // color = { true: blue, false: red }
}

/*
TODO:
------
1. Sort using dimension d
2. Recursively solve the problem by choosing 2 points as start and end range for dimension i
3. Filter points to have only those which are within range
4. Use the points to recursively perform the operation
5. Use 1 dimensional fair range algorithm to find all pairs of fair ranges
*/

// 5. Use 1 dimensional fair range algorithm to find all pairs of fair ranges
// Get fair range when we have reached last dimension
void get_fair_range_one_dim(vector<Point>& db, int diff, Range range) {
    int size_val = db.size(), total = 0, min = size_val, max = -size_val, dims_cur_range = range.start.size(), intersect=0;
    vector<int> cumultv;
    vector<int> count_intersection;
    unordered_map<int, vector<int>> locations_map;
    for(int i=0; i<size_val; ++i) {
        total += color_to_value(db.at(i).color);
        cumultv.push_back(total);

        intersect += int(db.at(i).intersect);
        count_intersection.push_back(intersect);
        
        if(total < min) {
            min = total;
        }
        if(total > max) {
            max = total;
        }
        if(locations_map.find(total) == locations_map.end()) {
            vector<int> index;
            index.push_back(i);
            locations_map.insert(make_pair(total, index));
        } else {
            locations_map.find(total)->second.push_back(i);
        }
    }

    // Find all pairs of fair ranges
    for(int i = min; i <= max; ++i) {
        int end_range = (i + diff) < max ? (i + diff) : max;
        if(locations_map.find(i) == locations_map.end()) {
            continue;
        }
        for(int j = i; j <= end_range; ++j) {
            if(locations_map.find(j) == locations_map.end()) {
                continue;
            }

            // All pairs are with cumulative sum i,j are fair
            for (auto it_i = locations_map.find(i)->second.begin(); it_i != locations_map.find(i)->second.end(); ++it_i) {
                int index_i = *it_i;
                for (auto it_j = locations_map.find(j)->second.begin(); it_j != locations_map.find(j)->second.end(); ++it_j) {
                    int index_j = *it_j;
                    if(index_i == index_j) {
                        continue;
                    }
                    // Enter entry into fair range
                    // Copy constructor
                    Range f_range = range;
                    double val_i = db.at(index_i).X.at(dims-1), val_j = db.at(index_j).X.at(dims-1);
                    int min_index, max_index;
                    double min_val, max_val;
                    if(index_i < index_j) {
                        min_index = index_i;
                        max_index = index_j;
                        min_val = val_i;
                        max_val = val_j;
                    } else {
                        min_index = index_j;
                        max_index = index_i;
                        min_val = val_j;
                        max_val = val_i;
                    }
                    f_range.start.push_back(db.at(min_index+1).X.at(dims-1));
                    f_range.end.push_back(max_val);
                    // Append the last coordinate to the range
                    // int intersection_size = count_intersection.at(max_index)-count_intersection.at(min_index) + int(db.at(min_index).color);
                    // int union_size = input_range_size + (max_index - min_index + 1)- intersection_size;
                    int intersection_size = count_intersection.at(max_index)-count_intersection.at(min_index);
                    int union_size = input_range_size + (max_index - min_index)- intersection_size;
                    double similarity = double(intersection_size)/double(union_size);

                    // TODO: check formula
                    if(similarity > best_similarity) {
                        best_similarity = similarity;
                        fair_range = f_range;
                        opt_intersect = intersection_size;
                        opt_union = union_size;
                        // cout << "Similarity is : " << best_similarity << endl;
                        // print_range(fair_range);
                    }
                }
            }
        }
    }
}

void get_fair_range_recursive(vector<Point> db, int diff, Range range, bool check=false){
    int dim = range.start.size(), db_dim = db.at(0).X.size();
    bool call_last = (dim == (db_dim - 2));
    cout << "Database size after filtering is " << db.size() << endl;
    for(int i=0; i<db.size(); ++i) {
        if(i%10 == 0) {
            cout << "Outer loop : " << i << endl;
        }

        double val_i = db.at(i).X.at(dim);
        for(int j=i+1; j<db.size(); ++j) {
            double val_j = db.at(j).X.at(dim);
            double min_val = val_i < val_j ? val_i : val_j, max_val = val_i > val_j ? val_i : val_j;
            vector<Point> new_db;
            for(int k=0; k<db.size(); ++k) {
                if(db.at(k).X.at(dim) >= min_val && db.at(k).X.at(dim) <= max_val) {
                    new_db.push_back(db.at(k));
                }
            }
            // Copy constructor
            Range new_range = range;
            new_range.start.push_back(min_val);
            new_range.end.push_back(max_val);
            if(call_last) {
                get_fair_range_one_dim(new_db, diff, new_range);
            } else {
                get_fair_range_recursive(new_db, diff, new_range);
            }
            if(check && similarity_expected <= best_similarity) {
                return;
            }
        }
    }
}

// Compare the X[dim] with elem to perform binary search
int binary_search(vector<Point>& db, int start, int end, int dim, double elem) {
    if(end < start) {
        // Should not happen
        return start;
    }
    if(start == end) {
        return start;
    }
    int mid = (start + end)/2;
    if(db.at(mid).X.at(dim) == elem) {
        return mid;
    } else if(db.at(mid).X.at(dim) < elem) {
        return binary_search(db, mid+1, end, dim, elem);
    } else{
        return binary_search(db, start, mid-1, dim, elem);
    }
}

void fair_ranges(Range input_range, bool check=false) {
    // Sort by first dimension
    sort(database.begin(), database.end()); // Sorted along last dimension

    vector<vector<Point>> database_sorted;
    for(int i=0; i < dims; ++i) { // Last dimension is already sorted in database variable
        vector<Point> db_copy = database;
        for(int j=0; j<db_copy.size(); ++j) {
            db_copy.at(j).sort_dim = i;// Set sort dim
        }
        sort(db_copy.begin(), db_copy.end());
        database_sorted.push_back(db_copy);
    }

    // Set intersection size
    int start_point = binary_search(database, 0, database.size(), dims-1, input_range.start.at(dims-1));
    int end_point = binary_search(database, 0, database.size(), dims-1, input_range.end.at(dims-1));
    input_range_size = 0;
    for(int i=start_point; i<end_point; ++i) {
        int j =0;
        for( ;j<dims; ++j) {
            if(database.at(i).X.at(j) < input_range.start.at(j) || database.at(i).X.at(j) > input_range.end.at(j)) {
                break;
            }
        }
        if(j == dims) {
            ++input_range_size;
            database.at(i).intersect = true;
        }
    }

    int expansion_size = ceil(input_range_size*(1-similarity_expected)/similarity_expected);
    cout << "Input range size is : " << input_range_size << endl;

    // Get new range
    vector<double> st, e;
    Range new_boundary=Range(st, e);
    for(int i=0; i< dims; i++) {
        int start_range_index, end_range_index;
        vector<Point> &ref = database_sorted.at(i);
        start_range_index = binary_search(ref, 0, database.size()-1, i, input_range.start.at(i));
        end_range_index = binary_search(ref, 0, database.size()-1, i, input_range.end.at(i));

        // Expand logic
        int count_expand = 0;
        while(end_range_index < ref.size() && count_expand <= expansion_size) {
            int j=0;
            for(; j<dims; ++j){
                if(i == j) {
                    continue;
                }
                if(ref.at(end_range_index).X.at(j) < input_range.start.at(j) || ref.at(end_range_index).X.at(j) > input_range.end.at(j)) {
                    break;
                }
            }
            count_expand += int(j == dims);
            ++end_range_index;
        }

        if(end_range_index == ref.size()) {
            --end_range_index;
        }
        new_boundary.end.push_back(ref.at(end_range_index).X.at(i));

        count_expand = 0;
        while(start_range_index >= 0 && count_expand <= expansion_size) {
            int j=0;
            for(; j<dims; ++j){
                if(i == j) {
                    continue;
                }
                if(ref.at(start_range_index).X.at(j) < input_range.start.at(j) || ref.at(start_range_index).X.at(j) > input_range.end.at(j)) {
                    break;
                }
            }
            count_expand += int(j == dims);
            --start_range_index;
        }

        if(start_range_index < 0) {
            ++start_range_index;
        }
        new_boundary.start.push_back(ref.at(start_range_index).X.at(i));
    }

    print_range(new_boundary);

    // Obtain new points
    start_point = binary_search(database, 0, database.size(), dims-1, new_boundary.start.at(dims-1));
    end_point = binary_search(database, 0, database.size(), dims-1, new_boundary.end.at(dims-1));
    vector<Point> new_db;
    for(int i=start_point; i <= end_point; ++i) {
        int j=0;
        for(; j<dims; ++j){
            if(database.at(i).X.at(j) < new_boundary.start.at(j) || database.at(i).X.at(j) > new_boundary.end.at(j)) {
                break;
            }
        }
        if(j == dims) {
            new_db.push_back(database.at(i));
        }
    }

    // cout << "New database" << endl;
    // for(int i=0; i< new_db.size(); ++i) {
    //     for(int j=0; j < dims; j++) {
    //         cout << new_db.at(i).X.at(j) << " ";
    //     }
    //     cout << endl;
    // }

    vector<double> dummy_start, dummy_end;
    Range r = Range(dummy_start, dummy_end);
    // Call fair range on that new boundary database to find all fair ranges within similarity value
    get_fair_range_recursive(new_db, diff_val, r, check);
}

bool check_fairness(vector<Point> db, Range r, int eps, bool print=true) {
    int total=0, total_dims = r.start.size();
    for(int i=0; i < db.size(); ++i) {
        int j=0;
        for(; j<total_dims; ++j) {
            if(db.at(i).X.at(j) < r.start.at(j) || db.at(i).X.at(j) > r.end.at(j)) {
                break;
            }
        }
        if(j == total_dims) {
            total += color_to_value(db.at(i).color);
        }
    }
    if(print) {
        cout << "Total dispariry is : " << total << endl;
    }
    return abs(total) <= eps;
}

double check_similarity(vector<Point> db, Range input, Range other, bool print=true) {
    if(input.start.size() ==0 || other.start.size() == 0) {
        return 0.;
    }
    double similarity;
    int intersect = 0, union_points = 0, total_dims = input.start.size();
    for(int i=0; i < db.size(); ++i) {
        int j=0, k=0;
        for(; j<total_dims; ++j) {
            if(db.at(i).X.at(j) < input.start.at(j) || db.at(i).X.at(j) > input.end.at(j)) {
                break;
            }
        }
        for(; k<total_dims; ++k) {
            if(db.at(i).X.at(k) < other.start.at(k) || db.at(i).X.at(k) > other.end.at(k)) {
                break;
            }
        }
        intersect += int((j == total_dims) && (k == total_dims));
        union_points += int((j == total_dims) || (k == total_dims));
        // if(k == total_dims) {
            // for(int x_val=0; x_val < total_dims; ++x_val) {
            //     cout << db.at(i).X.at(x_val) << " ";
            // }
            // cout << endl;
        // }
    }
    if(print) {
        cout << "Intersect : " << intersect << endl;
        cout << "Union : " << union_points << endl;
    }
    similarity = double(intersect)/double(union_points);
    return similarity;
}

vector<double> brute_s, brute_e;
Range brute_opt = Range(brute_s, brute_e);
double brute_best_similarity = 0; 

vector<Point> sampling(vector<Point> db, int total_points) {
    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(0.0,db.size());

    vector<Point> sample_points;
    for (int i=0; i<total_points; ++i) {
        sample_points.push_back(db.at(int(distribution(generator))));
    }
    return sample_points;
}

void recursive_brute_force(vector<Point> db, Range input, int eps, vector<Point> state) {
    int dims = input.start.size();
    if(state.size() == (dims*2)) {
        vector<double> start_new_range, end_new_range;
        double inf = std::numeric_limits<double>::infinity();
        for(int i=0;i<dims; ++i) {
            start_new_range.push_back(inf);
            end_new_range.push_back(-inf);
        }
        for(int i=0;i<state.size(); ++i) {
            for(int j=0;j<dims; ++j) {
                if(start_new_range.at(j) > state.at(i).X.at(j)) {
                    start_new_range.at(j) = state.at(i).X.at(j);
                }
                if(end_new_range.at(j) < state.at(i).X.at(j)) {
                    end_new_range.at(j) = state.at(i).X.at(j);
                }
            }
        }
        Range r(start_new_range, end_new_range);
        if(check_fairness(db, r, eps, false)) {
            double sim = check_similarity(db, input, r, false);
            if(sim > brute_best_similarity) {
                brute_best_similarity = sim;
                brute_opt = r;
            }
        }
        return;
    }
    for(int i=0; i<db.size(); ++i) {
        if(state.size() == 0 && i%10 == 0) {
            cout << "Outer : " << i << endl;
        }
        vector<Point> state_copy = state; // Copy constructor
        state_copy.push_back(db.at(i));
        recursive_brute_force(db, input, eps, state_copy);
    }
}

int main() {

    read_urban_gb_file_global_db("data/10k_urban_gb.csv");
    cout << "Total database size : " << points_in_db << endl;
    cout << "Blue weight : " << blue_weight << endl << "Red weight : " << red_weight << endl;
    cout << "Diff value set to : " << diff_val << endl;
    vector<double> start_range{-0.015, 51.45};
    vector<double> end_range{0, 51.51};
    Range r=Range(start_range, end_range);

    bool check = false;
    // Set precission to 17 digits
    cout.precision(17);
    similarity_expected = 0.75;
    print_range(r);
    cout << "=================================================" << endl;
    cout << "No early stop criteria" << endl;
    cout << "Similarity expected " << similarity_expected << endl;
    auto start = std::chrono::high_resolution_clock::now();
    fair_ranges(r, check);
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    cout << "Points in the input range : " << input_range_size << endl;
    cout << "Similarity is : " << best_similarity << endl;
    std::cout << "Duration MPQA:" <<  duration.count()/1000000.0 << std::endl;
    // cout << "Most similar range is : " << fair_range.start << " "  << fair_range.end << endl;
    print_range(fair_range);
    cout << "Opt_intersect : " << opt_intersect << endl; 
    cout << "Opt_Union : " << opt_union << endl; 
    cout << "Similarity : " << check_similarity(database, r, fair_range) << endl;


    vector<double> st1, e1;
    fair_range=Range(st1, e1); // Place holder
    best_similarity=0.;
    cout << "=================================================" << endl;
    cout << "Early stop criteria enabled" << endl;
    cout << "Similarity expected " << similarity_expected << endl;
    check = true;
    start = std::chrono::high_resolution_clock::now();
    fair_ranges(r, check);
    stop = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    cout << "Points in the input range : " << input_range_size << endl;
    cout << "Similarity is : " << best_similarity << endl;
    std::cout << "Duration MPQA:" <<  duration.count()/1000000.0 << std::endl;
    // cout << "Most similar range is : " << fair_range.start << " "  << fair_range.end << endl;
    print_range(fair_range);
    cout << "Opt_intersect : " << opt_intersect << endl;
    cout << "Opt_Union : " << opt_union << endl;
    cout << "Similarity : " << check_similarity(database, r, fair_range) << endl;

    cout << "=================================================" << endl;
    cout << "Sampling" << endl;
    // Sampling method
    // Reset global variables
    vector<double> st, e;
    fair_range=Range(st, e); // Place holder
    best_similarity=0.;
    start = std::chrono::high_resolution_clock::now();
    vector<Point> sampled_db = sampling(database, 2*max(200, input_range_size)), state;
    cout << "Size of sampled points is " << sampled_db.size() << endl;
    for(int j=0; j<sampled_db.size(); ++j) {
        sampled_db.at(j).sort_dim = dims-1;// Set sort dim
    }
    sort(sampled_db.begin(), sampled_db.end());
    vector<double> dummy_start, dummy_end;
    Range r1 = Range(dummy_start, dummy_end);
    get_fair_range_recursive(sampled_db, diff_val, r1);
    stop = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    cout << "Duration sampling :" <<  duration.count()/1000000.0 << std::endl;
    cout << "Similarity as per sampling is : " << best_similarity << endl;
    cout << "Check similarity : " << check_similarity(database, fair_range, r) << endl;
    cout << "Check fairness : " << check_fairness(database, fair_range, diff_val) << endl;
    print_range(fair_range);
    return 0;
}
