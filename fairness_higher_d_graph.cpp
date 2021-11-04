#include<iostream>
#include<vector>
#include <algorithm>
#include <random>
#include<unordered_map>
#include <set>

using namespace std;

// Number of dimensions
int dims=2, diff_val=1, points_in_db = 1000;

class Point {
    public:
    
    int id;
    vector<double> X;
    bool color; // color = { true: blue, false: red }
    bool intersect; // intersect with query rectangle

    Point(vector<double> coordinates, bool col, int identifier) {
        X = coordinates;
        color = col;
        id = identifier;
    }

    bool operator < (const Point& other) const
    {
        // Use last dimension for comparision, sorting
        return X[X.size()-1] < other.X[X.size()-1];
    }
};

vector<Point> generate_database(int count, int dim) {
    vector<Point> db;
    for(int i=0; i< count;++i) {
        vector<double> pnt;
        for(int j=0; j< dim;++j) {
            pnt.push_back(400.*double(rand())/double(RAND_MAX));
        }
        db.push_back(Point(pnt, rand() % 2 == 0, i));
    }
    return db;
}

vector<Point> database = generate_database(points_in_db, dims);

class Range {
    public:
    vector<double> start, end;

    Range(vector<double> st, vector<double> e) {
        start = st;
        end = e;
    }
};

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
    return color ? 1 : -1; // color = { true: blue, false: red }
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
void get_fair_range_one_dim(vector<Point>& db, int diff, set<Range, decltype(compare_range)*>& fair_ranges, Range range) {
    int size_val = db.size(), total = 0, min = size_val, max = -size_val, dims = range.start.size();
    vector<int> cumultv;
    unordered_map<int, vector<int>> locations_map;
    for(int i=0; i<size_val; ++i) {
        total += color_to_value(db.at(i).color);
        cumultv.push_back(total);
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
        int end_range = (min + diff) < max ? (min + diff) : max;
        for(int j = i; j <= end_range; ++j) {
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
                    Range fair_range = range;

                    // Append the last coordinate to the range
                    fair_range.start.push_back(db.at(index_i).X.at(dims));
                    fair_range.end.push_back(db.at(index_j).X.at(dims));
                    fair_ranges.insert(fair_range);
                }
            }
        }
    }
}

void get_fair_range_recursive(vector<Point> db, int diff, set<Range, decltype(compare_range)*>& fair_ranges, Range range){
    int dim = range.start.size(), db_dim = db.at(0).X.size();
    bool call_last = (dim == (db_dim - 2));
    for(int i=0; i<db.size(); ++i) {
        cout << "Outer loop : " << i << endl;
        double val_i = db.at(i).X.at(dim);
        for(int j=0; j<db.size(); ++j) {
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
                get_fair_range_one_dim(new_db, diff, fair_ranges, new_range);
            } else {
                get_fair_range_recursive(new_db, diff, fair_ranges, new_range);
            }
        }
    }
}

/*
1. Get unique set of ranges : Mark them with ids
2. For each id find the set of elements which are a delta from optimal 
3. Apply greedy hitting set algorithm
*/
void get_unique_range_one_dim(vector<Point> db, set<Range, decltype(compare_range)*>& unique_ranges, Range range) {
    int dim = db.at(0).X.size() - 1;
    for(unsigned int i=0; i<db.size(); ++i) {
        double val_i = db.at(i).X.at(dim);
        for(unsigned int j=i+1; j<db.size(); ++j) {
            double val_j = db.at(j).X.at(dim);
            double min_ij = val_i < val_j ? val_i : val_j, max_ij = val_i > val_j ? val_i : val_j;
            Range copy_range = range;
            copy_range.start.push_back(min_ij);
            copy_range.end.push_back(max_ij);
            unique_ranges.insert(copy_range);
        }
    }
}

void get_unique_recursive(vector<Point> db, set<Range, decltype(compare_range)*>& unique_ranges, Range range) {
    int dim = range.start.size(), db_dim = db.at(0).X.size();
    bool call_last = (dim == (db_dim - 2));
    for(int i=0; i<db.size(); ++i) {
        cout << "Outer loop : " << i << endl;
        double val_i = db.at(i).X.at(dim);
        for(int j=0; j<db.size(); ++j) {
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
                get_unique_range_one_dim(new_db, unique_ranges, new_range);
            } else {
                get_unique_recursive(new_db, unique_ranges, new_range);
            }
        }
    }
}

 void greedy(vector<Range>& fair_representative, vector<Point>& db, vector<Range>& fair_range, vector<Range>& unique_range, double diff) {
    int fair_range_size = fair_range.size(), unique_range_size = unique_range.size(), db_size = db.size(), dim=db.at(0).X.size(), id=0;
    vector<set<int>> hs_form;
    vector<int> fr_set_size, hs_form_count;
    for(int j=0; j < fair_range_size; ++j) {
        set<int> s;
        hs_form.push_back(s);
        hs_form_count.push_back(0);
    }
    cout << "Allocated empty hitting sets" << endl;
    for(auto range_iter=fair_range.begin(); range_iter != fair_range.end(); ++range_iter) {
        int total;
        for(int j=0; j < db_size; ++j) {
            int k=0;
            for(; k < dim; ++k) {
                double coordinate = db.at(j).X.at(k);
                if(!(coordinate >= range_iter->start.at(k) && coordinate <= range_iter->end.at(k))) {
                    break;
                }
            }
            if(k == dim) {
                ++total;
            }
        }
        fr_set_size.push_back(total);
    }
    cout << "Computed fair range set sizes" << endl;
    int max_hs_count = 0, hs_index = -1;
    bool complete = false;
    for(auto range_iter=unique_range.begin(); range_iter != unique_range.end(); ++range_iter, ++id) {
        cout << "Processing id : " << id << endl;
        vector<Point> filtered;
        for(int j=0; j < db_size; ++j) {
            int k=0;
            for(; k < dim; ++k) {
                double coordinate = db.at(j).X.at(k);
                if(!(coordinate >= range_iter->start.at(k) && coordinate <= range_iter->end.at(k))) {
                    break;
                }
            }
            if(k == dim) {
                filtered.push_back(db.at(j));
            }
        }
        // Check which all fair ranges does filtered set of points belong to
        int filtered_size = filtered.size();
        vector<double> similar;
        int fr_id = 0;
        double max_similar = 0;
        for(auto fr_iter = fair_range.begin(); fr_iter != fair_range.end(); ++fr_iter, ++fr_id) {
            int intersect = 0;
            for(int j=0; j < filtered_size; ++j) {
                int k=0;
                for(; k < dim; ++k) {
                    if(!(filtered.at(j).X.at(k) <= fr_iter->start.at(k)) && (filtered.at(j).X.at(k) >= fr_iter->end.at(k))) {
                        break;
                    }
                }
                if(k == dim) {
                    ++intersect;
                }
            }
            int set_size = fr_set_size.at(fr_id);
            double val = (1. * intersect)/(set_size + filtered_size - set_size );
            similar.push_back(val);
            if(val > max_similar) {
                max_similar = val;
            }
        }
        double similarity_accepted = max_similar - diff;
        for(int i=0; i < fair_range_size; ++i) {
            if(similar.at(i) >= similarity_accepted) {
                hs_form.at(i).insert(id);
                ++hs_form_count.at(i);
                if(hs_form_count.at(i) > max_hs_count) {
                    max_hs_count = hs_form_count.at(i);
                    hs_index = i;
                }
            }
        }
    }
    vector<int> hs;
    cout << "Starting hitting set computation" << endl;
    while(!complete) {
        hs.push_back(hs_index);
        fair_representative.push_back(fair_range.at(hs_index));
        set<int> copy_ideal = hs_form.at(hs_index);
        max_hs_count = 0;
        hs_index = -1;
        complete= true;
        for(int i=0; i < fair_range_size; ++i) {
            set<int> temp;
            set_difference(hs_form.at(i).begin(), hs_form.at(i).end(), copy_ideal.begin(), copy_ideal.end(), inserter(temp, temp.begin()));
            hs_form.at(i) = temp;
            if(temp.size() > 0) {
                complete = false;
            }
            if(temp.size() > max_hs_count) {
                max_hs_count = temp.size();
                hs_index = i;
            }
        }
    }
    cout << "Size of set is " << fair_representative.size() << endl;
}

void fair_ranges() {
    // Sort by first dimension
    sort(database.begin(), database.end());
    vector<double> start, end;
    Range range(start, end);
    set<Range, decltype(compare_range)*> fair_range(compare_range);
    get_fair_range_recursive(database, diff_val, fair_range, range);
    cout << "Fair ranges size : " << fair_range.size() << endl;
    return;
    set<Range, decltype(compare_range)*> unique_range(compare_range);
    vector<double> start1, end1;
    Range range1(start1, end1);
    get_unique_recursive(database, unique_range, range1);
    cout << "Unique ranges size : " << unique_range.size() << endl;
    vector<Range> unique_vector, fair_vector;
    for(auto iter=unique_range.begin(); iter != unique_range.end(); ++iter) {
        unique_vector.push_back(*iter);
    }
    for(auto iter=fair_range.begin(); iter != fair_range.end(); ++iter) {
        fair_vector.push_back(*iter);
    }
    vector<Range> fair_repr;
    double delta = 0.1;
    greedy(fair_repr, database, fair_vector, unique_vector, delta);
}

int main() {
    fair_ranges();
    return 0;
}