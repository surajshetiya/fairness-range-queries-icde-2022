#include<iostream>
#include<vector>
#include <algorithm>
#include <random>
#include<unordered_map>

using namespace std;

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
        return X[0] < other.X[0];
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

// Number of dimensions
int dims=2;
vector<Point> database = generate_database(100000, dims);

class Range {
    public:
    vector<double> start, end;

    Range(vector<double> st, vector<double> e) {
        start = st;
        end = e;
    }
};
vector<double> _start = {400.*double(rand())/double(RAND_MAX), 400.*double(rand())/double(RAND_MAX)};
vector<double> _end = {400.*double(rand())/double(RAND_MAX), 400.*double(rand())/double(RAND_MAX)};
Range query = Range(_start, _end);
Range opt = Range(vector<double>(), vector<double>());
double best_jaccard = 0;

int color_to_value(bool color) {
    return color ? 1 : -1; // color = { true: blue, false: red }
}

void mark_points(vector<Point> db) {
    int dim = db.at(0).X.size();
    for(int j=0; j < db.size(); ++j) {
        int i=0;
        for(; i < dim; ++i) {
            double coordinate_val = db.at(j).X.at(i);
            if(coordinate_val <= query.start.at(i) || coordinate_val >= query.start.at(i)) {
                break;
            }
        }
        if(i == dim) {
            db.at(j).intersect = true;
        } else {
            db.at(j).intersect = false;
        }
    }
}

void preprocess() {
    // Sort by first dimension
    sort(database.begin(), database.end());
    mark_points(database);
}

double get_jaccard_similarity(Range output_range) {
    int in_range = 0, in_query_range = 0, dim = database.at(0).X.size();
    for(int i=0; i < database.size(); ++i) {
        int j=0;
        for(; j < dim; ++j) {
            double coordinate_val = database.at(i).X.at(j);
            if(coordinate_val <= output_range.start.at(j) || coordinate_val >= output_range.end.at(j)) {
                break;
            }
            if(coordinate_val <= query.start.at(j) || coordinate_val >= query.end.at(j)) {
                break;
            }
        }
        if(j == dim) {
            in_range++;
            in_query_range += database.at(i).intersect ? 1 : 0;
        }
    }
    if(in_range == 0) {
        return 0;
    }
    cout << double(in_query_range)/double(in_range) << endl;
    return double(in_query_range)/double(in_range);
}

void get_fair_range_last_dim_new(vector<Point>& db, int eps, Range output_range) {
    vector<int> cumulative;
    int cumulative_total = 0;
    unordered_map<int, vector<int>> map = unordered_map<int, vector<int>>();
    int min_cumlative = 0, max_cumlative = 0;
    for(int i=0; i< db.size(); ++i) {
        cumulative_total += color_to_value(db.at(i).color);
        cumulative.push_back(cumulative_total);
        if(cumulative_total < min_cumlative) {
            min_cumlative = cumulative_total;
        }
        if(cumulative_total > max_cumlative) {
            max_cumlative = cumulative_total;
        }
        auto key = map.find(cumulative_total);
        if(key != map.end()) {
            key->second.push_back(i);
        } else {
            vector<int> val = vector<int>();
            val.push_back(i);
            map.insert(make_pair(cumulative_total, val));
        }
    }
    for(int i=min_cumlative; i <= max_cumlative; ++i) {
        vector<int> vec_i = map.find(i)->second;
        for(int j=min_cumlative; j <= min_cumlative + eps; ++j) {
            vector<int> vec_j = map.find(j)->second;
        }
    }
}


void get_fair_range_last_dim(vector<Point>& db, int eps, Range output_range) {
    // TODO: Make cumulative output sensitive - make it O(output sensitive)
    vector<int> cumulative;
    int cumulative_total = 0;
    for(int i=0; i< db.size(); ++i) {
        cumulative_total += color_to_value(db.at(i).color);
        cumulative.push_back(cumulative_total);
    }
    for(int i=0; i< db.size(); ++i) {
        for(int j=i+1; j< db.size(); ++j) {
            if(abs(cumulative.at(i)-cumulative.at(j)) < eps) {
                Range _output_range = output_range; // Copy constructor
                _output_range.start.insert(_output_range.start.begin(), i);
                _output_range.end.insert(_output_range.end.begin(), j);
                double similarity = get_jaccard_similarity(_output_range);
                if(similarity > best_jaccard) {
                    best_jaccard = similarity;
                    opt = _output_range; // Copy constructor
                    cout << best_jaccard << endl;
                }
            }
        }
    }
}

void get_fair_range_recursive(vector<Point> &db, int eps, int dim, Range output_range) {
    cout << dim << " dimensions" << endl;
    cout << db.size() << endl;
    for(int i=0; i < db.size(); ++i) {
        // i is starting point for dimension dim.
        cout << i << endl;
        double start_dim = db.at(i).X.at(dim);
        for(int j=0; j < db.size(); ++j) {
            cout << j << endl;
            if(i == j || db.at(j).X.at(dim) < start_dim) {
                continue;
            }
            vector<Point> db_copy;
            double end_dim = db.at(j).X.at(dim);
            for(int k=0; k < db.size(); ++k) {
                // Remove all objects whose value for dimension dim is 
                // less than that of object i or greater than object j
                if(db.at(k).X.at(dim) > end_dim || db.at(k).X.at(dim) < start_dim) {
                    continue;
                }
                db_copy.push_back(Point(db.at(k).X, db.at(k).color, db.at(k).id));
            }
            // See if any dimensions are left to explore.
            // If not call the last_dim
            if(dim == dims - 1) {
                Range _output_range = output_range; // Copy constructor
                _output_range.start.push_back(db.at(i).X.at(dim));
                _output_range.end.push_back(db.at(j).X.at(dim));
                get_fair_range_last_dim(db_copy, eps, _output_range);
            } else {
                Range _output_range = output_range; // Copy constructor
                _output_range.start.push_back(db.at(i).X.at(dim));
                _output_range.end.push_back(db.at(j).X.at(dim));
                get_fair_range_recursive(db_copy, eps, dim + 1, _output_range);
            }
        }
    }
}

vector<int> generate_samples(int sample_size, int db_size) {
    default_random_engine generator;
    uniform_int_distribution<int> distribution(0,db_size-1);
    vector<int> samples;
    for(int i=0; i< sample_size;++i) {
        samples.push_back(distribution(generator));
    }
    return samples;
}

void get_fair_range(int eps) {
    int total_samples = 150; // 500;
    vector<int> samples = generate_samples(total_samples, database.size());
    vector<Point> samples_database = vector<Point>();
    for(int i=0; i<total_samples; ++i) {
        samples_database.push_back(database.at(samples.at(i)));
    }
    get_fair_range_recursive(samples_database, eps, 1, Range(vector<double>(), vector<double>()));
}

int main() {    
    preprocess();
    get_fair_range(5);
    for(int i=0; i< opt.start.size(); ++i) {
        cout << opt.start.at(i) << " ";
    }
    cout << endl;
    for(int i=0; i< opt.end.size(); ++i) {
        cout << opt.end.at(i) << " ";
    }
    cout << endl;
    return 0;
}
