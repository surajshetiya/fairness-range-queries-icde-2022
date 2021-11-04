#include<iostream>
#include<vector>
#include <algorithm>
#include <limits>
#include <map>
#include <stack>
#include <sstream>
#include <string>
#include <fstream>
#include <chrono> 

using namespace std;

class Point_1D {
    public:
    
    double X;
    bool color; // color = { true: blue, false: red }
    int cumulative;
    int forward_ptr, backward_ptr;

    Point_1D(double coordinate, bool col) {
        X = coordinate;
        color = col;
        forward_ptr = backward_ptr = -1;
        cumulative = 0;
    }

    bool operator < (const Point_1D& other) const
    {
        return X < other.X;
    }
};

class Range {
    public:
    int start, end;

    Range(int st, int e) {
        start = st;
        end = e;
    }
};

int color_to_value(bool color) {
    return color ? 1 : -1; // color = { true: blue, false: red }
}

void preprocess(vector<Point_1D>* db) {
    // Sort by X
    sort(db->begin(), db->end());

    // Insert dummy start and end
    db->insert(db->begin(), Point_1D(-numeric_limits<double>::infinity(), 0));
    db->push_back(Point_1D(numeric_limits<double>::infinity(), 0));

    // Calculate backward_ptr
    map<int, vector<int>> ptr_map;
    int cumul = 0;
    for(int i = db->size()-1; i > 0 ; --i) {
        if(db->size() - 1 == i) {
            db->at(i).cumulative = cumul;
            continue;
        }
        int cur_col = color_to_value(db->at(i).color);
        cumul += cur_col;
        db->at(i).cumulative  = cumul;
        if(ptr_map.count(cumul-2*cur_col)) {
            ptr_map.find(cumul-2*cur_col)->second.push_back(i+1);
        } else {
            ptr_map.insert(make_pair(cumul-2*cur_col, vector<int>{i+1}));
        }
        if(ptr_map.count(cumul)) {
            for(int j: ptr_map.find(cumul)->second) {
                db->at(j).backward_ptr = i;
            }
            ptr_map.erase(cumul);
        }
    }

    ptr_map.clear();
    // Calculate forward_ptr
    // Calculate cumulative
    cumul = 0;
    for(int i = 0; i < db->size() - 1; ++i) {
        if(i == 0) {
            db->at(i).cumulative = cumul;
            continue;
        }
        int cur_col = color_to_value(db->at(i).color);
        cumul += cur_col;
        db->at(i).cumulative  = cumul;
        if(ptr_map.count(cumul-2*cur_col)) {
            ptr_map.find(cumul-2*cur_col)->second.push_back(i-1);
        } else {
            ptr_map.insert(make_pair(cumul-2*cur_col, vector<int>{i-1}));
        }
        if(ptr_map.count(cumul)) {
            for(int j: ptr_map.find(cumul)->second) {
                db->at(j).forward_ptr = i;
            }
            ptr_map.erase(cumul);
        }
    }
    ptr_map.clear();
}

Range get_fair_range(vector<Point_1D>* db, Range input_range, int eps, double similarity(Range, Range)) {
    int disparity = db->at(input_range.end).cumulative - db->at(input_range.start-1).cumulative;
    if(abs(disparity) <= eps) {
        return input_range;
    }
    double best_similarity = 0;
    Range fair_range = input_range;
    int disparity_difference = abs(disparity) - eps;
    int jumps_to_make = disparity_difference;
    stack<int> start_left, start_right, end_left, end_right;
    bool excess_color = disparity > 0;
    start_left.push(input_range.start);
    start_right.push(input_range.start);
    end_left.push(input_range.end);
    end_right.push(input_range.end);

    while(jumps_to_make) {
        int top_index = start_right.top();
        if(db->at(top_index).color == excess_color) {
            start_right.push(top_index + 1);
        } else {
            start_right.push(db->at(top_index-1).forward_ptr + 1);
        }
        --jumps_to_make;
    }

    jumps_to_make = disparity_difference;
    while(jumps_to_make) {
        int sr = start_right.top();
        int el = end_left.top();
        int er = end_right.top();
        Range sr_el = Range(sr, el);
        double similarity_sr_el = similarity(input_range, sr_el);
        if(similarity_sr_el > best_similarity) {
            best_similarity = similarity_sr_el;
            fair_range = sr_el;
        }
        if(end_right.size() + start_right.size() == disparity_difference + 2) {
            Range sr_er = Range(sr, er);
            double similarity_sr_er = similarity(input_range, sr_er);
            if(similarity_sr_er > best_similarity) {
                best_similarity = similarity_sr_er;
                fair_range = sr_er;
            }
        }
        if(er < db->size() - 1) {
            if(db->at(er + 1).color != excess_color) {
                end_right.push(er + 1);
            } else if(db->at(er).forward_ptr != -1) {
                end_right.push(db->at(er).forward_ptr);
            }
        }
        if(db->at(el).color == excess_color) {
            end_left.push(el - 1);
        } else {
            end_left.push(db->at(el + 1).backward_ptr - 1);
        }
        start_right.pop();
        --jumps_to_make;
    }

    jumps_to_make = disparity_difference + 1;
    while(jumps_to_make) {
        int sl = start_left.top();
        int el = end_left.top();
        int er = end_right.top();
        Range sl_el = Range(sl, el);
        double similarity_sl_el = similarity(input_range, sl_el);
        if(similarity_sl_el > best_similarity) {
            best_similarity = similarity_sl_el;
            fair_range = sl_el;
        }
        if(end_right.size() + start_left.size() == disparity_difference + 2) {
            Range sl_er = Range(sl, er);
            double similarity_sl_er = similarity(input_range, sl_er);
            if(similarity_sl_er > best_similarity) {
                best_similarity = similarity_sl_er;
                fair_range = sl_er;
            }
            end_right.pop();
        }
        end_left.pop();
        if(sl > 1) {
            if(db->at(sl-1).color != excess_color) {
                start_left.push(sl-1);
            } else if(db->at(sl).backward_ptr != -1){
                start_left.push(db->at(sl).backward_ptr);
            }
            else{
                break;
            }
        } else {
            break;
        }
        --jumps_to_make;
    }
    cout << "Similarity : " << best_similarity << endl;
    return fair_range;
}

Range brute_force(vector<Point_1D>* db, Range input_range, int eps, double similarity(Range, Range)) {
    int disparity = db->at(input_range.end).cumulative - db->at(input_range.start-1).cumulative;
    if(abs(disparity) < eps) {
        return input_range;
    }
    double best_similarity = 0;
    Range fair_range = input_range;
    for(int i=1; i < db->size() - 1; ++i) {
        for(int j=1; j <= i; ++j) {
            int disparity = db->at(i).cumulative - db->at(j-1).cumulative;
            Range r_ij = Range(j, i);
            if(abs(disparity) <= eps && similarity(r_ij, input_range) > best_similarity) {
                best_similarity = similarity(r_ij, input_range);
                fair_range = r_ij;
            }
        }
    }
    cout << "Similarity : " << best_similarity << endl;
    return fair_range;
}

double Jaccard_similarity(Range r1, Range r2){
    int union_r1_r2 = max(r1.end, r2.end) - min(r1.start, r2.start)  + 1;
    int int_r1_r2 = min(r1.end, r2.end) - max(r1.start, r2.start) + 1;
    double sim = double(int_r1_r2)/double(union_r1_r2);
    //cout << "[" << r1.start << "," << r1.end << "], [" << r2.start << "," << r2.end << "] : " << sim << endl;  
    return sim;
}

// Compare the X[dim] with elem to perform binary search
int binary_search(vector<Point_1D>* db, int start, int end, double elem) {
    if(end < start) {
        // Should not happen
        return start;
    }
    if(start == end) {
        return start;
    }
    int mid = (start + end)/2;
    if(db->at(mid).X == elem) {
        return mid;
    } else if(db->at(mid).X < elem) {
        return binary_search(db, mid+1, end, elem);
    } else{
        return binary_search(db, start, mid-1, elem);
    }
}

vector<Point_1D>* read_texas_tribune_db(string filename) {
    vector<Point_1D>* database = new vector<Point_1D>();
    ifstream file_handler(filename);
    if (!file_handler.is_open())  // check file is open, quit if not
    {
        std::cerr << "failed to open file\n";
        return database;
    }

    double salary;
    int color;
    string line;
    while (getline(file_handler, line))
    {
        istringstream iss(line);
        if (!(iss >> color >> salary)) {
            cout << "Issues parsing file!" << endl;
            break; 
        } // error
        Point_1D p(salary, color == 1); 
        database->push_back(p);
    }
    return database;
}

int main() {
    int eps=0;
    //vector<Point_1D>* db = read_texas_tribune_db("data/texas_tribune.csv");
    vector<Point_1D>* db = read_texas_tribune_db("data/uniform_1k.csv");
    double salary_start = 0.35, salary_end = 0.62;

    preprocess(db);

    int start_index = binary_search(db, 1, db->size()-1, salary_start);
    int end_index = binary_search(db, 1, db->size()-1, salary_end);

    Range out = get_fair_range(db, Range(start_index, end_index), eps, Jaccard_similarity);
    cout << "Value of epsilon is : " << eps << endl;
    cout << "Value of disparity is : " << abs(db->at(end_index).cumulative - db->at(start_index-1).cumulative) << endl;
    
    cout << "Fair range indices using our algorithm is : (" << out.start << " , " << out.end << ")" << endl;
    cout << "Fair range using our algorithm is : (" << db->at(out.start).X << " , " << db->at(out.end).X << ")" << endl;
    cout << "Value of disparity using our algorithm is : " << abs(db->at(out.end).cumulative - db->at(out.start).cumulative) << endl;
    cout << "Total points in output is : " << (out.end-out.start+1) << endl;
   /*
    double salary_start = 55000., salary_end=90000.;
    auto start = std::chrono::high_resolution_clock::now();
    preprocess(db);
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << "Duration SPQA Preprocessing: " <<  duration.count()/1000000.0 << std::endl;
    for(int i=5000; i < 100000; i+=5000) {
        salary_start = i;
        for(int j=i+5000; j < 100000; j+=5000) {
            cout << "-----------------------------------" << endl;
            salary_end = j;
            start = std::chrono::high_resolution_clock::now();
            int start_index = binary_search(db, 1, db->size()-1, salary_start);
            int end_index = binary_search(db, 1, db->size()-1, salary_end);
            Range out = get_fair_range(db, Range(start_index, end_index), eps, Jaccard_similarity);
            stop = std::chrono::high_resolution_clock::now();
            duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
            std::cout << "Duration SPQA : " <<  duration.count()/1000000.0 << std::endl;
            int disparity = abs(db->at(end_index).cumulative - db->at(start_index-1).cumulative);
            cout << "Value of epsilon is : " << eps << endl;
            cout << "Value of disparity is : " << disparity << endl;
            cout << "Input range is : (" << salary_start << " , " << salary_end << ")" << endl;
            cout << "Input range indices is : (" << start_index << " , " << end_index << ")" << endl;
            cout << "Total points in input is : " << (end_index-start_index+1) << endl;
            cout << "Fair range indices using our algorithm is : (" << out.start << " , " << out.end << ")" << endl;
            cout << "Fair range using our algorithm is : (" << db->at(out.start).X << " , " << db->at(out.end).X << ")" << endl;
            cout << "Total points in output is : " << (out.end-out.start+1) << endl;
            cout << "-----------------------------------" << endl;
        }
    }
    */
    delete db;
}

int main1() {
    vector<Point_1D>* db = new vector<Point_1D>();
    db->push_back(Point_1D(0, false));
    db->push_back(Point_1D(1, false));
    db->push_back(Point_1D(3, false));
    db->push_back(Point_1D(5, false));
    db->push_back(Point_1D(9, false));
    db->push_back(Point_1D(10, false));
    db->push_back(Point_1D(12, false));
    db->push_back(Point_1D(13, false));
    db->push_back(Point_1D(2, true));
    db->push_back(Point_1D(4, true));
    db->push_back(Point_1D(6, true));
    db->push_back(Point_1D(7, true));
    db->push_back(Point_1D(8, true));
    db->push_back(Point_1D(11, true));

    preprocess(db);
    for(int i=0; i<db->size(); ++i) {
        Point_1D p = db->at(i);
        cout << p.X << " " << (p.color ? "blue" : "red") << " " << p.cumulative << " " << 
            p.forward_ptr << " " << p.backward_ptr << endl;
    }
    
    for(int i=1; i < db->size() - 1; ++i) {
        for(int j=1; j <= i; ++j) {
            cout << j << " " << i << endl;
            Range out = get_fair_range(db, Range(j, i), 0, Jaccard_similarity);
            cout << out.start-1 << " " << out.end-1 << endl;
            Range out_brute = brute_force(db, Range(j, i), 0, Jaccard_similarity);
            cout << out_brute.start-1 << " " << out_brute.end-1 << endl;
        }
    }
    Range out_brute = get_fair_range(db, Range(9, 14), 0, Jaccard_similarity);
    delete db;
    return 0;
}
