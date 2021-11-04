#include<iostream>
#include<vector>
#include <algorithm>
#include <limits>
#include <map>
#include <stack>
#include <chrono> 
#include <sstream>
#include <string>
#include <fstream>

/*
TODO:
1. Update forward and backward pointers in preprocess
2. Update get_fair_range method to work for weighted case
*/

//int blue_weight = 5, red_weight = -4;

int blue_weight = 2, red_weight = -1;

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
    return color ? blue_weight : red_weight; // color = { true: blue, false: red }
}

void preprocess(vector<Point_1D>* db) {
    // Sort by X
    sort(db->begin(), db->end());

    // Insert dummy start and end
    db->insert(db->begin(), Point_1D(-numeric_limits<double>::infinity(), 0));
    db->push_back(Point_1D(numeric_limits<double>::infinity(), 0));

    // Calculate backward_ptr
    map<int, vector<int>> large_ptr, small_ptr;
    int cumul = 0;
    for(int i = db->size()-1; i > 0 ; --i) {
        if(db->size() - 1 == i) {
            db->at(i).cumulative = cumul;
            continue;
        }
        int cur_col = color_to_value(db->at(i).color);
        if(cur_col < 0) { // Color is RED
            if(small_ptr.count(cumul)) {
                small_ptr.find(cumul)->second.push_back(i+1);
            } else {
                small_ptr.insert(make_pair(cumul, vector<int>{i+1}));
            }
        } else { // Color is BLUE
            if(large_ptr.count(cumul)) {
                large_ptr.find(cumul)->second.push_back(i+1);
            } else {
                large_ptr.insert(make_pair(cumul, vector<int>{i+1}));
            }
        }
        cumul += cur_col;
        db->at(i).cumulative  = cumul;
        
        // Go over respective pointers and remove items
        if(cur_col > 0) { // Color is BLUE, search for small_ptr
            map<int, vector<int>>::iterator it = small_ptr.begin();
            while(it != small_ptr.end()) {
                if(it->first >= cumul) {
                    break;
                }
                for(auto iter=it->second.begin();iter != it->second.end(); ++iter){
                    db->at(*iter).backward_ptr = i;
                }
                auto temp = it++;
                small_ptr.erase(temp);
            }
        } else {
            map<int, vector<int>>::iterator it = large_ptr.upper_bound(cumul);
            while(it != large_ptr.end()) {
                if(it->first == cumul) {
                    continue;
                }
                for(auto iter=it->second.begin();iter != it->second.end(); ++iter){
                    db->at(*iter).backward_ptr = i;
                }
                auto temp = it++;
                large_ptr.erase(temp);
            }
        }
    }

    small_ptr.clear();
    large_ptr.clear();
    // Calculate forward_ptr
    // Calculate cumulative
    cumul = 0;
    for(int i = 0; i < db->size() - 1; ++i) {
        if(0 == i) {
            db->at(i).cumulative = cumul;
            continue;
        }
        int cur_col = color_to_value(db->at(i).color);
        if(cur_col < 0) { // Color is RED
            if(small_ptr.count(cumul)) {
                small_ptr.find(cumul)->second.push_back(i-1);
            } else {
                small_ptr.insert(make_pair(cumul, vector<int>{i-1}));
            }
        } else { // Color is BLUE
            if(large_ptr.count(cumul)) {
                large_ptr.find(cumul)->second.push_back(i-1);
            } else {
                large_ptr.insert(make_pair(cumul, vector<int>{i-1}));
            }
        }
        cumul += cur_col;
        db->at(i).cumulative  = cumul;
        
        // Go over respective pointers and remove items
        if(cur_col > 0) { // Color is BLUE, search for small_ptr
            map<int, vector<int>>::iterator it = small_ptr.begin();
            while(it != small_ptr.end()) {
                if(it->first >= cumul) {
                    break;
                }
                for(auto iter=it->second.begin();iter != it->second.end(); ++iter){
                    db->at(*iter).forward_ptr = i;
                }
                auto temp = it++;
                small_ptr.erase(temp);
            }
        } else {
            map<int, vector<int>>::iterator it = large_ptr.upper_bound(cumul);
            while(it != large_ptr.end()) {
                if(it->first == cumul) {
                    continue;
                }
                for(auto iter=it->second.begin();iter != it->second.end(); ++iter){
                    db->at(*iter).forward_ptr = i;
                }
                auto temp = it++;
                large_ptr.erase(temp);
            }
        }
    }
}

Range get_fair_range(vector<Point_1D>* db, Range input_range, int eps, double similarity(Range, Range)) {
    int disparity = db->at(input_range.end).cumulative - db->at(input_range.start-1).cumulative;
    if(abs(disparity) <= eps) {
        return input_range;
    }
    double best_similarity = 0;
    Range fair_range = input_range;
    stack<int> expand_left, expand_right, shrink_left, shrink_right;
    bool excess_color = disparity > 0; // If disparity > 0 then excess blue else red

    // First expand right
    expand_left.push(input_range.start);
    shrink_left.push(input_range.start);
    expand_right.push(input_range.end);
    shrink_right.push(input_range.end);
    
    // Expand right till either we explore entire array 
    while(expand_right.top() < db->size() - 2) { // Cannot expand beyond borders
        int top = expand_right.top();
        if(db->at(top + 1).color != excess_color) { // Check if the color is what we want
            expand_right.push(top + 1);
        } else {
            if(db->at(top).forward_ptr == -1) { // If no forward pointer found then break
                break;
            }
            expand_right.push(db->at(top).forward_ptr); // Follow forward pointer
        }
        disparity = db->at(expand_right.top()).cumulative - db->at(input_range.start-1).cumulative;
        // cout << "Expanded right to : " << expand_right.top() << endl;
        if(abs(disparity) <= eps) {
            break;
        }
    }

    while(expand_right.size() >= 1) {
        while(expand_left.top() > 1) { // Cannot expand beyond borders
            disparity = db->at(expand_right.top()).cumulative - db->at(expand_left.top()-1).cumulative;
            if(abs(disparity) <= eps) {
                break;
            }
            int top = expand_left.top();
            if(db->at(top - 1).color != excess_color) { // Check if the color is what we want
                expand_left.push(top - 1);
            } else {
                if(db->at(top).backward_ptr == -1) { // If no backward pointer found then break
                    break;
                }
                expand_left.push(db->at(top).backward_ptr); // Follow backward pointer
            }
            // cout << "Expanded left to : " << expand_left.top() << endl;
        }

        while(shrink_left.top() < input_range.end) { // Cannot expand beyond borders
            disparity = db->at(expand_right.top()).cumulative - db->at(shrink_left.top()-1).cumulative;
            if(abs(disparity) <= eps) {
                break;
            }
            int top = shrink_left.top();
            if(db->at(top).color == excess_color) { // Check if the color is what we want
                shrink_left.push(top + 1);
            } else {
                if(db->at(top - 1).forward_ptr == -1) { // If no backward pointer found then break
                    break;
                }
                // Follow forward pointer and add1 because we do not want to include the extra color
                shrink_left.push(db->at(top - 1).forward_ptr + 1);
            }
            // cout << "Shrinked left to : " << shrink_left.top() << endl;
        }

        disparity = db->at(expand_right.top()).cumulative - db->at(expand_left.top()-1).cumulative;
        if(abs(disparity) <= eps) {
            Range r = Range(expand_left.top(), expand_right.top());
            double sim = similarity(r, input_range);
            if(sim > best_similarity) {
                best_similarity = sim;
                fair_range = r;
            }
        }
        disparity = db->at(expand_right.top()).cumulative - db->at(shrink_left.top()-1).cumulative;
        if(abs(disparity) <= eps) {
            Range r = Range(shrink_left.top(), expand_right.top());
            double sim = similarity(r, input_range);
            if(sim > best_similarity) {
                best_similarity = sim;
                fair_range = r;
            }
        }
        expand_right.pop(); // Remove element from top of stack  
        /* if(expand_right.size()) {
            cout << "Expand right removed to : " << expand_right.top() << endl;
        }*/
    }

    // Look to shrink_right
    while(shrink_right.top() > input_range.start) {
        disparity = db->at(shrink_right.top()).cumulative - db->at(expand_left.top()-1).cumulative;
        if(abs(disparity) <= eps) {
            Range r = Range(expand_left.top(), shrink_right.top());
            double sim = similarity(r, input_range);
            if(sim > best_similarity) {
                best_similarity = sim;
                fair_range = r;
            }
        }
        disparity = db->at(shrink_right.top()).cumulative - db->at(shrink_left.top()-1).cumulative;
        if(abs(disparity) <= eps) {
            Range r = Range(shrink_left.top(), shrink_right.top());
            double sim = similarity(r, input_range);
            if(sim > best_similarity) {
                best_similarity = sim;
                fair_range = r;
            }
        }
        
        if(expand_left.size() == 1 && shrink_left.size() == 1) {
            break;
        }

        int top = shrink_right.top();
        if(db->at(top).color == excess_color) {
            shrink_right.push(top - 1);
        } else {
            if(db->at(top + 1).backward_ptr == -1) { // If no backward pointer found then break
                break;
            }
            shrink_right.push(db->at(top + 1).backward_ptr - 1); // Follow backward pointer
        }
        // cout << "Skrinked right to : " << shrink_right.top() << endl;

        while(shrink_left.size() > 1) { // Cannot expand beyond borders
            if(shrink_left.top() <= shrink_right.top()) {
                disparity = db->at(shrink_right.top()).cumulative - db->at(shrink_left.top()-1).cumulative;
                // Check if we have overrun because of some constraints
                if(abs(disparity) <= eps) {
                    Range r = Range(shrink_left.top(), shrink_right.top());
                    double sim = similarity(r, input_range);
                    if(sim > best_similarity) {
                        best_similarity = sim;
                        fair_range = r;
                    }
                }
                if(disparity >= 0 == excess_color && abs(disparity) > eps) {
                    break;
                }
            }
            shrink_left.pop();
            // cout << "Shrink left removed to : " << shrink_left.top() << endl;
        }

        while(expand_left.size() > 1) { // Cannot expand beyond borders
            disparity = db->at(shrink_right.top()).cumulative - db->at(expand_left.top()-1).cumulative;
            // Check if we have overrun because of some constraints
            if(abs(disparity) <= eps) {
                Range r = Range(expand_left.top(), shrink_right.top());
                double sim = similarity(r, input_range);
                if(sim > best_similarity) {
                    best_similarity = sim;
                    fair_range = r;
                }
            }
            if(disparity >= 0 == excess_color && abs(disparity) > eps) {
                break;
            }
            expand_left.pop();
            // cout << "Expand right removed to : " << expand_left.top() << endl;
        }
        if(expand_left.size() == 0 && shrink_left.size() == 0) {
            break;
        }
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


vector<Point_1D>* read_urbangb(string filename) {
    vector<Point_1D>* database = new vector<Point_1D>();
    ifstream file_handler(filename);
    if (!file_handler.is_open())  // check file is open, quit if not
    {
        std::cerr << "failed to open file\n";
        return database;
    }

    double x, y;
    int color;
    string line;
    while (getline(file_handler, line))
    {
        istringstream iss(line);
        if (!(iss >> x >> y >> color)) {
            cout << "Issues parsing file!" << endl;
            break;
        } // error
        Point_1D p(x, color == 1);
        database->push_back(p);
    }
    return database;
}

double Jaccard_similarity(Range r1, Range r2){
    int union_r1_r2 = max(r1.end, r2.end) - min(r1.start, r2.start)  + 1;
    int int_r1_r2 = min(r1.end, r2.end) - max(r1.start, r2.start) + 1;
    double sim = double(int_r1_r2)/double(union_r1_r2);
    //cout << "[" << r1.start << "," << r1.end << "], [" << r2.start << "," << r2.end << "] : " << sim << endl;  
    return max(0., sim); // Ensure that atleast 0 is returned as similarity
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

int main() {
    //int eps=500;
    //vector<Point_1D>* db = read_texas_tribune_db("data/texas_tribune_ethnicity_0s5_1s4.csv");
    //double salary_start = 20000., salary_end=45000.;
    int eps=350;
    vector<Point_1D>* db = read_urbangb("data/10k_new.csv");
    double salary_start = -0.1, salary_end=0.2;
    std::cout << "Epsilon value set to : " << eps << std::endl;
    std::cout << "Initial range value set to : " << salary_start << ":" << salary_end << std::endl;
    auto start = std::chrono::high_resolution_clock::now();
    preprocess(db);
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << "Duration weighted SPQA preprocessing: " <<  duration.count()/1000000.0 << std::endl;
    start = std::chrono::high_resolution_clock::now();
    int start_index = binary_search(db, 1, db->size()-1, salary_start);
    int end_index = binary_search(db, 1, db->size()-1, salary_end);
    Range out = get_fair_range(db, Range(start_index, end_index), eps, Jaccard_similarity);
    stop = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << "Duration weighted SPQA : " <<  duration.count()/1000000.0 << std::endl;
    std::cout << "Range indices are " << out.start << ":" << out.end << std::endl;
    std::cout << "Range is " << db->at(out.start).X << ":" << db->at(out.end).X << std::endl;
    std::cout << "Cumulative is " << db->at(out.start).cumulative << ":" << db->at(out.end).cumulative << std::endl;
    /*
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
            std::cout << "Duration weighted SPQA : " <<  duration.count()/1000000.0 << std::endl;
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


            /*
            start = std::chrono::high_resolution_clock::now();
            int start_index = binary_search(db, 1, db->size()-1, salary_start);
            int end_index = binary_search(db, 1, db->size()-1, salary_end);
            Range out = get_fair_range(db, Range(start_index, end_index), eps, Jaccard_similarity);
            stop = std::chrono::high_resolution_clock::now();
            duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
            std::cout << "Duration weighted SPQA : " <<  duration.count()/1000000.0 << std::endl;
            cout << "Fair range using our algorithm is : (" << out.start << " " << out.end << ")" << endl;
            cout << "Fair range using our algorithm is : (" << db->at(out.start).X << " " << db->at(out.end).X << ")" << endl;
            */
    delete db;
}

int main1() {
    int eps = 2;
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
            cout << "Input range is : (" << j << " " << i << ")" << endl << "===================================" << endl;
            Range out = get_fair_range(db, Range(j, i), eps, Jaccard_similarity);
            cout << "Fair range using our algorithm is : (" << out.start << " " << out.end << ")" << endl;
            Range out_brute = brute_force(db, Range(j, i), eps, Jaccard_similarity);
            cout << "Fair range using brute force algorithm is : (" << out_brute.start << " " << out_brute.end << ")" << endl;
            cout << "===================================" << endl;
        }
    }
    delete db;
    return 0;
}
