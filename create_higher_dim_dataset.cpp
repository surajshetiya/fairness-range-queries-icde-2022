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
int dims=3, points_in_db = 10000;


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

void write_file(vector<Point> db, string filename) {
    ofstream myfile;
    myfile.open(filename, ios::out);
    for(int i=0; i < db.size() ;++i) {
        for(int j=0; j < db.at(i).X.size();++j) {
            myfile << db.at(i).X.at(j) << " ";
        }
        myfile << (db.at(i).color ? 1 : 0) << endl;
    }
    myfile.close();

}


int main() {

    vector<Point> database = generate_database(points_in_db, dims);
    write_file(database, string("data/uniform_3d.csv"));
    return 0;
}
