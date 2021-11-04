#include<iostream>
#include<unordered_map>
#include<boost/heap/pairing_heap.hpp>
#include<cmath>
#include<vector>
#include"range_tree.h"
#include<chrono>


int dims = 3;

double get_double(double fMin, double fMax)
{
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

struct rangeHD {
    double *start;
    double *end;
    bool inclusiveXMin;
    bool inclusiveXMax;
    bool inclusiveYMin;
    bool inclusiveYMax;
    int dimns;

    rangeHD() {
        start = nullptr;
        end = nullptr;
    }

    rangeHD(double *s, double *e, bool ixmin, bool ixmax, bool iymin, bool iymax, int dimensions) {
        //start = s;
        //end = e;
        dimns = dimensions;
        start = new double[dims];
        end = new double[dims];
        for(int i=0; i<dims; ++i) {
            start[i] = s[i];
            end[i] = e[i];
        }
        /*
        start[0] = s[0];
        end[0] = e[0];
        start[1] = s[1];
        end[1] = e[1];
        */
        inclusiveXMin = ixmin;
        inclusiveXMax = ixmax;
        inclusiveYMin = iymin;
        inclusiveYMax = iymax; 
    }

    rangeHD& operator=(const rangeHD& other) {
        inclusiveXMin = other.inclusiveXMin;
        inclusiveXMax = other.inclusiveXMax;
        inclusiveYMin = other.inclusiveYMin;
        inclusiveYMax = other.inclusiveYMax;
        if (start != nullptr) {
            delete [] start;
        }
        if (end != nullptr) {
            delete [] end;
        }
        start = new double[dims];
        end = new double[dims];
        for(int i=0; i<dims; ++i) {
            start[i] = other.start[i];
            end[i] = other.end[i];
        }
        /*
        start = new double[2];
        end = new double[2];
        start[0] = other.start[0];
        end[0] = other.end[0];
        start[1] = other.start[1];
        end[1] = other.end[1];
        */
        return *this;
    }

    rangeHD(const range& other) {
        inclusiveXMin = other.inclusiveXMin;
        inclusiveXMax = other.inclusiveXMax;
        inclusiveYMin = other.inclusiveYMin;
        inclusiveYMax = other.inclusiveYMax;
        start = new double[dims];
        end = new double[dims];
        for(int i=0; i<dims; ++i) {
            start[i] = other.start[i];
            end[i] = other.end[i];
        }
        /*
        start = new double[2];
        end = new double[2];
        start[0] = other.start[0];
        end[0] = other.end[0];
        start[1] = other.start[1];
        end[1] = other.end[1];
        */
    }

    ~rangeHD() {
        if (start != nullptr) {
            delete [] start;
        }
        start = nullptr;
        if (end != nullptr) {
            delete [] end;
        }
        end = nullptr;    
    }
};


vector<Point> *findPointsInRangeLocal(RangeTree &rt, rangeHD& initialRange) {
    // Filter the first 2 dimensions
    range initialRange2D = range(initialRange.start, initialRange.end, initialRange.inclusiveXMin, initialRange.inclusiveXMax, initialRange.inclusiveYMin, initialRange.inclusiveYMax);
    vector<Point> *pointsInRange = rt.findPointsInRange(initialRange2D);
    if(dims <= 2) {
        return pointsInRange;
    }
    vector<Point> *output = new vector<Point>();
    if(!pointsInRange || pointsInRange->size() == 0) {
        free(pointsInRange);
        return output;
    }
    for(int i=0; i < pointsInRange->size(); ++i) {
        int j = 2;
        for(; j < dims; ++j) {
            double value = pointsInRange->at(i).coordinates[j];
            if(value > initialRange.end[j] || value < initialRange.start[j]) {
                break;
            }
        }
        if(j == dims) {
            output->push_back(pointsInRange->at(i));
        }
    }
    free(pointsInRange);
    return output;
}

using namespace std;

int main1(int argc, char *argv[]) {
    unordered_map<int, vector<pair<pair<double, double>, pair<double, double>>>> map_result;
    double startRangeDummy[2], endRangeDummy[2];
    //RangeTree rt("data/10k_new.csv", startRangeDummy, endRangeDummy);
    //double startX=-0.507015, endX=0.297345, startY=51.306584, endY=51.660974;
    /*
    startRangeDummy[0] = 0;
    startRangeDummy[1] = 1;
    endRangeDummy[0] = 0;
    endRangeDummy[1] = 1;
    */
    int start=2000, end=4000;
    RangeTree rt("data/uniform.csv", startRangeDummy, endRangeDummy);
    double startX=0., endX=1000., startY=0., endY=1000.;
    for(int i=start; i<= end; i += 200) {
        vector<pair<pair<double, double>, pair<double, double>>> empty;
        // Insert an empty vector for each value of i
        map_result.insert(make_pair(i, empty));
    }
    for(int i=0; i < 1000; ++i) {
        double t0 = get_double(startX, endX), t1 = get_double(startX, endX);
        double x0, x1, y0, y1;
        x0 = min(t0, t1);
        x1 = max(t0, t1);
        t0 = get_double(startY,  endY);
        t1 = get_double(startY, endY);
        y0 = min(t0, t1);
        y1 = max(t0, t1);
        double skylineStart[2], skylineEnd[2];
        skylineStart[0] = x0;
        skylineStart[1] = y0;
        skylineEnd[0] = x1;
        skylineEnd[1] = y1;
        range skylineRange(skylineStart, skylineEnd, true, true, true, true);
        vector<Point> *pointsInSkyline = rt.findPointsInRange(skylineRange);
        int length = pointsInSkyline->size();
        delete pointsInSkyline;
        int bucket = 200*ceil(length/200);
        if(bucket > end || bucket < start || map_result.find(bucket)->second.size() >= 10) {
            continue;
        }
        pair<double, double> startRange = make_pair(x0, y0), endRange = make_pair(x1, y1);
        map_result.find(bucket)->second.push_back(make_pair(startRange, endRange));
    }
    for(int i=start; i<= end; i += 200) {
        auto range_vectors = map_result.find(i)->second;
        cout << "Exploring bucket " << i << endl << endl;
        for(int j=0; j < range_vectors.size(); ++j) {
            cout << "[" << range_vectors.at(j).first.first << "," << range_vectors.at(j).first.second << "]";
            cout << "[" << range_vectors.at(j).second.first << "," << range_vectors.at(j).second.second << "]";
            cout << endl;
        }
    }
}


int main(int argc, char *argv[]) {
    unordered_map<int, vector<pair<vector<double>, vector<double>>>> map_result;
    double startRangeDummy[2], endRangeDummy[2];
    int start=200, end=1000;
    RangeTree rt("data/uniform_3d.csv", startRangeDummy, endRangeDummy, 3);
    double startVal=0., endVal=1000.;
    for(int i=start; i<= end; i += 200) {
        vector<pair<vector<double>, vector<double>>> empty;
        // Insert an empty vector for each value of i
        map_result.insert(make_pair(i, empty));
    }
    for(int i=0; i < 1000; ++i) {
        double skylineStart[dims], skylineEnd[dims];
        for(int j=0; j<dims; ++j) {
            double v1 = get_double(startVal, endVal), v2 = get_double(startVal, endVal);
            skylineStart[j] = min(v1, v2);
            skylineEnd[j] = max(v1, v2);
        }
        rangeHD skylineRange(skylineStart, skylineEnd, true, true, true, true, dims);
        vector<Point> *pointsInSkyline = findPointsInRangeLocal(rt, skylineRange);
        int length = pointsInSkyline->size();
        delete pointsInSkyline;
        int bucket = 200*ceil(length/200);
        if(bucket > end || bucket < start || map_result.find(bucket)->second.size() >= 10) {
            continue;
        }
        vector<double> startRange, endRange;
        for(int j=0; j<dims; ++j) {
            startRange.push_back(skylineStart[j]);
            endRange.push_back(skylineEnd[j]);
        }
        map_result.find(bucket)->second.push_back(make_pair(startRange, endRange));
    }
    for(int i=start; i<= end; i += 200) {
        auto range_vectors = map_result.find(i)->second;
        cout << "Exploring bucket " << i << endl << endl;
        for(int j=0; j < range_vectors.size(); ++j) {
            char c = '[';
            for(double startVal: range_vectors.at(j).first) {
                cout << c << startVal;
                c = ' ';
            }
            cout << "], ";
            c = '[';
            for(double endVal: range_vectors.at(j).second) {
                cout << c << endVal;
                c = ' ';
            }
            cout << "]";
            cout << endl;
        }
    }
    return 0;
}
