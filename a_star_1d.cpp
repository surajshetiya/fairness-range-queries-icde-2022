#include<iostream>
#include<set>
#include<boost/heap/pairing_heap.hpp>
#include<cmath>
#include<vector>
#include<string>
#include<chrono>
#include<fstream>
#include<sstream>

int blueWeight = 2, redWeight = 1, diff = 10;
int blueTotal = 0, redTotal = 0;
double percentageDiff = 0.1;
using namespace std;

struct VectorHash {
    size_t operator()(const std::vector<double>& v) const {
        std::hash<double> hasher;
        size_t seed = 0;
        for (double i : v) {
            seed ^= hasher(i) + 0x9e3779b9 + (seed<<6) + (seed>>2);
        }
        return seed;
    }
};

struct Point {
    double x;
    bool isBlue;
    bool inRange;

    bool operator<(Point const & point2) const {
        return x < point2.x;
    }
};

vector<Point> *readPoints(string filename, double minimum, double maximum) {
    vector<Point> *points = new vector<Point>();

    string line;

    ifstream file(filename);

    if (file.bad()) {
        cout << "BAD FILE!";
        exit;
    }

    getline(file, line, '\n');

    while (!file.eof()) {
        stringstream ss(line);
        Point point;
        double x;
        ss >> x;
        point.x = x;
        ss >> x;
        int isBlue;
        //1 is blue?
        ss >> isBlue;
        point.isBlue = isBlue == 1;
        point.inRange = point.x >= minimum && point.x <= maximum;
        if(point.inRange) {
            if(point.isBlue) {
                blueTotal++;
            } else {
                redTotal++;
            }
        }
        points->push_back(point);
        getline(file, line, '\n');
    }
    // Sort in ascending order for quick movement
    sort(points->begin(), points->end());
    return points;
}

struct sps {
    double start, end;
    int startIndex, endIndex;
    int blues, reds, intersection, union_, diffVal;

    sps(double s, double e, int sInd, int eInd, int blueCount, int redCount, int itemsIntersection, int itemsUnion) {
        start = s;
        end = e;
        startIndex = sInd;
        endIndex = eInd;
        blues= blueCount;
        reds = redCount;
        intersection = itemsIntersection;
        union_ = itemsUnion;
        diffVal = abs(blueWeight * blues - redWeight * reds);
    }

    double similarity() {
        return (double)intersection/((double)union_);
    }

    double distanceApprox() const{
        if(diffVal > diff) {
            // Dominated by blues
            double diff_ = diffVal;
            double J1 = (double(intersection) - ceil((double(diff_ - diff))/blueWeight))/union_;
            double J2 = (double(intersection))/double(union_ + ceil((double(diff_ - diff))/redWeight));
            return J1 > J2 ? (1-J1) : (1-J2); // J1 and J2 represent similarity, larger the better
        } else if(diffVal < -diff) {
            // Dominated by reds
            double diff_ = -diffVal; // Get the absolute value of diff
            double J1 = (double(intersection) - ceil((double(diff_ - diff))/redWeight))/union_;
            double J2 = (double(intersection))/double(union_ + ceil((double(diff_ - diff))/blueWeight));
            return J1 > J2 ? (1-J1) : (1-J2); // J1 and J2 represent similarity, larger the better
        } else {
            // Already fair - return distance of 0
            return 0;
        }
    }

    int numberPointsAway() const {
        if(diffVal > diff) {
            // Dominated by blues
            double diff_ = diffVal;
            double J1 = (double(intersection) - ceil((double(diff_ - diff))/blueWeight))/union_;
            double J2 = (double(intersection))/double(union_ + ceil((double(diff_ - diff))/redWeight));
            return J1 > J2 ? ceil((double(diff_ - diff))/blueWeight) : ceil((double(diff_ - diff))/redWeight);
        } else if(diffVal < -diff) {
            // Dominated by reds
            double diff_ = -diffVal; // Get the absolute value of diff
            double J1 = (double(intersection) - ceil((double(diff_ - diff))/redWeight))/union_;
            double J2 = (double(intersection))/double(union_ + ceil((double(diff_ - diff))/blueWeight));
            return J1 > J2 ? ceil((double(diff_ - diff))/redWeight) : ceil((double(diff_ - diff))/blueWeight);
        } else {
            // Already fair - return distance of 0
            return 0;
        }
    }

    bool operator<(sps const & sps2) const {
        if(distanceApprox() == sps2.distanceApprox()) {
            return numberPointsAway() > sps2.numberPointsAway();
        }
        return distanceApprox() > sps2.distanceApprox();
    }
};

int main(int argc, char *argv[]) {
    double inputRangeStart, inputRangeEnd;
    //unordered_set<vector<double>> uniqueRanges;
    inputRangeStart = -0.1;
    inputRangeEnd = 0.2;


    double threshold = 0.2; // Distance threshold

    //RangeTree rt(argv[1], inputRange[0], inputRange[1]);
    vector<Point>* points = readPoints("data/10k_new.csv", inputRangeStart, inputRangeEnd);
    cout << "Initital range is " << inputRangeStart << ":" << inputRangeEnd << endl;
    auto startTime = chrono::high_resolution_clock::now();
    vector<Point>::iterator low,up;
    Point p;
    p.x = inputRangeStart;
    low=std::lower_bound (points->begin(), points->end(), p);
    p.x = inputRangeEnd;
    up = std::upper_bound (points->begin(), points->end(), p);
    int startIndex = low - points->begin(), endIndex = (up - points->begin()) - 1;
    boost::heap::pairing_heap<sps> H;
    sps first(points->at(startIndex).x, points->at(endIndex).x, startIndex, endIndex, blueTotal, redTotal, blueTotal + redTotal, blueTotal + redTotal);
    int dummyCounter=0;
    H.push(first);
    cout << "Points in input range " << first.union_ << endl;
    cout << "Iniital disparity " << first.diffVal << endl;
    set<pair<int, int>> uniqueVals;
    pair<int, int> firstVal = make_pair(startIndex, endIndex);
    uniqueVals.insert(firstVal);
    diff = ceil( first.union_ * 6514. * percentageDiff/ 10000.);
    cout << "Allowed disparity is " << diff << endl;
    int minDisp = first.diffVal;
    while(!H.empty()) {
        sps topItem = H.top();
        H.pop();
        if(topItem.diffVal < minDisp) {
            minDisp = topItem.diffVal;
            cout << "Found better disparity of " << minDisp << endl;
        }
        if(++dummyCounter % 100000 == 0) {
            cout << "Disparity : " << topItem.diffVal << " " << topItem.intersection << ":" << topItem.union_ << " ; Distance approx : " << topItem.distanceApprox() << endl;
        }
        if(topItem.diffVal < diff) {
            // Fair range
            cout << "Found fair range with similarity " << topItem.similarity() << endl;
            auto endTime = chrono::high_resolution_clock::now();
            auto duration = chrono::duration_cast<chrono::microseconds>(endTime - startTime);
            cout << "Total time is : " << duration.count()/1000000. << " seconds" << endl;
            cout << "Fair range is " << topItem.start << ":" << topItem.end << std::endl;
            cout << "Fair range indices are " << topItem.startIndex << ":" << topItem.endIndex << std::endl;
            return 0;
        }
        if(points->at(topItem.startIndex).inRange) {
            pair<int, int> setCheck = make_pair(topItem.startIndex + 1, topItem.endIndex);
            if(uniqueVals.count(setCheck) <= 0) {
                uniqueVals.insert(setCheck);
                int blueCountVal = topItem.blues, redCountVal = topItem.reds;
                if(points->at(topItem.startIndex).isBlue) {
                    blueCountVal--;
                } else {
                    redCountVal--;
                }
                sps leftShrink(points->at(topItem.startIndex + 1).x, topItem.end, topItem.startIndex + 1, topItem.endIndex, blueCountVal, redCountVal, topItem.intersection-1, topItem.union_);
                H.push(leftShrink);
            }
        }
        if(points->at(topItem.endIndex).inRange) {
            pair<int, int> setCheck = make_pair(topItem.startIndex, topItem.endIndex-1);
            if(uniqueVals.count(setCheck) <= 0) {
                uniqueVals.insert(setCheck);
                int blueCountVal = topItem.blues, redCountVal = topItem.reds;
                if(points->at(topItem.endIndex).isBlue) {
                    blueCountVal--;
                } else {
                    redCountVal--;
                }
                sps rightShrink(topItem.start, points->at(topItem.endIndex-1).x, topItem.startIndex, topItem.endIndex-1, blueCountVal, redCountVal, topItem.intersection-1, topItem.union_);
                H.push(rightShrink);
            }
        }
        if(topItem.startIndex > 0 && (!points->at(topItem.startIndex-1).inRange)) {
            pair<int, int> setCheck = make_pair(topItem.startIndex - 1, topItem.endIndex);
            if(uniqueVals.count(setCheck) <= 0) {
                uniqueVals.insert(setCheck);
                int blueCountVal = topItem.blues, redCountVal = topItem.reds;
                if(points->at(topItem.startIndex-1).isBlue) {
                    blueCountVal++;
                } else {
                    redCountVal++;
                }
                sps leftExpand(points->at(topItem.startIndex - 1).x, topItem.end, topItem.startIndex - 1, topItem.endIndex, blueCountVal, redCountVal, topItem.intersection, topItem.union_+1);
                H.push(leftExpand);
            }
        }
        if(topItem.endIndex < (points->size()-1) && (!points->at(topItem.endIndex+1).inRange)) {
            pair<int, int> setCheck = make_pair(topItem.startIndex, topItem.endIndex+1);
            if(uniqueVals.count(setCheck) <= 0) {
                uniqueVals.insert(setCheck);
                int blueCountVal = topItem.blues, redCountVal = topItem.reds;
                if(points->at(topItem.endIndex + 1).isBlue) {
                    blueCountVal++;
                } else {
                    redCountVal++;
                }
                sps rightExpand(topItem.start, points->at(topItem.endIndex + 1).x, topItem.startIndex, topItem.endIndex + 1, blueCountVal, redCountVal, topItem.intersection, topItem.union_+1);
                H.push(rightExpand);
            }
        }

    }
    auto endTime = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::microseconds>(endTime - startTime);
    cout << "Total time is : " << duration.count()/1000000. << " seconds" << endl;
    delete points;
    return 0;
}
