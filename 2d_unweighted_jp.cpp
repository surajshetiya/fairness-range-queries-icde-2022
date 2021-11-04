#include <vector>
#include <map>
#include <fstream>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <stdlib.h>
#include <math.h>
#include <chrono>

using namespace std;

struct Point {
    bool isBlue;
    bool inRange;
    double *coordinates;

    Point() {
        coordinates = nullptr;
    }   

    Point(const Point &other) {
        isBlue = other.isBlue;
        inRange = other.inRange;
        coordinates = new double[2];
        coordinates[0] = other.coordinates[0];
        coordinates[1] = other.coordinates[1];
    }

    Point& operator=(const Point &other) {
        isBlue = other.isBlue;
        inRange = other.inRange;
        if (coordinates != nullptr) {
            delete [] coordinates;
        }
        coordinates = new double[2];
        coordinates[0] = other.coordinates[0];
        coordinates[1] = other.coordinates[1];
        return *this;
    }

    ~Point() {
        if (coordinates != nullptr) {
            delete [] coordinates;
        }
        coordinates = nullptr;
    }
};

struct RangePoint {
    Point point;
    int value;
    int pointsAddedOrRemoved;
    bool added;

    bool operator==(const RangePoint& other) {
        return this->point.coordinates[0] == other.point.coordinates[0] && this->point.coordinates[1] == other.point.coordinates[1];
    }
};

bool operator<(const RangePoint& p1, const RangePoint& p2) {
    return true;
}

struct JumpPointer {
    int pointsAddedOrRemoved;
    bool added;
    int value;
    int disparity;
    Point point;
    vector<Point> *intermediates;

    JumpPointer() {
        intermediates = new vector<Point>();
    }
};

struct JumpPointerTracker {
    vector<JumpPointer> *jumpPointers;
    int currentValue;
    int currentMax;
    int currentPointsAddedOrRemoved;
    vector<Point> *topIntermediates;
    vector<Point> *bottomIntermediates;
    int minusTopValue;
    int minusBottomValue;
    int centerCount;

    JumpPointerTracker() {
        currentValue = 0;
        currentMax = 0;
        currentPointsAddedOrRemoved = 0;
        minusTopValue = 0;
        minusBottomValue = 0;
        centerCount = 0;
        jumpPointers = new vector<JumpPointer>();
        topIntermediates = new vector<Point>();
        bottomIntermediates = new vector<Point>();
    }
};

struct Level {
    vector<JumpPointerTracker> *trackers;
    int disparity;
    JumpPointer top;
    JumpPointer bottom;
    double previousTop;
    double previousBottom;

    Level(int d) {
        disparity = d; 
        trackers = new vector<JumpPointerTracker>(4);
    }
};

class pointSorter {
    int d;
    bool asc;
    public:
        pointSorter(int d, bool asc) {
            this->d = d;
            this->asc = asc;
        }
        bool operator()(Point const o1, Point const o2) const {
            return asc ? o1.coordinates[d] < o2.coordinates[d] : o1.coordinates[d] > o2.coordinates[d];
        }

};

class pointSearcher {
    int d;
    public:
        pointSearcher(int d) {
            this->d = d;
        }
        bool operator()(Point const p, double const val) {
            return p.coordinates[d] < val;
        }
        bool operator()(double const val, Point const p) {
            return val < p.coordinates[d];
        }
};

vector<Point> *readPoints(string filename, double *minimum, double *maximum, int &redInRange, int &blueInRange) {
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
        point.coordinates = new double[2];
        ss >> point.coordinates[0];
        ss >> point.coordinates[1];
        int isBlue;
        //1 is blue?
        ss >> isBlue;
        point.isBlue = isBlue == 1;
        point.inRange = point.coordinates[0] >= minimum[0] && point.coordinates[1] >= minimum[1] &&
            point.coordinates[0] <= maximum[0] && point.coordinates[1] <= maximum[1];
        (*points).push_back(point);
        if (point.inRange) {
            point.isBlue ? blueInRange++ : redInRange++; 
        }
        getline(file, line, '\n');
    }
    return points;  
}

vector<JumpPointer> *createJumpPointersFromPoints(vector<Point> *points, Point initial, bool needsBlues, int threshold, bool added, bool stretching) {
    int currentVal = 0;
    int currentMax = 0;
    int totalPoints = 0;

    vector<JumpPointer> *jumpPointers = new vector<JumpPointer>();
    vector<Point> *intermediates = new vector<Point>();

    JumpPointer initialJp;
    initialJp.point = initial;
    initialJp.pointsAddedOrRemoved = 0;
    initialJp.added = added;
    initialJp.value = 0;
    jumpPointers->push_back(initialJp);

    for (Point p : *points) {
        totalPoints++;
        if (stretching ? p.isBlue == needsBlues : p.isBlue != needsBlues) {
            currentVal += 1;
        } else {
            currentVal -= 1;
        }
        if (currentVal > currentMax) {
            currentMax = currentVal;
            JumpPointer jp;
            jp.point = p;
            jp.pointsAddedOrRemoved = totalPoints;
            jp.added = added;
            jp.value = currentVal;
            copy(intermediates->begin(), intermediates->end(), back_inserter(*jp.intermediates));
            jumpPointers->push_back(jp);
            if (currentMax == threshold) {
                return jumpPointers;
            }
            intermediates->clear();
        } else {
            intermediates->push_back(p);
        }
    }
    return jumpPointers;
}

vector<JumpPointer> **generateInitialJumpPoints(vector<Point> *points, double *minimum, double *maximum, bool needsBlues, int threshold) {
    sort(points->begin(), points->end(), pointSorter(0, true));
    
    auto lb = lower_bound(points->begin(), points->end(), minimum[0], pointSearcher(0)) + 1;
    auto ub = upper_bound(points->begin(), points->end(), maximum[0], pointSearcher(0)); 

    vector<Point> *pointsInD1Range = new vector<Point>(); 
    copy(lb, ub, back_inserter(*pointsInD1Range));

    sort(pointsInD1Range->begin(), pointsInD1Range->end(), pointSorter(1, true));

    lb = pointsInD1Range->begin();
    ub = upper_bound(pointsInD1Range->begin(), pointsInD1Range->end(), minimum[1], pointSearcher(1)); 
    vector<Point> *pointsStretchDown = new vector<Point>();
    copy(lb, ub, back_inserter(*pointsStretchDown));

    lb = lower_bound(pointsInD1Range->begin(), pointsInD1Range->end(), maximum[1], pointSearcher(1)) + 1;
    ub = pointsInD1Range->end(); 
    vector<Point> *pointsStretchUp = new vector<Point>();
    copy(lb, ub, back_inserter(*pointsStretchUp));

    lb = lower_bound(pointsInD1Range->begin(), pointsInD1Range->end(), minimum[1], pointSearcher(1)) + 1;
    ub = upper_bound(pointsInD1Range->begin(), pointsInD1Range->end(), maximum[1], pointSearcher(1));

    vector<Point> *pointsShrinkDown = new vector<Point>();
    copy(lb, ub, back_inserter(*pointsShrinkDown));

    lb = pointsShrinkDown->begin();
    ub = pointsShrinkDown->end();
    vector<Point> *pointsShrinkUp = new vector<Point>();
    copy(lb, ub, back_inserter(*pointsShrinkUp));


    reverse(pointsStretchDown->begin(), pointsStretchDown->end());
    reverse(pointsShrinkDown->begin(), pointsShrinkDown->end());

    vector<Point> *pointsInInitialRange = new vector<Point>();
    copy_if(pointsInD1Range->begin(), pointsInD1Range->end(), back_inserter(*pointsInInitialRange), [minimum, maximum](Point pnt) {
        return pnt.coordinates[1] >= minimum[1] && pnt.coordinates[1] <= maximum[1];
    });

    auto maximumPoint = pointsInInitialRange->at(pointsInInitialRange->size()-1);
    auto minimumPoint = pointsInInitialRange->at(0);

    vector<JumpPointer> **jps = new vector<JumpPointer> *[4];
    jps[0] = createJumpPointersFromPoints(pointsStretchUp, maximumPoint, needsBlues, threshold, true, true);
    jps[1] = createJumpPointersFromPoints(pointsShrinkUp, minimumPoint, needsBlues, threshold, false, false);
    jps[2] = createJumpPointersFromPoints(pointsStretchDown, minimumPoint, needsBlues, threshold, true, true);
    jps[3] = createJumpPointersFromPoints(pointsShrinkDown, maximumPoint, needsBlues, threshold, false, false);

    return jps;
}

void addOrShrinkPoints(vector<Point>::iterator initial, Level level, bool needsBlues, int tracker, vector<Point>::iterator lowerBound, vector<Point>::iterator upperBound, bool adding, int epsilon) {
    JumpPointer initialJp;
    initialJp.point = *initial;
    initialJp.pointsAddedOrRemoved = 0;
    initialJp.value = 0;
    initialJp.added = adding;
    initialJp.disparity = level.disparity;
    (*level.trackers)[tracker].jumpPointers->push_back(initialJp);
    for (auto iter = lowerBound; iter < upperBound; iter++) {
        JumpPointer rp1 = level.top;
        JumpPointer rp2 = level.bottom;
        if (iter->coordinates[1] < rp1.point.coordinates[1] && iter->coordinates[1] > rp2.point.coordinates[1]) {
            if (level.disparity - (*level.trackers)[tracker].jumpPointers->size() + 1 <= epsilon) {
                continue;
            }
            (*level.trackers)[tracker].currentPointsAddedOrRemoved += 1;
            if (adding ? iter->isBlue == needsBlues : iter->isBlue != needsBlues) {
                (*level.trackers)[tracker].currentValue += 1;
            } else {
                (*level.trackers)[tracker].currentValue -= 1;
            }
            if ((*level.trackers)[tracker].currentValue > (*level.trackers)[tracker].currentMax) {
                (*level.trackers)[tracker].currentMax = (*level.trackers)[tracker].currentValue;
                JumpPointer jp;
                jp.point = *iter;
                jp.pointsAddedOrRemoved = (*level.trackers)[tracker].currentPointsAddedOrRemoved;
                jp.value = (*level.trackers)[tracker].currentMax;
                jp.added = adding;
                jp.disparity = level.disparity - (*level.trackers)[tracker].currentMax;
                (*level.trackers)[tracker].jumpPointers->push_back(jp); 
            }
            if (iter->coordinates[1] > level.previousTop) {
                (*level.trackers)[tracker].topIntermediates->push_back(*iter);
            } else {
                (*level.trackers)[tracker].minusTopValue += (adding ? iter->isBlue == needsBlues : iter->isBlue != needsBlues) ? 1 : -1;
            }
            if (iter->coordinates[1] < level.previousBottom) {
                (*level.trackers)[tracker].bottomIntermediates->push_back(*iter);
            } else {
                (*level.trackers)[tracker].minusBottomValue += (adding ? iter->isBlue == needsBlues : iter->isBlue != needsBlues) ? 1 : -1;
            }
            if (iter->coordinates[1] < level.previousTop && iter->coordinates[1] > level.previousBottom) {
                (*level.trackers)[tracker].centerCount += 1;
            }
        }
    }
}


void findMaxInRange(int& count, double &position, int currentValue, int desiredDesparity, vector<Point> *centerIntermediates, vector<Point> *extremeIntermediates, bool ascending, bool needsBlues) {
    bool setToAny = false;
    int i = 0;
    int j = 0;
    while (i < centerIntermediates->size() && j < extremeIntermediates->size()) {
        count++;
        if (ascending ? centerIntermediates->at(i).coordinates[1] < extremeIntermediates->at(j).coordinates[1] : centerIntermediates->at(i).coordinates[1] > extremeIntermediates->at(j).coordinates[1]) {
            currentValue += centerIntermediates->at(i).isBlue == needsBlues ? 1 : -1;
            if (setToAny) {
                position = centerIntermediates->at(i).coordinates[1];
            }
            i++; 
        } else {
            currentValue += extremeIntermediates->at(j).isBlue == needsBlues ? 1 : -1;
            position = extremeIntermediates->at(j).coordinates[1];
            setToAny = true;
            j++;
        }
        if (currentValue == desiredDesparity) {
            return;
        }
    }
    while (j < extremeIntermediates->size()) {
        currentValue += extremeIntermediates->at(j).isBlue == needsBlues ? 1 : -1;
        position = extremeIntermediates->at(j).coordinates[1];
        j++;
        if (currentValue == desiredDesparity) {
            return;
        }
    }
}

void findMostSimilarFairRange(Level level, vector<Point> *pointsInRange, double &mostSimilarFairRange, double* minimum, double* maximum, int initialPoints, int epsilon, bool needsBlues) {
    int threshold = level.disparity - epsilon;
    int offset0 = threshold < (*level.trackers)[0].jumpPointers->size() ? threshold : (*level.trackers)[0].jumpPointers->size() - 1;
    if (offset0 >= 0 && (*level.trackers)[0].jumpPointers->at(offset0).disparity - epsilon < (*level.trackers)[2].jumpPointers->size()) {
        int offset1 = (*level.trackers)[0].jumpPointers->at(offset0).disparity - epsilon;
        while (offset0 >= 0 && offset1 < (*level.trackers)[2].jumpPointers->size()) {
            vector<Point> *pointsStrictlyInRange = new vector<Point>();
            copy_if(pointsInRange->begin(), pointsInRange->end(), back_inserter(*pointsStrictlyInRange), [level, offset0, offset1](Point p) { return p.coordinates[0] < (*level.trackers)[0].jumpPointers->at(offset0).point.coordinates[0] && p.coordinates[0] > (*level.trackers)[2].jumpPointers->at(offset1).point.coordinates[0];});
            vector<Point>::iterator minElement = min_element(pointsStrictlyInRange->begin(), pointsStrictlyInRange->end(), pointSorter(1, true));
            vector<Point>::iterator maxElement = max_element(pointsStrictlyInRange->begin(), pointsStrictlyInRange->end(), pointSorter(1, true));
            double minY = (*minElement).coordinates[1];
            double maxY = (*maxElement).coordinates[1];
            double maxX = (*level.trackers)[0].jumpPointers->at(offset0).point.coordinates[0];
            double minX = (*level.trackers)[2].jumpPointers->at(offset1).point.coordinates[0];
            int pointsInTopRange = 0;
            int pointsInBottomRange = 0;
            vector<Point> *topIntermediates = new vector<Point>();
            vector<Point> *bottomIntermediates = new vector<Point>();
            copy((*level.trackers)[0].topIntermediates->begin(), upper_bound((*level.trackers)[0].topIntermediates->begin(), (*level.trackers)[0].topIntermediates->end(), maxX, pointSearcher(0)), back_inserter(*topIntermediates));
            copy((*level.trackers)[2].topIntermediates->begin(), upper_bound((*level.trackers)[2].topIntermediates->begin(), (*level.trackers)[2].topIntermediates->end(), minX, pointSearcher(0)), back_inserter(*topIntermediates));
            sort(topIntermediates->begin(), topIntermediates->end(), pointSorter(1, true));
            copy((*level.trackers)[0].bottomIntermediates->begin(), upper_bound((*level.trackers)[0].bottomIntermediates->begin(), (*level.trackers)[0].bottomIntermediates->end(), maxX, pointSearcher(0)), back_inserter(*bottomIntermediates));
            copy((*level.trackers)[2].bottomIntermediates->begin(), upper_bound((*level.trackers)[2].bottomIntermediates->begin(), (*level.trackers)[2].bottomIntermediates->end(), minX, pointSearcher(0)), back_inserter(*bottomIntermediates));
            sort(bottomIntermediates->begin(), bottomIntermediates->end(), pointSorter(1, false));
            findMaxInRange(pointsInTopRange, maxY, (*level.trackers)[0].minusTopValue + (*level.trackers)[2].minusTopValue, (*level.trackers)[2].jumpPointers->at((*level.trackers)[2].jumpPointers->size() - 1).value + (*level.trackers)[0].jumpPointers->at((*level.trackers)[0].jumpPointers->size() - 1).value, (*level.trackers)[0].jumpPointers->at(offset0).intermediates, topIntermediates, true, needsBlues);
            findMaxInRange(pointsInBottomRange, minY, (*level.trackers)[0].minusBottomValue + (*level.trackers)[2].minusBottomValue, (*level.trackers)[2].jumpPointers->at((*level.trackers)[2].jumpPointers->size() - 1).value + (*level.trackers)[0].jumpPointers->at((*level.trackers)[0].jumpPointers->size() - 1).value, (*level.trackers)[0].jumpPointers->at(offset0).intermediates, bottomIntermediates, false, needsBlues);
            int newPoints = (*level.trackers)[0].centerCount + pointsInBottomRange + pointsInTopRange + (level.top.added ? level.top.pointsAddedOrRemoved : 0) + (level.bottom.added ? level.bottom.pointsAddedOrRemoved : 0);
            int removedPoints = (!level.top.added ? level.top.pointsAddedOrRemoved : 0) + (!level.bottom.added ? level.bottom.pointsAddedOrRemoved : 0);
            double similarity = (double) (initialPoints - removedPoints) / (double) (newPoints + initialPoints);
            if (similarity > mostSimilarFairRange) {
                
                mostSimilarFairRange = similarity;
                minimum[0] = minX;
                minimum[1] = minY;
                maximum[0] = maxX;
                maximum[1] = maxY;
            }
            offset0--;
            offset1++;
        }
    }

    offset0 = threshold < (*level.trackers)[0].jumpPointers->size() ? threshold : (*level.trackers)[0].jumpPointers->size() -1;
    if (offset0 >= 0 && (*level.trackers)[0].jumpPointers->at(offset0).disparity- epsilon < (*level.trackers)[1].jumpPointers->size()) {
        int offset1 = (*level.trackers)[0].jumpPointers->at(offset0).disparity- epsilon;
        while (offset0 >= 0 && offset1 < (*level.trackers)[1].jumpPointers->size()) {
            vector<Point> *pointsStrictlyInRange = new vector<Point>();
            copy_if(pointsInRange->begin(), pointsInRange->end(), back_inserter(*pointsStrictlyInRange), [level, offset0, offset1](Point p) { return p.coordinates[0] < (*level.trackers)[0].jumpPointers->at(offset0).point.coordinates[0] && p.coordinates[0] > (*level.trackers)[1].jumpPointers->at(offset1).point.coordinates[0];});
            vector<Point>::iterator minElement = min_element(pointsStrictlyInRange->begin(), pointsStrictlyInRange->end(), pointSorter(1, true));
            vector<Point>::iterator maxElement = max_element(pointsStrictlyInRange->begin(), pointsStrictlyInRange->end(), pointSorter(1, true));
            double minY = (*minElement).coordinates[1];
            double maxY = (*maxElement).coordinates[1];
            double maxX = (*level.trackers)[0].jumpPointers->at(offset0).point.coordinates[0];
            auto lb = lower_bound(pointsStrictlyInRange->begin(), pointsStrictlyInRange->end(), (*level.trackers)[1].jumpPointers->at(offset1).point.coordinates[0], pointSearcher(0)) + 1;
            double minX = (lb == pointsStrictlyInRange->end() ? pointsStrictlyInRange->end()-1 : lb)->coordinates[0];
            int pointsInTopRange = 0;
            int pointsInBottomRange = 0;
            vector<Point> *topIntermediates = new vector<Point>();
            vector<Point> *bottomIntermediates = new vector<Point>();
            copy((*level.trackers)[0].topIntermediates->begin(), upper_bound((*level.trackers)[0].topIntermediates->begin(), (*level.trackers)[0].topIntermediates->end(), maxX, pointSearcher(0)), back_inserter(*topIntermediates));
            sort(topIntermediates->begin(), topIntermediates->end(), pointSorter(1, true));
            copy((*level.trackers)[0].bottomIntermediates->begin(), upper_bound((*level.trackers)[0].bottomIntermediates->begin(), (*level.trackers)[0].bottomIntermediates->end(), maxX, pointSearcher(0)), back_inserter(*bottomIntermediates));
            sort(bottomIntermediates->begin(), bottomIntermediates->end(), pointSorter(1, false));
            findMaxInRange(pointsInTopRange, maxY, (*level.trackers)[0].minusTopValue, (*level.trackers)[0].jumpPointers->at((*level.trackers)[0].jumpPointers->size() - 1).value, (*level.trackers)[0].jumpPointers->at(offset0).intermediates, topIntermediates, true, needsBlues);
            findMaxInRange(pointsInBottomRange, minY, (*level.trackers)[0].minusBottomValue, (*level.trackers)[0].jumpPointers->at((*level.trackers)[0].jumpPointers->size() - 1).value, (*level.trackers)[0].jumpPointers->at(offset0).intermediates, bottomIntermediates, false, needsBlues);
            int newPoints = (*level.trackers)[0].centerCount + pointsInBottomRange + pointsInTopRange + (level.top.added ? level.top.pointsAddedOrRemoved : 0) + (level.bottom.added ? level.bottom.pointsAddedOrRemoved : 0);
            int removedPoints = (*level.trackers)[1].currentPointsAddedOrRemoved + (!level.top.added ? level.top.pointsAddedOrRemoved : 0) + (!level.bottom.added ? level.bottom.pointsAddedOrRemoved : 0);
            double similarity = (double) (initialPoints - removedPoints) / (double) (newPoints + initialPoints);
            if (similarity > mostSimilarFairRange) {
                mostSimilarFairRange = similarity;
                minimum[0] = minX;
                minimum[1] = minY;
                maximum[0] = maxX;
                maximum[1] = maxY;
            }
            offset0--;
            offset1++;
        }
    }

    offset0 = threshold < (*level.trackers)[3].jumpPointers->size() ? threshold: (*level.trackers)[3].jumpPointers->size() - 1;
    if (offset0 >= 0 && (*level.trackers)[3].jumpPointers->at(offset0).disparity- epsilon < (*level.trackers)[2].jumpPointers->size()) {
        int offset1 = (*level.trackers)[3].jumpPointers->at(offset0).disparity- epsilon;
        while (offset0 >= 0 && offset1 < (*level.trackers)[2].jumpPointers->size()) {
            vector<Point> *pointsStrictlyInRange = new vector<Point>();
            copy_if(pointsInRange->begin(), pointsInRange->end(), back_inserter(*pointsStrictlyInRange), [level, offset0, offset1](Point p) { return p.coordinates[0] < (*level.trackers)[3].jumpPointers->at(offset0).point.coordinates[0] && p.coordinates[0] > (*level.trackers)[2].jumpPointers->at(offset1).point.coordinates[0];});
            vector<Point>::iterator minElement = min_element(pointsStrictlyInRange->begin(), pointsStrictlyInRange->end(), pointSorter(1, true));
            vector<Point>::iterator maxElement = max_element(pointsStrictlyInRange->begin(), pointsStrictlyInRange->end(), pointSorter(1, true));
            double minY = (*minElement).coordinates[1];
            double maxY = (*maxElement).coordinates[1];
            auto ub = upper_bound(pointsStrictlyInRange->begin(), pointsStrictlyInRange->end(), (*level.trackers)[3].jumpPointers->at(offset0).point.coordinates[0], pointSearcher(0));
            double maxX = (ub == pointsStrictlyInRange->end() ? pointsStrictlyInRange->end()-1 : ub)->coordinates[0];
            double minX = (*level.trackers)[2].jumpPointers->at(offset1).point.coordinates[0];
            int pointsInTopRange = 0;
            int pointsInBottomRange = 0;
            vector<Point> *topIntermediates = new vector<Point>();
            vector<Point> *bottomIntermediates = new vector<Point>();
            copy((*level.trackers)[2].topIntermediates->begin(), upper_bound((*level.trackers)[2].topIntermediates->begin(), (*level.trackers)[2].topIntermediates->end(), maxX, pointSearcher(0)), back_inserter(*topIntermediates));
            sort(topIntermediates->begin(), topIntermediates->end(), pointSorter(1, true));
            copy((*level.trackers)[2].bottomIntermediates->begin(), upper_bound((*level.trackers)[2].bottomIntermediates->begin(), (*level.trackers)[2].bottomIntermediates->end(), maxX, pointSearcher(0)), back_inserter(*bottomIntermediates));
            sort(bottomIntermediates->begin(), bottomIntermediates->end(), pointSorter(1, false));
            findMaxInRange(pointsInTopRange, maxY, (*level.trackers)[2].minusTopValue, (*level.trackers)[2].jumpPointers->at((*level.trackers)[2].jumpPointers->size() - 1).value, (*level.trackers)[2].jumpPointers->at(offset1).intermediates, topIntermediates, true, needsBlues);
            findMaxInRange(pointsInBottomRange, minY, (*level.trackers)[2].minusBottomValue, (*level.trackers)[2].jumpPointers->at((*level.trackers)[2].jumpPointers->size() - 1).value, (*level.trackers)[2].jumpPointers->at(offset1).intermediates, bottomIntermediates, false, needsBlues);
            int newPoints = (*level.trackers)[2].centerCount + pointsInBottomRange + pointsInTopRange + (level.top.added ? level.top.pointsAddedOrRemoved : 0) + (level.bottom.added ? level.bottom.pointsAddedOrRemoved : 0);
            int removedPoints = (*level.trackers)[3].currentPointsAddedOrRemoved + (!level.top.added ? level.top.pointsAddedOrRemoved : 0) + (!level.bottom.added ? level.bottom.pointsAddedOrRemoved : 0);
            double similarity = (double) (initialPoints - removedPoints) / (double) (newPoints + initialPoints);
            if (similarity > mostSimilarFairRange) {
                mostSimilarFairRange = similarity;
                minimum[0] = minX;
                minimum[1] = minY;
                maximum[0] = maxX;
                maximum[1] = maxY;
            }
            offset0--;
            offset1++;
        }
    }

    offset0 = threshold < (*level.trackers)[3].jumpPointers->size() ? threshold : (*level.trackers)[3].jumpPointers->size() - 1;
    if (offset0 >= 0 && (*level.trackers)[3].jumpPointers->at(offset0).disparity - epsilon < (*level.trackers)[1].jumpPointers->size()) {
        int offset1 = (*level.trackers)[3].jumpPointers->at(offset0).disparity - epsilon ;
        while (offset0 >= 0 && offset1 < (*level.trackers)[1].jumpPointers->size()) {
            vector<Point> *pointsStrictlyInRange = new vector<Point>();
            copy_if(pointsInRange->begin(), pointsInRange->end(), back_inserter(*pointsStrictlyInRange), [level, offset0, offset1](Point p) { return p.coordinates[0] < (*level.trackers)[3].jumpPointers->at(offset0).point.coordinates[0] && p.coordinates[0] > (*level.trackers)[1].jumpPointers->at(offset1).point.coordinates[0];});
            vector<Point>::iterator minElement = min_element(pointsStrictlyInRange->begin(), pointsStrictlyInRange->end(), pointSorter(1, true));
            vector<Point>::iterator maxElement = max_element(pointsStrictlyInRange->begin(), pointsStrictlyInRange->end(), pointSorter(1, true));
            double minY = (*minElement).coordinates[1];
            double maxY = (*maxElement).coordinates[1];
            auto ub = upper_bound(pointsStrictlyInRange->begin(), pointsStrictlyInRange->end(), (*level.trackers)[3].jumpPointers->at(offset0).point.coordinates[0], pointSearcher(0));
            double maxX = (ub == pointsStrictlyInRange->end() ? pointsStrictlyInRange->end()-1 : ub)->coordinates[0];
            auto lb = lower_bound(pointsStrictlyInRange->begin(), pointsStrictlyInRange->end(), (*level.trackers)[1].jumpPointers->at(offset1).point.coordinates[0], pointSearcher(0)) + 1;
            double minX = (lb == pointsStrictlyInRange->end() ? pointsStrictlyInRange->end()-1 : lb)->coordinates[0];
            int newPoints = (level.top.added ? level.top.pointsAddedOrRemoved : 0) + (level.bottom.added ? level.bottom.pointsAddedOrRemoved : 0);
            int removedPoints = (*level.trackers)[3].currentPointsAddedOrRemoved + (*level.trackers)[1].currentPointsAddedOrRemoved + (!level.top.added ? level.top.pointsAddedOrRemoved : 0) + (!level.bottom.added ? level.bottom.pointsAddedOrRemoved : 0);
            double similarity = (double) (initialPoints - removedPoints) / (double) (newPoints + initialPoints);
            if (similarity > mostSimilarFairRange) {
                mostSimilarFairRange = similarity;
                minimum[0] = minX;
                minimum[1] = minY;
                maximum[0] = maxX;
                maximum[1] = maxY;
            }
            offset0--;
            offset1++;
        }
    }
}

vector<Level> *generateLevels(vector<Point> *points, vector<JumpPointer> **jumpPointers, int disparity, double* minimum, double* maximum, bool needsBlues, int epsilon, int initialPoints, chrono::time_point<chrono::high_resolution_clock> start) {
    sort(points->begin(), points->end(), pointSorter(0, true));
    vector<Point> *pointsReversed = new vector<Point>();
    auto lowerBoundRight = lower_bound(points->begin(), points->end(), maximum[0], pointSearcher(0)) + 1;
    auto rightInitial = lowerBoundRight;
    auto lowerBoundLeft = lower_bound(points->begin(), points->end(), minimum[0], pointSearcher(0)) + 1;
    vector<Point> *upperBoundLeftReversed = new vector<Point>();
    auto upperBoundLeft = upper_bound(points->begin(), points->end(), minimum[0], pointSearcher(0));
    reverse_copy(points->begin(), upperBoundLeft, back_inserter(*upperBoundLeftReversed));
    auto leftInitial = lowerBoundLeft;
    auto upperBoundRight = upper_bound(points->begin(), points->end(), maximum[0], pointSearcher(0));

    double mostSimilarFairRange = 0;
    double* minimumResult = new double[2];
    double* maximumResult = new double[2];


    for (int i = 0; i < jumpPointers[0]->size(); i++) {
        for (int j = 0; j < jumpPointers[1]->size(); j++) {
            JumpPointer p1 = (*jumpPointers[0])[i];
            JumpPointer p2 = (*jumpPointers[1])[j];
            if (p1.value + p2.value > disparity - epsilon) {
                continue;
            }
            vector<Point> *pointsInRange = new vector<Point>();
            copy_if(lowerBoundLeft, upperBoundRight, back_inserter(*pointsInRange), [jumpPointers, i, j](Point pnt) { 
                return pnt.coordinates[1] <= (*jumpPointers[0])[i].point.coordinates[1] && pnt.coordinates[1] >= (*jumpPointers[1])[j].point.coordinates[1]; 
            });
            auto lowerInRangeLeft = lower_bound(pointsInRange->begin(), pointsInRange->end(), minimum[0], pointSearcher(0)) + 1;
            auto upperInRangeRight = upper_bound(pointsInRange->begin(), pointsInRange->end(), maximum[0], pointSearcher(0));
            vector<Point> *centerReversed = new vector<Point>();
            reverse_copy(lowerInRangeLeft, upperInRangeRight, back_inserter(*centerReversed));


            Level level(disparity - p1.value - p2.value);
            level.disparity = disparity - p1.value - p2.value;
            level.top = p1;
            level.bottom = p2;
            level.previousTop = i > 1 ? (*jumpPointers[0])[i-1].point.coordinates[1] : maximum[1];;
            level.previousBottom = j < jumpPointers[1]->size() - 1 ? (*jumpPointers[1])[j+1].point.coordinates[1] : minimum[1];

            addOrShrinkPoints(rightInitial, level, needsBlues, 0, lowerBoundRight, points->end(), true, epsilon);
            addOrShrinkPoints(leftInitial, level, needsBlues, 1, lowerInRangeLeft, upperInRangeRight, false, epsilon);
            addOrShrinkPoints(leftInitial, level, needsBlues, 2, upperBoundLeftReversed->begin(), upperBoundLeftReversed->end(), true, epsilon);
            addOrShrinkPoints(rightInitial, level, needsBlues, 3, centerReversed->begin(), centerReversed->end(), false, epsilon);
            
            findMostSimilarFairRange(level, pointsInRange, mostSimilarFairRange, minimumResult, maximumResult, initialPoints, epsilon, needsBlues);
        }   
    }
    for (int i = 0; i < jumpPointers[0]->size(); i++) {
        for (int j = 0; j < jumpPointers[2]->size(); j++) {
            JumpPointer p1 = (*jumpPointers[0])[i];
            JumpPointer p2 = (*jumpPointers[2])[j];
            
            if (p1.value + p2.value > disparity - epsilon) {
                continue;
            }

            vector<Point> *pointsInRange = new vector<Point>();
            copy_if(lowerBoundLeft, upperBoundRight, back_inserter(*pointsInRange), [jumpPointers, i, j](Point pnt) { 
                return pnt.coordinates[1] <= (*jumpPointers[0])[i].point.coordinates[1] && pnt.coordinates[1] >= (*jumpPointers[2])[j].point.coordinates[1]; 
            });
            auto lowerInRangeLeft = lower_bound(pointsInRange->begin(), pointsInRange->end(), minimum[0], pointSearcher(0)) + 1;
            auto upperInRangeRight = upper_bound(pointsInRange->begin(), pointsInRange->end(), maximum[0], pointSearcher(0));
            vector<Point> *centerReversed = new vector<Point>();
            reverse_copy(lowerInRangeLeft, upperInRangeRight, back_inserter(*centerReversed));

            Level level(disparity - p1.value - p2.value);
            level.disparity = disparity - p1.value - p2.value;
            level.top = p1;
            level.bottom = p2;
            level.previousTop = i > 1 ? (*jumpPointers[0])[i-1].point.coordinates[1] : maximum[1];;
            level.previousBottom = j > 1 ? (*jumpPointers[2])[j-1].point.coordinates[1] : minimum[1];

            addOrShrinkPoints(rightInitial, level, needsBlues, 0, lowerBoundRight, points->end(), true, epsilon);
            addOrShrinkPoints(leftInitial, level, needsBlues, 1, lowerInRangeLeft, upperInRangeRight, false, epsilon);
            addOrShrinkPoints(leftInitial, level, needsBlues, 2, upperBoundLeftReversed->begin(), upperBoundLeftReversed->end(), true, epsilon);
            addOrShrinkPoints(rightInitial, level, needsBlues, 3, centerReversed->begin(), centerReversed->end(), false, epsilon);
            
            findMostSimilarFairRange(level, pointsInRange, mostSimilarFairRange, minimumResult, maximumResult, initialPoints, epsilon, needsBlues);
        }   
    }
    for (int i = 0; i < jumpPointers[3]->size(); i++) {
        for (int j = 0; j < jumpPointers[1]->size(); j++) {
            JumpPointer p1 = (*jumpPointers[3])[i];
            JumpPointer p2 = (*jumpPointers[1])[j];
            
            if (p1.value + p2.value > disparity - epsilon) {
                continue;
            }

            if (p1.point.coordinates[1] > p2.point.coordinates[1]) {
                vector<Point> *pointsInRange = new vector<Point>();
                copy_if(lowerBoundLeft, upperBoundRight, back_inserter(*pointsInRange), [jumpPointers, i, j](Point pnt) { 
                    return pnt.coordinates[1] <= (*jumpPointers[3])[i].point.coordinates[1] && pnt.coordinates[1] >= (*jumpPointers[1])[j].point.coordinates[1]; 
                });
                auto lowerInRangeLeft = lower_bound(pointsInRange->begin(), pointsInRange->end(), minimum[0], pointSearcher(0)) + 1;
                auto upperInRangeRight = upper_bound(pointsInRange->begin(), pointsInRange->end(), maximum[0], pointSearcher(0));
                vector<Point> *centerReversed = new vector<Point>();
                reverse_copy(lowerInRangeLeft, upperInRangeRight, back_inserter(*centerReversed));

                Level level(disparity - p1.value - p2.value);
                level.disparity = disparity - p1.value - p2.value;
                level.top = p1;
                level.bottom = p2;
                level.previousTop = i < jumpPointers[3]->size() - 1 ? (*jumpPointers[3])[i+1].point.coordinates[1] : maximum[1];
                level.previousBottom = j < jumpPointers[1]->size() - 1 ? (*jumpPointers[1])[j+1].point.coordinates[1] : minimum[1];
           
                addOrShrinkPoints(rightInitial, level, needsBlues, 0, lowerBoundRight, points->end(), true, epsilon);
                addOrShrinkPoints(leftInitial, level, needsBlues, 1, lowerInRangeLeft, upperInRangeRight, false, epsilon);
                addOrShrinkPoints(leftInitial, level, needsBlues, 2, upperBoundLeftReversed->begin(), upperBoundLeftReversed->end(), true, epsilon);
                addOrShrinkPoints(rightInitial, level, needsBlues, 3, centerReversed->begin(), centerReversed->end(), false, epsilon);
                
                findMostSimilarFairRange(level, pointsInRange, mostSimilarFairRange, minimumResult, maximumResult, initialPoints, epsilon, needsBlues);
            }
        }   
    }
    for (int i = 0; i < jumpPointers[3]->size(); i++) {
        for (int j = 0; j < jumpPointers[2]->size(); j++) {
            JumpPointer p1 = (*jumpPointers[3])[i];
            JumpPointer p2 = (*jumpPointers[2])[j];
            
            if (p1.value + p2.value > disparity - epsilon) {
                continue;
            }

            vector<Point> *pointsInRange = new vector<Point>();
            copy_if(lowerBoundLeft, upperBoundRight, back_inserter(*pointsInRange), [jumpPointers, i, j](Point pnt) {
                return pnt.coordinates[1] <= (*jumpPointers[3])[i].point.coordinates[1] && pnt.coordinates[1] >= (*jumpPointers[2])[j].point.coordinates[1];
            });
            auto lowerInRangeLeft = lower_bound(pointsInRange->begin(), pointsInRange->end(), minimum[0], pointSearcher(0)) + 1;
            auto upperInRangeRight = upper_bound(pointsInRange->begin(), pointsInRange->end(), maximum[0], pointSearcher(0));
            vector<Point> *centerReversed = new vector<Point>();
            reverse_copy(lowerInRangeLeft, upperInRangeRight, back_inserter(*centerReversed));

            Level level(disparity - p1.value - p2.value);
            level.disparity = disparity - p1.value - p2.value;
            level.top = p1;
            level.bottom = p2;
            level.previousTop = i < jumpPointers[3]->size() - 1 ? (*jumpPointers[3])[i+1].point.coordinates[1] : maximum[1];
            level.previousBottom = j > 1 ? (*jumpPointers[2])[j-1].point.coordinates[1] : minimum[1];

            addOrShrinkPoints(rightInitial, level, needsBlues, 0, lowerBoundRight, points->end(), true, epsilon);
            addOrShrinkPoints(leftInitial, level, needsBlues, 1, lowerInRangeLeft, upperInRangeRight, false, epsilon);
            addOrShrinkPoints(leftInitial, level, needsBlues, 2, upperBoundLeftReversed->begin(), upperBoundLeftReversed->end(), true, epsilon);
            addOrShrinkPoints(rightInitial, level, needsBlues, 3, centerReversed->begin(), centerReversed->end(), false, epsilon);
            
            findMostSimilarFairRange(level, pointsInRange, mostSimilarFairRange, minimumResult, maximumResult, initialPoints, epsilon, needsBlues);
        }   
    }
    chrono::time_point<chrono::high_resolution_clock> end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(end- start).count();
    cout << "Duration: " << duration << endl;

    cout << "Similarity: " << mostSimilarFairRange << endl;
    cout << "{" << minimumResult[0] << ", " << minimumResult[1] << "}, {" << maximumResult[0] << ", " << maximumResult[1] << "}" << endl;
    vector<Point> *pointsInRange = new vector<Point>();
    copy_if(points->begin(), points->end(), back_inserter(*pointsInRange), 
        [minimumResult, maximumResult](Point pnt) { 
            return pnt.coordinates[0] >= minimumResult[0] && pnt.coordinates[1] >= minimumResult[1] && pnt.coordinates[0] <= maximumResult[0] && pnt.coordinates[1] <= maximumResult[1];
        });
    int finalBlueCount = 0;
    int finalRedCount = 0;
    for (Point p : *pointsInRange) {
        if (p.isBlue) {
            finalBlueCount++;
        } else {
            finalRedCount++;
        }
    }
    cout << "Final blue count: " << finalBlueCount << endl;
    cout << "Final red count: " << finalRedCount << endl; 
    cout << "Epsilon: " << epsilon << endl;
}

int main() {
    // cout << "enter filename: " << endl;
    string filename = "./data/uniform.csv";
    // cin >> filename;
    double *minimum = new double[2];
    double *maximum = new double[2];
    // cout << "enter min[0], min[1], max[0], max[1]: " << endl;
    // cin >> minimum[0];
    // cin >> minimum[1];
    // cin >> maximum[0];
    // cin >> maximum[1];
/*[,][,]
[,][,]
[,][,]
[368.664,232.262][749.771,294.16]
[,], [,]
[656.564,639.458][967.405,759.735]*/

    minimum[0] = 667.341;
    minimum[1] = 313.23;
    maximum[0] = 879.009;
    maximum[1] = 653.305;
    int redInRange = 0;
    int blueInRange = 0;
    chrono::time_point<chrono::high_resolution_clock> start = chrono::high_resolution_clock::now();

    vector<Point> *points = readPoints(filename, minimum, maximum, redInRange, blueInRange);
    int epsilon = ceil(0.025*(redInRange+blueInRange));
    int threshold = abs(redInRange-blueInRange) - epsilon;
    if (threshold <= 0) {
        cout << "Initial range is fair" << endl;
        return 0;
    }
    vector<JumpPointer> **jps = generateInitialJumpPoints(points, minimum, maximum, redInRange > blueInRange, threshold);
    vector<Level> *levels = generateLevels(points, jps, abs(redInRange - blueInRange), minimum, maximum, redInRange > blueInRange, epsilon, redInRange + blueInRange, start);
}
