#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <math.h>
#include <iomanip>

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

struct D1IP {
    vector<Point> *inRange; // [previousImprovingPoint, improvingPoint)
    int improvement;
    Point previousImprovingPoint;
    Point improvingPoint;
    int pointsBelowRange;
};

struct D2IP {
    Point improvingPoint;
    vector<Point> *inMaximumRange;
    vector<Point> *inMinimumRange;
    int improvementBelowMaxRange;
    int improvementAboveMinRange;
    int improvement;
    int improvementInCenter;
    bool shrinking;
    int pointsInStartingRange;
    int pointsInCenter;
    Point next;
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


//This method could be optimized, but made more complex. It uses a filtered array and a (potentially merged) array. It could
//filter the points as it adds them, and keep the arrays unmerged. This would reduce the upper bound for finding the maximum
//from the #points in final improving point the region to the points necessary to reach fairness. 
double maxD1(vector<Point> *inRangeD1, vector<Point> *inRangeHorizontalExpansion, bool ascending, int d1, int improvementSoFar, int desiredImprovement, bool needsBlues, bool &willNotCompute, Point previousImprovingPoint, int &pointsAdded) {
    int i = 0;
    int j = 0;
    if (improvementSoFar > desiredImprovement) {
        willNotCompute = true;
        return -1;
    }
    if (inRangeD1->size() > 0 && inRangeHorizontalExpansion->size() > 0) {
        if (ascending ? (inRangeD1->at(i).coordinates[d1] < inRangeHorizontalExpansion->at(j).coordinates[d1]) : (inRangeD1->at(i).coordinates[d1] > inRangeHorizontalExpansion->at(j).coordinates[d1])) {
            if (improvementSoFar == desiredImprovement) {
                if (inRangeD1->at(i).isBlue != needsBlues && inRangeD1->at(i).coordinates[d1] != previousImprovingPoint.coordinates[d1]) {
                    willNotCompute = true;
                    return -1;
                } else {
                    return inRangeD1->at(i).coordinates[d1];
                }                
            }
            if (inRangeD1->at(0).coordinates[d1] == previousImprovingPoint.coordinates[d1]) {
                i++;
            }
        } else {
            if (improvementSoFar == desiredImprovement) {
                if (inRangeHorizontalExpansion->at(j).isBlue != needsBlues) {
                    willNotCompute = true;
                    return -1;
                } else {
                    return inRangeHorizontalExpansion->at(j).coordinates[d1];
                }
                
            }
        }
    } else if (inRangeD1->size() > 0) {
        if (improvementSoFar == desiredImprovement) {
            if (inRangeD1->at(i).isBlue != needsBlues && inRangeD1->at(i).coordinates[d1] != previousImprovingPoint.coordinates[d1]) {
                willNotCompute = true;
                return -1;
            } else {
                return inRangeD1->at(i).coordinates[d1];
            }           
        }
        if (inRangeD1->at(0).coordinates[d1] == previousImprovingPoint.coordinates[d1]) {
            i++;
        }
    } else if (inRangeHorizontalExpansion->size() > 0) {
        if (improvementSoFar == desiredImprovement) {
            if (inRangeHorizontalExpansion->at(j).isBlue != needsBlues) {
                willNotCompute = true;
                return -1;
            } else {
                return inRangeHorizontalExpansion->at(j).coordinates[d1];
            }
        }
    } else {
        willNotCompute = true;
        return -1;
    }


    while (i < inRangeD1->size() && j < inRangeHorizontalExpansion->size()) {
        pointsAdded++;
        if (ascending ? (inRangeD1->at(i).coordinates[d1] < inRangeHorizontalExpansion->at(j).coordinates[d1]) : (inRangeD1->at(i).coordinates[d1] > inRangeHorizontalExpansion->at(j).coordinates[d1])) {
            if (needsBlues == inRangeD1->at(i).isBlue) {
                improvementSoFar++;
            } else {
                improvementSoFar--;
            }
            if (improvementSoFar >= desiredImprovement) {
                return inRangeD1->at(i).coordinates[d1];
            }
            i++;
        } else {
            if (needsBlues == inRangeHorizontalExpansion->at(j).isBlue) {
                improvementSoFar++;
            } else {
                improvementSoFar--;
            } if (improvementSoFar >= desiredImprovement) {
                return inRangeHorizontalExpansion->at(j).coordinates[d1];
            }
            j++;
        }
    }
    while (i < inRangeD1->size()) {
        pointsAdded++;
        if (needsBlues == inRangeD1->at(i).isBlue) {
            improvementSoFar++;
        } else {
            improvementSoFar--;
        }
        if (improvementSoFar >= desiredImprovement) {
            return inRangeD1->at(i).coordinates[d1];
        }
        i++;
    }
    while (j < inRangeHorizontalExpansion->size()) {
        pointsAdded++;
        if (needsBlues == inRangeHorizontalExpansion->at(j).isBlue) {
            improvementSoFar++;
        } else {
            improvementSoFar--;
        } if (improvementSoFar >= desiredImprovement) {
            return inRangeHorizontalExpansion->at(j).coordinates[d1];
        }
        j++;
    }
    throw(runtime_error("Couldn't find fair in range"));
}

vector<Point> *mergePoints(vector<Point> *v1, vector<Point> *v2, int d, bool ascending) {
    vector<Point> *result = new vector<Point>();
    int i = 0;
    int j = 0;

    while (i < v1->size() && j < v2->size()) {
        if (ascending ? v1->at(i).coordinates[d] < v2->at(j).coordinates[d] : v1->at(i).coordinates[d] > v2->at(j).coordinates[d]) {
            result->push_back(v1->at(i));
            i++;
        } else {
            result->push_back(v2->at(j));
            j++;
        }
    }
    while (i < v1->size()) {
        result->push_back(v1->at(i));
        i++;
    }
    while (j < v2->size()) {
        result->push_back(v2->at(j));
        j++;
    }
    return result;

}

double** findD1D2Range(D1IP d1minimum, D1IP d1maximum, D2IP d2minimum, D2IP d2maximum, bool needsBlues, int d1, int d2, int improvementNeeded, vector<Point> *allPoints, int improvementAboveMaxRange, int improvementBelowMinRange, bool &willNotCompute, int &pointsAdded) {
    double **range = new double*[2];
    range[0] = new double[2]; //minimum
    range[1] = new double[2]; //maximum

    pointsAdded = 0;

    vector<Point> *inRangeD1Increasing = new vector<Point>();
    copy_if(d1maximum.inRange->begin(), d1maximum.inRange->end(), back_inserter(*inRangeD1Increasing), [d2minimum, d2maximum, d2] (Point pnt) {
        if (d2minimum.shrinking && pnt.coordinates[d2] <= d2minimum.improvingPoint.coordinates[d2]) {
            return false;
        }
        if (d2maximum.shrinking && pnt.coordinates[d2] >= d2maximum.improvingPoint.coordinates[d2]) {
            return false;
        }
        return true;
    });
    vector<Point> *inRangeHorizontalExpansionIncreasing;
    if (d2minimum.shrinking && d2maximum.shrinking) {
        inRangeHorizontalExpansionIncreasing = new vector<Point>();   
    } else if (d2minimum.shrinking) {
        inRangeHorizontalExpansionIncreasing = new vector<Point>();
        copy(d2maximum.inMaximumRange->begin(), d2maximum.inMaximumRange->end(), back_inserter(*inRangeHorizontalExpansionIncreasing));
    } else if (d2maximum.shrinking) {
        inRangeHorizontalExpansionIncreasing = new vector<Point>();
        copy(d2minimum.inMaximumRange->begin(), d2minimum.inMaximumRange->end(), back_inserter(*inRangeHorizontalExpansionIncreasing));
    } else {
        inRangeHorizontalExpansionIncreasing = mergePoints(d2minimum.inMaximumRange, d2maximum.inMaximumRange, d1, true);
    }
    sort(inRangeHorizontalExpansionIncreasing->begin(), inRangeHorizontalExpansionIncreasing->end(), pointSorter(d1, true));
    
    range[1][d1] = maxD1(inRangeD1Increasing, inRangeHorizontalExpansionIncreasing, true, d1, 0, improvementAboveMaxRange, needsBlues, willNotCompute, d1maximum.previousImprovingPoint, pointsAdded);
    
    if (willNotCompute) {
        return range;
    }

    vector<Point> *inRangeD1Decreasing = new vector<Point>();
    copy_if(d1minimum.inRange->begin(), d1minimum.inRange->end(), back_inserter(*inRangeD1Decreasing), [d2minimum, d2maximum, d2] (Point pnt) {
        if (d2minimum.shrinking && pnt.coordinates[d2] <= d2minimum.improvingPoint.coordinates[d2]) {
            return false;
        }
        if (d2maximum.shrinking && pnt.coordinates[d2] >= d2maximum.improvingPoint.coordinates[d2]) {
            return false;
        }
        return true;
    });

    vector<Point> *inRangeHorizontalExpansionDecreasing;
    if (d2minimum.shrinking && d2maximum.shrinking) {
        inRangeHorizontalExpansionDecreasing = new vector<Point>();   
    } else if (d2minimum.shrinking) {
        inRangeHorizontalExpansionDecreasing = new vector<Point>();
        copy(d2maximum.inMinimumRange->begin(), d2maximum.inMinimumRange->end(), back_inserter(*inRangeHorizontalExpansionDecreasing));
    } else if (d2maximum.shrinking) {
        inRangeHorizontalExpansionDecreasing = new vector<Point>();
        copy(d2minimum.inMinimumRange->begin(), d2minimum.inMinimumRange->end(), back_inserter(*inRangeHorizontalExpansionDecreasing));
    } else {
        inRangeHorizontalExpansionDecreasing = mergePoints(d2minimum.inMinimumRange, d2maximum.inMinimumRange, d1, false);
    }
    sort(inRangeHorizontalExpansionDecreasing->begin(), inRangeHorizontalExpansionDecreasing->end(), pointSorter(d1, false));

    range[0][d1] = maxD1(inRangeD1Decreasing, inRangeHorizontalExpansionDecreasing, false, d1, 0, improvementBelowMinRange, needsBlues, willNotCompute, d1minimum.previousImprovingPoint, pointsAdded);

    if (willNotCompute) {
        return range;
    }

    if (d2minimum.shrinking) {
        range[0][d2] = d2minimum.next.coordinates[d2];
    } else {
        range[0][d2] = d2minimum.improvingPoint.coordinates[d2];
    }

    if (d2maximum.shrinking) {
        range[1][d2] = d2maximum.next.coordinates[d2];
    } else {
        range[1][d2] = d2maximum.improvingPoint.coordinates[d2];
    }

    delete inRangeD1Increasing;
    delete inRangeHorizontalExpansionIncreasing;
    delete inRangeD1Decreasing;
    delete inRangeHorizontalExpansionDecreasing;

    return range;    
}

double** findD1D2RangeUpperBound(D1IP d1minimum, D1IP d1maximum, D2IP d2minimum, D2IP d2maximum, bool needsBlues, int d1, int d2, int improvementNeeded, vector<Point> *allPoints, int improvementNeededOnBottom, bool &willNotCompute, int &pointsAdded) {
    double **range = new double*[2];
    range[0] = new double[2]; //minimum
    range[1] = new double[2]; //maximum

    pointsAdded = 0;

    if ((d2minimum.shrinking && d1maximum.improvingPoint.coordinates[d2] < d2minimum.next.coordinates[d2]) || (d2maximum.shrinking && d1maximum.improvingPoint.coordinates[d2] > d2maximum.next.coordinates[d2])) {
        willNotCompute = true;
        return range;
    }

    range[1][d1] = d1maximum.improvingPoint.coordinates[d1];

    pointsAdded += d1maximum.inRange->size();

    vector<Point> *inRangeD1Decreasing = new vector<Point>();
    copy_if(d1minimum.inRange->begin(), d1minimum.inRange->end(), back_inserter(*inRangeD1Decreasing), [d2minimum, d2maximum, d2] (Point pnt) {
        if (d2minimum.shrinking && pnt.coordinates[d2] <= d2minimum.improvingPoint.coordinates[d2]) {
            return false;
        }
        if (d2maximum.shrinking && pnt.coordinates[d2] >= d2maximum.improvingPoint.coordinates[d2]) {
            return false;
        }
        return true;
    });

    vector<Point> *inRangeHorizontalExpansionDecreasing;
    if (d2minimum.shrinking && d2maximum.shrinking) {
        inRangeHorizontalExpansionDecreasing = new vector<Point>();   
    } else if (d2minimum.shrinking) {
        inRangeHorizontalExpansionDecreasing = new vector<Point>();
        copy(d2maximum.inMinimumRange->begin(), d2maximum.inMinimumRange->end(), back_inserter(*inRangeHorizontalExpansionDecreasing));
    } else if (d2maximum.shrinking) {
        inRangeHorizontalExpansionDecreasing = new vector<Point>();
        copy(d2minimum.inMinimumRange->begin(), d2minimum.inMinimumRange->end(), back_inserter(*inRangeHorizontalExpansionDecreasing));
    } else {
        inRangeHorizontalExpansionDecreasing = mergePoints(d2minimum.inMinimumRange, d2maximum.inMinimumRange, d1, false);
    }
    sort(inRangeHorizontalExpansionDecreasing->begin(), inRangeHorizontalExpansionDecreasing->end(), pointSorter(d1, false));

    range[0][d1] = maxD1(inRangeD1Decreasing, inRangeHorizontalExpansionDecreasing, false, d1, 0, improvementNeededOnBottom, needsBlues, willNotCompute, d1minimum.previousImprovingPoint, pointsAdded);

    if (willNotCompute) {
        return range;
    }

    if (d2minimum.shrinking) {
        range[0][d2] = d2minimum.next.coordinates[d2];
    } else {
        range[0][d2] = d2minimum.improvingPoint.coordinates[d2];
    }

    if (d2maximum.shrinking) {
        range[1][d2] = d2maximum.next.coordinates[d2];
    } else {
        range[1][d2] = d2maximum.improvingPoint.coordinates[d2];
    }
    delete inRangeD1Decreasing;
    delete inRangeHorizontalExpansionDecreasing;

    return range;    
}

double** findD1D2RangeLowerBound(D1IP d1minimum, D1IP d1maximum, D2IP d2minimum, D2IP d2maximum, bool needsBlues, int d1, int d2, int improvementNeeded, vector<Point> *allPoints, int improvementNeededOnTop, bool &willNotCompute, int &pointsAdded) {
    double **range = new double*[2];
    range[0] = new double[2]; //minimum
    range[1] = new double[2]; //maximum

    pointsAdded = 0;

    vector<Point> *inRangeD1Increasing = new vector<Point>();
    copy_if(d1maximum.inRange->begin(), d1maximum.inRange->end(), back_inserter(*inRangeD1Increasing), [d2minimum, d2maximum, d2] (Point pnt) {
        if (d2minimum.shrinking && pnt.coordinates[d2] <= d2minimum.improvingPoint.coordinates[d2]) {
            return false;
        }
        if (d2maximum.shrinking && pnt.coordinates[d2] >= d2maximum.improvingPoint.coordinates[d2]) {
            return false;
        }
        return true;
    });
    vector<Point> *inRangeHorizontalExpansionIncreasing;
    if (d2minimum.shrinking && d2maximum.shrinking) {
        inRangeHorizontalExpansionIncreasing = new vector<Point>();   
    } else if (d2minimum.shrinking) {
        inRangeHorizontalExpansionIncreasing = new vector<Point>();
        copy(d2maximum.inMaximumRange->begin(), d2maximum.inMaximumRange->end(), back_inserter(*inRangeHorizontalExpansionIncreasing));
    } else if (d2maximum.shrinking) {
        inRangeHorizontalExpansionIncreasing = new vector<Point>();
        copy(d2minimum.inMaximumRange->begin(), d2minimum.inMaximumRange->end(), back_inserter(*inRangeHorizontalExpansionIncreasing));
    } else {
        inRangeHorizontalExpansionIncreasing = mergePoints(d2minimum.inMaximumRange, d2maximum.inMaximumRange, d1, true);
    }
    sort(inRangeHorizontalExpansionIncreasing->begin(), inRangeHorizontalExpansionIncreasing->end(), pointSorter(d1, true));

    range[1][d1] = maxD1(inRangeD1Increasing, inRangeHorizontalExpansionIncreasing, true, d1, 0, improvementNeededOnTop, needsBlues, willNotCompute, d1maximum.previousImprovingPoint, pointsAdded);

    if (willNotCompute) {
        return range;
    }

    range[0][d1] = d1minimum.improvingPoint.coordinates[d1];

    pointsAdded += d1minimum.inRange->size();

    if (d2minimum.shrinking) {
        range[0][d2] = d2minimum.next.coordinates[d2];
    } else {
        range[0][d2] = d2minimum.improvingPoint.coordinates[d2];
    }

    if (d2maximum.shrinking) {
        range[1][d2] = d2maximum.next.coordinates[d2];
    } else {
        range[1][d2] = d2maximum.improvingPoint.coordinates[d2];
    }

    delete inRangeD1Increasing;
    delete inRangeHorizontalExpansionIncreasing;

    return range;    
}

//end of minor inefficencies mentioned above

//new inefficiency! It checks all points to find the points in the range. Instead, it should keep track of
//the points in range as it goes and use that. This is a factor of #points total! difference
void checkD1D2Similarity(double** range, int pointsInInitialRange, vector<Point> *allPoints, double &maxSimilarity, double **maxRange, int epsilon, int pointsAdded, int pointsRemoved) {
    // vector<Point> *pointsInTopIPRange = new vector<Point>();
    // copy_if(allPoints->begin(), allPoints->end(), back_inserter(*pointsInTopIPRange), [d1minimum, d1maximum, range, d1, d2] (Point pnt) {
    //     return pnt.coordinates[d2] >= range[0][d2] && pnt.coordinates[d2] <= range[1][d2] && pnt.coordinates[d1] > d1maximum.previousImprovingPoint.coordinates[d1] && pnt.coordinates[d1] <= range[1][d1];
    // });
    // vector<Point> *pointsInRange = new vector<Point>();
    // copy_if(allPoints->begin(), allPoints->end(), back_inserter(*pointsInRange), [range] (Point pnt) {
    //     return pnt.coordinates[0] >= range[0][0] && pnt.coordinates[0] <= range[1][0] && pnt.coordinates[1] >= range[0][1] && pnt.coordinates[1] <= range[1][1];
    // });

    // int bluesInRange = 0;
    // int redsInRange = 0;
    // for (Point pnt : *pointsInRange) {
    //     if (pnt.isBlue) {
    //         bluesInRange++;
    //     } else {
    //         redsInRange++;
    //     }
    // }

    // if (abs(redsInRange-bluesInRange) > epsilon) {
    //     throw(runtime_error("Unfair range!"));
    // }
    // int countTop = 0;
    // for (Point p : *pointsInTopIPRange) {
    //     countTop += p.isBlue ? 1 : -1;
    // }
    // int countBottom = 0;
    // for (Point p : *pointsInBottomIPRange) {
    //     countBottom += p.isBlue ? 1 : -1;
    // }
    // cout << "Improvement in top and below top range: " << countTop + d2minimum.improvementBelowMaxRange + d2maximum.improvementBelowMaxRange + d1maximum.improvement + d1minimum.improvement - 2 << endl;
    // cout << "Improvement in bottom and above bottom range: " << countBottom + d2minimum.improvementAboveMinRange + d2maximum.improvementAboveMinRange + d1maximum.improvement + d1minimum.improvement- 2 << endl;
    // cout << "------------" << endl;

    

    double similarity = (double) (pointsInInitialRange - pointsRemoved) / (double) (pointsAdded + pointsInInitialRange);

    if (similarity > maxSimilarity) {
        maxSimilarity = similarity;
        maxRange[0][0] = range[0][0];
        maxRange[0][1] = range[0][1];
        maxRange[1][0] = range[1][0];
        maxRange[1][1] = range[1][1];
    }

    delete [] range[0];
    delete [] range[1];
    delete [] range;
}

vector<D2IP> *stretchD2(vector<Point> *points, D1IP d1minimum, D1IP d1maximum, bool needsBlues, int d1, int improvementNeeded) {
    vector<D2IP> *d2ips = new vector<D2IP>();
    vector<Point> *currentInMinRange = new vector<Point>();
    vector<Point> *currentInMaxRange = new vector<Point>();

    int currentImprovement = 0;
    int currentImprovementAboveMinRange = 0;
    int currentImprovementBelowMaxRange = 0;
    int maxImprovement = 0;
    int pointsInCenter = 0;
    int improvementInCenter = 0;

    D2IP first;
    first.improvement = 0;
    first.improvementAboveMinRange = 0;
    first.improvementBelowMaxRange = 0;
    first.shrinking = false;
    first.inMaximumRange = new vector<Point>();
    first.inMinimumRange = new vector<Point>();
    first.improvingPoint = points->at(0);
    first.next.coordinates = new double[2];
    first.next.coordinates[0] = 0;
    first.next.coordinates[1] = 0;
    first.pointsInCenter = 0;
    first.pointsInStartingRange = 0;
    first.improvementInCenter = 0;
    d2ips->push_back(first);
    
    for (vector<Point>::iterator pnt = points->begin()+1; pnt != points->end(); pnt++) {
        if (pnt->isBlue == needsBlues) {
            currentImprovement++;
            if (pnt->coordinates[d1] > d1minimum.previousImprovingPoint.coordinates[d1]) {
                currentImprovementAboveMinRange++;
            }
            if (pnt->coordinates[d1] < d1maximum.previousImprovingPoint.coordinates[d1]) {
                currentImprovementBelowMaxRange++;
            }
            if (pnt->coordinates[d1] >= d1minimum.previousImprovingPoint.coordinates[d1] && pnt->coordinates[d1] <= d1maximum.previousImprovingPoint.coordinates[d1]) {
                improvementInCenter++;
            }
        } else {
            currentImprovement--;
            if (pnt->coordinates[d1] > d1minimum.previousImprovingPoint.coordinates[d1]) {
                currentImprovementAboveMinRange--;
            }
            if (pnt->coordinates[d1] < d1maximum.previousImprovingPoint.coordinates[d1]) {
                currentImprovementBelowMaxRange--;
            }
            if (pnt->coordinates[d1] >= d1minimum.previousImprovingPoint.coordinates[d1] && pnt->coordinates[d1] <= d1maximum.previousImprovingPoint.coordinates[d1]) {
                improvementInCenter--;
            }
        }
        if (pnt->coordinates[d1] <= d1minimum.previousImprovingPoint.coordinates[d1]) {
            currentInMinRange->push_back(*pnt);            
        }
        if (pnt->coordinates[d1] >= d1maximum.previousImprovingPoint.coordinates[d1]) {
            currentInMaxRange->push_back(*pnt);
        }
        if (pnt->coordinates[d1] > d1minimum.previousImprovingPoint.coordinates[d1] && pnt->coordinates[d1] < d1maximum.previousImprovingPoint.coordinates[d1]) {
            pointsInCenter++;
        }

        if (currentImprovement > maxImprovement) {
            maxImprovement = currentImprovement;
            D2IP d2ip;
            d2ip.improvement = maxImprovement;
            d2ip.improvementAboveMinRange = currentImprovementAboveMinRange;
            d2ip.improvementBelowMaxRange = currentImprovementBelowMaxRange;
            d2ip.inMaximumRange = new vector<Point>();
            copy(currentInMaxRange->begin(), currentInMaxRange->end(), back_inserter(*d2ip.inMaximumRange));
            d2ip.inMinimumRange = new vector<Point>();
            copy(currentInMinRange->begin(), currentInMinRange->end(), back_inserter(*d2ip.inMinimumRange));
            d2ip.shrinking = false;
            d2ip.improvingPoint = *pnt;
            d2ip.next.coordinates = new double [2];
            d2ip.next.coordinates[0] = 0;
            d2ip.next.coordinates[1] = 0;
            d2ip.pointsInCenter = pointsInCenter;
            d2ip.improvementInCenter = improvementInCenter;
            d2ips->push_back(d2ip);

            if (maxImprovement >= improvementNeeded) {
                delete currentInMinRange;
                delete currentInMaxRange;
                return d2ips;
            }
        }
    }

    delete currentInMinRange;
    delete currentInMaxRange;
    return d2ips;
}

vector<D2IP> *shrinkD2(vector<Point> *points, D1IP d1minimum, D1IP d1maximum, bool needsBlues, int d1, int improvementNeeded, Point pointPreceeding) {
    vector<D2IP> *d2ips = new vector<D2IP>();
    vector<Point> *currentInMinRange = new vector<Point>();
    vector<Point> *currentInMaxRange = new vector<Point>();

    int currentImprovement = 0;
    int currentImprovementAboveMinRange = 0;
    int currentImprovementBelowMaxRange = 0;
    int maxImprovement = 0;
    int pointsInStartingRange = 0;
    int pointsInCenter = 0;
    int improvementInCenter = 0;

    D2IP first;
    first.improvement = 0;
    first.improvementAboveMinRange = 0;
    first.improvementBelowMaxRange = 0;
    first.shrinking = true;
    first.inMaximumRange = new vector<Point>();
    first.inMinimumRange = new vector<Point>();
    first.improvingPoint = pointPreceeding;
    first.next = points->at(0);
    first.pointsInStartingRange = 0;
    first.pointsInCenter = 0;
    first.improvementInCenter = 0;
    d2ips->push_back(first);


    for (vector<Point>::iterator pnt = points->begin(); pnt != points->end()-1; pnt++) {
        if (pnt->isBlue != needsBlues) {
            currentImprovement++;
            if (pnt->coordinates[d1] >= d1minimum.previousImprovingPoint.coordinates[d1]) {
                currentImprovementAboveMinRange++;
            }
            if (pnt->coordinates[d1] <= d1maximum.previousImprovingPoint.coordinates[d1]) {
                currentImprovementBelowMaxRange++;
            }
            if (pnt->coordinates[d1] >= d1minimum.previousImprovingPoint.coordinates[d1] && pnt->coordinates[d1] <= d1maximum.previousImprovingPoint.coordinates[d1]) {
                improvementInCenter++;
            }

        } else {
            currentImprovement--;
            if (pnt->coordinates[d1] >= d1minimum.previousImprovingPoint.coordinates[d1]) {
                currentImprovementAboveMinRange--;
            }
            if (pnt->coordinates[d1] <= d1maximum.previousImprovingPoint.coordinates[d1]) {
                currentImprovementBelowMaxRange--;
            }
            if (pnt->coordinates[d1] >= d1minimum.previousImprovingPoint.coordinates[d1] && pnt->coordinates[d1] <= d1maximum.previousImprovingPoint.coordinates[d1]) {
                improvementInCenter--;
            }

        }
        if (pnt->coordinates[d1] < d1minimum.previousImprovingPoint.coordinates[d1]) {
            currentInMinRange->push_back(*pnt);            
        }
        if (pnt->coordinates[d1] > d1maximum.previousImprovingPoint.coordinates[d1]) {
            currentInMaxRange->push_back(*pnt);
        }
        if (pnt->inRange) {
            pointsInStartingRange++;
        }
        if (pnt->coordinates[d1] > d1minimum.previousImprovingPoint.coordinates[d1] && pnt->coordinates[d1] < d1maximum.previousImprovingPoint.coordinates[d1]) {
            pointsInCenter++;
        }
        if (currentImprovement > maxImprovement) {
            maxImprovement = currentImprovement;
            D2IP d2ip;
            d2ip.improvement = maxImprovement;
            d2ip.improvementAboveMinRange = currentImprovementAboveMinRange;
            d2ip.improvementBelowMaxRange = currentImprovementBelowMaxRange;
            d2ip.inMaximumRange = new vector<Point>();
            copy(currentInMaxRange->begin(), currentInMaxRange->end(), back_inserter(*d2ip.inMaximumRange));
            d2ip.inMinimumRange = new vector<Point>();
            copy(currentInMinRange->begin(), currentInMinRange->end(), back_inserter(*d2ip.inMinimumRange));
            d2ip.shrinking = true;
            d2ip.improvingPoint = *pnt;
            d2ip.next = *(pnt+1);
            d2ip.pointsInCenter = pointsInCenter;
            d2ip.pointsInStartingRange = pointsInStartingRange;
            d2ip.improvementInCenter = improvementInCenter;
            d2ips->push_back(d2ip);

            if (maxImprovement >= improvementNeeded) {
                delete currentInMinRange;
                delete currentInMaxRange;
                return d2ips;
            }
        }
    }

    delete currentInMinRange;
    delete currentInMaxRange;
    return d2ips;
}


void computeMinimumFairnessForD1Combination(D1IP d1minimum, D1IP d1maximum, vector<Point> *pointsAscendingD2, vector<Point> *pointsDescendingD2, int d1, int d2, double d2InitialMax, double d2InitialMin, bool needsBlues, int improvementNeeded, int pointsInInitialRange, double &maxSimilarity, double **maxRange, int epsilon) {
    vector<Point> *pointsAscendingGreater = new vector<Point>();
    copy_if(pointsAscendingD2->begin(), pointsAscendingD2->end(), back_inserter(*pointsAscendingGreater), [d1minimum, d1maximum, d1, d2, d2InitialMax](Point pnt) {
        return pnt.coordinates[d1] <= d1maximum.improvingPoint.coordinates[d1] && pnt.coordinates[d1] >= d1minimum.improvingPoint.coordinates[d1] && pnt.coordinates[d2] >= d2InitialMax;
    });
    vector<D2IP> *d2ipAscendingGreater = stretchD2(pointsAscendingGreater, d1minimum, d1maximum, needsBlues, d1, improvementNeeded - d1minimum.improvement - d1maximum.improvement);

    vector<Point> *pointsDescendingLess = new vector<Point>();
    copy_if(pointsDescendingD2->begin(), pointsDescendingD2->end(), back_inserter(*pointsDescendingLess), [d1minimum, d1maximum, d1, d2, d2InitialMin](Point pnt) {
        return pnt.coordinates[d1] <= d1maximum.improvingPoint.coordinates[d1] && pnt.coordinates[d1] >= d1minimum.improvingPoint.coordinates[d1] && pnt.coordinates[d2] <= d2InitialMin;
    });
    vector<D2IP> *d2ipDescendingLess = stretchD2(pointsDescendingLess, d1minimum, d1maximum, needsBlues, d1, improvementNeeded - d1minimum.improvement - d1maximum.improvement);

    vector<Point> *pointsAscendingInside = new vector<Point>();
    copy_if(pointsAscendingD2->begin(), pointsAscendingD2->end(), back_inserter(*pointsAscendingInside), [d1minimum, d1maximum, d1, d2, d2InitialMax, d2InitialMin](Point pnt) {
        return pnt.coordinates[d1] <= d1maximum.improvingPoint.coordinates[d1] && pnt.coordinates[d1] >= d1minimum.improvingPoint.coordinates[d1] && pnt.coordinates[d2] >= d2InitialMin && pnt.coordinates[d2] <= d2InitialMax;
    });
    Point pointDirectlyMin = pointsAscendingD2->at(0);
    for (Point p : *pointsAscendingD2) {
        if (p.coordinates[d2] < pointsAscendingInside->at(0).coordinates[d2]) {
            pointDirectlyMin = p;
        }
    }
    vector<D2IP> *d2ipAscendingInside = shrinkD2(pointsAscendingInside, d1minimum, d1maximum, needsBlues, d1, improvementNeeded - d1minimum.improvement - d1maximum.improvement, pointDirectlyMin);

    vector<Point> *pointsDescendingInside = new vector<Point>();
    copy_if(pointsDescendingD2->begin(), pointsDescendingD2->end(), back_inserter(*pointsDescendingInside), [d1minimum, d1maximum, d1, d2, d2InitialMax, d2InitialMin](Point pnt) {
        return pnt.coordinates[d1] <= d1maximum.improvingPoint.coordinates[d1] && pnt.coordinates[d1] >= d1minimum.improvingPoint.coordinates[d1] && pnt.coordinates[d2] >= d2InitialMin && pnt.coordinates[d2] <= d2InitialMax;
    });
    Point pointDirectlyMax = pointsDescendingD2->at(0);
    for (Point p : *pointsDescendingD2) {
        if (p.coordinates[d2] > pointsDescendingInside->at(0).coordinates[d2]) {
            pointDirectlyMax = p;
        }
    }
    vector<D2IP> *d2ipDescendingInside = shrinkD2(pointsDescendingInside, d1minimum, d1maximum, needsBlues, d1, improvementNeeded - d1minimum.improvement - d1maximum.improvement, pointDirectlyMax);

    int pointsAdded;

    for (D2IP ascendingGreater : *d2ipAscendingGreater) {
        for (D2IP descendingLess : *d2ipDescendingLess) {
            int improvementNeededInRegions = improvementNeeded - ascendingGreater.improvementInCenter - descendingLess.improvementInCenter - d1maximum.improvement - d1minimum.improvement + 2;
            int improvementAboveMaxRange = ascendingGreater.improvementAboveMinRange - ascendingGreater.improvementInCenter + descendingLess.improvementAboveMinRange - descendingLess.improvementInCenter;
            int improvementBelowMinRange = ascendingGreater.improvementBelowMaxRange - ascendingGreater.improvementInCenter + descendingLess.improvementBelowMaxRange - descendingLess.improvementInCenter;
            if (ascendingGreater.improvement + descendingLess.improvement + d1minimum.improvement + d1maximum.improvement - 2 >= improvementNeeded) {
                int improvementNeededOnTop = improvementNeededInRegions; 
                while (improvementNeededOnTop >= 0) {
                    int improvementNeededOnBottom = improvementNeededInRegions - improvementNeededOnTop;
                    if (improvementNeededOnTop <= improvementAboveMaxRange && improvementNeededOnBottom <= improvementBelowMinRange) {
                        bool willNotCompute = false;
                        double **range = findD1D2Range(d1minimum, d1maximum, descendingLess, ascendingGreater, needsBlues, d1, d2, improvementNeeded, pointsAscendingD2, improvementNeededOnTop, improvementNeededOnBottom, willNotCompute, pointsAdded);
                        if (!willNotCompute) {
                            checkD1D2Similarity(range, pointsInInitialRange, pointsAscendingD2, maxSimilarity, maxRange, epsilon, pointsAdded + d1minimum.pointsBelowRange + d1maximum.pointsBelowRange + ascendingGreater.pointsInCenter + descendingLess.pointsInCenter, 0);
                        }
                    }
                    improvementNeededOnTop--;
                } 
            }
            if (ascendingGreater.improvement + descendingLess.improvement + d1maximum.improvement + d1minimum.improvement - 1 >= improvementNeeded) {
                int improvementNeededOnBottom = improvementNeededInRegions - improvementAboveMaxRange - 1;
                if (improvementNeededOnBottom <= improvementBelowMinRange) {
                    bool willNotCompute = false;
                    double **range = findD1D2RangeUpperBound(d1minimum, d1maximum, descendingLess, ascendingGreater, needsBlues, d1, d2, improvementNeeded, pointsAscendingD2, improvementNeededOnBottom, willNotCompute, pointsAdded);
                    if (!willNotCompute) {
                        checkD1D2Similarity(range, pointsInInitialRange, pointsAscendingD2, maxSimilarity, maxRange, epsilon, pointsAdded + d1minimum.pointsBelowRange + d1maximum.pointsBelowRange + ascendingGreater.pointsInCenter + descendingLess.pointsInCenter, 0);
                    }
                }
                int improvementNeededOnTop = improvementNeededInRegions - improvementBelowMinRange - 1;
                if (improvementNeededOnTop <= improvementAboveMaxRange) {
                    bool willNotCompute = false;
                    double **range = findD1D2RangeLowerBound(d1minimum, d1maximum, descendingLess, ascendingGreater, needsBlues, d1, d2, improvementNeeded, pointsAscendingD2, improvementNeededOnTop, willNotCompute, pointsAdded);
                    if (!willNotCompute) {
                        checkD1D2Similarity(range, pointsInInitialRange, pointsAscendingD2, maxSimilarity, maxRange, epsilon, pointsAdded + d1minimum.pointsBelowRange + d1maximum.pointsBelowRange + ascendingGreater.pointsInCenter + descendingLess.pointsInCenter, 0);
                    }
                }
            }
            if (d1maximum.improvement + d1minimum.improvement + ascendingGreater.improvement + descendingLess.improvement >= improvementNeeded) {
                double **range = new double*[2];
                range[0] = new double[2];
                range[1] = new double[2];
                range[1][d1] = d1maximum.improvingPoint.coordinates[d1];
                range[0][d1] = d1minimum.improvingPoint.coordinates[d1];
                range[1][d2] = ascendingGreater.improvingPoint.coordinates[d2];
                range[0][d2] = descendingLess.improvingPoint.coordinates[d2];
                checkD1D2Similarity(range, pointsInInitialRange, pointsAscendingD2, maxSimilarity, maxRange, epsilon, d1maximum.pointsBelowRange + d1minimum.pointsBelowRange + d1maximum.inRange->size() + d1minimum.inRange->size() + ascendingGreater.pointsInCenter + descendingLess.pointsInCenter, 0);
            }
        }
    }


    for (D2IP ascendingGreater : *d2ipAscendingGreater) {
        for (D2IP ascendingInside : *d2ipAscendingInside) {
            int improvementNeededInRegions = improvementNeeded - ascendingGreater.improvementInCenter - ascendingInside.improvementInCenter - d1maximum.improvement - d1minimum.improvement + 2;
            int improvementAboveMaxRange = ascendingGreater.improvementAboveMinRange - ascendingGreater.improvementInCenter + ascendingInside.improvementAboveMinRange - ascendingInside.improvementInCenter;
            int improvementBelowMinRange = ascendingGreater.improvementBelowMaxRange - ascendingGreater.improvementInCenter + ascendingInside.improvementBelowMaxRange - ascendingInside.improvementInCenter;
            if (ascendingGreater.improvement + ascendingInside.improvement + d1minimum.improvement + d1maximum.improvement - 2 >= improvementNeeded) {
                int improvementNeededOnTop = improvementNeededInRegions; 
                while (improvementNeededOnTop >= 0) {
                    int improvementNeededOnBottom = improvementNeededInRegions - improvementNeededOnTop;
                    if (improvementNeededOnTop <= improvementAboveMaxRange && improvementNeededOnBottom <= improvementBelowMinRange) {    
                        bool willNotCompute = false;
                        double **range = findD1D2Range(d1minimum, d1maximum, ascendingInside, ascendingGreater, needsBlues, d1, d2, improvementNeeded, pointsAscendingD2, improvementNeededOnTop, improvementNeededOnBottom, willNotCompute, pointsAdded);
                        if (!willNotCompute) {
                            checkD1D2Similarity(range, pointsInInitialRange, pointsAscendingD2, maxSimilarity, maxRange, epsilon, pointsAdded + d1minimum.pointsBelowRange + d1maximum.pointsBelowRange + ascendingGreater.pointsInCenter, ascendingInside.pointsInStartingRange);
                        }
                    }
                    improvementNeededOnTop--;
                } 
            }
            if (ascendingGreater.improvement + ascendingInside.improvement + d1maximum.improvement + d1minimum.improvement - 1 >= improvementNeeded) {
                int improvementNeededOnBottom = improvementNeededInRegions - improvementAboveMaxRange - 1;
                if (improvementNeededOnBottom <= improvementBelowMinRange) {
                    bool willNotCompute = false;
                    double **range = findD1D2RangeUpperBound(d1minimum, d1maximum, ascendingInside, ascendingGreater, needsBlues, d1, d2, improvementNeeded, pointsAscendingD2, improvementNeededOnBottom, willNotCompute, pointsAdded);
                    if (!willNotCompute) {
                        checkD1D2Similarity(range, pointsInInitialRange, pointsAscendingD2, maxSimilarity, maxRange, epsilon, pointsAdded + d1minimum.pointsBelowRange + d1maximum.pointsBelowRange + ascendingGreater.pointsInCenter, ascendingInside.pointsInStartingRange);
                    }
                }
                int improvementNeededOnTop = improvementNeededInRegions - improvementBelowMinRange - 1;
                if (improvementNeededOnTop <= improvementAboveMaxRange) {
                    bool willNotCompute = false;
                    double **range = findD1D2RangeLowerBound(d1minimum, d1maximum, ascendingInside, ascendingGreater, needsBlues, d1, d2, improvementNeeded, pointsAscendingD2, improvementNeededOnTop, willNotCompute, pointsAdded);
                    if (!willNotCompute) {
                        checkD1D2Similarity(range, pointsInInitialRange, pointsAscendingD2, maxSimilarity, maxRange, epsilon, pointsAdded + d1minimum.pointsBelowRange + d1maximum.pointsBelowRange + ascendingGreater.pointsInCenter, ascendingInside.pointsInStartingRange);
                    }
                }
            }
            if (d1maximum.improvement + d1minimum.improvement + ascendingGreater.improvement + ascendingInside.improvement >= improvementNeeded) {
                double **range = new double*[2];
                range[0] = new double[2];
                range[1] = new double[2];
                range[1][d1] = d1maximum.improvingPoint.coordinates[d1];
                range[0][d1] = d1minimum.improvingPoint.coordinates[d1];
                range[1][d2] = ascendingGreater.improvingPoint.coordinates[d2];
                range[0][d2] = ascendingInside.next.coordinates[d2];
                checkD1D2Similarity(range, pointsInInitialRange, pointsAscendingD2, maxSimilarity, maxRange, epsilon, d1maximum.pointsBelowRange + d1minimum.pointsBelowRange + d1maximum.inRange->size() - ascendingInside.inMaximumRange->size() + d1minimum.inRange->size() - ascendingInside.inMinimumRange->size() - ascendingInside.pointsInCenter + ascendingInside.pointsInStartingRange + ascendingGreater.pointsInCenter, ascendingInside.pointsInStartingRange);
            }
        }
    }

    for (D2IP descendingInside : *d2ipDescendingInside) {
        for (D2IP descendingLess : *d2ipDescendingLess) {
            int improvementNeededInRegions = improvementNeeded - descendingLess.improvementInCenter - descendingInside.improvementInCenter - d1maximum.improvement - d1minimum.improvement + 2;
            int improvementAboveMaxRange = descendingInside.improvementAboveMinRange - descendingInside.improvementInCenter + descendingLess.improvementAboveMinRange - descendingLess.improvementInCenter;
            int improvementBelowMinRange = descendingInside.improvementBelowMaxRange - descendingInside.improvementInCenter + descendingLess.improvementBelowMaxRange - descendingLess.improvementInCenter;
            if (descendingInside.improvement + descendingLess.improvement + d1minimum.improvement + d1maximum.improvement - 2 >= improvementNeeded) {
                int improvementNeededOnTop = improvementNeededInRegions; 
                while (improvementNeededOnTop >= 0) {
                    int improvementNeededOnBottom = improvementNeededInRegions - improvementNeededOnTop;
                    if (improvementNeededOnTop <= improvementAboveMaxRange && improvementNeededOnBottom <= improvementBelowMinRange) {    
                        bool willNotCompute = false;
                        double **range = findD1D2Range(d1minimum, d1maximum, descendingLess, descendingInside, needsBlues, d1, d2, improvementNeeded, pointsAscendingD2, improvementNeededOnTop, improvementNeededOnBottom, willNotCompute, pointsAdded);
                        if (!willNotCompute) {
                            checkD1D2Similarity(range, pointsInInitialRange, pointsAscendingD2, maxSimilarity, maxRange, epsilon, pointsAdded + d1minimum.pointsBelowRange + d1maximum.pointsBelowRange + descendingLess.pointsInCenter, descendingInside.pointsInStartingRange);
                        }
                    }
                    improvementNeededOnTop--;
                }
            }
            if (descendingLess.improvement + descendingInside.improvement + d1maximum.improvement + d1minimum.improvement - 1 >= improvementNeeded) {
                int improvementNeededOnBottom = improvementNeededInRegions - improvementAboveMaxRange - 1;
                if (improvementNeededOnBottom <= improvementBelowMinRange) {
                    bool willNotCompute = false;
                    double **range = findD1D2RangeUpperBound(d1minimum, d1maximum, descendingLess, descendingInside, needsBlues, d1, d2, improvementNeeded, pointsAscendingD2, improvementNeededOnBottom, willNotCompute, pointsAdded);
                    if (!willNotCompute) {
                        checkD1D2Similarity(range, pointsInInitialRange, pointsAscendingD2, maxSimilarity, maxRange, epsilon, pointsAdded + d1minimum.pointsBelowRange + d1maximum.pointsBelowRange + descendingLess.pointsInCenter, descendingInside.pointsInStartingRange);
                    }
                }
                int improvementNeededOnTop = improvementNeededInRegions - improvementBelowMinRange - 1;
                if (improvementNeededOnTop <= improvementAboveMaxRange) {
                    bool willNotCompute = false;
                    double **range = findD1D2RangeLowerBound(d1minimum, d1maximum, descendingLess, descendingInside, needsBlues, d1, d2, improvementNeeded, pointsAscendingD2, improvementNeededOnTop, willNotCompute, pointsAdded);
                    if (!willNotCompute) {
                        checkD1D2Similarity(range, pointsInInitialRange, pointsAscendingD2, maxSimilarity, maxRange, epsilon, pointsAdded + d1minimum.pointsBelowRange + d1maximum.pointsBelowRange + descendingLess.pointsInCenter, descendingInside.pointsInStartingRange);
                    }
                }
            }
            if (d1maximum.improvement + d1minimum.improvement + descendingLess.improvement + descendingInside.improvement >= improvementNeeded) {
                double **range = new double*[2];
                range[0] = new double[2];
                range[1] = new double[2];
                range[1][d1] = d1maximum.improvingPoint.coordinates[d1];
                range[0][d1] = d1minimum.improvingPoint.coordinates[d1];
                range[1][d2] = descendingInside.next.coordinates[d2];
                range[0][d2] = descendingLess.improvingPoint.coordinates[d2];
                checkD1D2Similarity(range, pointsInInitialRange, pointsAscendingD2, maxSimilarity, maxRange, epsilon, d1maximum.pointsBelowRange - descendingInside.inMaximumRange->size() + d1minimum.pointsBelowRange - descendingInside.inMinimumRange->size() + d1maximum.inRange->size() + d1minimum.inRange->size() + descendingLess.pointsInCenter - descendingInside.pointsInCenter + descendingInside.pointsInStartingRange, descendingInside.pointsInStartingRange);
            }
        }
    }

    for (D2IP descendingInside : *d2ipDescendingInside) {
        for (D2IP ascendingInside : *d2ipAscendingInside) {
            int improvementNeededInRegions = improvementNeeded - descendingInside.improvementInCenter - ascendingInside.improvementInCenter - d1maximum.improvement - d1minimum.improvement + 2;
            int improvementAboveMaxRange = descendingInside.improvementAboveMinRange - descendingInside.improvementInCenter + ascendingInside.improvementAboveMinRange - ascendingInside.improvementInCenter;
            int improvementBelowMinRange = descendingInside.improvementBelowMaxRange - descendingInside.improvementInCenter + ascendingInside.improvementBelowMaxRange - ascendingInside.improvementInCenter;
            if (descendingInside.improvement + ascendingInside.improvement + d1minimum.improvement + d1maximum.improvement - 2 >= improvementNeeded) {
                int improvementNeededOnTop = improvementNeededInRegions; 
                while (improvementNeededOnTop >= 0) {
                    int improvementNeededOnBottom = improvementNeededInRegions - improvementNeededOnTop;
                    if (improvementNeededOnTop <= improvementAboveMaxRange && improvementNeededOnBottom <= improvementBelowMinRange) {    
                        bool willNotCompute = false;
                        double **range = findD1D2Range(d1minimum, d1maximum, ascendingInside, descendingInside, needsBlues, d1, d2, improvementNeeded, pointsAscendingD2, improvementNeededOnTop, improvementNeededOnBottom, willNotCompute, pointsAdded);
                        if (!willNotCompute) {
                            checkD1D2Similarity(range, pointsInInitialRange, pointsAscendingD2, maxSimilarity, maxRange, epsilon, pointsAdded + d1minimum.pointsBelowRange + d1maximum.pointsBelowRange, descendingInside.pointsInStartingRange + ascendingInside.pointsInStartingRange);
                        }
                    }
                    improvementNeededOnTop--;
                }
            }
            if (ascendingInside.improvement + descendingInside.improvement + d1maximum.improvement + d1minimum.improvement - 1 >= improvementNeeded) {
                int improvementNeededOnBottom = improvementNeededInRegions - improvementAboveMaxRange - 1;
                if (improvementNeededOnBottom <= improvementBelowMinRange) {                
                    bool willNotCompute = false;
                    double **range = findD1D2RangeUpperBound(d1minimum, d1maximum, ascendingInside, descendingInside, needsBlues, d1, d2, improvementNeeded, pointsAscendingD2, improvementNeededOnBottom, willNotCompute, pointsAdded);
                    if (!willNotCompute) {
                        checkD1D2Similarity(range, pointsInInitialRange, pointsAscendingD2, maxSimilarity, maxRange, epsilon, pointsAdded + d1minimum.pointsBelowRange + d1maximum.pointsBelowRange, descendingInside.pointsInStartingRange + ascendingInside.pointsInStartingRange);
                    }
                }
                int improvementNeededOnTop = improvementNeededInRegions - improvementBelowMinRange - 1;
                if (improvementNeededOnTop <= improvementAboveMaxRange) {
                    bool willNotCompute = false;
                    double **range = findD1D2RangeLowerBound(d1minimum, d1maximum, ascendingInside, descendingInside, needsBlues, d1, d2, improvementNeeded, pointsAscendingD2, improvementNeededOnTop, willNotCompute, pointsAdded);
                    if (!willNotCompute) {
                        checkD1D2Similarity(range, pointsInInitialRange, pointsAscendingD2, maxSimilarity, maxRange, epsilon, pointsAdded + d1minimum.pointsBelowRange + d1maximum.pointsBelowRange, descendingInside.pointsInStartingRange + ascendingInside.pointsInStartingRange);
                    }
                }
            }
            if (d1maximum.improvement + d1minimum.improvement + ascendingInside.improvement + descendingInside.improvement >= improvementNeeded) {
                double **range = new double*[2];
                range[0] = new double[2];
                range[1] = new double[2];
                range[1][d1] = d1maximum.improvingPoint.coordinates[d1];
                range[0][d1] = d1minimum.improvingPoint.coordinates[d1];
                range[1][d2] = descendingInside.next.coordinates[d2];
                range[0][d2] = ascendingInside.next.coordinates[d2];
                if (d1maximum.improvingPoint.coordinates[d2] >= range[0][d2]  && d1maximum.improvingPoint.coordinates[d2] <= range[1][d2] && d1minimum.improvingPoint.coordinates[d2] >= range[0][d2] && d1minimum.improvingPoint.coordinates[d2] <= range[1][d2] ) {
                    checkD1D2Similarity(range, pointsInInitialRange, pointsAscendingD2, maxSimilarity, maxRange, epsilon, d1maximum.pointsBelowRange + d1minimum.pointsBelowRange + d1maximum.inRange->size() - ascendingInside.inMaximumRange->size() - descendingInside.inMaximumRange->size() + d1minimum.inRange->size() - ascendingInside.inMinimumRange->size() - descendingInside.inMinimumRange->size() - descendingInside.pointsInCenter + descendingInside.pointsInStartingRange - ascendingInside.pointsInCenter + ascendingInside.pointsInStartingRange, ascendingInside.pointsInCenter + descendingInside.pointsInCenter);
                }
            }
        }
    }

    delete pointsAscendingGreater;
    delete pointsDescendingLess;
    delete pointsDescendingInside;
    delete pointsAscendingInside;
    delete d2ipAscendingGreater;
    delete d2ipDescendingLess;
    delete d2ipDescendingInside;
    delete d2ipAscendingInside;
}

void computeExpansionFairness(vector<D1IP> *d1ipGreater, vector<D1IP> *d1ipLess, int improvementNeeded, vector<Point> *pointsAscendingD2, vector<Point> *pointsDescendingD2, int d1, int d2, double d2InitialMax, double d2InitialMin, bool needsBlues, int pointsInInitialRange, double &maxSimilarity, double **maxRange, int epsilon) {
    for (D1IP greater : *d1ipGreater) {
        for (D1IP less : *d1ipLess) {
            if (greater.improvement + less.improvement > improvementNeeded) {
                continue;
            }
            computeMinimumFairnessForD1Combination(less, greater, pointsAscendingD2, pointsDescendingD2, d1, d2, d2InitialMax, d2InitialMin, needsBlues, improvementNeeded, pointsInInitialRange, maxSimilarity, maxRange, epsilon);
        }
    }
}

//Bug here! We're not considering the final range s.t. extending into it could potentially have an improvement, but
//the points run out before there's an improving point. This is a really rare edge case, so I'm not going to handle
//it right now, but it needs to be fixed in order to find all potential ranges.
//more...
//Additionally, I'm not 100% sure about the improvementNeeded + 1, but it is there to allow for the point s.t. the
//improving point range is == improvement needed.
vector<D1IP> *buildD1Expansion(vector<Point> *points, bool needsBlues, int improvementNeeded) {
    vector<D1IP> *d1ips = new vector<D1IP>();

    vector<Point> *inRange = new vector<Point>();
    Point previousPoint = points->at(0);
    inRange->push_back(previousPoint);

    int currentImprovement = 0;
    int maxImprovement = 0;
    int pointsBelowRange = 0;
    int pointsAdded = 0;

    for (vector<Point>::iterator point = points->begin()+1; point != points->end(); point++) {
        pointsAdded++;
        if (point->isBlue == needsBlues) {
            currentImprovement++;
        } else {
            currentImprovement--;
        }
        if (currentImprovement > maxImprovement) {
            maxImprovement = currentImprovement;
            D1IP d1ip;
            d1ip.improvement = maxImprovement;
            d1ip.inRange = inRange;
            d1ip.previousImprovingPoint = previousPoint;
            d1ip.improvingPoint = *point;
            d1ip.pointsBelowRange = pointsBelowRange;

            d1ips->push_back(d1ip);
            
            if (maxImprovement == improvementNeeded + 1) {
                return d1ips;
            }

            previousPoint = *point;
            inRange = new vector<Point>();
            inRange->push_back(previousPoint);
            pointsBelowRange = pointsAdded;
        } else {
            inRange->push_back(*point);
        }
    }

    //we're currently deleting the final range, but to fix the above mentioned likely bug, we need to do something with it.
    delete inRange;
    return d1ips;
}

void buildAndComputeExpansionFairness(vector<Point> *pointsAscendingD1, vector<Point> *pointsDescendingD1, bool needsBlues, int improvementNeeded, double *minimumPointStartingRange, double *maximumPointStartingRange, int d1, int d2, vector<Point> *pointsAscendingD2, vector<Point> *pointsDescendingD2, int pointsInInitialRange, double &maxSimilarity, double **maxRange, int epsilon) {
    vector<Point> *ascendingGreaterD1 = new vector<Point>();
    copy_if(pointsAscendingD1->begin(), pointsAscendingD1->end(), back_inserter(*ascendingGreaterD1), [minimumPointStartingRange, maximumPointStartingRange, d1, d2] (Point pnt) {
        return pnt.coordinates[d2] >= minimumPointStartingRange[d2] && pnt.coordinates[d2] <= maximumPointStartingRange[d2] && pnt.coordinates[d1] >= maximumPointStartingRange[d1];
    });
    vector<D1IP> *ascendingD1IP = buildD1Expansion(ascendingGreaterD1, needsBlues, improvementNeeded);

    vector<Point> *descendingLessD1 = new vector<Point>();
    copy_if(pointsDescendingD1->begin(), pointsDescendingD1->end(), back_inserter(*descendingLessD1), [minimumPointStartingRange, maximumPointStartingRange, d1, d2] (Point pnt) {
        return pnt.coordinates[d2] >= minimumPointStartingRange[d2] && pnt.coordinates[d2] <= maximumPointStartingRange[d2] && pnt.coordinates[d1] <= minimumPointStartingRange[d1];
    });
    vector<D1IP> *descendingD1IP = buildD1Expansion(descendingLessD1, needsBlues, improvementNeeded);

    computeExpansionFairness(ascendingD1IP, descendingD1IP, improvementNeeded, pointsAscendingD2, pointsDescendingD2, d1, d2, maximumPointStartingRange[d2], minimumPointStartingRange[d2], needsBlues, pointsInInitialRange, maxSimilarity, maxRange, epsilon);

    delete ascendingGreaterD1;
    delete ascendingD1IP;
    delete descendingLessD1;
    delete descendingD1IP;

}

struct IPSD1 {
    Point improvingPoint;
    Point next;
    vector<Point> *intermediates;
    int improvement;
    int pointsRemoved;
};

struct IPSD2 {
    Point next;
    int improvement;
    int pointsRemoved;
};

double **findD1D2ShrinkRange(IPSD1 minD1, IPSD1 maxD1, IPSD2 minD2, IPSD2 maxD2, int d1, int d2, vector<Point> *allPoints, bool needsBlues) {
    double **range = new double*[2];
    range[0] = new double[2];
    range[1] = new double[2];

    range[0][d1] = minD1.next.coordinates[d1];
    range[1][d1] = maxD1.next.coordinates[d1];
    range[0][d2] = minD2.next.coordinates[d2];
    range[1][d2] = maxD2.next.coordinates[d2];

    return range;
}

void checkD1D2ShrinkingSimilarity(IPSD1 minD1, IPSD1 maxD1, IPSD2 minD2, IPSD2 maxD2, int d1, int d2, vector<Point> *allPoints, double &maxSimilarity, double **maxRange, int pointsInInitialRange, int epsilon, int improvementNeeded, bool needsBlues) {
    double **range = findD1D2ShrinkRange(minD1, maxD1, minD2, maxD2, d1, d2, allPoints, needsBlues);

    double similarity = (double) (pointsInInitialRange - minD1.pointsRemoved - maxD1.pointsRemoved - minD2.pointsRemoved - maxD2.pointsRemoved) / (double) pointsInInitialRange;

    if (similarity > maxSimilarity) {
        maxSimilarity = similarity;
        maxRange[0][0] = range[0][0];
        maxRange[0][1] = range[0][1];
        maxRange[1][0] = range[1][0];
        maxRange[1][1] = range[1][1];
    }

    delete [] range[0];
    delete [] range[1];
    delete [] range;
}

vector<IPSD2> *shrinkShrinkD2(vector<Point> *points, bool needsBlues, int improvementNeeded) {
    int currentMaximum = 0;
    int currentImprovement = 0; 
    int pointsRemoved = 0;

    vector<IPSD2> *ipsd2s = new vector<IPSD2>();

    
    IPSD2 first;
    first.improvement = 0;
    first.next = points->at(0);
    ipsd2s->push_back(first);

    for (vector<Point>::iterator point = points->begin(); point != points->end()-1; point++) {
        pointsRemoved++;
        if (point->isBlue != needsBlues) {
            currentImprovement++;
        } else {
            currentImprovement--;
        }
        if (currentImprovement > currentMaximum) {
            currentMaximum = currentImprovement;
            IPSD2 ipsd2;
            ipsd2.improvement = currentMaximum;
            ipsd2.next = *(point+1);
            ipsd2.pointsRemoved = pointsRemoved;
            ipsd2s->push_back(ipsd2);

            if (currentMaximum >= improvementNeeded) {
                return ipsd2s;
            }
        }
    }
    return ipsd2s;
}

void computeMinimumFairnessForD1ShrinkingCombination(IPSD1 d1minimum, IPSD1 d1maximum, vector<Point> *pointsAscendingD2, vector<Point> *pointsDescendingD2, int d1, int d2, double d2InitialMax, double d2InitialMin, bool needsBlues, int improvementNeeded, int pointsInInitialRange, double &maxSimilarity, double **maxRange, int epsilon) {
    vector<Point> *pointsDescendingInside = new vector<Point>();
    copy_if(pointsDescendingD2->begin(), pointsDescendingD2->end(), back_inserter(*pointsDescendingInside), [d1minimum, d1maximum, d1, d2, d2InitialMax, d2InitialMin](Point pnt) {
        return pnt.coordinates[d1] < d1maximum.improvingPoint.coordinates[d1] && pnt.coordinates[d1] > d1minimum.improvingPoint.coordinates[d1] && pnt.coordinates[d2] <= d2InitialMax && pnt.coordinates[d2] >= d2InitialMin;
    });
    if (pointsDescendingInside->size() == 0) {
        return;
    }

    vector<IPSD2> *ipsd2DescendingInside = shrinkShrinkD2(pointsDescendingInside, needsBlues, improvementNeeded - d1minimum.improvement - d1maximum.improvement);

    vector<Point> *pointsAscendingInside = new vector<Point>();
    copy_if(pointsAscendingD2->begin(), pointsAscendingD2->end(), back_inserter(*pointsAscendingInside), [d1minimum, d1maximum, d1, d2, d2InitialMax, d2InitialMin](Point pnt) {
        return pnt.coordinates[d1] < d1maximum.improvingPoint.coordinates[d1] && pnt.coordinates[d1] > d1minimum.improvingPoint.coordinates[d1] && pnt.coordinates[d2] <= d2InitialMax && pnt.coordinates[d2] >= d2InitialMin;
    });
    if (pointsAscendingInside->size() == 0) {
        return;
    }

    vector<IPSD2> *ipsd2AscendingInside = shrinkShrinkD2(pointsAscendingInside, needsBlues, improvementNeeded - d1minimum.improvement - d1maximum.improvement);

    int i = 0;
    int j = 0;
    while (i < (int) ipsd2DescendingInside->size() - 1  && d1minimum.improvement + d1maximum.improvement + ipsd2DescendingInside->at(i).improvement < improvementNeeded) {
        i++;
    }
        
    while (j < ipsd2AscendingInside->size() && d1minimum.improvement + d1maximum.improvement + ipsd2DescendingInside->at(i).improvement + ipsd2AscendingInside->at(j).improvement < improvementNeeded) {
        j++;
    }
    while (i >= 0 && j < ipsd2AscendingInside->size()) {
        if (ipsd2DescendingInside->at(i).next.coordinates[d2] >= ipsd2AscendingInside->at(j).next.coordinates[d2]) {
            checkD1D2ShrinkingSimilarity(d1minimum, d1maximum, ipsd2AscendingInside->at(j), ipsd2DescendingInside->at(i), d1, d2, pointsDescendingD2, maxSimilarity, maxRange, pointsInInitialRange, epsilon, improvementNeeded, needsBlues);
        } 
        i--;
        j++;
    }
    delete pointsDescendingInside;
    delete pointsAscendingInside;
    delete ipsd2AscendingInside;
    delete ipsd2DescendingInside;
}

void computeShrinkingFairness(vector<IPSD1> *descending, vector<IPSD1> *ascending, int d1, int d2, vector<Point> *pointsAscendingD2, vector<Point> *pointsDescendingD2, double d2InitialMax, double d2InitialMin, bool needsBlues, int improvementNeeded, int pointsInInitialRange, double &maxSimilarity, double **maxRange, int epsilon) {
    for (IPSD1 ipsMax : *descending) {
        for (IPSD1 ipsMin: *ascending) {
            if (ipsMax.next.coordinates[d1] >= ipsMin.next.coordinates[d1]) {
                computeMinimumFairnessForD1ShrinkingCombination(ipsMin, ipsMax, pointsAscendingD2, pointsDescendingD2, d1, d2, d2InitialMax, d2InitialMin, needsBlues, improvementNeeded, pointsInInitialRange, maxSimilarity, maxRange, epsilon);
            }
        }
    }
}

vector<IPSD1> *shrinkD1(vector<Point> *points, bool needsBlues, int improvementNeeded) {
    int maxImprovement = 0;
    int currentImprovement = 0;
    int pointsRemoved = 0;

    vector<IPSD1> *ipsd1s = new vector<IPSD1>();

    for (vector<Point>::iterator point = points->begin(); point < points->end()-1; point++) {
        pointsRemoved++;
        if (point->isBlue != needsBlues) {
            currentImprovement++;
        } else {
            currentImprovement--;
        }

        if (currentImprovement > maxImprovement) {
            maxImprovement = currentImprovement;
            IPSD1 ipsd1;
            ipsd1.next = *(point+1);
            ipsd1.improvement = maxImprovement;
            ipsd1.improvingPoint = *point;
            ipsd1.pointsRemoved = pointsRemoved;
            ipsd1s->push_back(ipsd1);

            if (maxImprovement >= improvementNeeded) {
                return ipsd1s;
            }
        }
    }
    return ipsd1s;
}

void buildAndComputeShrinkingFairness(vector<Point> *pointsD1Increasing, vector<Point> *pointsD1Decreasing, vector<Point> *pointsD2Increasing, vector<Point> *pointsD2Decreasing, bool needsBlues, int improvementNeeded, double *minimumPointStartingRange, double *maximumPointStartingRange, int d1, int d2, int pointsInInitialRange, double &maxSimilarity, double **maxRange, int epsilon) {    
    vector<Point> *pointsDescendingInsideD1 = new vector<Point>();
    copy_if(pointsD1Decreasing->begin(), pointsD1Decreasing->end(), back_inserter(*pointsDescendingInsideD1), [minimumPointStartingRange, maximumPointStartingRange, d1, d2](Point pnt) {
        return pnt.coordinates[d2] >= minimumPointStartingRange[d2] && pnt.coordinates[d2] <= maximumPointStartingRange[d2] && pnt.coordinates[d1] <= maximumPointStartingRange[d1] && pnt.coordinates[d1] >= minimumPointStartingRange[d1];
    });

    vector<Point> *pointsAscendingInsideD1 = new vector<Point>();
    copy_if(pointsD1Increasing->begin(), pointsD1Increasing->end(), back_inserter(*pointsAscendingInsideD1), [minimumPointStartingRange, maximumPointStartingRange, d1, d2](Point pnt) {
        return pnt.coordinates[d2] >= minimumPointStartingRange[d2] && pnt.coordinates[d2] <= maximumPointStartingRange[d2] && pnt.coordinates[d1] <= maximumPointStartingRange[d1] && pnt.coordinates[d1] >= minimumPointStartingRange[d1];
    });

    vector<IPSD1> *ipsd1Descending = shrinkD1(pointsDescendingInsideD1, needsBlues, improvementNeeded);
    vector<IPSD1> *ipsd1Ascending = shrinkD1(pointsAscendingInsideD1, needsBlues, improvementNeeded);

    computeShrinkingFairness(ipsd1Descending, ipsd1Ascending, d1, d2, pointsD2Increasing, pointsD2Decreasing, maximumPointStartingRange[d2], minimumPointStartingRange[d2], needsBlues, improvementNeeded, pointsInInitialRange, maxSimilarity, maxRange, epsilon);
}

void runSearch(vector<Point> *points, double *minimumPointStartingRange, double *maximumPointStartingRange) {
    vector<Point> *ascending0 = new vector<Point>();
    vector<Point> *descending0 = new vector<Point>();
    vector<Point> *ascending1 = new vector<Point>();
    vector<Point> *descending1 = new vector<Point>();

    copy(points->begin(), points->end(), back_inserter(*ascending0));
    copy(points->begin(), points->end(), back_inserter(*ascending1));
    copy(points->begin(), points->end(), back_inserter(*descending0));
    copy(points->begin(), points->end(), back_inserter(*descending1));

    sort(ascending0->begin(), ascending0->end(), pointSorter(0, true));
    sort(ascending1->begin(), ascending1->end(), pointSorter(1, true));
    sort(descending0->begin(), descending0->end(), pointSorter(0, false));
    sort(descending1->begin(), descending1->end(), pointSorter(1, false));

    double similarity = 0;
    double **range = new double*[2];
    range[0] = new double[2];
    range[1] = new double[2];

    //todo::
    int redsInRange = 0;
    int bluesInRange = 0;

    for (Point p : *points) {
        if (p.inRange && p.isBlue) {
            bluesInRange++;
        } else if (p.inRange && !p.isBlue) {
            redsInRange++;
        }
    }

    bool needsBlues = bluesInRange < redsInRange;
    int epsilon = ceil(0.025*(redsInRange+bluesInRange));
    int improvementNeeded = abs(redsInRange-bluesInRange) - epsilon;
    cout << "Epsilon: " << epsilon << endl;
    cout << "Initial disparity " << abs(redsInRange - bluesInRange) << endl;
    cout << "Initial range, {" << minimumPointStartingRange[0] << ", "<< minimumPointStartingRange[1] << "}, {" << maximumPointStartingRange[0] << ", " << maximumPointStartingRange[1] <<  "}" << endl;

    if (improvementNeeded <= 0) {
        cout << "Initial range is already fair" << endl;
        return;
    }
    int pointsInInitialRange = bluesInRange + redsInRange;

    buildAndComputeExpansionFairness(ascending0, descending0, needsBlues, improvementNeeded, minimumPointStartingRange, maximumPointStartingRange, 0, 1, ascending1, descending1, pointsInInitialRange, similarity, range, epsilon);
    buildAndComputeExpansionFairness(ascending1, descending1, needsBlues, improvementNeeded, minimumPointStartingRange, maximumPointStartingRange, 1, 0, ascending0, descending0, pointsInInitialRange, similarity, range, epsilon);
    buildAndComputeShrinkingFairness(ascending0, descending0, ascending1, descending1, needsBlues, improvementNeeded, minimumPointStartingRange, maximumPointStartingRange, 0, 1, pointsInInitialRange, similarity, range, epsilon);

    int blueInRange = 0;
    int redInRange = 0;

    for (Point p : *points) {
        if (p.coordinates[0] >= range[0][0] && p.coordinates[1] >= range[0][1] && p.coordinates[0] <= range[1][0] && p.coordinates[1] <= range[1][1]) {
            if (p.isBlue) {
                blueInRange++;
            } else {
                redInRange++;
            }
        }
    }

    cout << "Final disparity " << abs(blueInRange - redInRange) << endl;
    cout << "Similarity: " << similarity << endl;
    cout << "{" << range[0][0] << ", " << range[0][1] << "}, {" << range[1][0] << ", " << range[1][1] << "}" << endl;
    cout << endl;
}

vector<Point> *readPoints(string filename, double *minimum, double *maximum) {
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
        getline(file, line, '\n');
    }
    return points;  
}

int main() {
    string filename = "./data/uniform.csv";
    double *minimum = new double[2];
    double *maximum = new double[2];

    cout << setprecision(std::numeric_limits<double>::digits10 + 2);

    minimum[0] = 427.328;
    minimum[1] = 265.461;
    maximum[0] = 934.495;
    maximum[1] = 583.57;

    vector<Point> *points = readPoints(filename, minimum, maximum);

    double *trueMinimum = new double[2];
    double *trueMaximum = new double[2];
    trueMinimum[0] = maximum[0];
    trueMinimum[1] = maximum[1];
    trueMaximum[0] = minimum[0];
    trueMaximum[1] = minimum[1];
    for (Point p : *points) {
        if (p.inRange) {
            if (p.coordinates[0] < trueMinimum[0]) {
                trueMinimum[0] = p.coordinates[0];
            }
            if (p.coordinates[1] < trueMinimum[1]) {
                trueMinimum[1] = p.coordinates[1];
            }
            if (p.coordinates[0] > trueMaximum[0]) {
                trueMaximum[0] = p.coordinates[0];
            }
            if (p.coordinates[1] > trueMaximum[1]) {
                trueMaximum[1] = p.coordinates[1];
            }
        }
    }

    runSearch(points, trueMinimum, trueMaximum);
}