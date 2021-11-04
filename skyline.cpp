#include<iostream>
#include<boost/heap/pairing_heap.hpp>
#include<cmath>
#include<vector>
#include"range_tree.h"
#include<set>

const int dims = 2;
int blueWeight = 3, redWeight = 1, diff = 10;

using namespace std;

struct sps {
    double start[dims], end[dims];
    int blues, reds, intersection, union_, diffVal;

    sps(double s[dims], double e[dims], int blueCount, int redCount, int itemsIntersection, int itemsUnion) {
        for(int i=0; i < dims; ++i) {
            start[i] = s[i];
            end[i] = e[i];
        }
        // start = s;
        // end = e;
        blues= blueCount;
        reds = redCount;
        intersection = itemsIntersection;
        union_ = itemsUnion;
        diffVal = blueWeight * blues + redWeight * reds;
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

    double distanceJaccard() const {
        return 1-((double)intersection)/ ((double)union_);
    }

    bool operator<(sps const & sps2) const {
        return distanceJaccard() > sps2.distanceJaccard();
    }
};

int main(int argc, char *argv[]) {
    double inputRange[2][dims];
    set<vector<double>> uniqueRanges;
    inputRange[0][0] = -0.1;
    inputRange[0][1] = 51.31;
    inputRange[1][0] = 0;
    inputRange[1][1] = 51.44;

    double inputRangeCopy[2][dims];

    RangeTree rt(argv[1], inputRange[0], inputRange[1]);

    // Fix the input range
    double fixInput[2][dims];
    for(int i=0; i<dims; ++i) {
        fixInput[0][i] = inputRange[0][i];
        fixInput[1][i] = inputRange[1][i];
        inputRangeCopy[0][i] = inputRange[0][i];
        inputRangeCopy[1][i] = inputRange[1][i];
    }
    for(int i=0; i<dims; ++i) {
        for(int j=0; j<=1; ++j) {
            // Fix the range for dimension j
            fixInput[(j+1)%2][i] = inputRange[j][i];
            range fixEndPointsRange;
            fixEndPointsRange.start = fixInput[0];
            fixEndPointsRange.end = fixInput[1];
            fixEndPointsRange.inclusiveXMax = true;
            fixEndPointsRange.inclusiveXMin = true;
            fixEndPointsRange.inclusiveYMax = true;
            fixEndPointsRange.inclusiveYMin = true;
            vector<Point>* pointsFixRange = rt.findPointsInRange(fixEndPointsRange);
            if(pointsFixRange->size() == 0) {
                if(j == 0){
                    fixInput[0][i] = inputRange[0][i];
                    fixInput[1][i] = inputRange[1][i];
                } else {
                    fixInput[0][i] = inputRange[1][i];
                    fixInput[1][i] = inputRange[0][i];
                }
                Point* p = rt.findNearest(fixInput[0], fixInput[1], i);
                if(p) {
                    inputRange[j][i] = p->coordinates[i];
                }
                fixInput[0][i] = inputRange[0][i];
                fixInput[1][i] = inputRange[1][i];
                delete p;
            }
            delete pointsFixRange;
            fixInput[0][i] = inputRange[0][i];
            fixInput[1][i] = inputRange[1][i];
        }
    }

    // Create fake small and large boundary
    // Ideally use RangeTree to create the boundary
    double boundary[2][dims];
    boost::heap::pairing_heap<sps> H;
    // The two boundaries values are set to -1000. and 1000.
    // Works for the current application, may not work for other cases.
    for(int i=0; i<dims; ++i) {
        boundary[0][i] = -1000.;
        boundary[1][i] = 1000.;
    }

    int blueCount=0, redCount=0, itemsIntersection=0, itemsUnion=0;
    // Go through the data structure Ian has created and 
    // get blueCount, redCount, itemsIntersection, itemsUnion
    range initialRange;
    initialRange.start = inputRange[0];
    initialRange.end = inputRange[1];
    initialRange.inclusiveYMax = true;
    initialRange.inclusiveYMin = true;
    initialRange.inclusiveXMax = true;
    initialRange.inclusiveXMin = true;
    vector<Point> *pointsInRange = rt.findPointsInRange(initialRange);
    itemsIntersection = itemsUnion = pointsInRange->size();
    initialRange.start = inputRangeCopy[0];
    initialRange.end = inputRangeCopy[1];
    vector<Point> *pointsInRangeCopy = rt.findPointsInRange(initialRange);
    if(pointsInRange->size() != pointsInRangeCopy->size()) {
        cout << "Something fishy going on here!" << endl;
    }
    delete pointsInRangeCopy;
    for(int i=0; i<pointsInRange->size(); ++i) {
        if(pointsInRange->at(i).isBlue) {
            blueCount++;
        } else {
            redCount++;
        }
    }
    cout << "Number of points in range " << pointsInRange->size() << endl;
    cout << "Input range " << blueCount << " blues and " << redCount << "reds with a total of " << itemsUnion << "points." << endl;
    delete pointsInRange;
    struct sps first(inputRange[0], inputRange[1], blueCount, redCount, itemsIntersection, itemsUnion);
    set<vector<double>> ranges;
    vector<double> rangeFirst;
    for(int i=0; i<dims; ++i) {
        rangeFirst.push_back(first.start[i]);
    }
    for(int i=0; i<dims; ++i) {
        rangeFirst.push_back(first.end[i]);
    }
    ranges.insert(rangeFirst);
    H.push(first);
    int dummyCounter=0;
    while(!H.empty()) {

        struct sps topItem = H.top();
        if(++dummyCounter % 100 == 0) {
            cout << "Disparity : " << topItem.intersection << ":" << topItem.union_ << endl << "Heap size : " << H.size() << endl << "Set size" << ranges.size() << endl;
        }
        vector<double> rangeVectorRemove;
        for(int i=0; i<dims; ++i) {
            rangeVectorRemove.push_back(topItem.start[i]);
        }
        for(int i=0; i<dims; ++i) {
            rangeVectorRemove.push_back(topItem.end[i]);
        }
        ranges.erase(rangeVectorRemove);
        H.pop();
        int curDiffVal = topItem.blues * blueWeight - topItem.reds * redWeight;
        if(abs(curDiffVal) <= diff) {
            cout << "Range definition : " << endl;
            for(int i=0; i < dims; ++i) {
                cout << topItem.start[i] << "-" << topItem.end[i] << endl;
            }
            cout << "Disparity : " << topItem.intersection << ":" << topItem.union_ << ":" << double(topItem.intersection)/double(topItem.union_) << endl;
            return 0;
        }
        Point *boundaryPoints[2][dims];
        double boundaryVals[2][dims];
        for(int i=0; i<dims; ++i) {
            boundaryVals[0][i] = topItem.start[i];
            boundaryVals[1][i] = topItem.end[i];
        }
        // Find neighbours using Ian's method and also shrink range
        for(int dim=0; dim < dims; ++dim) {
            // Find the boundary
            boundaryVals[1][dim] = -1000.;
            boundaryPoints[0][dim] = rt.findNearest(boundaryVals[0], boundaryVals[1], dim);

            boundaryVals[0][dim] = topItem.end[dim];
            boundaryVals[1][dim] = 1000.;
            boundaryPoints[1][dim] = rt.findNearest(boundaryVals[0], boundaryVals[1], dim);

            boundaryVals[0][dim] = topItem.start[dim];
            boundaryVals[1][dim] = topItem.end[dim];

            // Skip points if any of the points boundaries is on the axis
            for(int boundaryType=0; boundaryType<2; ++boundaryType) {
                bool addRange = true;
                // If there is no point to expand the rectangle to or if the point being added is in range, skip
                if(boundaryPoints[boundaryType][dim] == nullptr || boundaryPoints[boundaryType][dim]->inRange) {
                    continue;
                }
                for(int dimVal=0; dimVal < dims; ++dimVal) {
                    if(boundaryPoints[boundaryType][dim]->coordinates[dimVal] == topItem.start[dim] || 
                            boundaryPoints[boundaryType][dim]->coordinates[dimVal] == topItem.end[dim]) {
                        addRange = false;
                        break;
                    }
                }
                if(addRange) {
                    // Generate new sps and add to heap
                    double rangeStart[dims], rangeEnd[dims];
                    for(int dim_=0; dim_ < dims; ++dim_) {
                        rangeStart[dim_] = topItem.start[dim_];
                        rangeEnd[dim_] = topItem.end[dim_];
                    }
                    if(boundaryType) {
                        rangeEnd[dim] = boundaryPoints[1][dim]->coordinates[dim];
                    } 
                    else {
                        rangeStart[dim] = boundaryPoints[0][dim]->coordinates[dim];
                    }
                    int blues = topItem.blues, reds = topItem.reds;
                    if(boundaryPoints[1][dim]->isBlue) {
                        blues++;
                    } else {
                        reds++;
                    }
                    vector<double> rangeVector;
                    for(int i=0; i<dims; ++i) {
                        rangeVector.push_back(rangeStart[i]);
                    }
                    for(int i=0; i<dims; ++i) {
                        rangeVector.push_back(rangeEnd[i]);
                    }
                    if(ranges.count(rangeVector) <= 0) {
                        ranges.insert(rangeVector);
                        H.push(sps(rangeStart, rangeEnd, blues, reds, topItem.intersection, topItem.union_ + 1));
                    }
                }
            }
        }
        bool corner[dims] = {false};
        int i=dims;
        while(i>=0) {
            if(i == dims) {
                // Success case process and update location of i
                // Generate range and call findSkyline
                double skylineRangeStart[dims], skylineRangeEnd[dims];
                for(int j=0; j<dims; j++) {
                    // 0 - start side, 1 is end side
                    if(corner[j]) {
                        // Swapping items here for Ian's method to know that its in the wrong direction
                        skylineRangeStart[j] = topItem.end[j];
                        if(boundaryPoints[1][j] != nullptr) {
                            skylineRangeEnd[j] = boundaryPoints[1][j]->coordinates[j];
                        } else {
                            skylineRangeEnd[j] = 1000.;// boundary[1][j];
                        }
                    } else {
                        skylineRangeStart[j] = topItem.start[j];
                        if(boundaryPoints[0][j] != nullptr) {
                            skylineRangeEnd[j] = boundaryPoints[0][j]->coordinates[j];
                        } else {
                            skylineRangeEnd[j] = -1000.;// boundary[1][j];
                        }
                        //skylineRangeEnd[j] = boundaryPoints[1][j]->coordinates[j];
                    }
                }
                std::vector<Point> *skylineVec;
                skylineVec = rt.findSkyline(skylineRangeStart, skylineRangeEnd);
                if(skylineVec && skylineVec->size() == 0) {
                    range skylineRange;
                    double skylineStart[dims], skylineEnd[dims];
                    for(int i_=0; i_ < dims; ++i_) {
                        skylineStart[i_] = min(skylineRangeStart[i_], skylineRangeEnd[i_]);
                        skylineEnd[i_] = max(skylineRangeStart[i_], skylineRangeEnd[i_]);
                    }
                    skylineRange.start = skylineStart;
                    skylineRange.end = skylineEnd;
                    skylineRange.inclusiveYMax = true;
                    skylineRange.inclusiveYMin = true;
                    skylineRange.inclusiveXMax = true;
                    skylineRange.inclusiveXMin = true;
                    vector<Point> *pointsInSkyline = rt.findPointsInRange(skylineRange);
                    if(pointsInSkyline->size() > 0) {
                        bool error = true;
                        if(pointsInSkyline->size() == 1) {
                            error = false;
                            for(int i_=0; i_ < dims; ++i_) {
                                if(corner[i_]) {
                                    error = error || (pointsInSkyline->at(0).coordinates[i_] != skylineStart[i_]);
                                } else {
                                    error = error || (pointsInSkyline->at(0).coordinates[i_] != skylineEnd[i_]);
                                }
                             }
                            // (pointsInSkyline->at(0).coordinates[0] == skylineStart[0] || pointsInSkyline->at(0).coordinates[0] == skylineEnd[0]) &&
                            // (pointsInSkyline->at(0).coordinates[1] == skylineStart[1] || pointsInSkyline->at(0).coordinates[1] == skylineEnd[1])
                            // Corner point - very loost check
                            // Should make use of corner variable to check exact corner
                        } 
                        if(error) {
                            cout << "Error occured" << endl;
                        }
                    }
                    delete pointsInSkyline;
                }
                for(int j=0; skylineVec && j < skylineVec->size(); ++j) {
                    double rangeStart[dims], rangeEnd[dims];
                    Point point = skylineVec->at(j);
                    if(point.inRange) {
                        continue;
                    }
                    int blues = topItem.blues, reds = topItem.reds;
                    if(point.isBlue) {
                        blues++;
                    } else {
                        reds++;
                    }
                    for(int k=0; k<dims; ++k) {
                        if(corner[k]) {
                            rangeStart[k] =  topItem.start[k];
                            rangeEnd[k] = point.coordinates[k];
                        } else {
                            rangeStart[k] = point.coordinates[k];
                            rangeEnd[k] = topItem.end[k];
                        }
                    }
                    vector<double> rangeVector;
                    for(int i=0; i<dims; ++i) {
                        rangeVector.push_back(rangeStart[i]);
                    }
                    for(int i=0; i<dims; ++i) {
                        rangeVector.push_back(rangeEnd[i]);
                    }
                    if(ranges.count(rangeVector) <= 0) {
                        ranges.insert(rangeVector);
                        H.push(sps(rangeStart, rangeEnd, blues, reds, topItem.intersection, topItem.union_ + 1));
                    }
                }
                delete skylineVec;
                i--;
                while(i >= 0 && corner[i]) i--;
            } else {
                corner[i] = (! corner[i]);
                i++;
            }
        }

        // Freeing allocated boundary points
        for(int i=0; i<1; ++i) {
            for(int j=0; j<dims; ++j) {
                delete boundaryPoints[i][j];
            }
        }

        // Add code for removal of point from boundary
        double start_[dims], end_[dims];
        bool completed[2][dims] = {false};
        for(int dim=0; dim < dims; ++dim) {
            start_[dim] = topItem.start[dim];
            end_[dim] = topItem.end[dim];
        }
        for(int dim=0; dim < dims; ++dim) {
            for(int direction=0; direction < 2; ++direction) {
                if(completed[direction][dim]) {
                    continue;
                }
                // Check if the point lies in multiple dimensions
                if(direction == 0) {
                    // Set the start dimension on both ends
                    start_[dim] = topItem.start[dim];
                    end_[dim] = topItem.start[dim];
                }
                else {
                    // Set the end dimension on both ends
                    start_[dim] = topItem.end[dim];
                    end_[dim] = topItem.end[dim];
                }
                range edgeRange;
                edgeRange.start = start_;
                edgeRange.end = end_;
                edgeRange.inclusiveXMax = true;
                edgeRange.inclusiveXMin = true;
                edgeRange.inclusiveYMax = true;
                edgeRange.inclusiveYMin = true;
                vector<Point>* pointsRange = rt.findPointsInRange(edgeRange);
                if(pointsRange->size() != 1) {
                    // Something is wrong here!
                    cout << "Points in range " << pointsRange->size() << endl;
                } else {
                    int blues = topItem.blues, reds = topItem.reds;
                    if(pointsRange->at(0).isBlue) {
                        blues--;
                    } else {
                        reds--;
                    }
                    // Reset the current dimension
                    start_[dim] = topItem.start[dim];
                    end_[dim] = topItem.end[dim];

                    if(!pointsRange->at(0).inRange) {
                        // The point that is trying to be delted was recently added
                        delete pointsRange;
                        continue;
                    }

                    double newRangeStart[dims], newRangeEnd[dims];

                    for(int i=0; i<=1; ++i) {
                        for(int j=0; j<dims; ++j) {
                            // TODO: Write code that shrinks the range along all dimemnsions that match and 
                            // update the new points into heap H
                            // Remove a point from dim dimension. The point is in range.
                            if(i == 0) {
                                // Should be equal to the start point's coordinates
                                if(pointsRange->at(0).coordinates[j] == topItem.start[j]) {
                                    completed[i][j] = true;
                                    // Push a point in from start direction
                                    start_[j] = topItem.start[j];
                                    end_[j] = topItem.end[j];
                                    Point* p = rt.findNearest(start_, end_, j);
                                    newRangeStart[j] = p->coordinates[j];
                                    delete p;
                                    start_[j] = topItem.start[j];
                                    end_[j] = topItem.end[j];
                                } else {
                                    newRangeStart[j] = topItem.start[j];
                                }
                            }
                            else {
                                // Should be equal to the end point's coordinates
                                if(pointsRange->at(0).coordinates[j] == topItem.end[j]) {
                                    completed[i][j] = true;

                                    // Push a point in from start direction
                                    start_[j] = topItem.end[j];
                                    end_[j] = topItem.start[j];
                                    Point* p = rt.findNearest(start_, end_, j);
                                    newRangeEnd[j] = p->coordinates[j];
                                    delete p;
                                    start_[j] = topItem.start[j];
                                    end_[j] = topItem.end[j];
                                } else {
                                    newRangeEnd[j] = topItem.end[j];
                                } 
                            }
                        }
                    }

                    vector<double> rangeVector;
                    for(int i=0; i<dims; ++i) {
                        rangeVector.push_back(newRangeStart[i]);
                    }
                    for(int i=0; i<dims; ++i) {
                        rangeVector.push_back(newRangeEnd[i]);
                    }
                    if(ranges.count(rangeVector) <= 0) {
                        //  cout << "shrinking range" << endl;
                        ranges.insert(rangeVector);
                        H.push(sps(newRangeStart, newRangeEnd, blues, reds, topItem.intersection - 1, topItem.union_));
                    }
                    // Set completed flag for all those entires which lie on the point
                    delete pointsRange;
                }
            }
        }
    }
    return 0;
}
