#include<iostream>
#include<boost/heap/pairing_heap.hpp>
#include<cmath>
#include<vector>
#include"range_tree.h"
//#include<unordered_set>
#include<chrono>
#include <fstream>

const int dims = 2;
int blueWeight = 1, redWeight = 1, diff = 10;
double percentageDiff = 0.05;
ofstream exploreTree("tree_exploration.txt");
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

struct sps {
    double start[dims], end[dims], distApproxToFair;
    int blues, reds, intersection, union_, diffVal;

    sps(double s[dims], double e[dims], int blueCount, int redCount, int itemsIntersection, int itemsUnion, int parentSimilarity) {
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
        diffVal = abs(blueWeight * blues - redWeight * reds);
        distApproxToFair = 1 - distanceApproxNew();
        if(blues == 1080 and reds > 976) {
            cout << "GDB" << endl;
        }
        if(intersection == 1900) {
            if(blues == 1080) {
                cout << "Reached the rectangle" << distanceApprox() << " " << (1- distanceApproxNew()) <<endl;
                cout << "Diff val is " << diffVal << endl;
            }
        }
        if(distApproxToFair < parentSimilarity) {
            cout << "Bug found" << endl;
            distanceApprox();
            distanceApproxNew();
        }
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

    double distanceApproxNew() {
        int delta = blueWeight * blues - redWeight * reds;
        if(abs(delta) < diff) {
            // Already fair - return a distance of 0
            return intersection/union_;
        }
        if(delta > diff) {
            // Case when dominanat color is blue
            if(blueWeight > redWeight) {
                // Save computation when the weight for blue is greater than that of red
                double maxSimilarity = -1.;
                for(int i=0 ;i <= ceil(double(delta - diff)/blueWeight); ++i) {
                    double C_r = max(ceil(double(delta-diff-i*blueWeight)/redWeight), 0.);
                    maxSimilarity = max(maxSimilarity, double(intersection - i)/double(union_ + C_r));
                }
                return maxSimilarity;
            } else {
                // Save computation when the weight for red is greater than that of blue
                double maxSimilarity = -1.;
                for(int i=0 ;i <= ceil(double(delta - diff)/redWeight); ++i) {
                    double C_b = max(ceil(double(delta-diff-i*redWeight)/blueWeight), 0.);
                    maxSimilarity = max(maxSimilarity, double(intersection - C_b)/double(union_ + i));
                }
                return maxSimilarity;
            }
        } else {
            // Case when dominanat color is red
            if(blueWeight > redWeight) {
                // Save computation when the weight for blue is greater than that of red
                double maxSimilarity = -1.;
                for(int i=0 ;i <= ceil(double(-delta - diff)/blueWeight); ++i) {
                    double C_r = max(ceil(double(-delta - diff - i*blueWeight)/redWeight), 0.);
                    maxSimilarity = max(maxSimilarity, double(intersection - C_r)/double(union_ + i));
                }
                return maxSimilarity;
            } else {
                // Save computation when the weight for red is greater than that of blue
                double maxSimilarity = -1.;
                for(int i=0; i <= ceil(double(-delta - diff)/redWeight); ++i) {
                    double C_b = max(ceil(double(-delta - diff - i*redWeight)/blueWeight), 0.);
                    maxSimilarity = max(maxSimilarity, double(intersection - i)/double(union_ + C_b));
                }
                return maxSimilarity;
            }
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

    /*
    bool operator<(sps const & sps2) const {
        if(distanceApprox() == sps2.distanceApprox()) {
            return numberPointsAway() > sps2.numberPointsAway();
        }
        return distanceApprox() > sps2.distanceApprox();
    }
    */

    bool operator<(sps const & sps2) const {
        return distApproxToFair > sps2.distApproxToFair;
    }

};

void fileWriter(sps node, double startParent[dims],double endParent[dims]){
    for(int i=0; i< dims; ++i) {
        exploreTree << startParent[i] << " ";
    }
    for(int i=0; i< dims; ++i) {
        exploreTree << endParent[i] << " ";
    }
    for(int i=0; i< dims; ++i) {
        exploreTree << node.start[i] << " ";
    }
    for(int i=0; i< dims; ++i) {
        exploreTree << node.end[i] << " ";
    }
    exploreTree << node.distApproxToFair << " ";
    exploreTree << node.blues << " ";
    exploreTree << node.reds << " ";
    exploreTree << node.intersection << " ";
    exploreTree << node.union_ << endl;
}

int main(int argc, char *argv[]) {
    double inputRange[2][dims];
    //unordered_set<vector<double>> uniqueRanges;
    /*
    inputRange[0][0] = -0.414781;
    inputRange[0][1] = 51.5621;
    inputRange[1][0] = -0.270792;
    inputRange[1][1] = 51.6133;
    */
    inputRange[0][0] = stod(argv[1]);
    inputRange[0][1] = stod(argv[2]);
    inputRange[1][0] = stod(argv[3]);
    inputRange[1][1] = stod(argv[4]);
    auto dummyTime = chrono::high_resolution_clock::now();
    auto skylineDuration = chrono::duration_cast<chrono::microseconds>(dummyTime - dummyTime);
    auto nearestDuration = chrono::duration_cast<chrono::microseconds>(dummyTime - dummyTime);
    auto setDuration = chrono::duration_cast<chrono::microseconds>(dummyTime - dummyTime);

    cout << "[" << inputRange[0][0] << "," << inputRange[0][1] << "], [" << inputRange[1][0] << "," << inputRange[1][1] << "]" << endl;

    double threshold = 0.1; // Distance threshold

    double inputRangeCopy[2][dims];

    //RangeTree rt(argv[1], inputRange[0], inputRange[1]);
    //RangeTree rt("data/10k_new.csv", inputRange[0], inputRange[1]);
    RangeTree rt("temp_test/data.csv", inputRange[0], inputRange[1]);

    // Fix the input range
    auto startTime = chrono::high_resolution_clock::now();
    double *fixInputStart = new double[dims];
    double *fixInputEnd = new double[dims];
    for(int i=0; i<dims; ++i) {
        fixInputStart[i] = inputRange[0][i];
        fixInputEnd[i] = inputRange[1][i];
        inputRangeCopy[0][i] = inputRange[0][i];
        inputRangeCopy[1][i] = inputRange[1][i];
    }
    for(int i=0; i<dims; ++i) {
        for(int j=0; j<=1; ++j) {
            // Fix the range for dimension j
            if(j == 0) {
                fixInputEnd[i] = inputRange[0][i];
            } else {
                fixInputStart[i] = inputRange[1][i];
            }
            // fixInput[(j+1)%2][i] = inputRange[j][i];
            range fixEndPointsRange(fixInputStart, fixInputEnd, true, true, true, true);
            /*
            fixEndPointsRange.start = fixInputStart;
            fixEndPointsRange.end = fixInputEnd;
            fixEndPointsRange.inclusiveXMax = true;
            fixEndPointsRange.inclusiveXMin = true;
            fixEndPointsRange.inclusiveYMax = true;
            fixEndPointsRange.inclusiveYMin = true;
            */
            vector<Point>* pointsFixRange = rt.findPointsInRange(fixEndPointsRange);
            if(pointsFixRange->size() == 0) {
                if(j == 0){
                    fixInputStart[i] = inputRange[0][i];
                    fixInputEnd[i] = inputRange[1][i];
                } else {
                    fixInputStart[i] = inputRange[1][i];
                    fixInputEnd[i] = inputRange[0][i];
                }
                double *fixInputSC = new double[dims];
                double *fixInputEC = new double[dims];
                for(int i=0; i<dims; ++i) {
                    fixInputSC[i] = fixInputStart[i];
                    fixInputEC[i] = fixInputEnd[i];
                }
                Point* p = rt.findNearest(fixInputSC, fixInputEC, i);
                if(p) {
                    inputRange[j][i] = p->coordinates[i];
                }
                fixInputStart[i] = inputRange[0][i];
                fixInputEnd[i] = inputRange[1][i];
                delete p;
            }
            delete pointsFixRange;
            fixInputStart[i] = inputRange[0][i];
            fixInputEnd[i] = inputRange[1][i];
        }
    }
    delete fixInputStart;
    delete fixInputEnd;

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
    range initialRange(inputRange[0], inputRange[1], true, true, true, true);
    /*
    initialRange.start = inputRange[0];
    initialRange.end = inputRange[1];
    initialRange.inclusiveYMax = true;
    initialRange.inclusiveYMin = true;
    initialRange.inclusiveXMax = true;
    initialRange.inclusiveXMin = true;
    */
    vector<Point> *pointsInRange = rt.findPointsInRange(initialRange);
    itemsIntersection = itemsUnion = pointsInRange->size();
    range initialRange1(inputRangeCopy[0], inputRangeCopy[1], true, true, true, true);
    /*
    initialRange.start = inputRangeCopy[0];
    initialRange.end = inputRangeCopy[1];
    */
    vector<Point> *pointsInRangeCopy = rt.findPointsInRange(initialRange1);
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
    //diff = ceil( itemsUnion * 6514. * percentageDiff/ 10000.);
    // diff = ceil( 0.025 * itemsUnion);
    diff = 99;
    cout << "Diff value chosen as " << diff << endl;
    cout << "Number of points in range " << pointsInRange->size() << endl;
    cout << "Input range " << blueCount << " blues and " << redCount << "reds with a total of " << itemsUnion << "points." << endl;
    delete pointsInRange;
    struct sps first(inputRange[0], inputRange[1], blueCount, redCount, itemsIntersection, itemsUnion, 0.);
    cout << "Initial disparity " << first.diffVal << endl;
    unordered_set<vector<double>, VectorHash> ranges;
    vector<double> rangeFirst;
    for(int i=0; i<dims; ++i) {
        rangeFirst.push_back(first.start[i]);
    }
    for(int i=0; i<dims; ++i) {
        rangeFirst.push_back(first.end[i]);
    }
    ranges.insert(rangeFirst);
    int minDiff = first.diffVal;
    H.push(first);
    int dummyCounter=0;
    while(!H.empty()) {

        struct sps topItem = H.top();
        if(minDiff > topItem.diffVal) {
            minDiff = topItem.diffVal;
            cout << "Found a better range with disparity of " << minDiff << endl;
        }
        if(topItem.distanceApprox() > threshold) {
            cout << "Cannot find fair range within threshold of " << threshold << endl;
            auto endTime = chrono::high_resolution_clock::now();
            auto duration = chrono::duration_cast<chrono::microseconds>(endTime - startTime); 
            cout << "Total time is : " << duration.count()/1000000. << " seconds" << endl;
            exploreTree.close();
            return 0;
        }
        if(topItem.blues == 1080 && topItem.reds == 978) {
            cout << "GDB" << endl;
        }
        if(++dummyCounter % 50000 == 0) {
            cout << "Disparity : " << topItem.intersection << ":" << topItem.union_ << " ; Distance approx : " << topItem.distApproxToFair << endl << "Skyline duration : " << skylineDuration.count() << "microsectonds" << endl;
        }
        H.pop();
        int curDiffVal = topItem.blues * blueWeight - topItem.reds * redWeight;
        if(abs(curDiffVal) <= diff) {
            cout << "Range definition : " << endl;
            for(int i=0; i < dims; ++i) {
                cout << topItem.start[i] << "-" << topItem.end[i] << endl;
            }
            cout << "Disparity : " << topItem.intersection << ":" << topItem.union_ << ":" << double(topItem.intersection)/double(topItem.union_) << endl;
            auto endTime = chrono::high_resolution_clock::now();
            auto duration = chrono::duration_cast<chrono::microseconds>(endTime - startTime); 
            cout << "Total time is : " << duration.count()/1000000. << " seconds" << endl;
            exploreTree.close();
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
                    if(boundaryPoints[boundaryType][dim]->isBlue) {
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
                        fileWriter(sps(rangeStart, rangeEnd, blues, reds, topItem.intersection, topItem.union_ + 1, topItem.distApproxToFair), topItem.start, topItem.end);
                        H.push(sps(rangeStart, rangeEnd, blues, reds, topItem.intersection, topItem.union_ + 1, topItem.distApproxToFair));
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
                double *skylineRangeStartCopy = new double [dims], *skylineRangeEndCopy = new double[dims];
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
                    skylineRangeEndCopy[j] = skylineRangeStart[j];
                    skylineRangeStartCopy[j] = skylineRangeEnd[j];
                }
                std::vector<Point> *skylineVec;
                auto skylineStartTime = chrono::high_resolution_clock::now();
                skylineVec = rt.findSkyline(skylineRangeStartCopy, skylineRangeEndCopy);
                auto skylineEndTime = chrono::high_resolution_clock::now();
                skylineDuration += chrono::duration_cast<chrono::microseconds>(skylineEndTime - skylineStartTime); 
                if(skylineVec && skylineVec->size() == 0) {
                    double skylineStart[dims], skylineEnd[dims];
                    for(int i_=0; i_ < dims; ++i_) {
                        skylineStart[i_] = min(skylineRangeStart[i_], skylineRangeEnd[i_]);
                        skylineEnd[i_] = max(skylineRangeStart[i_], skylineRangeEnd[i_]);
                    }
                    range skylineRange(skylineStart, skylineEnd, true, true, true, true);
                    /*
                    skylineRange.start = skylineStart;
                    skylineRange.end = skylineEnd;
                    skylineRange.inclusiveYMax = true;
                    skylineRange.inclusiveYMin = true;
                    skylineRange.inclusiveXMax = true;
                    skylineRange.inclusiveXMin = true;
                    */
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
                        fileWriter(sps(rangeStart, rangeEnd, blues, reds, topItem.intersection, topItem.union_ + 1, topItem.distApproxToFair), topItem.start, topItem.end);
                        H.push(sps(rangeStart, rangeEnd, blues, reds, topItem.intersection, topItem.union_ + 1, topItem.distApproxToFair));
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
                range edgeRange(start_, end_, true, true, true, true);
                /*
                edgeRange.start = start_;
                edgeRange.end = end_;
                edgeRange.inclusiveXMax = true;
                edgeRange.inclusiveXMin = true;
                edgeRange.inclusiveYMax = true;
                edgeRange.inclusiveYMin = true;
                */
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
                        fileWriter(sps(newRangeStart, newRangeEnd, blues, reds, topItem.intersection - 1, topItem.union_, topItem.distApproxToFair), topItem.start, topItem.end);
                        H.push(sps(newRangeStart, newRangeEnd, blues, reds, topItem.intersection - 1, topItem.union_, topItem.distApproxToFair));
                    }
                    // Set completed flag for all those entires which lie on the point
                    delete pointsRange;
                }
            }
        }
    }
    auto endTime = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::microseconds>(endTime - startTime); 
    //cout << "Total time is : " << duration.count() << "microseconds" << endl;
    cout << "Total time is : " << duration.count()/1000000. << " seconds" << endl;
    exploreTree.close();
    return 0;
}
