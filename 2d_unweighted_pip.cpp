#include <string>
#include <vector>
#include <chrono>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>

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

struct ExpandingVerticalLevel {
    Point furthest;
    Point closest;
    vector<Point> *intermediates;
    bool finished;
    int improvementUnderClosest;
    int pointsUnderClosest;
};

struct ShrinkingVerticalLevel {
    Point furthest;
    Point greatest;
    int pointsRemovedAtFurthest;
    int improvementAtFurthest;
};

struct ExpandingHorizontalPoint {
    Point expandedPoint;
    int improvement;
    int pointsInCenter;
    int improvementUnderTop;
    int improvementAboveBottom;
    vector<Point> *upperSubrange;
    vector<Point> *lowerSubrange;
};

struct ShrinkingHorizontalPoint {
    Point shrunkPoint;
    int improvement;
    int pointsRemoved;
    int pointsRemovedFromCenter;
    int pointsRemovedFromBottom;
    int pointsRemovedFromTop;
    int improveUnderTop;
    int improvementAboveBottom;
    int pointsRemovedFromStartingRange;
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

vector<ExpandingVerticalLevel>  *expandVerticalLevel(vector<Point> *points, Point boundaryPoint, bool needsBlues, int threshold) {
    int pointsUnderClosest = 0;
    int improvement = 0;
    int maxImprovement = 0;
    int totalPoints = 0;
    ExpandingVerticalLevel first;
    first.furthest = boundaryPoint;
    first.closest = boundaryPoint;
    first.improvementUnderClosest = 0;
    first.pointsUnderClosest = 0;
    first.intermediates = new vector<Point>();
    first.finished = true;
    Point previousPoint = boundaryPoint;
    vector<ExpandingVerticalLevel> *levels = new vector<ExpandingVerticalLevel>();
    levels->push_back(first);
    vector<Point> *intermediates = new vector<Point>();
    Point lastPoint = boundaryPoint;
    for (Point p : *points) {
        lastPoint = p;
        totalPoints++;
        if (p.isBlue == needsBlues) {
            improvement++;
        } else {
            improvement--;
        }
        if (improvement > maxImprovement) {
            ExpandingVerticalLevel level;
            level.furthest = p;
            level.closest = previousPoint;
            level.intermediates = intermediates;
            level.improvementUnderClosest = maxImprovement;
            level.pointsUnderClosest = pointsUnderClosest;
            level.finished = true;
            maxImprovement = improvement;
            pointsUnderClosest = totalPoints;
            levels->push_back(level);
            intermediates= new vector<Point>();
            if (threshold == maxImprovement) {
                ExpandingVerticalLevel level;
                level.furthest = p;
                level.closest = p;
                level.intermediates = new vector<Point>();
                level.improvementUnderClosest = maxImprovement;
                level.pointsUnderClosest = pointsUnderClosest;
                level.finished = true;
                levels->push_back(level);
                return levels;
            }
        } else {
            intermediates->push_back(p);
        }
    }
    ExpandingVerticalLevel last;
    last.furthest = lastPoint;
    last.closest = previousPoint;
    last.intermediates = intermediates;
    last.improvementUnderClosest = maxImprovement;
    last.pointsUnderClosest = pointsUnderClosest;
    last.finished = false;
    levels->push_back(last);
    return levels;
 }

vector<ShrinkingVerticalLevel> *shrinkVerticalLevel(vector<Point> *points, Point startingPoint, bool needsRemoveBlues, int threshold) {
    int maxImprovement = 0;
    int currentImprovement = 0;
    int pointsRemoved = 0;
    Point lastPoint = startingPoint;
    bool previouslyFinished = false;
    vector<ShrinkingVerticalLevel> *levels = new vector<ShrinkingVerticalLevel>();
    
    for (Point p : *points) {
        if (previouslyFinished) {
            levels->at(levels->size()-1).greatest = p;
            if (threshold == maxImprovement) {
                return levels;
            }
            previouslyFinished = false;
        }
        pointsRemoved++;
        if (p.isBlue == needsRemoveBlues) {
            currentImprovement++;
        } else {
            currentImprovement--;
        }
        if (currentImprovement > maxImprovement) {
            ShrinkingVerticalLevel level;
            level.improvementAtFurthest = currentImprovement;
            level.furthest = lastPoint;
            level.pointsRemovedAtFurthest = pointsRemoved;
            levels->push_back(level);
            previouslyFinished = true;
            maxImprovement = currentImprovement;
            lastPoint = p;
        }
    }
    return levels; 
}

vector<Point> *replaceWithCopy(vector<Point> *original) {
    vector<Point> *c = new vector<Point>();
    copy(original->begin(), original->end(), back_inserter(c));
    return c;
}

vector<ExpandingHorizontalPoint> *expandHorizontallyBothExpanding(vector<Point> *points, Point boundaryPoint, ExpandingVerticalLevel upperLevel, ExpandingVerticalLevel lowerLevel, bool needsBlues, int threshold) {
    int currentImprovement = 0;
    int pointsInCenter = 0;
    int maxImprovement = 0;
    int improvementUnderTop = 0;
    int improvementAboveBottom = 0;

    vector<Point> *upperSubrange = new vector<Point>();
    vector<Point> *lowerSubrange = new vector<Point>();

    vector<ExpandingHorizontalPoint> *horizontalPoints = new vector<ExpandingHorizontalPoint>();

    ExpandingHorizontalPoint start;
    start.improvement = 0;
    start.pointsInCenter = 0;
    start.expandedPoint = boundaryPoint;
    start.upperSubrange = new vector<Point>();
    start.lowerSubrange = new vector<Point>();
    start.improvementAboveBottom = 0;
    start.improvementUnderTop = 0;
    horizontalPoints->push_back(start);
    if (threshold - upperLevel.improvementUnderClosest - lowerLevel.improvementUnderClosest == 0) {
        return;
    }

    for (Point p : *points) {
        if (p.coordinates[1] <= upperLevel.closest.coordinates[1] && p.coordinates[1] >= lowerLevel.closest.coordinates[1]) {
            pointsInCenter++;
            improvementAboveBottom++;
            improvementUnderTop++;
        } else if (p.coordinates[1] > upperLevel.closest.coordinates[1]) {
            upperSubrange->push_back(p);
            improvementAboveBottom++;
        } else {
            lowerSubrange->push_back(p);
            improvementUnderTop++;
        }
        if (p.isBlue == needsBlues) {
            currentImprovement++;
        } else {
            currentImprovement--;
        }
        if (currentImprovement > maxImprovement) {
            maxImprovement = currentImprovement;
            ExpandingHorizontalPoint horizontalPoint;
            horizontalPoint.improvement = maxImprovement;
            horizontalPoint.pointsInCenter = pointsInCenter;
            horizontalPoint.expandedPoint = p;
            horizontalPoint.upperSubrange = replaceWithCopy(upperSubrange);
            horizontalPoint.lowerSubrange = replaceWithCopy(lowerSubrange);
            horizontalPoint.improvementUnderTop = improvementUnderTop;
            horizontalPoint.improvementAboveBottom = improvementAboveBottom;
            horizontalPoints->push_back(horizontalPoint);
            if (currentImprovement >= threshold - upperLevel.improvementUnderClosest - lowerLevel.improvementUnderClosest) {
                return horizontalPoints;
            }
        }
    }
    return horizontalPoints;
}

vector<ExpandingHorizontalPoint> *expandHorizontallyTopExpanding(vector<Point> *points, Point boundaryPoint, ExpandingVerticalLevel upperLevel, ShrinkingVerticalLevel lowerLevel, bool needsBlues, int threshold) {
    int currentImprovement = 0;
    int pointsInCenter = 0;
    int maxImprovement = 0;
    int improvementUnderTop = 0;
    int improvementAboveBottom = 0;

    vector<Point> *upperSubrange = new vector<Point>();
    vector<Point> *lowerSubrange = new vector<Point>();

    vector<ExpandingHorizontalPoint> *horizontalPoints = new vector<ExpandingHorizontalPoint>();

    ExpandingHorizontalPoint start;
    start.improvement = 0;
    start.pointsInCenter = 0;
    start.expandedPoint = boundaryPoint;
    start.upperSubrange = new vector<Point>();
    start.lowerSubrange = new vector<Point>();
    start.improvementAboveBottom = 0;
    start.improvementUnderTop = 0;
    horizontalPoints->push_back(start);
    if (threshold - upperLevel.improvementUnderClosest - lowerLevel.improvementAtFurthest == 0) {
        return;
    }

    for (Point p : *points) {
        if (p.coordinates[1] <= upperLevel.closest.coordinates[1] && p.coordinates[1] >= lowerLevel.furthest.coordinates[1]) {
            pointsInCenter++;
            improvementAboveBottom++;
            improvementUnderTop++;
        } else if (p.coordinates[1] > upperLevel.closest.coordinates[1]) {
            upperSubrange->push_back(p);
            improvementAboveBottom++;
        } else {
            lowerSubrange->push_back(p);
            improvementUnderTop++;
        }
        if (p.isBlue == needsBlues) {
            currentImprovement++;
        } else {
            currentImprovement--;
        }
        if (currentImprovement > maxImprovement) {
            maxImprovement = currentImprovement;
            ExpandingHorizontalPoint horizontalPoint;
            horizontalPoint.improvement = maxImprovement;
            horizontalPoint.pointsInCenter = pointsInCenter;
            horizontalPoint.expandedPoint = p;
            horizontalPoint.upperSubrange = replaceWithCopy(upperSubrange);
            horizontalPoint.lowerSubrange = replaceWithCopy(lowerSubrange);
            horizontalPoint.improvementUnderTop = improvementUnderTop;
            horizontalPoint.improvementAboveBottom = improvementAboveBottom;
            horizontalPoints->push_back(horizontalPoint);
            if (currentImprovement >= threshold - upperLevel.improvementUnderClosest - lowerLevel.improvementAtFurthest) {
                return horizontalPoints;
            }
        }
    }
    return horizontalPoints;
}

vector<ExpandingHorizontalPoint> *expandHorizontallyBottomExpanding(vector<Point> *points, Point boundaryPoint, ExpandingVerticalLevel lowerLevel, ShrinkingVerticalLevel upperLevel, bool needsBlues, int threshold) {
    int currentImprovement = 0;
    int pointsInCenter = 0;
    int maxImprovement = 0;
    int improvementUnderTop = 0;
    int improvementAboveBottom = 0;

    vector<Point> *upperSubrange = new vector<Point>();
    vector<Point> *lowerSubrange = new vector<Point>();

    vector<ExpandingHorizontalPoint> *horizontalPoints = new vector<ExpandingHorizontalPoint>();

    ExpandingHorizontalPoint start;
    start.improvement = 0;
    start.pointsInCenter = 0;
    start.expandedPoint = boundaryPoint;
    start.upperSubrange = new vector<Point>();
    start.lowerSubrange = new vector<Point>();
    start.improvementAboveBottom = 0;
    start.improvementUnderTop = 0;
    horizontalPoints->push_back(start);
    if (threshold - upperLevel.improvementAtFurthest - lowerLevel.improvementUnderClosest == 0) {
        return;
    }

    for (Point p : *points) {
        if (p.coordinates[1] <= upperLevel.furthest.coordinates[1] && p.coordinates[1] >= lowerLevel.closest.coordinates[1]) {
            pointsInCenter++;
            improvementAboveBottom++;
            improvementUnderTop++;
        } else if (p.coordinates[1] > upperLevel.furthest.coordinates[1]) {
            upperSubrange->push_back(p);
            improvementAboveBottom++;
        } else {
            lowerSubrange->push_back(p);
            improvementUnderTop++;
        }
        if (p.isBlue == needsBlues) {
            currentImprovement++;
        } else {
            currentImprovement--;
        }
        if (currentImprovement > maxImprovement) {
            maxImprovement = currentImprovement;
            ExpandingHorizontalPoint horizontalPoint;
            horizontalPoint.improvement = maxImprovement;
            horizontalPoint.pointsInCenter = pointsInCenter;
            horizontalPoint.expandedPoint = p;
            horizontalPoint.upperSubrange = replaceWithCopy(upperSubrange);
            horizontalPoint.lowerSubrange = replaceWithCopy(lowerSubrange);
            horizontalPoint.improvementUnderTop = improvementUnderTop;
            horizontalPoint.improvementAboveBottom = improvementAboveBottom;
            horizontalPoints->push_back(horizontalPoint);
            if (currentImprovement >= threshold - upperLevel.improvementAtFurthest - lowerLevel.improvementUnderClosest) {
                return horizontalPoints;
            }
        }
    }
    return horizontalPoints;
}

vector<ShrinkingHorizontalPoint> *shrinkHorizontallyBothExpanding(vector<Point> *points, ExpandingVerticalLevel upperLevel, ExpandingVerticalLevel lowerLevel, bool needsRemoveBlues, int threshold, double startingMinimumY, double startingMaximumY) {
    int maxImprovement = 0;
    int currentImprovement = 0;
    int pointsRemoved = 0;
    int pointsRemovedFromCenter = 0;
    int pointsRemovedFromTop = 0;
    int pointsRemovedFromBottom = 0;
    int improvementUnderTop = 0;
    int improvementAboveBottom = 0;
    int pointsRemovedFromStartingRange = 0;
    bool lastPointImproved = false;
    vector<ShrinkingHorizontalPoint> *horizontalPoints = new vector<ShrinkingHorizontalPoint>();
    ShrinkingHorizontalPoint start;
    start.improvement = 0;
    start.pointsRemoved = 0;
    start.pointsRemovedFromCenter = 0;
    start.pointsRemovedFromBottom = 0;
    start.pointsRemovedFromStartingRange = 0;
    start.pointsRemovedFromTop = 0;
    start.improveUnderTop = 0;
    start.improvementAboveBottom = 0;
    start.shrunkPoint = points->at(0);
    start.pointsRemovedFromStartingRange = 0;
    horizontalPoints->push_back(start);

    for (Point p : *points) {
        if (lastPointImproved == true) {
            ShrinkingHorizontalPoint horizontalPoint;
            horizontalPoint.pointsRemoved = pointsRemoved;
            horizontalPoint.improvement = maxImprovement;
            horizontalPoint.shrunkPoint = p;
            horizontalPoint.pointsRemovedFromCenter = pointsRemovedFromCenter;
            horizontalPoint.pointsRemovedFromBottom = pointsRemovedFromBottom;
            horizontalPoint.pointsRemovedFromTop = pointsRemovedFromTop;
            horizontalPoint.improveUnderTop = improvementUnderTop;
            horizontalPoint.improvementAboveBottom = improvementAboveBottom;
            horizontalPoints->push_back(horizontalPoint);
            lastPointImproved = false;
            if (maxImprovement == threshold - upperLevel.improvementUnderClosest - lowerLevel.improvementUnderClosest) {
                return horizontalPoints;
            }
        }
        //increase points after saving improving point so that we can use the number of points previously removed
        pointsRemoved++;
        if (p.coordinates[1] <= upperLevel.closest.coordinates[1] && p.coordinates[1] >= lowerLevel.closest.coordinates[1]) {
            pointsRemovedFromCenter++;
        } else if (p.coordinates[1] <= upperLevel.closest.coordinates[1]) {
            pointsRemovedFromBottom++;
        } else {
            pointsRemovedFromTop++;
        }
        if (p.coordinates[1] >= startingMinimumY && p.coordinates[1] <= startingMaximumY) {
            pointsRemovedFromStartingRange++;
        }
        if (p.isBlue == needsRemoveBlues) {
            currentImprovement++;
            if (p.coordinates[1] <= upperLevel.closest.coordinates[1]) {
                improvementUnderTop++;
            }
            if (p.coordinates[1] >= lowerLevel.closest.coordinates[1]) {
                improvementAboveBottom++;
            }
        } else {
            currentImprovement--;
            if (p.coordinates[1] <= upperLevel.closest.coordinates[1]) {
                improvementUnderTop--;
            }
            if (p.coordinates[1] >= lowerLevel.closest.coordinates[1]) {
                improvementAboveBottom--;
            }
        }
        if (currentImprovement > maxImprovement) {
            lastPointImproved = true;
            maxImprovement = currentImprovement;
        }
    }
    return horizontalPoints;
}


vector<ShrinkingHorizontalPoint> *shrinkHorizontallyTopExpanding(vector<Point> *points, ExpandingVerticalLevel upperLevel, ShrinkingVerticalLevel lowerLevel, bool needsRemoveBlues, int threshold, double startingMinimumY, double startingMaximumY) {
    int maxImprovement = 0;
    int currentImprovement = 0;
    int pointsRemoved = 0;
    int pointsRemovedFromCenter = 0;
    int pointsRemovedFromTop = 0;
    int pointsRemovedFromBottom = 0;
    int improvementUnderTop = 0;
    int improvementAboveBottom = 0;
    int pointsRemovedFromStartingRange = 0;
    bool lastPointImproved = false;
    vector<ShrinkingHorizontalPoint> *horizontalPoints = new vector<ShrinkingHorizontalPoint>();
    ShrinkingHorizontalPoint start;
    start.improvement = 0;
    start.pointsRemoved = 0;
    start.pointsRemovedFromCenter = 0;
    start.pointsRemovedFromBottom = 0;
    start.pointsRemovedFromStartingRange = 0;
    start.pointsRemovedFromTop = 0;
    start.improveUnderTop = 0;
    start.improvementAboveBottom = 0;
    start.shrunkPoint = points->at(0);
    start.pointsRemovedFromStartingRange = 0;
    horizontalPoints->push_back(start);

    for (Point p : *points) {
        if (lastPointImproved == true) {
            ShrinkingHorizontalPoint horizontalPoint;
            horizontalPoint.pointsRemoved = pointsRemoved;
            horizontalPoint.improvement = maxImprovement;
            horizontalPoint.shrunkPoint = p;
            horizontalPoint.pointsRemovedFromCenter = pointsRemovedFromCenter;
            horizontalPoint.pointsRemovedFromBottom = pointsRemovedFromBottom;
            horizontalPoint.pointsRemovedFromTop = pointsRemovedFromTop;
            horizontalPoint.improveUnderTop = improvementUnderTop;
            horizontalPoint.improvementAboveBottom = improvementAboveBottom;
            horizontalPoints->push_back(horizontalPoint);
            lastPointImproved = false;
            if (maxImprovement == threshold - upperLevel.improvementUnderClosest - lowerLevel.improvementAtFurthest) {
                return horizontalPoints;
            }
        }
        //increase points after saving improving point so that we can use the number of points previously removed
        pointsRemoved++;
        if (p.coordinates[1] <= upperLevel.closest.coordinates[1] && p.coordinates[1] >= lowerLevel.furthest.coordinates[1]) {
            pointsRemovedFromCenter++;
        } else if (p.coordinates[1] <= upperLevel.closest.coordinates[1]) {
            pointsRemovedFromBottom++;
        } else {
            pointsRemovedFromTop++;
        }
        if (p.coordinates[1] >= startingMinimumY && p.coordinates[1] <= startingMaximumY) {
            pointsRemovedFromStartingRange++;
        }
        if (p.isBlue == needsRemoveBlues) {
            currentImprovement++;
            if (p.coordinates[1] <= upperLevel.closest.coordinates[1]) {
                improvementUnderTop++;
            }
            if (p.coordinates[1] >= lowerLevel.furthest.coordinates[1]) {
                improvementAboveBottom++;
            }
        } else {
            currentImprovement--;
            if (p.coordinates[1] <= upperLevel.closest.coordinates[1]) {
                improvementUnderTop--;
            }
            if (p.coordinates[1] >= lowerLevel.furthest.coordinates[1]) {
                improvementAboveBottom--;
            }
        }
        if (currentImprovement > maxImprovement) {
            lastPointImproved = true;
            maxImprovement = currentImprovement;
        }
    }
    return horizontalPoints;
}

vector<ShrinkingHorizontalPoint> *shrinkHorizontallyBottomExpanding(vector<Point> *points, ExpandingVerticalLevel lowerLevel, ShrinkingVerticalLevel upperLevel, bool needsRemoveBlues, int threshold, double startingMinimumY, double startingMaximumY) {
    int maxImprovement = 0;
    int currentImprovement = 0;
    int pointsRemoved = 0;
    int pointsRemovedFromCenter = 0;
    int pointsRemovedFromTop = 0;
    int pointsRemovedFromBottom = 0;
    int improvementUnderTop = 0;
    int improvementAboveBottom = 0;
    int pointsRemovedFromStartingRange = 0;
    bool lastPointImproved = false;
    vector<ShrinkingHorizontalPoint> *horizontalPoints = new vector<ShrinkingHorizontalPoint>();
    ShrinkingHorizontalPoint start;
    start.improvement = 0;
    start.pointsRemoved = 0;
    start.pointsRemovedFromCenter = 0;
    start.pointsRemovedFromBottom = 0;
    start.pointsRemovedFromStartingRange = 0;
    start.pointsRemovedFromTop = 0;
    start.improveUnderTop = 0;
    start.improvementAboveBottom = 0;
    start.shrunkPoint = points->at(0);
    start.pointsRemovedFromStartingRange = 0;
    horizontalPoints->push_back(start);

    for (Point p : *points) {
        if (lastPointImproved == true) {
            ShrinkingHorizontalPoint horizontalPoint;
            horizontalPoint.pointsRemoved = pointsRemoved;
            horizontalPoint.improvement = maxImprovement;
            horizontalPoint.shrunkPoint = p;
            horizontalPoint.pointsRemovedFromCenter = pointsRemovedFromCenter;
            horizontalPoint.pointsRemovedFromBottom = pointsRemovedFromBottom;
            horizontalPoint.pointsRemovedFromTop = pointsRemovedFromTop;
            horizontalPoint.improveUnderTop = improvementUnderTop;
            horizontalPoint.improvementAboveBottom = improvementAboveBottom;
            horizontalPoints->push_back(horizontalPoint);
            lastPointImproved = false;
            if (maxImprovement == threshold - upperLevel.improvementAtFurthest - lowerLevel.improvementUnderClosest) {
                return horizontalPoints;
            }
        }
        //increase points after saving improving point so that we can use the number of points previously removed
        pointsRemoved++;
        if (p.coordinates[1] <= upperLevel.furthest.coordinates[1] && p.coordinates[1] >= lowerLevel.closest.coordinates[1]) {
            pointsRemovedFromCenter++;
        } else if (p.coordinates[1] <= upperLevel.furthest.coordinates[1]) {
            pointsRemovedFromBottom++;
        } else {
            pointsRemovedFromTop++;
        }
        if (p.coordinates[1] >= startingMinimumY && p.coordinates[1] <= startingMaximumY) {
            pointsRemovedFromStartingRange++;
        }
        if (p.isBlue == needsRemoveBlues) {
            currentImprovement++;
            if (p.coordinates[1] <= upperLevel.furthest.coordinates[1]) {
                improvementUnderTop++;
            }
            if (p.coordinates[1] >= lowerLevel.closest.coordinates[1]) {
                improvementAboveBottom++;
            }
        } else {
            currentImprovement--;
            if (p.coordinates[1] <= upperLevel.furthest.coordinates[1]) {
                improvementUnderTop--;
            }
            if (p.coordinates[1] >= lowerLevel.closest.coordinates[1]) {
                improvementAboveBottom--;
            }
        }
        if (currentImprovement > maxImprovement) {
            lastPointImproved = true;
            maxImprovement = currentImprovement;
        }
    }
    return horizontalPoints;
}

void findY(double &Y, int &topPointsAdded, vector<Point> *subrange, vector<Point> *intermediates, int startingImprovement, int goalImprovement, bool needsBlue, bool ascending) {
    int i = 0;
    int j = 0;
    topPointsAdded = 0;
    int currentValue = startingImprovement;
    while (i < subrange->size() && j < intermediates->size()) {
        if (ascending ? subrange->at(i).coordinates[1] < intermediates->at(j).coordinates[1] : subrange->at(i).coordinates[1] > intermediates->at(j).coordinates[1]) {
            currentValue += (subrange->at(i).isBlue == needsBlue) ? 1 : -1;
            Y = subrange->at(i).coordinates[1];
            i++;
        } else {
            currentValue += (subrange->at(j).isBlue == needsBlue) ? 1 : -1;
            Y = intermediates->at(j).coordinates[1];
            j++;
        }
        if (currentValue = goalImprovement) {
            return;
        }
    }
    while (i < subrange->size()) {
        currentValue += (subrange->at(i).isBlue == needsBlue) ? 1 : -1;
        Y = subrange->at(i).coordinates[1];
        i++;
        if (currentValue = goalImprovement) {
            return;
        }
    }
    while (j < intermediates->size()) {
        currentValue += (subrange->at(j).isBlue == needsBlue) ? 1 : -1;
        Y = intermediates->at(j).coordinates[1];
        j++;
        if (currentValue = goalImprovement) {
            return;
        }
    }
    throw(runtime_error("Couldn't find max Y!"));
}

int main() {
    string filename = "./data/uniform.csv";
    double *minimum = new double[2];
    double *maximum = new double[2];

    minimum[0] = 667.341;
    minimum[1] = 313.23;
    maximum[0] = 879.009;
    maximum[1] = 653.305;
    int redInRange = 0;
    int blueInRange = 0;
    int initialPoints = redInRange + blueInRange;
    chrono::time_point<chrono::high_resolution_clock> start = chrono::high_resolution_clock::now();

    vector<Point> *points = readPoints(filename, minimum, maximum, redInRange, blueInRange);
    int epsilon = ceil(0.025*(redInRange+blueInRange));
    int threshold = abs(redInRange-blueInRange) - epsilon;
    if (threshold <= 0) {
        cout << "Initial range is fair" << endl;
        return 0;
    }

    Point topPoint;
    bool pointInRangeFound = false;
    for (Point p : *points) {
        if (p.coordinates[0] >= minimum[0] && p.coordinates[0] <= maximum[0] && p.coordinates[1] >= minimum[1] && p.coordinates[1] <= minimum[1]) {
            if (!pointInRangeFound) {
                topPoint = p;
                pointInRangeFound = true;
            } else {
                if (p.coordinates[1] > topPoint.coordinates[1]) {
                    topPoint = p;
                }
            }
        }
    }
    Point bottomPoint;
    pointInRangeFound = false;
    for (Point p : *points) {
        if (p.coordinates[0] >= minimum[0] && p.coordinates[0] <= maximum[0] && p.coordinates[1] >= minimum[1] && p.coordinates[1] <= minimum[1]) {
            if (!pointInRangeFound) {
                bottomPoint = p;
                pointInRangeFound = true;
            } else {
                if (p.coordinates[1] < topPoint.coordinates[1]) {
                    bottomPoint = p;
                }
            }
        }
    }

    Point rightPoint;
    pointInRangeFound = false;
    for (Point p : *points) {
        if (p.coordinates[0] >= minimum[0] && p.coordinates[0] <= maximum[0] && p.coordinates[1] >= minimum[1] && p.coordinates[1] <= minimum[1]) {
            if (!pointInRangeFound) {
                rightPoint = p;
                pointInRangeFound = true;
            } else {
                if (p.coordinates[0] > rightPoint.coordinates[0]) {
                    rightPoint = p;
                }
            }
        }
    }

    Point leftPoint;
    pointInRangeFound = false;
    for (Point p : *points) {
        if (p.coordinates[0] >= minimum[0] && p.coordinates[0] <= maximum[0] && p.coordinates[1] >= minimum[1] && p.coordinates[1] <= minimum[1]) {
            if (!pointInRangeFound) {
                leftPoint = p;
                pointInRangeFound = true;
            } else {
                if (p.coordinates[0] < rightPoint.coordinates[0]) {
                    leftPoint = p;
                }
            }
        }
    }


    vector<Point> *pointsUpward = new vector<Point>();
    copy_if(points->begin(), points->end(), back_inserter(*pointsUpward), [minimum, maximum](Point p)
        {return p.coordinates[0] >= minimum[0] && p.coordinates[0] <= maximum[0] && p.coordinates[1] > maximum[1];}
    );
    sort(pointsUpward->begin(), pointsUpward->end(), pointSorter(1, true));

    vector<Point> *pointsDownward = new vector<Point>();
    copy_if(points->begin(), points->end(), back_inserter(*pointsDownward), [minimum, maximum](Point p)
        {return p.coordinates[0] >= minimum[0] && p.coordinates[0] <= maximum[0] && p.coordinates[1] < minimum[1];}
    );
    sort(pointsDownward->begin(), pointsDownward->end(), pointSorter(1, false));

    vector<Point> *pointsInsideDownward = new vector<Point>();
    copy_if(points->begin(), points->end(), back_inserter(*pointsInsideDownward), [minimum, maximum](Point p)
        {return p.coordinates[0] >= minimum[0] && p.coordinates[0] <= maximum[0] && p.coordinates[1] >= minimum[1] && p.coordinates[1] <= maximum[1];}
    );
    sort(pointsInsideDownward->begin(), pointsInsideDownward->end(), pointSorter(1, false));
    vector<Point> *pointsInsideUpward = new vector<Point>();
    copy(pointsInsideDownward->begin(), pointsInsideDownward->end(), back_inserter(*pointsInsideUpward));
    reverse(pointsInsideDownward->begin(), pointsInsideDownward->end());

    vector<ExpandingVerticalLevel> *expandingUpward = expandVerticalLevel(pointsUpward, topPoint, redInRange > blueInRange, threshold);
    vector<ExpandingVerticalLevel> *expandingDownward = expandVerticalLevel(pointsDownward, bottomPoint, redInRange > blueInRange, threshold);
    vector<ShrinkingVerticalLevel> *shrinkingDownward = shrinkVerticalLevel(pointsInsideDownward, topPoint, blueInRange > redInRange, threshold);
    vector<ShrinkingVerticalLevel> *shrinkingUpward = shrinkVerticalLevel(pointsInsideUpward, bottomPoint, blueInRange > redInRange, threshold);

    double similarity;
    double *fairMinimum = new double[2];
    double *fairMaximum = new double[2];

    for (ExpandingVerticalLevel upperLevel : *expandingUpward) {
        for (ExpandingVerticalLevel lowerLevel : *expandingDownward) {
            vector<Point> *pointsRight = new vector<Point>();
            copy_if(points->begin(), points->end(), back_inserter(*pointsRight), [minimum, maximum, upperLevel, lowerLevel](Point p)
                {return p.coordinates[0] > maximum[0] && p.coordinates[1] <= upperLevel.furthest.coordinates[1] && p.coordinates[1] >= lowerLevel.furthest.coordinates[1];}
            );
            sort(pointsRight->begin(), pointsRight->end(), pointSorter(0, true));
            vector<Point> *pointsLeft = new vector<Point>();
            copy_if(points->begin(), points->end(), back_inserter(*pointsLeft), [minimum, maximum, upperLevel, lowerLevel](Point p)
                {return p.coordinates[0] < minimum[0] && p.coordinates[1] <= upperLevel.furthest.coordinates[1] && p.coordinates[1] >= lowerLevel.furthest.coordinates[1];}
            );
            sort(pointsLeft->begin(), pointsLeft->end(), pointSorter(0, false));
            vector<Point> *pointsInside = new vector<Point>();
            copy_if(points->begin(), points->end(), back_inserter(*pointsInside), [minimum, maximum, upperLevel, lowerLevel](Point p) 
                {return p.coordinates[0] <= maximum[0] && p.coordinates[0] >= minimum[0] && p.coordinates[1] <= upperLevel.furthest.coordinates[1] && p.coordinates[1] >= lowerLevel.furthest.coordinates[1];}
            );
            sort(pointsInside->begin(), pointsInside->end(), pointSorter(0, false));
            vector<Point> *pointsInsideReversed = new vector<Point>();
            copy(pointsInside->begin(), pointsInside->end(), back_inserter(*pointsInsideReversed));
            reverse(pointsInsideReversed->begin(), pointsInsideReversed->end());

            vector<ExpandingHorizontalPoint> *expandRight = expandHorizontallyBothExpanding(pointsRight, rightPoint, upperLevel, lowerLevel, redInRange > blueInRange, threshold);
            vector<ExpandingHorizontalPoint> *expandLeft = expandHorizontallyBothExpanding(pointsLeft, leftPoint, upperLevel, lowerLevel, redInRange > blueInRange, threshold);
            vector<ShrinkingHorizontalPoint> *shrinkRight = shrinkHorizontallyBothExpanding(pointsInside, upperLevel, lowerLevel, blueInRange > redInRange, threshold, minimum[1], maximum[1]);
            vector<ShrinkingHorizontalPoint> *shrinkLeft = shrinkHorizontallyBothExpanding(pointsInsideReversed, upperLevel, lowerLevel, blueInRange > redInRange, threshold, minimum[1], maximum[1]); 

            int i = 0;
            while (i < expandRight->size() - 1 && expandRight->at(i).improvement + upperLevel.improvementUnderClosest + lowerLevel.improvementUnderClosest < threshold) {
                i++;
            }
            int j = 0;
            while (j < expandLeft->size() && expandLeft->at(j).improvement + expandRight->at(i).improvement + upperLevel.improvementUnderClosest + lowerLevel.improvementUnderClosest < threshold) {
                j++;
            }

            while (i >= 0 && j < expandLeft->size()) {

                //expansions with 0 shrinks

                vector<Point> *upperSubrange = new vector<Point>();
                copy(expandRight->at(i).upperSubrange->begin(), expandRight->at(i).upperSubrange->end(), back_inserter(*upperSubrange));
                copy(expandLeft->at(j).upperSubrange->begin(), expandLeft->at(i).upperSubrange->end(), back_inserter(*upperSubrange));
                double maxY;
                int topPointsAdded = 0;
                if (upperSubrange->empty()) {
                    maxY = upperLevel.closest.coordinates[1];
                } else {
                    sort(upperLevel.intermediates->begin(), upperLevel.intermediates->end(), pointSorter(1, true));
                    sort(upperSubrange->begin(), upperSubrange->end(), pointSorter(1, true));
                    findY(maxY, topPointsAdded, upperSubrange, upperLevel.intermediates, expandRight->at(i).improvementUnderTop + expandLeft->at(j).improvementUnderTop, expandRight->at(i).improvement + expandLeft->at(j).improvement, blueInRange < redInRange, true);
                }
                vector<Point> *lowerSubrange = new vector<Point>();
                copy(expandRight->at(i).lowerSubrange->begin(), expandRight->at(i).lowerSubrange->end(), back_inserter(*lowerSubrange));
                copy(expandLeft->at(j).lowerSubrange->begin(), expandLeft->at(j).lowerSubrange->end(), back_inserter(*lowerSubrange));
                double minY;
                int bottomPointsAdded = 0;
                if (lowerSubrange->empty()) {
                    minY = lowerLevel.closest.coordinates[1];
                } else {
                    sort(lowerLevel.intermediates->begin(), lowerLevel.intermediates->end(), pointSorter(1, false));
                    sort(lowerSubrange->begin(), lowerSubrange->end(), pointSorter(1, false));
                    findY(minY, bottomPointsAdded, lowerSubrange, lowerLevel.intermediates, expandRight->at(i).improvementAboveBottom + expandLeft->at(j).improvementAboveBottom, expandRight->at(i).improvement + expandLeft->at(j).improvement, blueInRange < redInRange, false);
                }

                int pointsAdded = upperLevel.pointsUnderClosest + lowerLevel.pointsUnderClosest + expandRight->at(i).pointsInCenter + expandLeft->at(j).pointsInCenter + topPointsAdded + bottomPointsAdded;
                double currentSimilarity = (double) initialPoints / (double) (initialPoints + pointsAdded) ;

                if (currentSimilarity > similarity) {
                    currentSimilarity = similarity;
                    fairMinimum[0] = expandLeft->at(j).expandedPoint.coordinates[0];
                    fairMinimum[1] = minY;
                    fairMaximum[0] = expandRight->at(i).expandedPoint.coordinates[0];
                    fairMaximum[1] = maxY;
                }

                //expansions with one shrink on the right side
                if (i > 0) {
                    vector<Point> *upperSubrange = new vector<Point>();
                    copy(expandRight->at(i-1).upperSubrange->begin(), expandRight->at(i-1).upperSubrange->end(), back_inserter(*upperSubrange));
                    copy(expandLeft->at(j).upperSubrange->begin(), expandLeft->at(j).upperSubrange->end(), back_inserter(*upperSubrange));
                    double maxY;
                    int topPointsAdded = 0;
                    if (upperSubrange->empty()) {
                        maxY = upperLevel.closest.coordinates[1];
                    } else {
                        sort(upperLevel.intermediates->begin(), upperLevel.intermediates->end(), pointSorter(1, true));
                        sort(upperSubrange->begin(), upperSubrange->end(), pointSorter(1, true));
                        findY(maxY, topPointsAdded, upperSubrange, upperLevel.intermediates, expandRight->at(i-1).improvementUnderTop + expandLeft->at(j).improvementUnderTop, expandRight->at(i-1).improvement + expandLeft->at(j).improvement, blueInRange < redInRange, true);
                    }
                    vector<Point> *lowerSubrange = new vector<Point>();
                    copy(expandRight->at(i-1).lowerSubrange->begin(), expandRight->at(i-1).lowerSubrange->end(), back_inserter(*lowerSubrange));
                    copy(expandLeft->at(j).lowerSubrange->begin(), expandLeft->at(j).lowerSubrange->end(), back_inserter(*lowerSubrange));
                    double minY;
                    int bottomPointsAdded = 0;
                    if (lowerSubrange->empty()) {
                        minY = lowerLevel.closest.coordinates[1];
                    } else {
                        sort(lowerLevel.intermediates->begin(), lowerLevel.intermediates->end(), pointSorter(1, false));
                        sort(lowerSubrange->begin(), lowerSubrange->end(), pointSorter(1, false));
                        findY(minY, bottomPointsAdded, lowerSubrange, lowerLevel.intermediates, expandRight->at(i-1).improvementAboveBottom + expandLeft->at(j).improvementAboveBottom, expandRight->at(i-1).improvement + expandLeft->at(j).improvement, blueInRange < redInRange, false);
                    }

                    if (lowerLevel.finished) {
                        int pointsAdded = upperLevel.pointsUnderClosest + lowerLevel.pointsUnderClosest + lowerLevel.intermediates->size() + 1 + expandRight->at(i-1).pointsInCenter + expandRight->at(i-1).lowerSubrange->size() + expandLeft->at(j).pointsInCenter + topPointsAdded;
                        double currentSimilarity = (double) initialPoints / (double) (initialPoints + pointsAdded);

                        if (currentSimilarity > similarity) {
                            currentSimilarity = similarity;
                            fairMinimum[0] = expandLeft->at(j).expandedPoint.coordinates[0];
                            fairMinimum[1] = lowerLevel.furthest.coordinates[1];
                            fairMaximum[0] = expandRight->at(i-1).expandedPoint.coordinates[0];
                            fairMaximum[1] = maxY;
                        }
                    }

                    if (upperLevel.finished) {
                        int pointsAdded = upperLevel.pointsUnderClosest + upperLevel.intermediates->size() + 1 + lowerLevel.pointsUnderClosest + expandRight->at(i-1).pointsInCenter + expandRight->at(i-1).upperSubrange->size() + expandLeft->at(j).pointsInCenter + bottomPointsAdded;
                        double currentSimilarity = (double) initialPoints / (double) (initialPoints + pointsAdded);

                        if (currentSimilarity > similarity) {
                            currentSimilarity = similarity;
                            fairMinimum[0] = expandLeft->at(j).expandedPoint.coordinates[0];
                            fairMinimum[1] = minY;
                            fairMaximum[0] = expandRight->at(i-1).expandedPoint.coordinates[0];
                            fairMaximum[1] = upperLevel.furthest.coordinates[1];
                        }
                    }
                }

                //expansions with one shrink on the left side
                if (j > 0) {
                    vector<Point> *upperSubrange = new vector<Point>();
                    copy(expandRight->at(i).upperSubrange->begin(), expandRight->at(i).upperSubrange->end(), back_inserter(*upperSubrange));
                    copy(expandLeft->at(j-1).upperSubrange->begin(), expandLeft->at(j-1).upperSubrange->end(), back_inserter(*upperSubrange));
                    double maxY;
                    int topPointsAdded = 0;
                    if (upperSubrange->empty()) {
                        maxY = upperLevel.closest.coordinates[1];
                    } else {
                        sort(upperLevel.intermediates->begin(), upperLevel.intermediates->end(), pointSorter(1, true));
                        sort(upperSubrange->begin(), upperSubrange->end(), pointSorter(1, true));
                        findY(maxY, topPointsAdded, upperSubrange, upperLevel.intermediates, expandRight->at(i).improvementUnderTop + expandLeft->at(j-1).improvementUnderTop, expandRight->at(i).improvement + expandLeft->at(j-1).improvement, blueInRange < redInRange, true);
                    }
                    vector<Point> *lowerSubrange = new vector<Point>();
                    copy(expandRight->at(i).lowerSubrange->begin(), expandRight->at(i).lowerSubrange->end(), back_inserter(*lowerSubrange));
                    copy(expandLeft->at(j-1).lowerSubrange->begin(), expandLeft->at(j-1).lowerSubrange->end(), back_inserter(*lowerSubrange));
                    double minY;
                    int bottomPointsAdded = 0;
                    if (lowerSubrange->empty()) {
                        minY = lowerLevel.closest.coordinates[1];
                    } else {
                        sort(lowerLevel.intermediates->begin(), lowerLevel.intermediates->end(), pointSorter(1, false));
                        sort(lowerSubrange->begin(), lowerSubrange->end(), pointSorter(1, false));
                        findY(minY, bottomPointsAdded, lowerSubrange, lowerLevel.intermediates, expandRight->at(i).improvementAboveBottom + expandLeft->at(j-1).improvementAboveBottom, expandRight->at(i).improvement + expandLeft->at(j-1).improvement, blueInRange < redInRange, false);
                    }

                    if (lowerLevel.finished) {
                        int pointsAdded = upperLevel.pointsUnderClosest + lowerLevel.pointsUnderClosest + lowerLevel.intermediates->size() + 1 + expandRight->at(i).pointsInCenter + expandLeft->at(j-1).pointsInCenter + expandLeft->at(j-1).lowerSubrange->size() + topPointsAdded;
                        double currentSimilarity = (double) initialPoints / (double) (initialPoints + pointsAdded);

                        if (currentSimilarity > similarity) {
                            currentSimilarity = similarity;
                            fairMinimum[0] = expandLeft->at(j-1).expandedPoint.coordinates[0];
                            fairMinimum[1] = lowerLevel.furthest.coordinates[1];
                            fairMaximum[0] = expandRight->at(i).expandedPoint.coordinates[0];
                            fairMaximum[1] = maxY;
                        }
                    }
                    if (upperLevel.finished) {
                        int pointsAdded = upperLevel.pointsUnderClosest + upperLevel.intermediates->size() + 1 + lowerLevel.pointsUnderClosest + expandRight->at(i).pointsInCenter + expandLeft->at(j-1).pointsInCenter + expandLeft->at(j-1).upperSubrange->size() + bottomPointsAdded;
                        double currentSimilarity = (double) initialPoints / (double) (initialPoints + pointsAdded);

                        if (currentSimilarity > similarity) {
                            currentSimilarity = similarity;
                            fairMinimum[0] = expandLeft->at(j-1).expandedPoint.coordinates[0];
                            fairMinimum[1] = minY;
                            fairMaximum[0] = expandRight->at(i).expandedPoint.coordinates[0];
                            fairMaximum[1] = upperLevel.furthest.coordinates[1];
                        }
                    }
                }

                //expansions with one shrink on both sides

                if (i > 0 && j > 0 && lowerLevel.finished && upperLevel.finished) {
                    int pointsAdded = upperLevel.pointsUnderClosest + upperLevel.intermediates->size() + 1 + lowerLevel.pointsUnderClosest + lowerLevel.intermediates->size() + 1 + expandRight->at(i-1).pointsInCenter + expandLeft->at(j-1).pointsInCenter + expandLeft->at(j-1).lowerSubrange->size() + expandRight->at(i-1).lowerSubrange->size() + expandLeft->at(j-1).upperSubrange->size() + expandRight->at(i-1).upperSubrange->size();
                    double currentSimilarity = (double) initialPoints / (double) (initialPoints + pointsAdded);

                    if (currentSimilarity > similarity) {
                        currentSimilarity = similarity;
                        fairMinimum[0] = expandLeft->at(j-1).expandedPoint.coordinates[0];
                        fairMinimum[1] = lowerLevel.furthest.coordinates[1];
                        fairMaximum[0] = expandRight->at(i-1).expandedPoint.coordinates[0];
                        fairMaximum[1] = upperLevel.furthest.coordinates[1];
                    }
                }

                //expansions with 2 shrinks on right side

                if (i > 1 && lowerLevel.finished && upperLevel.finished) {
                    int pointsAdded = upperLevel.pointsUnderClosest + upperLevel.intermediates->size() + 1 + lowerLevel.pointsUnderClosest + lowerLevel.intermediates->size() + 1 + expandRight->at(i-2).pointsInCenter + expandLeft->at(j).pointsInCenter + expandLeft->at(j).lowerSubrange->size() + expandRight->at(i-2).lowerSubrange->size() + expandLeft->at(j).upperSubrange->size() + expandRight->at(i-2).upperSubrange->size();
                    double currentSimilarity = (double) initialPoints / (double) (initialPoints + pointsAdded);

                    if (currentSimilarity > similarity) {
                        currentSimilarity = similarity;
                        fairMinimum[0] = expandLeft->at(j).expandedPoint.coordinates[0];
                        fairMinimum[1] = lowerLevel.furthest.coordinates[1];
                        fairMaximum[0] = expandRight->at(i-2).expandedPoint.coordinates[0];
                        fairMaximum[1] = upperLevel.furthest.coordinates[1];
                    }
                }
                
                //expansions with 2 shrinks on left side

                if (j > 1 && lowerLevel.finished && upperLevel.finished) {
                    int pointsAdded = upperLevel.pointsUnderClosest + upperLevel.intermediates->size() + 1 + lowerLevel.pointsUnderClosest + lowerLevel.intermediates->size() + 1 + expandRight->at(i).pointsInCenter + expandLeft->at(j-2).pointsInCenter + expandLeft->at(j-2).lowerSubrange->size() + expandRight->at(i).lowerSubrange->size() + expandLeft->at(j-2).upperSubrange->size() + expandRight->at(i).upperSubrange->size();
                    double currentSimilarity = (double) initialPoints / (double) (initialPoints + pointsAdded);

                    if (currentSimilarity > similarity) {
                        currentSimilarity = similarity;
                        fairMinimum[0] = expandLeft->at(j-2).expandedPoint.coordinates[0];
                        fairMinimum[1] = lowerLevel.furthest.coordinates[1];
                        fairMaximum[0] = expandRight->at(i).expandedPoint.coordinates[0];
                        fairMaximum[1] = upperLevel.furthest.coordinates[1];
                    }
                }
                i--;
                j++;
            }

            int i = 0;
            while (i < expandRight->size() - 1 && expandRight->at(i).improvement + upperLevel.improvementUnderClosest + lowerLevel.improvementUnderClosest < threshold) {
                i++;
            }
            int j = 0;
            while (j < shrinkRight->size() && shrinkRight->at(j).improvement + expandRight->at(i).improvement + upperLevel.improvementUnderClosest + lowerLevel.improvementUnderClosest < threshold) {
                j++;
            }

            while (i >= 0 && j < shrinkRight->size()) {
                 //expansions with 0 shrinks
                double maxY;
                int topPointsAdded = 0;
                if (expandRight->at(i).upperSubrange->empty() && shrinkRight->at(j).pointsRemovedFromTop == 0) {
                    maxY = upperLevel.closest.coordinates[1];
                } else {
                    sort(upperLevel.intermediates->begin(), upperLevel.intermediates->end(), pointSorter(1, true));
                    vector<Point> *intermediates = new vector<Point>();
                    copy_if(upperLevel.intermediates->begin(), upperLevel.intermediates->end(), back_inserter(intermediates), [shrinkRight, j](Point pnt) {
                        return pnt.coordinates[0] >= shrinkRight->at(j).shrunkPoint.coordinates[0];
                    });
                    sort(expandRight->at(i).upperSubrange->begin(), expandRight->at(i).upperSubrange->end(), pointSorter(1, true));
                    findY(maxY, topPointsAdded, expandRight->at(i).upperSubrange, intermediates, expandRight->at(i).improvementUnderTop + shrinkRight->at(j).improveUnderTop, expandRight->at(i).improvement + shrinkRight->at(j).improvement, blueInRange < redInRange, true);
                }
                double minY;
                int bottomPointsAdded = 0;
                if (expandRight->at(i).lowerSubrange->empty() && shrinkRight->at(j).pointsRemovedFromBottom == 0) {
                    minY = lowerLevel.closest.coordinates[1];
                } else {
                    sort(lowerLevel.intermediates->begin(), lowerLevel.intermediates->end(), pointSorter(1, false));
                    vector<Point> *intermediates = new vector<Point>();
                    copy_if(lowerLevel.intermediates->begin(), lowerLevel.intermediates->end(), back_inserter(intermediates), [shrinkRight, j](Point pnt) {
                        return pnt.coordinates[0] >= shrinkRight->at(j).shrunkPoint.coordinates[0];                    });

                    sort(expandRight->at(i).lowerSubrange->begin(), expandRight->at(i).lowerSubrange->end(), pointSorter(1, false));
                    findY(minY, bottomPointsAdded, expandRight->at(i).lowerSubrange, intermediates, expandRight->at(i).improvementAboveBottom + shrinkRight->at(j).improvementAboveBottom, expandRight->at(i).improvement + shrinkRight->at(j).improvement, blueInRange < redInRange, false);
                }

                int pointsAdded = upperLevel.pointsUnderClosest + lowerLevel.pointsUnderClosest + expandRight->at(i).pointsInCenter - shrinkRight->at(j).pointsRemovedFromCenter + shrinkRight->at(j).pointsRemovedFromStartingRange + topPointsAdded + bottomPointsAdded;
                int pointsRemoved = shrinkRight->at(j).pointsRemovedFromStartingRange;
                double currentSimilarity = (double) (initialPoints - pointsRemoved) / (double) (initialPoints + pointsAdded);

                if (currentSimilarity > similarity) {
                    currentSimilarity = similarity;
                    fairMinimum[0] = shrinkRight->at(j).shrunkPoint.coordinates[0];
                    fairMinimum[1] = minY;
                    fairMaximum[0] = expandRight->at(i).expandedPoint.coordinates[0];
                    fairMaximum[1] = maxY;
                }

                 //expansions with one shrink on the right side
                if (i > 0) {
                    double maxY;
                    int topPointsAdded = 0;
                    if (expandRight->at(i-1).upperSubrange->empty() && shrinkRight->at(i-1).pointsRemovedFromTop == 0) {
                        maxY = upperLevel.closest.coordinates[1];
                    } else {
                        sort(upperLevel.intermediates->begin(), upperLevel.intermediates->end(), pointSorter(1, true));
                        vector<Point> *intermediates = new vector<Point>();
                        copy_if(upperLevel.intermediates->begin(), upperLevel.intermediates->end(), back_inserter(intermediates), [shrinkRight, j](Point pnt) {
                            return pnt.coordinates[0] >= shrinkRight->at(j).shrunkPoint.coordinates[0];
                        });
                        sort(expandRight->at(i-1).upperSubrange->begin(), expandRight->at(i-1).upperSubrange->end(), pointSorter(1, true));
                        findY(maxY, topPointsAdded, expandRight->at(i-1).upperSubrange, intermediates, expandRight->at(i-1).improvementUnderTop + shrinkLeft->at(j).improveUnderTop, expandRight->at(i-1).improvement + shrinkRight->at(j).improvement, blueInRange < redInRange, true);
                    }
                    double minY;
                    int bottomPointsAdded = 0;
                    if (expandRight->at(i-1).lowerSubrange->empty() && shrinkRight->at(i-1).pointsRemovedFromBottom == 0) {
                        minY = lowerLevel.closest.coordinates[1];
                    } else {
                        sort(lowerLevel.intermediates->begin(), lowerLevel.intermediates->end(), pointSorter(1, false));
                        vector<Point> *intermediates = new vector<Point>();
                        copy_if(upperLevel.intermediates->begin(), upperLevel.intermediates->end(), back_inserter(intermediates), [shrinkRight, j](Point pnt) {
                            return pnt.coordinates[0] >= shrinkRight->at(j).shrunkPoint.coordinates[0];
                        });
                        sort(expandRight->at(i-1).lowerSubrange->begin(), expandRight->at(i-1).lowerSubrange->end(), pointSorter(1, false));
                        findY(minY, bottomPointsAdded, expandRight->at(i-1).lowerSubrange, intermediates, expandRight->at(i-1).improvementAboveBottom + shrinkLeft->at(j).improvementAboveBottom, expandRight->at(i-1).improvement + shrinkLeft->at(j).improvement, blueInRange < redInRange, false);
                    }

                    if (lowerLevel.finished) {
                        int pointsAdded = upperLevel.pointsUnderClosest + lowerLevel.pointsUnderClosest + lowerLevel.intermediates->size() + 1 + expandRight->at(i-1).pointsInCenter - shrinkRight->at(j).pointsRemovedFromCenter + shrinkRight->at(j).pointsRemovedFromStartingRange + expandRight->at(i-1).lowerSubrange->size() + topPointsAdded;
                        int pointsRemoved = shrinkRight->at(j).pointsRemovedFromStartingRange;
                        double currentSimilarity = (double) (initialPoints - pointsRemoved) / (double) (initialPoints + pointsAdded);

                        if (currentSimilarity > similarity) {
                            currentSimilarity = similarity;
                            fairMinimum[0] = shrinkRight->at(j).shrunkPoint.coordinates[0];
                            fairMinimum[1] = lowerLevel.furthest.coordinates[1];
                            fairMaximum[0] = expandRight->at(i-1).expandedPoint.coordinates[0];
                            fairMaximum[1] = maxY;
                        }
                    }

                    if (upperLevel.finished) {
                        int pointsAdded = upperLevel.pointsUnderClosest + upperLevel.intermediates->size() + 1 + lowerLevel.pointsUnderClosest + expandRight->at(i-1).pointsInCenter - shrinkRight->at(j).pointsRemovedFromCenter + shrinkRight->at(j).pointsRemovedFromStartingRange + expandRight->at(i-1).upperSubrange->size() + bottomPointsAdded;
                        int pointsRemoved = shrinkRight->at(j).pointsRemovedFromStartingRange;
                        double currentSimilarity = (double) (initialPoints - pointsRemoved) / (double) (initialPoints + pointsAdded);
                        if (currentSimilarity > similarity) {
                            currentSimilarity = similarity;
                            fairMinimum[0] = shrinkRight->at(j).shrunkPoint.coordinates[0];
                            fairMinimum[1] = minY;
                            fairMaximum[0] = expandRight->at(i-1).expandedPoint.coordinates[0];
                            fairMaximum[1] = upperLevel.furthest.coordinates[1];;
                        }
                    }
                }

                //expansions with one shrink on the left side
                if (j > 0) {
                    double maxY;
                    int topPointsAdded = 0;
                    if (expandRight->at(i).upperSubrange->empty() && shrinkRight->at(j-1).pointsRemovedFromTop == 0) {
                        maxY = upperLevel.closest.coordinates[1];
                    } else {
                        sort(upperLevel.intermediates->begin(), upperLevel.intermediates->end(), pointSorter(1, true));
                        vector<Point> *intermediates = new vector<Point>();
                        copy_if(upperLevel.intermediates->begin(), upperLevel.intermediates->end(), back_inserter(intermediates), [shrinkRight, j](Point pnt) {
                            return pnt.coordinates[0] >= shrinkRight->at(j-1).shrunkPoint.coordinates[0];
                        });
                        sort(expandRight->at(i).upperSubrange->begin(), expandRight->at(i).upperSubrange->end(), pointSorter(1, true));
                        findY(maxY, topPointsAdded, expandRight->at(i).upperSubrange, intermediates, expandRight->at(i).improvementUnderTop + shrinkRight->at(j-1).improveUnderTop, expandRight->at(i).improvement + shrinkRight->at(j-1).improvement, blueInRange < redInRange, true);
                    }
                    double minY;
                    int bottomPointsAdded = 0;
                    if (expandRight->at(i).lowerSubrange->empty() && shrinkRight->at(j-1).pointsRemovedFromBottom == 0) {
                        minY = lowerLevel.closest.coordinates[1];
                    } else {
                        sort(lowerLevel.intermediates->begin(), lowerLevel.intermediates->end(), pointSorter(1, false));
                        vector<Point> *intermediates = new vector<Point>();
                        copy_if(upperLevel.intermediates->begin(), upperLevel.intermediates->end(), back_inserter(intermediates), [shrinkRight, j](Point pnt) {
                            return pnt.coordinates[0] >= shrinkRight->at(j-1).shrunkPoint.coordinates[0];
                        });

                        sort(expandRight->at(i).lowerSubrange->begin(), expandRight->at(i).lowerSubrange->end(), pointSorter(1, false));
                        findY(minY, bottomPointsAdded, expandRight->at(i).lowerSubrange, intermediates, expandRight->at(i).improvementAboveBottom + shrinkRight->at(j-1).improvementAboveBottom, expandRight->at(i).improvement + shrinkRight->at(j-1).improvement, blueInRange < redInRange, false);
                    }

                    if (lowerLevel.finished) {
                        int pointsAdded = upperLevel.pointsUnderClosest +lowerLevel.pointsUnderClosest + lowerLevel.intermediates->size() + 1 + expandRight->at(i).pointsInCenter - shrinkRight->at(j-1).pointsRemovedFromCenter + shrinkRight->at(j-1).pointsRemovedFromStartingRange + expandRight->at(i).lowerSubrange->size() + topPointsAdded;
                        int pointsRemoved = shrinkRight->at(j-1).pointsRemovedFromStartingRange;
                        double currentSimilarity = (double) (initialPoints - pointsRemoved) / (double) (initialPoints + pointsAdded);

                        if (currentSimilarity > similarity) {
                            currentSimilarity = similarity;
                            fairMinimum[0] = shrinkRight->at(j-1).shrunkPoint.coordinates[0];
                            fairMinimum[1] = lowerLevel.furthest.coordinates[1];
                            fairMaximum[0] = expandRight->at(i).expandedPoint.coordinates[0];
                            fairMaximum[1] = maxY;
                        }
                    }

                    if (upperLevel.finished) {
                        int pointsAdded = upperLevel.pointsUnderClosest + upperLevel.intermediates->size() + 1 + lowerLevel.pointsUnderClosest + expandRight->at(i).pointsInCenter - shrinkRight->at(j-1).pointsRemovedFromCenter + shrinkRight->at(j-1).pointsRemovedFromStartingRange + expandRight->at(i).upperSubrange->size() + bottomPointsAdded;
                        int pointsRemoved = shrinkRight->at(j-1).pointsRemovedFromStartingRange;
                        double currentSimilarity = (double) (initialPoints - pointsRemoved) / (double) (initialPoints + pointsAdded);

                        if (currentSimilarity > similarity) {
                            currentSimilarity = similarity;
                            fairMinimum[0] = shrinkRight->at(j-1).shrunkPoint.coordinates[0];
                            fairMinimum[1] = minY;
                            fairMaximum[0] = expandRight->at(i).expandedPoint.coordinates[0];
                            fairMaximum[1] = upperLevel.furthest.coordinates[1];
                        }
                    }
                }
                
                //expansions with one shrink on both sides
                if (i > 0 && j > 0 && lowerLevel.finished && upperLevel.finished) {
                    int pointsAdded = upperLevel.pointsUnderClosest + upperLevel.intermediates->size() + 1 + lowerLevel.pointsUnderClosest + lowerLevel.intermediates->size() + 1 + expandRight->at(i-1).pointsInCenter - shrinkRight->at(j-1).pointsRemovedFromCenter + shrinkRight->at(j-1).pointsRemovedFromStartingRange + expandRight->at(i-1).upperSubrange->size() + expandRight->at(i-1).upperSubrange->size();
                    int pointsRemoved = shrinkRight->at(j-1).pointsRemovedFromStartingRange;
                    double currentSimilarity = (double) (initialPoints - pointsRemoved) / (double) (initialPoints + pointsAdded);

                    if (currentSimilarity > similarity) {
                        currentSimilarity = similarity;
                        fairMinimum[0] = shrinkRight->at(j-1).shrunkPoint.coordinates[0];
                        fairMinimum[1] = lowerLevel.furthest.coordinates[1];
                        fairMaximum[0] = expandRight->at(i-1).expandedPoint.coordinates[0];
                        fairMaximum[1] = upperLevel.furthest.coordinates[1];
                    }
                }

                //expansions with 2 shrinks on right side

                if (i > 1 && lowerLevel.finished && upperLevel.finished) {
                    int pointsAdded = upperLevel.pointsUnderClosest + upperLevel.intermediates->size() + 1 + lowerLevel.pointsUnderClosest + lowerLevel.intermediates->size() + 1 + expandRight->at(i-2).pointsInCenter - shrinkRight->at(j).pointsRemovedFromCenter + shrinkRight->at(j).pointsRemovedFromStartingRange + expandRight->at(i-2).upperSubrange->size() + expandRight->at(i-2).upperSubrange->size();
                    int pointsRemoved = shrinkRight->at(j).pointsRemovedFromStartingRange;
                    double currentSimilarity = (double) (initialPoints - pointsRemoved) / (double) (initialPoints + pointsAdded);

                    if (currentSimilarity > similarity) {
                        currentSimilarity = similarity;
                        fairMinimum[0] = shrinkRight->at(j).shrunkPoint.coordinates[0];
                        fairMinimum[1] = lowerLevel.furthest.coordinates[1];
                        fairMaximum[0] = expandRight->at(i-2).expandedPoint.coordinates[0];
                        fairMaximum[1] = upperLevel.furthest.coordinates[1];
                    }
                }
                
                //expansions with 2 shrinks on left side

                if (j > 1 && lowerLevel.finished && upperLevel.finished) {
                    int pointsAdded = upperLevel.pointsUnderClosest + upperLevel.intermediates->size() + 1 + lowerLevel.pointsUnderClosest + lowerLevel.intermediates->size() + 1 + expandRight->at(i).pointsInCenter - shrinkRight->at(j-2).pointsRemovedFromCenter + shrinkRight->at(j-2).pointsRemovedFromStartingRange + expandRight->at(i).lowerSubrange->size() + expandRight->at(i).upperSubrange->size();
                    int pointsRemoved = shrinkRight->at(j-2).pointsRemovedFromStartingRange;
                    double currentSimilarity = (double) (initialPoints - pointsRemoved) / (double) (initialPoints + pointsAdded);

                    if (currentSimilarity > similarity) {
                        currentSimilarity = similarity;
                        fairMinimum[0] = shrinkRight->at(j-2).shrunkPoint.coordinates[0];
                        fairMinimum[1] = lowerLevel.furthest.coordinates[1];
                        fairMaximum[0] = expandRight->at(i).expandedPoint.coordinates[0];
                        fairMaximum[1] = upperLevel.furthest.coordinates[1];
                    }
                }
                i--;
                j++;
            }

            int i = 0;
            while (i < expandLeft->size() - 1 && expandLeft->at(i).improvement + upperLevel.improvementUnderClosest + lowerLevel.improvementUnderClosest < threshold) {
                i++;
            }
            int j = 0;
            while (j < shrinkLeft->size() && shrinkLeft->at(j).improvement + expandLeft->at(i).improvement + upperLevel.improvementUnderClosest + lowerLevel.improvementUnderClosest < threshold) {
                j++;
            }

            while (i >= 0 && j < shrinkLeft->size()) {
                 //expansions with 0 shrinks
                double maxY;
                int topPointsAdded = 0;
                if (expandLeft->at(i).upperSubrange->empty() && shrinkLeft->at(j).pointsRemovedFromTop == 0) {
                    maxY = upperLevel.closest.coordinates[1];
                } else {
                    sort(upperLevel.intermediates->begin(), upperLevel.intermediates->end(), pointSorter(1, true));
                    vector<Point> *intermediates = new vector<Point>();
                    copy_if(upperLevel.intermediates->begin(), upperLevel.intermediates->end(), back_inserter(intermediates), [shrinkLeft, j](Point pnt) {
                        return pnt.coordinates[0] <= shrinkLeft->at(j).shrunkPoint.coordinates[0];
                    });
                    sort(expandLeft->at(i).upperSubrange->begin(), expandLeft->at(i).upperSubrange->end(), pointSorter(1, true));
                    findY(maxY, topPointsAdded, expandLeft->at(i).upperSubrange, intermediates, expandLeft->at(i).improvementUnderTop + shrinkLeft->at(j).improveUnderTop, expandLeft->at(i).improvement + shrinkLeft->at(j).improvement, blueInRange < redInRange, true);
                }
                double minY;
                int bottomPointsAdded = 0;
                if (expandLeft->at(i).lowerSubrange->empty() && shrinkLeft->at(j).pointsRemovedFromBottom == 0) {
                    minY = lowerLevel.closest.coordinates[1];
                } else {
                    sort(lowerLevel.intermediates->begin(), lowerLevel.intermediates->end(), pointSorter(1, false));
                    vector<Point> *intermediates = new vector<Point>();
                    copy_if(lowerLevel.intermediates->begin(), lowerLevel.intermediates->end(), back_inserter(intermediates), [shrinkLeft, j](Point pnt) {
                        return pnt.coordinates[0] <= shrinkLeft->at(j).shrunkPoint.coordinates[0];                    
                    });

                    sort(expandLeft->at(i).lowerSubrange->begin(), expandLeft->at(i).lowerSubrange->end(), pointSorter(1, false));
                    findY(minY, bottomPointsAdded, expandLeft->at(i).lowerSubrange, intermediates, expandLeft->at(i).improvementAboveBottom + shrinkLeft->at(j).improvementAboveBottom, expandLeft->at(i).improvement + shrinkLeft->at(j).improvement, blueInRange < redInRange, false);
                }

                int pointsAdded = upperLevel.pointsUnderClosest + lowerLevel.pointsUnderClosest + expandLeft->at(i).pointsInCenter - shrinkLeft->at(j).pointsRemovedFromCenter + shrinkLeft->at(j).pointsRemovedFromStartingRange + topPointsAdded + bottomPointsAdded;
                int pointsRemoved = shrinkLeft->at(j).pointsRemovedFromStartingRange;
                double currentSimilarity = (double) (initialPoints - pointsRemoved) / (double) (initialPoints + pointsAdded);

                if (currentSimilarity > similarity) {
                    currentSimilarity = similarity;
                    fairMinimum[0] = expandLeft->at(i).expandedPoint.coordinates[0];
                    fairMinimum[1] = minY;
                    fairMaximum[0] = shrinkLeft->at(j).shrunkPoint.coordinates[0];
                    fairMaximum[1] = maxY;
                }

                 //expansions with one shrink on the left side
                if (i > 0) {
                    double maxY;
                    int topPointsAdded = 0;
                    if (expandLeft->at(i-1).upperSubrange->empty() && shrinkLeft->at(i-1).pointsRemovedFromTop == 0) {
                        maxY = upperLevel.closest.coordinates[1];
                    } else {
                        sort(upperLevel.intermediates->begin(), upperLevel.intermediates->end(), pointSorter(1, true));
                        vector<Point> *intermediates = new vector<Point>();
                        copy_if(upperLevel.intermediates->begin(), upperLevel.intermediates->end(), back_inserter(intermediates), [shrinkLeft, j](Point pnt) {
                            return pnt.coordinates[0] <= shrinkLeft->at(j).shrunkPoint.coordinates[0];
                        });
                        sort(expandLeft->at(i-1).upperSubrange->begin(), expandLeft->at(i-1).upperSubrange->end(), pointSorter(1, true));
                        findY(maxY, topPointsAdded, expandLeft->at(i-1).upperSubrange, intermediates, expandLeft->at(i-1).improvementUnderTop + shrinkLeft->at(j).improveUnderTop, expandLeft->at(i-1).improvement + shrinkLeft->at(j).improvement, blueInRange < redInRange, true);
                    }
                    double minY;
                    int bottomPointsAdded = 0;
                    if (expandLeft->at(i-1).lowerSubrange->empty() && shrinkLeft->at(i-1).pointsRemovedFromBottom == 0) {
                        minY = lowerLevel.closest.coordinates[1];
                    } else {
                        sort(lowerLevel.intermediates->begin(), lowerLevel.intermediates->end(), pointSorter(1, false));
                        vector<Point> *intermediates = new vector<Point>();
                        copy_if(upperLevel.intermediates->begin(), upperLevel.intermediates->end(), back_inserter(intermediates), [shrinkLeft, j](Point pnt) {
                            return pnt.coordinates[0] <= shrinkLeft->at(j).shrunkPoint.coordinates[0];
                        });
                        sort(expandLeft->at(i-1).lowerSubrange->begin(), expandLeft->at(i-1).lowerSubrange->end(), pointSorter(1, false));
                        findY(minY, bottomPointsAdded, expandLeft->at(i-1).lowerSubrange, intermediates, expandLeft->at(i-1).improvementAboveBottom + shrinkLeft->at(j).improvementAboveBottom, expandLeft->at(i-1).improvement + shrinkLeft->at(j).improvement, blueInRange < redInRange, false);
                    }

                    if (lowerLevel.finished) {
                        int pointsAdded = upperLevel.pointsUnderClosest + lowerLevel.pointsUnderClosest + lowerLevel.intermediates->size() + 1 + expandLeft->at(i-1).pointsInCenter - shrinkLeft->at(j).pointsRemovedFromCenter + shrinkLeft->at(j).pointsRemovedFromStartingRange + expandLeft->at(i-1).lowerSubrange->size() + topPointsAdded;
                        int pointsRemoved = shrinkLeft->at(j).pointsRemovedFromStartingRange;
                        double currentSimilarity = (double) (initialPoints - pointsRemoved) / (double) (initialPoints + pointsAdded);

                        if (currentSimilarity > similarity) {
                            currentSimilarity = similarity;
                            fairMinimum[0] = expandLeft->at(i-1).expandedPoint.coordinates[0];
                            fairMinimum[1] = lowerLevel.furthest.coordinates[1];
                            fairMaximum[0] = shrinkLeft->at(j).shrunkPoint.coordinates[0];
                            fairMaximum[1] = maxY;
                        }
                    }

                    if (upperLevel.finished) {
                        int pointsAdded = upperLevel.pointsUnderClosest + upperLevel.intermediates->size() + 1 + lowerLevel.pointsUnderClosest + expandLeft->at(i-1).pointsInCenter - shrinkLeft->at(j).pointsRemovedFromCenter + shrinkLeft->at(j).pointsRemovedFromStartingRange + expandLeft->at(i-1).upperSubrange->size() + bottomPointsAdded;
                        int pointsRemoved = shrinkLeft->at(j).pointsRemovedFromStartingRange;
                        double currentSimilarity = (double) (initialPoints - pointsRemoved) / (double) (initialPoints + pointsAdded);
                        if (currentSimilarity > similarity) {
                            currentSimilarity = similarity;
                            fairMinimum[0] = expandLeft->at(i-1).expandedPoint.coordinates[0];
                            fairMinimum[1] = minY;
                            fairMaximum[0] = shrinkLeft->at(j).shrunkPoint.coordinates[0];
                            fairMaximum[1] = upperLevel.furthest.coordinates[1];
                        }
                    }
                }

                //expansions with one shrink on the right side
                if (j > 0) {
                    double maxY;
                    int topPointsAdded = 0;
                    if (expandLeft->at(i).upperSubrange->empty() && shrinkLeft->at(j-1).pointsRemovedFromTop == 0) {
                        maxY = upperLevel.closest.coordinates[1];
                    } else {
                        sort(upperLevel.intermediates->begin(), upperLevel.intermediates->end(), pointSorter(1, true));
                        vector<Point> *intermediates = new vector<Point>();
                        copy_if(upperLevel.intermediates->begin(), upperLevel.intermediates->end(), back_inserter(intermediates), [shrinkLeft, j](Point pnt) {
                            return pnt.coordinates[0] <= shrinkLeft->at(j-1).shrunkPoint.coordinates[0];
                        });
                        sort(expandLeft->at(i).upperSubrange->begin(), expandLeft->at(i).upperSubrange->end(), pointSorter(1, true));
                        findY(maxY, topPointsAdded, expandLeft->at(i).upperSubrange, intermediates, expandLeft->at(i).improvementUnderTop + shrinkLeft->at(j-1).improveUnderTop, expandLeft->at(i).improvement + shrinkLeft->at(j-1).improvement, blueInRange < redInRange, true);
                    }
                    double minY;
                    int bottomPointsAdded = 0;
                    if (expandLeft->at(i).lowerSubrange->empty() && shrinkLeft->at(j-1).pointsRemovedFromBottom == 0) {
                        minY = lowerLevel.closest.coordinates[1];
                    } else {
                        sort(lowerLevel.intermediates->begin(), lowerLevel.intermediates->end(), pointSorter(1, false));
                        vector<Point> *intermediates = new vector<Point>();
                        copy_if(upperLevel.intermediates->begin(), upperLevel.intermediates->end(), back_inserter(intermediates), [shrinkLeft, j](Point pnt) {
                            return pnt.coordinates[0] <= shrinkLeft->at(j-1).shrunkPoint.coordinates[0];
                        });

                        sort(expandLeft->at(i).lowerSubrange->begin(), expandLeft->at(i).lowerSubrange->end(), pointSorter(1, false));
                        findY(minY, bottomPointsAdded, expandLeft->at(i).lowerSubrange, intermediates, expandLeft->at(i).improvementAboveBottom + shrinkLeft->at(j-1).improvementAboveBottom, expandLeft->at(i).improvement + shrinkLeft->at(j-1).improvement, blueInRange < redInRange, false);
                    }

                    if (lowerLevel.finished) {
                        int pointsAdded = upperLevel.pointsUnderClosest +lowerLevel.pointsUnderClosest + lowerLevel.intermediates->size() + 1 + expandLeft->at(i).pointsInCenter - shrinkLeft->at(j-1).pointsRemovedFromCenter + shrinkLeft->at(j-1).pointsRemovedFromStartingRange + expandLeft->at(i).lowerSubrange->size() + topPointsAdded;
                        int pointsRemoved = shrinkLeft->at(j-1).pointsRemovedFromStartingRange;
                        double currentSimilarity = (double) (initialPoints - pointsRemoved) / (double) (initialPoints + pointsAdded);

                        if (currentSimilarity > similarity) {
                            currentSimilarity = similarity;
                            fairMinimum[0] = expandLeft->at(i).expandedPoint.coordinates[0];
                            fairMinimum[1] = lowerLevel.furthest.coordinates[1];
                            fairMaximum[0] = shrinkLeft->at(j-1).shrunkPoint.coordinates[0];
                            fairMaximum[1] = maxY;
                        }
                    }

                    if (upperLevel.finished) {
                        int pointsAdded = upperLevel.pointsUnderClosest + upperLevel.intermediates->size() + 1 + lowerLevel.pointsUnderClosest + expandLeft->at(i).pointsInCenter - shrinkLeft->at(j-1).pointsRemovedFromCenter + shrinkLeft->at(j-1).pointsRemovedFromStartingRange + expandLeft->at(i).upperSubrange->size() + bottomPointsAdded;
                        int pointsRemoved = shrinkLeft->at(j-1).pointsRemovedFromStartingRange;
                        double currentSimilarity = (double) (initialPoints - pointsRemoved) / (double) (initialPoints + pointsAdded);

                        if (currentSimilarity > similarity) {
                            currentSimilarity = similarity;
                            fairMinimum[0] = expandLeft->at(i).expandedPoint.coordinates[0];
                            fairMinimum[1] = minY;
                            fairMaximum[0] = shrinkLeft->at(j-1).shrunkPoint.coordinates[0];
                            fairMaximum[1] = upperLevel.furthest.coordinates[1];

                        }
                    }
                }
                
                //expansions with one shrink on both sides
                if (i > 0 && j > 0 && lowerLevel.finished && upperLevel.finished) {
                    int pointsAdded = upperLevel.pointsUnderClosest + upperLevel.intermediates->size() + 1 + lowerLevel.pointsUnderClosest + lowerLevel.intermediates->size() + 1 + expandLeft->at(i-1).pointsInCenter - shrinkLeft->at(j-1).pointsRemovedFromCenter + shrinkLeft->at(j-1).pointsRemovedFromStartingRange + expandLeft->at(i-1).upperSubrange->size() + expandLeft->at(i-1).lowerSubrange->size();
                    int pointsRemoved = shrinkLeft->at(j-1).pointsRemovedFromStartingRange;
                    double currentSimilarity = (double) (initialPoints - pointsRemoved) / (double) (initialPoints + pointsAdded);

                    if (currentSimilarity > similarity) {
                        currentSimilarity = similarity;
                        fairMinimum[0] = expandLeft->at(i-1).expandedPoint.coordinates[0];
                        fairMinimum[1] = lowerLevel.furthest.coordinates[1];
                        fairMaximum[0] = shrinkLeft->at(j-1).shrunkPoint.coordinates[0];
                        fairMaximum[1] = upperLevel.furthest.coordinates[1];
                    }
                }

                //expansions with 2 shrinks on left side

                if (i > 1 && lowerLevel.finished && upperLevel.finished) {
                    int pointsAdded = upperLevel.pointsUnderClosest + upperLevel.intermediates->size() + 1 + lowerLevel.pointsUnderClosest + lowerLevel.intermediates->size() + 1 + expandLeft->at(i-2).pointsInCenter - shrinkLeft->at(j).pointsRemovedFromCenter + shrinkLeft->at(j).pointsRemovedFromStartingRange + expandLeft->at(i-2).upperSubrange->size() + expandLeft->at(i-2).lowerSubrange->size();
                    int pointsRemoved = shrinkLeft->at(j).pointsRemovedFromStartingRange;
                    double currentSimilarity = (double) (initialPoints - pointsRemoved) / (double) (initialPoints + pointsAdded);

                    if (currentSimilarity > similarity) {
                        currentSimilarity = similarity;
                        fairMinimum[0] = expandLeft->at(i-2).expandedPoint.coordinates[0];
                        fairMinimum[1] = lowerLevel.furthest.coordinates[1];
                        fairMaximum[0] = shrinkLeft->at(j).shrunkPoint.coordinates[0];
                        fairMaximum[1] = upperLevel.furthest.coordinates[1];
                    }
                }
                
                //expansions with 2 shrinks on right side

                if (j > 1 && lowerLevel.finished && upperLevel.finished) {
                    int pointsAdded = upperLevel.pointsUnderClosest + upperLevel.intermediates->size() + 1 + lowerLevel.pointsUnderClosest + lowerLevel.intermediates->size() + 1 + expandLeft->at(i).pointsInCenter - shrinkLeft->at(j-2).pointsRemovedFromCenter + shrinkLeft->at(j-2).pointsRemovedFromStartingRange + expandLeft->at(i).lowerSubrange->size() + expandLeft->at(i).upperSubrange->size();
                    int pointsRemoved = shrinkLeft->at(j-2).pointsRemovedFromStartingRange;
                    double currentSimilarity = (double) (initialPoints - pointsRemoved) / (double) (initialPoints + pointsAdded);

                    if (currentSimilarity > similarity) {
                        currentSimilarity = similarity;
                        fairMinimum[0] = expandLeft->at(i).expandedPoint.coordinates[0];
                        fairMinimum[1] = lowerLevel.furthest.coordinates[1];
                        fairMaximum[0] = shrinkLeft->at(j-2).shrunkPoint.coordinates[0];
                        fairMaximum[1] = upperLevel.furthest.coordinates[1];
                    }
                }
                i--;
                j++;
            }

            int i = 0;
            while (i < shrinkRight->size() - 1 && shrinkRight->at(i).improvement + upperLevel.improvementUnderClosest + lowerLevel.improvementUnderClosest < threshold) {
                i++;
            }
            int j = 0;
            while (j < shrinkLeft->size() && shrinkLeft->at(j).improvement + expandLeft->at(i).improvement + upperLevel.improvementUnderClosest + lowerLevel.improvementUnderClosest < threshold) {
                j++;
            }

            if (shrinkRight->at(i).shrunkPoint.coordinates[0] <= shrinkLeft->at(j).shrunkPoint.coordinates[0]) {
                while (i >= 0 && j < shrinkLeft->size()) {
                    //expansions with 0 shrinks
                    double maxY;
                    int topPointsAdded = 0;
                    if (shrinkRight->at(i).pointsRemovedFromTop == 0 && shrinkLeft->at(j).pointsRemovedFromTop == 0) {
                        maxY = upperLevel.closest.coordinates[1];
                    } else {
                        sort(upperLevel.intermediates->begin(), upperLevel.intermediates->end(), pointSorter(1, true));
                        vector<Point> *intermediates = new vector<Point>();
                        copy_if(upperLevel.intermediates->begin(), upperLevel.intermediates->end(), back_inserter(intermediates), [shrinkLeft, shrinkRight, i, j](Point pnt) {
                            return pnt.coordinates[0] <= shrinkLeft->at(j).shrunkPoint.coordinates[0] && pnt.coordinates[0] >= shrinkRight->at(i).shrunkPoint.coordinates[0];
                        });
                        findY(maxY, topPointsAdded, new vector<Point>(), intermediates, shrinkRight->at(i).improveUnderTop + shrinkLeft->at(j).improveUnderTop, shrinkRight->at(i).improvement + shrinkLeft->at(j).improvement, blueInRange < redInRange, true);
                    }
                    double minY;
                    int bottomPointsAdded = 0;
                    if (shrinkRight->at(i).pointsRemovedFromBottom == 0 && shrinkLeft->at(j).pointsRemovedFromBottom == 0) {
                        minY = lowerLevel.closest.coordinates[1];
                    } else {
                        sort(lowerLevel.intermediates->begin(), lowerLevel.intermediates->end(), pointSorter(1, false));
                        vector<Point> *intermediates = new vector<Point>();
                        copy_if(lowerLevel.intermediates->begin(), lowerLevel.intermediates->end(), back_inserter(intermediates), [shrinkLeft, shrinkRight, i, j](Point pnt) {
                            return pnt.coordinates[0] <= shrinkLeft->at(j).shrunkPoint.coordinates[0] && pnt.coordinates[0] >= shrinkRight->at(i).shrunkPoint.coordinates[0];
                        });
                        findY(minY, bottomPointsAdded, new vector<Point>(), intermediates, shrinkRight->at(i).improvementAboveBottom + shrinkLeft->at(j).improvementAboveBottom, shrinkRight->at(i).improvement + shrinkLeft->at(j).improvement, blueInRange < redInRange, false);
                    }

                    int pointsAdded = upperLevel.pointsUnderClosest + lowerLevel.pointsUnderClosest - shrinkLeft->at(j).pointsRemovedFromCenter + shrinkLeft->at(j).pointsRemovedFromStartingRange - shrinkRight->at(i).pointsRemovedFromCenter + shrinkRight->at(i).pointsRemovedFromStartingRange + topPointsAdded + bottomPointsAdded;
                    int pointsRemoved = shrinkLeft->at(j).pointsRemovedFromStartingRange + shrinkRight->at(i).pointsRemovedFromStartingRange;
                    double currentSimilarity = (double) (initialPoints - pointsRemoved) / (double) (initialPoints + pointsAdded);

                    if (currentSimilarity > similarity) {
                        currentSimilarity = similarity;
                        fairMinimum[0] = shrinkRight->at(i).shrunkPoint.coordinates[0];
                        fairMinimum[1] = minY;
                        fairMaximum[0] = shrinkLeft->at(j).shrunkPoint.coordinates[0];
                        fairMaximum[1] = maxY;
                    }

                    //expansions with one shrink on the left side
                    if (i > 0) {
                        double maxY;
                        int topPointsAdded = 0;
                        if (shrinkRight->at(i-1).pointsRemovedFromTop == 0 && shrinkLeft->at(i-1).pointsRemovedFromTop == 0) {
                            maxY = upperLevel.closest.coordinates[1];
                        } else {
                            sort(upperLevel.intermediates->begin(), upperLevel.intermediates->end(), pointSorter(1, true));
                            vector<Point> *intermediates = new vector<Point>();
                            copy_if(lowerLevel.intermediates->begin(), lowerLevel.intermediates->end(), back_inserter(intermediates), [shrinkLeft, shrinkRight, i, j](Point pnt) {
                                return pnt.coordinates[0] <= shrinkLeft->at(j).shrunkPoint.coordinates[0] && pnt.coordinates[0] >= shrinkRight->at(i-1).shrunkPoint.coordinates[0];
                            });
                            findY(maxY, topPointsAdded, new vector<Point>(), intermediates, shrinkRight->at(i-1).improveUnderTop + shrinkLeft->at(j).improveUnderTop, shrinkRight->at(i-1).improvement + shrinkLeft->at(j).improvement, blueInRange < redInRange, true);
                        }
                        double minY;
                        int bottomPointsAdded = 0;
                        if (shrinkRight->at(i-1).pointsRemovedFromBottom == 0 && shrinkLeft->at(i-1).pointsRemovedFromBottom == 0) {
                            minY = lowerLevel.closest.coordinates[1];
                        } else {
                            sort(lowerLevel.intermediates->begin(), lowerLevel.intermediates->end(), pointSorter(1, false));
                            vector<Point> *intermediates = new vector<Point>();
                            copy_if(lowerLevel.intermediates->begin(), lowerLevel.intermediates->end(), back_inserter(intermediates), [shrinkLeft, shrinkRight, i, j](Point pnt) {
                                return pnt.coordinates[0] <= shrinkLeft->at(j).shrunkPoint.coordinates[0] && pnt.coordinates[0] >= shrinkRight->at(i-1).shrunkPoint.coordinates[0];
                            });
                            findY(minY, bottomPointsAdded, new vector<Point>(), intermediates, shrinkRight->at(i-1).improvementAboveBottom + shrinkLeft->at(j).improvementAboveBottom, shrinkRight->at(i-1).improvement + shrinkLeft->at(j).improvement, blueInRange < redInRange, false);
                        }

                        if (lowerLevel.finished) {
                            int pointsAdded = upperLevel.pointsUnderClosest + lowerLevel.pointsUnderClosest + lowerLevel.intermediates->size() + 1 - shrinkRight->at(i-1).pointsRemovedFromCenter + shrinkRight->at(i-1).pointsRemovedFromStartingRange - shrinkLeft->at(j).pointsRemovedFromCenter + shrinkLeft->at(j).pointsRemovedFromStartingRange + topPointsAdded;
                            int pointsRemoved = shrinkRight->at(i-1).pointsRemovedFromStartingRange + shrinkLeft->at(j).pointsRemovedFromStartingRange;
                            double currentSimilarity = (double) (initialPoints - pointsRemoved) / (double) (initialPoints + pointsAdded);

                            if (currentSimilarity > similarity) {
                                currentSimilarity = similarity;
                                fairMinimum[0] = shrinkRight->at(i-1).shrunkPoint.coordinates[0];
                                fairMinimum[1] = lowerLevel.furthest.coordinates[1];
                                fairMaximum[0] = shrinkLeft->at(j).shrunkPoint.coordinates[0];
                                fairMaximum[1] = maxY;
                            }
                        }

                        if (upperLevel.finished) {
                            int pointsAdded = upperLevel.pointsUnderClosest + upperLevel.intermediates->size() + 1 + lowerLevel.pointsUnderClosest - shrinkRight->at(i-1).pointsRemovedFromCenter + shrinkRight->at(i-1).pointsRemovedFromStartingRange - shrinkLeft->at(j).pointsRemovedFromCenter + shrinkLeft->at(j).pointsRemovedFromStartingRange + bottomPointsAdded;
                            int pointsRemoved = shrinkRight->at(i-1).pointsRemovedFromStartingRange + shrinkLeft->at(j).pointsRemovedFromStartingRange;
                            double currentSimilarity = (double) (initialPoints - pointsRemoved) / (double) (initialPoints + pointsAdded);
                            if (currentSimilarity > similarity) {
                                currentSimilarity = similarity;
                                fairMinimum[0] = shrinkRight->at(i-1).shrunkPoint.coordinates[0];
                                fairMinimum[1] = minY;
                                fairMaximum[0] = shrinkLeft->at(j).shrunkPoint.coordinates[0];
                                fairMaximum[1] = upperLevel.furthest.coordinates[1];
                            }
                        }
                    }

                    //expansions with one shrink on the right side
                    if (j > 0) {
                        double maxY;
                        int topPointsAdded = 0;
                        if (shrinkRight->at(i).pointsRemovedFromTop == 0 && shrinkLeft->at(j-1).pointsRemovedFromTop == 0) {
                            maxY = upperLevel.closest.coordinates[1];
                        } else {
                            sort(upperLevel.intermediates->begin(), upperLevel.intermediates->end(), pointSorter(1, true));
                            vector<Point> *intermediates = new vector<Point>();
                            copy_if(lowerLevel.intermediates->begin(), lowerLevel.intermediates->end(), back_inserter(intermediates), [shrinkLeft, shrinkRight, i, j](Point pnt) {
                                return pnt.coordinates[0] <= shrinkLeft->at(j-1).shrunkPoint.coordinates[0] && pnt.coordinates[0] >= shrinkRight->at(i).shrunkPoint.coordinates[0];
                            });
                            findY(maxY, topPointsAdded, new vector<Point>(), intermediates, shrinkRight->at(i).improveUnderTop + shrinkLeft->at(j-1).improveUnderTop, expandLeft->at(i).improvement + shrinkLeft->at(j-1).improvement, blueInRange < redInRange, true);
                        }
                        double minY;
                        int bottomPointsAdded = 0;
                        if (shrinkRight->at(i).pointsRemovedFromBottom == 0 && shrinkLeft->at(j-1).pointsRemovedFromBottom == 0) {
                            minY = lowerLevel.closest.coordinates[1];
                        } else {
                            sort(lowerLevel.intermediates->begin(), lowerLevel.intermediates->end(), pointSorter(1, false));
                            vector<Point> *intermediates = new vector<Point>();
                            copy_if(lowerLevel.intermediates->begin(), lowerLevel.intermediates->end(), back_inserter(intermediates), [shrinkLeft, shrinkRight, i, j](Point pnt) {
                                return pnt.coordinates[0] <= shrinkLeft->at(j-1).shrunkPoint.coordinates[0] && pnt.coordinates[0] >= shrinkRight->at(i).shrunkPoint.coordinates[0];
                            });
                            findY(minY, bottomPointsAdded, new vector<Point>(), intermediates, shrinkRight->at(i).improvementAboveBottom + shrinkLeft->at(j-1).improvementAboveBottom, shrinkRight->at(i).improvement + shrinkLeft->at(j-1).improvement, blueInRange < redInRange, false);
                        }

                        if (lowerLevel.finished) {
                            int pointsAdded = upperLevel.pointsUnderClosest + lowerLevel.pointsUnderClosest + lowerLevel.intermediates->size() + 1 - shrinkRight->at(i).pointsRemovedFromCenter + shrinkRight->at(i).pointsRemovedFromStartingRange - shrinkLeft->at(j-1).pointsRemovedFromCenter + shrinkLeft->at(j-1).pointsRemovedFromStartingRange + topPointsAdded;
                            int pointsRemoved = shrinkRight->at(i).pointsRemovedFromStartingRange + shrinkLeft->at(j-1).pointsRemovedFromStartingRange;
                            double currentSimilarity = (double) (initialPoints - pointsRemoved) / (double) (initialPoints + pointsAdded);

                            if (currentSimilarity > similarity) {
                                currentSimilarity = similarity;
                                fairMinimum[0] = shrinkRight->at(i).shrunkPoint.coordinates[0];
                                fairMinimum[1] = lowerLevel.furthest.coordinates[1];
                                fairMaximum[0] = shrinkLeft->at(j-1).shrunkPoint.coordinates[0];
                                fairMaximum[1] = maxY;
                            }
                        }

                        if (upperLevel.finished) {
                            int pointsAdded = upperLevel.pointsUnderClosest + upperLevel.intermediates->size() + 1 + lowerLevel.pointsUnderClosest - shrinkRight->at(i).pointsRemovedFromCenter + shrinkRight->at(i).pointsRemovedFromStartingRange - shrinkLeft->at(j-1).pointsRemovedFromCenter + shrinkLeft->at(j-1).pointsRemovedFromStartingRange + bottomPointsAdded;
                            int pointsRemoved = shrinkRight->at(i).pointsRemovedFromStartingRange + shrinkLeft->at(j-1).pointsRemovedFromStartingRange;
                            double currentSimilarity = (double) (initialPoints - pointsRemoved) / (double) (initialPoints + pointsAdded);

                            if (currentSimilarity > similarity) {
                                currentSimilarity = similarity;
                                fairMinimum[0] = shrinkRight->at(i).shrunkPoint.coordinates[0];
                                fairMinimum[1] = minY;
                                fairMaximum[0] = shrinkLeft->at(j-1).shrunkPoint.coordinates[0];
                                fairMaximum[1] = upperLevel.furthest.coordinates[1];
                            }
                        }
                    }
                    
                    //expansions with one shrink on both sides
                    if (i > 0 && j > 0 && lowerLevel.finished && upperLevel.finished) {
                        int pointsAdded = upperLevel.pointsUnderClosest + upperLevel.intermediates->size() + 1 + lowerLevel.pointsUnderClosest + lowerLevel.intermediates->size() + 1 - shrinkRight->at(i-1).pointsRemovedFromCenter + shrinkRight->at(i-1).pointsRemovedFromStartingRange - shrinkLeft->at(j-1).pointsRemovedFromCenter + shrinkLeft->at(j-1).pointsRemovedFromStartingRange;
                        int pointsRemoved = shrinkRight->at(i-1).pointsRemovedFromStartingRange + shrinkLeft->at(j-1).pointsRemovedFromStartingRange;
                        double currentSimilarity = (double) (initialPoints - pointsRemoved) / (double) (initialPoints + pointsAdded);

                        if (currentSimilarity > similarity) {
                            currentSimilarity = similarity;
                            fairMinimum[0] = shrinkRight->at(i-1).shrunkPoint.coordinates[0];
                            fairMinimum[1] = lowerLevel.furthest.coordinates[1];
                            fairMaximum[0] = shrinkLeft->at(j-1).shrunkPoint.coordinates[0];
                            fairMaximum[1] = upperLevel.furthest.coordinates[1];
                        }
                    }

                    //expansions with 2 shrinks on left side

                    if (i > 1 && lowerLevel.finished && upperLevel.finished) {
                        int pointsAdded = upperLevel.pointsUnderClosest + upperLevel.intermediates->size() + 1 + lowerLevel.pointsUnderClosest + lowerLevel.intermediates->size() + 1 - shrinkRight->at(i-2).pointsRemovedFromCenter + shrinkRight->at(i-2).pointsRemovedFromStartingRange - shrinkLeft->at(j).pointsRemovedFromCenter + shrinkLeft->at(j).pointsRemovedFromStartingRange;
                        int pointsRemoved = shrinkRight->at(i-2).pointsRemovedFromStartingRange + shrinkLeft->at(j).pointsRemovedFromStartingRange;
                        double currentSimilarity = (double) (initialPoints - pointsRemoved) / (double) (initialPoints + pointsAdded);

                        if (currentSimilarity > similarity) {
                            currentSimilarity = similarity;
                            fairMinimum[0] = shrinkRight->at(i-2).shrunkPoint.coordinates[0];
                            fairMinimum[1] = lowerLevel.furthest.coordinates[1];
                            fairMaximum[0] = shrinkLeft->at(j).shrunkPoint.coordinates[0];
                            fairMaximum[1] = upperLevel.furthest.coordinates[1];
                        }
                    }
                    
                    //expansions with 2 shrinks on right side

                    if (j > 1 && lowerLevel.finished && upperLevel.finished) {
                        int pointsAdded = upperLevel.pointsUnderClosest + upperLevel.intermediates->size() + 1 + lowerLevel.pointsUnderClosest + lowerLevel.intermediates->size() + 1 - shrinkRight->at(i).pointsRemovedFromCenter + shrinkRight->at(i).pointsRemovedFromStartingRange - shrinkLeft->at(j-2).pointsRemovedFromCenter + shrinkLeft->at(j-2).pointsRemovedFromStartingRange;
                        int pointsRemoved = shrinkRight->at(i).pointsRemovedFromStartingRange + shrinkLeft->at(j-2).pointsRemovedFromStartingRange;
                        double currentSimilarity = (double) (initialPoints - pointsRemoved) / (double) (initialPoints + pointsAdded);

                        if (currentSimilarity > similarity) {
                            currentSimilarity = similarity;
                            fairMinimum[0] = shrinkRight->at(i).shrunkPoint.coordinates[0];
                            fairMinimum[1] = lowerLevel.furthest.coordinates[1];
                            fairMaximum[0] = shrinkLeft->at(j-2).shrunkPoint.coordinates[0];
                            fairMaximum[1] = upperLevel.furthest.coordinates[1];
                        }
                    }
                    i--;
                    j++;
                }
            }
        } 
    }

    for (ExpandingVerticalLevel upperLevel : *expandingUpward) {
        for (ShrinkingVerticalLevel lowerLevel : *shrinkingUpward) {
            vector<Point> *pointsRight = new vector<Point>();
            copy_if(points->begin(), points->end(), back_inserter(*pointsRight), [minimum, maximum, upperLevel, lowerLevel](Point p)
                {return p.coordinates[0] > maximum[0] && p.coordinates[1] <= upperLevel.furthest.coordinates[1] && p.coordinates[1] >= lowerLevel.furthest.coordinates[1];}
            );
            sort(pointsRight->begin(), pointsRight->end(), pointSorter(0, true));
            vector<Point> *pointsLeft = new vector<Point>();
            copy_if(points->begin(), points->end(), back_inserter(*pointsLeft), [minimum, maximum, upperLevel, lowerLevel](Point p)
                {return p.coordinates[0] < minimum[0] && p.coordinates[1] <= upperLevel.furthest.coordinates[1] && p.coordinates[1] >= lowerLevel.furthest.coordinates[1];}
            );
            sort(pointsLeft->begin(), pointsLeft->end(), pointSorter(0, false));
            vector<Point> *pointsInside = new vector<Point>();
            copy_if(points->begin(), points->end(), back_inserter(*pointsInside), [minimum, maximum, upperLevel, lowerLevel](Point p) 
                {return p.coordinates[0] <= maximum[0] && p.coordinates[0] >= minimum[0] && p.coordinates[1] <= upperLevel.furthest.coordinates[1] && p.coordinates[1] >= lowerLevel.furthest.coordinates[1];}
            );
            sort(pointsInside->begin(), pointsInside->end(), pointSorter(0, false));
            vector<Point> *pointsInsideReversed = new vector<Point>();
            copy(pointsInside->begin(), pointsInside->end(), back_inserter(*pointsInsideReversed));
            reverse(pointsInsideReversed->begin(), pointsInsideReversed->end());

            vector<ExpandingHorizontalPoint> *expandRight = expandHorizontallyTopExpanding(pointsRight, rightPoint, upperLevel, lowerLevel, redInRange > blueInRange, threshold);
            vector<ExpandingHorizontalPoint> *expandLeft = expandHorizontallyTopExpanding(pointsLeft, leftPoint, upperLevel, lowerLevel, redInRange > blueInRange, threshold);
            vector<ShrinkingHorizontalPoint> *shrinkRight = shrinkHorizontallyTopExpanding(pointsInside, upperLevel, lowerLevel, blueInRange > redInRange, threshold, minimum[1], maximum[1]);
            vector<ShrinkingHorizontalPoint> *shrinkLeft = shrinkHorizontallyTopExpanding(pointsInsideReversed, upperLevel, lowerLevel, blueInRange > redInRange, threshold, minimum[1], maximum[1]);

            int i = 0;
            while (i < expandRight->size() - 1 && expandRight->at(i).improvement + upperLevel.improvementUnderClosest + lowerLevel.improvementAtFurthest < threshold) {
                i++;
            }
            int j = 0;
            while (j < expandLeft->size() && expandLeft->at(j).improvement + expandRight->at(i).improvement + upperLevel.improvementUnderClosest + lowerLevel.improvementAtFurthest < threshold) {
                j++;
            }
           
            while (i >= 0 && j < expandLeft->size()) {

                //expansions with 0 shrinks

                vector<Point> *upperSubrange = new vector<Point>();
                copy(expandRight->at(i).upperSubrange->begin(), expandRight->at(i).upperSubrange->end(), back_inserter(*upperSubrange));
                copy(expandLeft->at(j).upperSubrange->begin(), expandLeft->at(i).upperSubrange->end(), back_inserter(*upperSubrange));
                double maxY;
                int topPointsAdded = 0;
                if (upperSubrange->empty()) {
                    maxY = upperLevel.closest.coordinates[1];
                } else {
                    sort(upperLevel.intermediates->begin(), upperLevel.intermediates->end(), pointSorter(1, true));
                    sort(upperSubrange->begin(), upperSubrange->end(), pointSorter(1, true));
                    findY(maxY, topPointsAdded, upperSubrange, upperLevel.intermediates, expandRight->at(i).improvementUnderTop + expandLeft->at(j).improvementUnderTop, expandRight->at(i).improvement + expandLeft->at(j).improvement, blueInRange < redInRange, true);
                }

                double minY;
                vector<Point> *pointsInRange = new vector<Point>();
                copy_if(points->begin(), points->end(), back_inserter(pointsInRange), [expandLeft, expandRight, i, j, lowerLevel, maxY](Point pnt){
                    return (pnt.coordinates[0] >= expandLeft->at(j).expandedPoint.coordinates[0] && pnt.coordinates[1] > lowerLevel.furthest.coordinates[1] && pnt.coordinates[0] <= expandRight->at(i).expandedPoint.coordinates[0] && pnt.coordinates[1] <= maxY);
                });
                minY = pointsInRange->at(0).coordinates[1];
                for (Point pnt : *pointsInRange) {
                    if (pnt.coordinates[1] <= minY) {
                        minY = pnt.coordinates[1];
                    }
                }

                int pointsAdded = upperLevel.pointsUnderClosest + expandRight->at(i).pointsInCenter + expandLeft->at(j).pointsInCenter + topPointsAdded + expandRight->at(i).lowerSubrange->size() + expandLeft->at(j).lowerSubrange->size();
                int pointsRemoved = lowerLevel.pointsRemovedAtFurthest;
                double currentSimilarity = (double) (pointsAdded - pointsRemoved) / (double) (initialPoints + pointsAdded) ;

                if (currentSimilarity > similarity) {
                    currentSimilarity = similarity;
                    fairMinimum[0] = expandLeft->at(j).expandedPoint.coordinates[0];
                    fairMinimum[1] = minY;
                    fairMaximum[0] = expandRight->at(i).expandedPoint.coordinates[0];
                    fairMaximum[1] = maxY;
                }

                //expansions with 1 shrink on right

                if (i > 0 && upperLevel.finished) {
                    double maxY = upperLevel.furthest.coordinates[1];
                    int bottomPointsAdded = 0;
                    int topPointsAdded = upperLevel.intermediates->size() + 1 + expandRight->at(i-1).upperSubrange->size() + expandLeft->at(j).upperSubrange->size();
                    
                    double minY;
                    vector<Point> *pointsInRange = new vector<Point>();
                    copy_if(points->begin(), points->end(), back_inserter(pointsInRange), [expandLeft, expandRight, i, j, lowerLevel, maxY](Point pnt){
                        return (pnt.coordinates[0] >= expandLeft->at(j).expandedPoint.coordinates[0] && pnt.coordinates[1] > lowerLevel.furthest.coordinates[1] && pnt.coordinates[0] <= expandRight->at(i-1).expandedPoint.coordinates[0] && pnt.coordinates[1] <= maxY);
                    });
                    minY = pointsInRange->at(0).coordinates[1];
                    for (Point pnt : *pointsInRange) {
                        if (pnt.coordinates[1] <= minY) {
                            minY = pnt.coordinates[1];
                        }
                    }
                    int pointsAdded = upperLevel.pointsUnderClosest + expandRight->at(i-1).pointsInCenter + expandLeft->at(j).pointsInCenter + topPointsAdded + expandRight->at(i-1).lowerSubrange->size() + expandLeft->at(j).lowerSubrange->size();
                    int pointsRemoved = lowerLevel.pointsRemovedAtFurthest;
                    double currentSimilarity = (double) (pointsAdded - pointsRemoved) / (double) (initialPoints + pointsAdded) ;

                    if (currentSimilarity > similarity) {
                        currentSimilarity = similarity;
                        fairMinimum[0] = expandLeft->at(j).expandedPoint.coordinates[0];
                        fairMinimum[1] = minY;
                        fairMaximum[0] = expandRight->at(i-1).expandedPoint.coordinates[0];
                        fairMaximum[1] = maxY;
                    }
                }                

                //expansions with 1 shrink on left

                if (j > 0 && upperLevel.finished) {
                    double maxY = upperLevel.furthest.coordinates[1];
                    int bottomPointsAdded = 0;
                    int topPointsAdded = upperLevel.intermediates->size() + 1 + expandRight->at(i).upperSubrange->size() + expandLeft->at(j-1).upperSubrange->size();
               
                    double minY;
                    vector<Point> *pointsInRange = new vector<Point>();
                    copy_if(points->begin(), points->end(), back_inserter(pointsInRange), [expandLeft, expandRight, i, j, lowerLevel, maxY](Point pnt){
                        return (pnt.coordinates[0] >= expandLeft->at(j-1).expandedPoint.coordinates[0] && pnt.coordinates[1] > lowerLevel.furthest.coordinates[1] && pnt.coordinates[0] <= expandRight->at(i).expandedPoint.coordinates[0] && pnt.coordinates[1] <= maxY);
                    });
                    minY = pointsInRange->at(0).coordinates[1];
                    for (Point pnt : *pointsInRange) {
                        if (pnt.coordinates[1] <= minY) {
                            minY = pnt.coordinates[1];
                        }
                    }
                    int pointsAdded = upperLevel.pointsUnderClosest + expandRight->at(i).pointsInCenter + expandLeft->at(j-1).pointsInCenter + topPointsAdded + expandRight->at(i).lowerSubrange->size() + expandLeft->at(j-1).lowerSubrange->size();
                    int pointsRemoved = lowerLevel.pointsRemovedAtFurthest;
                    double currentSimilarity = (double) (pointsAdded - pointsRemoved) / (double) (initialPoints + pointsAdded) ;

                    if (currentSimilarity > similarity) {
                        currentSimilarity = similarity;
                        fairMinimum[0] = expandLeft->at(j-1).expandedPoint.coordinates[0];
                        fairMinimum[1] = minY;
                        fairMaximum[0] = expandRight->at(i).expandedPoint.coordinates[0];
                        fairMaximum[1] = maxY;
                    }               
                }
                i--;
                j++;
            }

            int i = 0;
            while (i < expandRight->size() - 1 && expandRight->at(i).improvement + upperLevel.improvementUnderClosest + lowerLevel.improvementAtFurthest - 1 < threshold) {
                i++;
            }
            int j = 0;
            while (j < shrinkRight->size() && shrinkRight->at(j).improvement + expandRight->at(i).improvement + upperLevel.improvementUnderClosest + lowerLevel.improvementAtFurthest - 1 < threshold) {
                j++;
            }

            while (i >= 0 && j < shrinkRight->size()) {
                //expansions with 0 shrinks

                double maxY;
                int topPointsAdded = 0;
                if (expandRight->at(i).upperSubrange->empty() && shrinkRight->at(j).pointsRemovedFromTop == 0) {
                    maxY = upperLevel.closest.coordinates[1];
                } else {
                    sort(upperLevel.intermediates->begin(), upperLevel.intermediates->end(), pointSorter(1, true));
                    vector<Point> *intermediates = new vector<Point>();
                    copy_if(upperLevel.intermediates->begin(), upperLevel.intermediates->end(), back_inserter(intermediates), [shrinkRight, j](Point pnt) {
                        return pnt.coordinates[0] >= shrinkRight->at(j).shrunkPoint.coordinates[0];
                    });

                    sort(expandRight->at(i).upperSubrange->begin(), expandRight->at(i).upperSubrange->end(), pointSorter(1, true));
                    findY(maxY, topPointsAdded, expandRight->at(i).upperSubrange, intermediates, expandRight->at(i).improvementUnderTop + shrinkRight->at(j).improveUnderTop, expandRight->at(i).improvement + shrinkRight->at(j).improvement, blueInRange < redInRange, true);
                }

                double minY;
                vector<Point> *pointsInRange = new vector<Point>();
                copy_if(points->begin(), points->end(), back_inserter(pointsInRange), [expandRight, shrinkRight, i, j, lowerLevel, maxY](Point pnt){
                    return (pnt.coordinates[0] >= shrinkRight->at(j).shrunkPoint.coordinates[0] && pnt.coordinates[1] > lowerLevel.furthest.coordinates[1] && pnt.coordinates[0] <= expandRight->at(i).expandedPoint.coordinates[0] && pnt.coordinates[1] <= maxY);
                });
                minY = pointsInRange->at(0).coordinates[1];
                for (Point pnt : *pointsInRange) {
                    if (pnt.coordinates[1] <= minY) {
                        minY = pnt.coordinates[1];
                    }
                }

                int pointsAdded = upperLevel.pointsUnderClosest - shrinkRight->at(j).pointsRemovedFromCenter + shrinkRight->at(j).pointsRemovedFromStartingRange  + expandRight->at(i).pointsInCenter + topPointsAdded + expandRight->at(i).lowerSubrange->size();            int pointsRemoved = lowerLevel.pointsRemovedAtFurthest + shrinkRight->at(j).pointsRemovedFromStartingRange;
                double currentSimilarity = (double) (pointsAdded - pointsRemoved) / (double) (initialPoints + pointsAdded) ;

                if (currentSimilarity > similarity) {
                    currentSimilarity = similarity;
                    fairMinimum[0] = shrinkRight->at(j).shrunkPoint.coordinates[0];
                    fairMinimum[1] = minY;
                    fairMaximum[0] = expandRight->at(i).expandedPoint.coordinates[0];
                    fairMaximum[1] = maxY;
                }

                //expansions with 1 shrink on right

                if (i > 0 && upperLevel.finished) {
                    double maxY = upperLevel.furthest.coordinates[1];
                    int topPointsAdded = upperLevel.intermediates->size() + 1 + expandRight->at(i-1).upperSubrange->size() - shrinkRight->at(j).pointsRemovedFromTop;

                    double minY;
                    vector<Point> *pointsInRange = new vector<Point>();
                    copy_if(points->begin(), points->end(), back_inserter(pointsInRange), [expandRight, shrinkRight, i, j, lowerLevel, maxY](Point pnt){
                        return (pnt.coordinates[0] >= shrinkRight->at(j).shrunkPoint.coordinates[0] && pnt.coordinates[1] > lowerLevel.furthest.coordinates[1] && pnt.coordinates[0] <= expandRight->at(i-1).expandedPoint.coordinates[0] && pnt.coordinates[1] <= maxY);
                    });
                    minY = pointsInRange->at(0).coordinates[1];
                    for (Point pnt : *pointsInRange) {
                        if (pnt.coordinates[1] <= minY) {
                            minY = pnt.coordinates[1];
                        }
                    }

                    int pointsAdded = upperLevel.pointsUnderClosest - shrinkRight->at(j).pointsRemovedFromCenter + shrinkRight->at(j).pointsRemovedFromStartingRange  + expandRight->at(i-1).pointsInCenter + topPointsAdded + expandRight->at(i-1).lowerSubrange->size();         
                    int pointsRemoved = lowerLevel.pointsRemovedAtFurthest + shrinkRight->at(j).pointsRemovedFromStartingRange;
                    double currentSimilarity = (double) (pointsAdded - pointsRemoved) / (double) (initialPoints + pointsAdded) ;

                    if (currentSimilarity > similarity) {
                        currentSimilarity = similarity;
                        fairMinimum[0] = shrinkRight->at(j).shrunkPoint.coordinates[0];
                        fairMinimum[1] = minY;
                        fairMaximum[0] = expandRight->at(i-1).expandedPoint.coordinates[0];
                        fairMaximum[1] = maxY;
                    }
                }                

                //expansions with 1 shrink on left

                if (j > 0 && upperLevel.finished) {
                    double maxY = upperLevel.furthest.coordinates[1];
                    int topPointsAdded = upperLevel.intermediates->size() + 1 + expandRight->at(i).upperSubrange->size() - shrinkRight->at(j-1).pointsRemovedFromTop;

                    double minY;
                    vector<Point> *pointsInRange = new vector<Point>();
                    copy_if(points->begin(), points->end(), back_inserter(pointsInRange), [expandRight, shrinkRight, i, j, lowerLevel, maxY](Point pnt){
                        return (pnt.coordinates[0] >= shrinkRight->at(j-1).shrunkPoint.coordinates[0] && pnt.coordinates[1] > lowerLevel.furthest.coordinates[1] && pnt.coordinates[0] <= expandRight->at(i).expandedPoint.coordinates[0] && pnt.coordinates[1] <= maxY);
                    });
                    minY = pointsInRange->at(0).coordinates[1];
                    for (Point pnt : *pointsInRange) {
                        if (pnt.coordinates[1] <= minY) {
                            minY = pnt.coordinates[1];
                        }
                    }

                    int pointsAdded = upperLevel.pointsUnderClosest - shrinkRight->at(j-1).pointsRemovedFromCenter + shrinkRight->at(j-1).pointsRemovedFromStartingRange  + expandRight->at(i).pointsInCenter + topPointsAdded + expandRight->at(i).lowerSubrange->size();
                    int pointsRemoved = lowerLevel.pointsRemovedAtFurthest + shrinkRight->at(j-1).pointsRemovedFromStartingRange;
                    double currentSimilarity = (double) (pointsAdded - pointsRemoved) / (double) (initialPoints + pointsAdded) ;

                    if (currentSimilarity > similarity) {
                        currentSimilarity = similarity;
                        fairMinimum[0] = shrinkRight->at(j-1).shrunkPoint.coordinates[0];
                        fairMinimum[1] = minY;
                        fairMaximum[0] = expandRight->at(i).expandedPoint.coordinates[0];
                        fairMaximum[1] = maxY;
                    }
                }
                i--;
                j++;
           }

            int i = 0;
            while (i < expandLeft->size() - 1 && expandLeft->at(i).improvement + upperLevel.improvementUnderClosest + lowerLevel.improvementAtFurthest - 1 < threshold) {
                i++;
            }
            int j = 0;
            while (j < shrinkLeft->size() && shrinkLeft->at(j).improvement + expandLeft->at(i).improvement + upperLevel.improvementUnderClosest + lowerLevel.improvementAtFurthest - 1 < threshold) {
                j++;
            }

            while (i >= 0 && j < shrinkLeft->size()) {
                //expansions with 0 shrinks

                double maxY;
                int topPointsAdded = 0;
                if (expandLeft->at(i).upperSubrange->empty() && shrinkLeft->at(j).pointsRemovedFromTop == 0) {
                    maxY = upperLevel.closest.coordinates[1];
                } else {
                    sort(upperLevel.intermediates->begin(), upperLevel.intermediates->end(), pointSorter(1, true));
                    vector<Point> *intermediates = new vector<Point>();
                    copy_if(upperLevel.intermediates->begin(), upperLevel.intermediates->end(), back_inserter(intermediates), [shrinkLeft, j](Point pnt) {
                        return pnt.coordinates[0] <= shrinkLeft->at(j).shrunkPoint.coordinates[0];
                    });
                    sort(expandLeft->at(i).upperSubrange->begin(), expandLeft->at(i).upperSubrange->end(), pointSorter(1, true));
                    findY(maxY, topPointsAdded, expandLeft->at(i).upperSubrange, intermediates, expandLeft->at(i).improvementUnderTop + shrinkLeft->at(j).improveUnderTop, expandLeft->at(i).improvement + shrinkLeft->at(j).improvement, blueInRange < redInRange, true);
                }

                double minY;
                vector<Point> *pointsInRange = new vector<Point>();
                copy_if(points->begin(), points->end(), back_inserter(pointsInRange), [expandLeft, shrinkLeft, i, j, lowerLevel, maxY](Point pnt){
                    return (pnt.coordinates[0] >= expandLeft->at(i).expandedPoint.coordinates[0] && pnt.coordinates[1] > lowerLevel.furthest.coordinates[1] && pnt.coordinates[0] <= shrinkLeft->at(j).shrunkPoint.coordinates[0] && pnt.coordinates[1] <= maxY);
                });
                minY = pointsInRange->at(0).coordinates[1];
                for (Point pnt : *pointsInRange) {
                    if (pnt.coordinates[1] < minY) {
                        minY = pnt.coordinates[1];
                    }
                }

                int pointsAdded = upperLevel.pointsUnderClosest - shrinkLeft->at(j).pointsRemovedFromCenter + shrinkLeft->at(j).pointsRemovedFromStartingRange  + expandLeft->at(i).pointsInCenter + topPointsAdded + expandLeft->at(i).lowerSubrange->size();
                int pointsRemoved = lowerLevel.pointsRemovedAtFurthest + shrinkLeft->at(j).pointsRemovedFromStartingRange;
                double currentSimilarity = (double) (pointsAdded - pointsRemoved) / (double) (initialPoints + pointsAdded) ;

                if (currentSimilarity > similarity) {
                    currentSimilarity = similarity;
                    fairMinimum[0] = expandLeft->at(i).expandedPoint.coordinates[0];
                    fairMinimum[1] = minY;
                    fairMaximum[0] = shrinkLeft->at(j).shrunkPoint.coordinates[0];
                    fairMaximum[1] = maxY;
                }

                //expansions with 1 shrink on right

                if (i > 0 && upperLevel.finished) {
                    double maxY = upperLevel.furthest.coordinates[1];
                    int topPointsAdded = upperLevel.intermediates->size() + 1 + expandLeft->at(i-1).upperSubrange->size() - shrinkLeft->at(j).pointsRemovedFromTop;

                    double minY;
                    vector<Point> *pointsInRange = new vector<Point>();
                    copy_if(points->begin(), points->end(), back_inserter(pointsInRange), [expandLeft, shrinkLeft, i, j, lowerLevel, maxY](Point pnt){
                        return (pnt.coordinates[0] >= expandLeft->at(i-1).expandedPoint.coordinates[0] && pnt.coordinates[1] > lowerLevel.furthest.coordinates[1] && pnt.coordinates[0] <= shrinkLeft->at(j).shrunkPoint.coordinates[0] && pnt.coordinates[1] <= maxY);
                    });
                    minY = pointsInRange->at(0).coordinates[1];
                    for (Point pnt : *pointsInRange) {
                        if (pnt.coordinates[1] < minY) {
                            minY = pnt.coordinates[1];
                        }
                    }

                    int pointsAdded = upperLevel.pointsUnderClosest - shrinkLeft->at(j).pointsRemovedFromCenter + shrinkLeft->at(j).pointsRemovedFromStartingRange + expandLeft->at(i-1).pointsInCenter + topPointsAdded + expandLeft->at(i-1).lowerSubrange->size();
                    int pointsRemoved = lowerLevel.pointsRemovedAtFurthest + shrinkLeft->at(j).pointsRemovedFromStartingRange;
                    double currentSimilarity = (double) (pointsAdded - pointsRemoved) / (double) (initialPoints + pointsAdded) ;

                    if (currentSimilarity > similarity) {
                        currentSimilarity = similarity;
                        fairMinimum[0] = expandLeft->at(i-1).expandedPoint.coordinates[0];
                        fairMinimum[1] = minY;
                        fairMaximum[0] = shrinkLeft->at(j).shrunkPoint.coordinates[0];
                        fairMaximum[1] = maxY;
                    }

                }                

                //expansions with 1 shrink on left

                if (j > 0 && upperLevel.finished) {
                    double maxY = upperLevel.furthest.coordinates[1];
                    int topPointsAdded = upperLevel.intermediates->size() + 1 + expandLeft->at(i).upperSubrange->size() - shrinkLeft->at(j-1).pointsRemovedFromTop;

                    double minY;
                    vector<Point> *pointsInRange = new vector<Point>();
                    copy_if(points->begin(), points->end(), back_inserter(pointsInRange), [expandLeft, shrinkLeft, i, j, lowerLevel, maxY](Point pnt){
                        return (pnt.coordinates[0] >= expandLeft->at(i).expandedPoint.coordinates[0] && pnt.coordinates[1] > lowerLevel.furthest.coordinates[1] && pnt.coordinates[0] <= shrinkLeft->at(j-1).shrunkPoint.coordinates[0] && pnt.coordinates[1] <= maxY);
                    });
                    minY = pointsInRange->at(0).coordinates[1];
                    for (Point pnt : *pointsInRange) {
                        if (pnt.coordinates[1] < minY) {
                            minY = pnt.coordinates[1];
                        }
                    }

                    int pointsAdded = upperLevel.pointsUnderClosest - shrinkLeft->at(j-1).pointsRemovedFromCenter + shrinkLeft->at(j-1).pointsRemovedFromStartingRange  + expandLeft->at(i).pointsInCenter + topPointsAdded + expandLeft->at(i).lowerSubrange->size();
                    int pointsRemoved = lowerLevel.pointsRemovedAtFurthest + shrinkLeft->at(j-1).pointsRemovedFromStartingRange;
                    double currentSimilarity = (double) (pointsAdded - pointsRemoved) / (double) (initialPoints + pointsAdded) ;

                    if (currentSimilarity > similarity) {
                        currentSimilarity = similarity;
                        fairMinimum[0] = expandLeft->at(i).expandedPoint.coordinates[0];
                        fairMinimum[1] = minY;
                        fairMaximum[0] = shrinkLeft->at(j-1).shrunkPoint.coordinates[0];
                        fairMaximum[1] = maxY;
                    }
                }
                i--;
                j++;
            }

            int i = 0;
            while (i < shrinkRight->size() - 1 && shrinkRight->at(i).improvement + upperLevel.improvementUnderClosest + lowerLevel.improvementAtFurthest - 1 < threshold) {
                i++;
            }
            int j = 0;
            while (j < shrinkLeft->size() && shrinkLeft->at(j).improvement + expandLeft->at(i).improvement + upperLevel.improvementUnderClosest + lowerLevel.improvementAtFurthest - 1 < threshold) {
                j++;
            }

            while (i >= 0 && j < shrinkLeft->size()) {
                //expansions with 0 shrinks

                double maxY;
                int topPointsAdded = 0;
                if (shrinkRight->at(i).pointsRemovedFromTop && shrinkLeft->at(j).pointsRemovedFromTop == 0) {
                    maxY = upperLevel.closest.coordinates[1];
                } else {
                    sort(upperLevel.intermediates->begin(), upperLevel.intermediates->end(), pointSorter(1, true));
                    vector<Point> *intermediates = new vector<Point>();
                    copy_if(upperLevel.intermediates->begin(), upperLevel.intermediates->end(), back_inserter(intermediates), [shrinkLeft, j](Point pnt) {
                        return pnt.coordinates[0] <= shrinkLeft->at(j).shrunkPoint.coordinates[0];
                    });
                    findY(maxY, topPointsAdded, new vector<Point>(), intermediates, shrinkRight->at(i).improveUnderTop + shrinkLeft->at(j).improveUnderTop, shrinkRight->at(i).improvement + shrinkLeft->at(j).improvement, blueInRange < redInRange, true);
                }

                double minY;
                vector<Point> *pointsInRange = new vector<Point>();
                copy_if(points->begin(), points->end(), back_inserter(pointsInRange), [shrinkRight, shrinkLeft, i, j, lowerLevel, maxY](Point pnt){
                    return (pnt.coordinates[0] >= shrinkRight->at(i).shrunkPoint.coordinates[0] && pnt.coordinates[1] > lowerLevel.furthest.coordinates[1] && pnt.coordinates[0] <= shrinkLeft->at(j).shrunkPoint.coordinates[0] && pnt.coordinates[1] <= maxY);
                });
                minY = pointsInRange->at(0).coordinates[1];
                for (Point pnt : *pointsInRange) {
                    if (pnt.coordinates[1] < minY) {
                        minY = pnt.coordinates[1];
                    }
                }

                int pointsAdded = upperLevel.pointsUnderClosest - shrinkLeft->at(j).pointsRemovedFromCenter + shrinkLeft->at(j).pointsRemovedFromStartingRange - shrinkRight->at(i).pointsRemovedFromCenter + shrinkRight->at(i).pointsRemovedFromStartingRange + topPointsAdded;
                int pointsRemoved = lowerLevel.pointsRemovedAtFurthest + shrinkLeft->at(j).pointsRemovedFromStartingRange + shrinkRight->at(i).pointsRemovedFromStartingRange;
                double currentSimilarity = (double) (pointsAdded - pointsRemoved) / (double) (initialPoints + pointsAdded) ;

                if (currentSimilarity > similarity) {
                    currentSimilarity = similarity;
                    fairMinimum[0] = shrinkRight->at(i).shrunkPoint.coordinates[0];
                    fairMinimum[1] = minY;
                    fairMaximum[0] = shrinkLeft->at(j).shrunkPoint.coordinates[0];
                    fairMaximum[1] = maxY;
                }


                //expansions with 1 shrink on left

                if (i > 0 && upperLevel.finished) {
                    double maxY = upperLevel.furthest.coordinates[1];
                    int topPointsAdded = upperLevel.intermediates->size() + 1 - shrinkRight->at(i-1).pointsRemovedFromTop - shrinkLeft->at(j).pointsRemovedFromTop;

                    double minY;
                    vector<Point> *pointsInRange = new vector<Point>();
                    copy_if(points->begin(), points->end(), back_inserter(pointsInRange), [shrinkRight, shrinkLeft, i, j, lowerLevel, maxY](Point pnt){
                        return (pnt.coordinates[0] >= shrinkRight->at(i-1).shrunkPoint.coordinates[0] && pnt.coordinates[1] > lowerLevel.furthest.coordinates[1] && pnt.coordinates[0] <= shrinkLeft->at(j).shrunkPoint.coordinates[0] && pnt.coordinates[1] <= maxY);
                    });
                    minY = pointsInRange->at(0).coordinates[1];
                    for (Point pnt : *pointsInRange) {
                        if (pnt.coordinates[1] < minY) {
                            minY = pnt.coordinates[1];
                        }
                    }

                    int pointsAdded = upperLevel.pointsUnderClosest - shrinkLeft->at(j).pointsRemovedFromCenter + shrinkLeft->at(j).pointsRemovedFromStartingRange - shrinkRight->at(i-1).pointsRemovedFromCenter + shrinkRight->at(i-1).pointsRemovedFromStartingRange + topPointsAdded;
                    int pointsRemoved = lowerLevel.pointsRemovedAtFurthest + shrinkLeft->at(j).pointsRemovedFromStartingRange + shrinkRight->at(i-1).pointsRemovedFromStartingRange;
                    double currentSimilarity = (double) (pointsAdded - pointsRemoved) / (double) (initialPoints + pointsAdded) ;

                    if (currentSimilarity > similarity) {
                        currentSimilarity = similarity;
                        fairMinimum[0] = shrinkRight->at(i-1).shrunkPoint.coordinates[0];
                        fairMinimum[1] = minY;
                        fairMaximum[0] = shrinkLeft->at(j).shrunkPoint.coordinates[0];
                        fairMaximum[1] = maxY;
                    }
                }                

                //expansions with 1 shrink on right

                if (j > 0 && upperLevel.finished) {
                    double maxY = upperLevel.furthest.coordinates[1];
                    int topPointsAdded = upperLevel.intermediates->size() + 1 - shrinkRight->at(i).pointsRemovedFromTop - shrinkLeft->at(j).pointsRemovedFromTop;

                    double minY;
                    vector<Point> *pointsInRange = new vector<Point>();
                    copy_if(points->begin(), points->end(), back_inserter(pointsInRange), [shrinkRight, shrinkLeft, i, j, lowerLevel, maxY](Point pnt){
                        return (pnt.coordinates[0] >= shrinkRight->at(i).shrunkPoint.coordinates[0] && pnt.coordinates[1] > lowerLevel.furthest.coordinates[1] && pnt.coordinates[0] <= shrinkLeft->at(j-1).shrunkPoint.coordinates[0] && pnt.coordinates[1] <= maxY);
                    });
                    minY = pointsInRange->at(0).coordinates[1];
                    for (Point pnt : *pointsInRange) {
                        if (pnt.coordinates[1] < minY) {
                            minY = pnt.coordinates[1];
                        }
                    }

                    int pointsAdded = upperLevel.pointsUnderClosest - shrinkLeft->at(j-1).pointsRemovedFromCenter + shrinkLeft->at(j-1).pointsRemovedFromStartingRange - shrinkRight->at(i).pointsRemovedFromCenter + shrinkRight->at(i).pointsRemovedFromStartingRange + topPointsAdded;
                    int pointsRemoved = lowerLevel.pointsRemovedAtFurthest + shrinkLeft->at(j-1).pointsRemovedFromStartingRange + shrinkRight->at(i).pointsRemovedFromStartingRange;
                    double currentSimilarity = (double) (pointsAdded - pointsRemoved) / (double) (initialPoints + pointsAdded) ;

                    if (currentSimilarity > similarity) {
                        currentSimilarity = similarity;
                        fairMinimum[0] = shrinkRight->at(i).shrunkPoint.coordinates[0];
                        fairMinimum[1] = minY;
                        fairMaximum[0] = shrinkLeft->at(j-1).shrunkPoint.coordinates[0];
                        fairMaximum[1] = maxY;
                    }
                }
                i--;
                j++;
            }
        }
    }

    for (ExpandingVerticalLevel lowerLevel : *expandingDownward) {
        for (ShrinkingVerticalLevel upperLevel : *shrinkingDownward) {
            vector<Point> *pointsRight = new vector<Point>();
            copy_if(points->begin(), points->end(), back_inserter(*pointsRight), [minimum, maximum, upperLevel, lowerLevel](Point p)
                {return p.coordinates[0] > maximum[0] && p.coordinates[1] <= upperLevel.furthest.coordinates[1] && p.coordinates[1] >= lowerLevel.furthest.coordinates[1];}
            );
            sort(pointsRight->begin(), pointsRight->end(), pointSorter(0, true));
            vector<Point> *pointsLeft = new vector<Point>();
            copy_if(points->begin(), points->end(), back_inserter(*pointsLeft), [minimum, maximum, upperLevel, lowerLevel](Point p)
                {return p.coordinates[0] < minimum[0] && p.coordinates[1] <= upperLevel.furthest.coordinates[1] && p.coordinates[1] >= lowerLevel.furthest.coordinates[1];}
            );
            sort(pointsLeft->begin(), pointsLeft->end(), pointSorter(0, false));
            vector<Point> *pointsInside = new vector<Point>();
            copy_if(points->begin(), points->end(), back_inserter(*pointsInside), [minimum, maximum, upperLevel, lowerLevel](Point p) 
                {return p.coordinates[0] <= maximum[0] && p.coordinates[0] >= minimum[0] && p.coordinates[1] <= upperLevel.furthest.coordinates[1] && p.coordinates[1] >= lowerLevel.furthest.coordinates[1];}
            );
            sort(pointsInside->begin(), pointsInside->end(), pointSorter(0, false));
            vector<Point> *pointsInsideReversed = new vector<Point>();
            copy(pointsInside->begin(), pointsInside->end(), back_inserter(*pointsInsideReversed));
            reverse(pointsInsideReversed->begin(), pointsInsideReversed->end());

            vector<ExpandingHorizontalPoint> *expandRight = expandHorizontallyBottomExpanding(pointsRight, rightPoint, lowerLevel, upperLevel, redInRange > blueInRange, threshold);
            vector<ExpandingHorizontalPoint> *expandLeft = expandHorizontallyBottomExpanding(pointsLeft, leftPoint, lowerLevel,  upperLevel, redInRange > blueInRange, threshold);
            vector<ShrinkingHorizontalPoint> *shrinkRight = shrinkHorizontallyBottomExpanding(pointsInside, lowerLevel, upperLevel, blueInRange > redInRange, threshold, minimum[1], maximum[1]);
            vector<ShrinkingHorizontalPoint> *shrinkLeft = shrinkHorizontallyBottomExpanding(pointsInsideReversed, lowerLevel, upperLevel, blueInRange > redInRange, threshold, minimum[1], maximum[1]);

            int i = 0;
            while (i < expandRight->size() - 1 && expandRight->at(i).improvement + upperLevel.improvementAtFurthest + lowerLevel.improvementUnderClosest -1 < threshold) {
                i++;
            }
            int j = 0;
            while (j < expandLeft->size() && expandLeft->at(j).improvement + expandRight->at(i).improvement + upperLevel.improvementAtFurthest + lowerLevel.improvementUnderClosest - 1 < threshold) {
                j++;
            }
           
            while (i >= 0 && j < expandLeft->size()) {
                //expansions with 0 shrinks

                vector<Point> *lowerSubrange = new vector<Point>();
                copy(expandRight->at(i).lowerSubrange->begin(), expandRight->at(i).lowerSubrange->end(), back_inserter(*lowerSubrange));
                copy(expandLeft->at(j).lowerSubrange->begin(), expandLeft->at(i).lowerSubrange->end(), back_inserter(*lowerSubrange));
                double minY;
                int bottomPointsAdded = 0;
                if (lowerSubrange->empty()) {
                    minY = lowerLevel.closest.coordinates[1];
                } else {
                    sort(lowerLevel.intermediates->begin(), lowerLevel.intermediates->end(), pointSorter(1, false));
                    sort(lowerSubrange->begin(),     lowerSubrange->end(), pointSorter(1, false));
                    findY(minY, bottomPointsAdded, lowerSubrange, lowerLevel.intermediates, expandRight->at(i).improvementAboveBottom + expandLeft->at(j).improvementAboveBottom, expandRight->at(i).improvement + expandLeft->at(j).improvement, blueInRange < redInRange, false);
                }

                double maxY;
                vector<Point> *pointsInRange = new vector<Point>();
                copy_if(points->begin(), points->end(), back_inserter(pointsInRange), [expandLeft, expandRight, i, j, upperLevel, minY](Point pnt){
                    return (pnt.coordinates[0] >= expandLeft->at(j).expandedPoint.coordinates[0] && pnt.coordinates[1] >= minY && pnt.coordinates[0] <= expandRight->at(i).expandedPoint.coordinates[0] && pnt.coordinates[1] <upperLevel.furthest.coordinates[1]);
                });
                maxY = pointsInRange->at(0).coordinates[1];
                for (Point pnt : *pointsInRange) {
                    if (pnt.coordinates[1] >= maxY) {
                        maxY = pnt.coordinates[1];
                    }
                }

                int pointsAdded = expandRight->at(i).pointsInCenter + expandLeft->at(j).pointsInCenter + bottomPointsAdded + expandRight->at(i).upperSubrange->size() + expandLeft->at(j).upperSubrange->size();
                int pointsRemoved = upperLevel.pointsRemovedAtFurthest;
                double currentSimilarity = (double) (pointsAdded - pointsRemoved) / (double) (initialPoints + pointsAdded);

                if (currentSimilarity > similarity) {
                    currentSimilarity = similarity;
                    fairMinimum[0] = expandLeft->at(j).expandedPoint.coordinates[0];
                    fairMinimum[1] = minY;
                    fairMaximum[0] = expandRight->at(i).expandedPoint.coordinates[0];
                    fairMaximum[1] = maxY;
                }

                //expansions with 1 shrink on right

                if (i > 0 && lowerLevel.finished) {
                    double minY = lowerLevel.furthest.coordinates[1];
                    int bottomPointsAdded = lowerLevel.intermediates->size() + 1 + expandRight->at(i-1).lowerSubrange->size() + expandLeft->at(j).lowerSubrange->size();
                    double maxY;
                    vector<Point> *pointsInRange = new vector<Point>();
                    copy_if(points->begin(), points->end(), back_inserter(pointsInRange), [expandLeft, expandRight, i, j, upperLevel, minY](Point pnt){
                        return (pnt.coordinates[0] >= expandLeft->at(j).expandedPoint.coordinates[0] && pnt.coordinates[1] >= minY && pnt.coordinates[0] <= expandRight->at(i-1).expandedPoint.coordinates[0] && pnt.coordinates[1] < upperLevel.furthest.coordinates[1]);
                    });
                    maxY = pointsInRange->at(0).coordinates[1];
                    for (Point pnt : *pointsInRange) {
                        if (pnt.coordinates[1] >= maxY) {
                            maxY = pnt.coordinates[1];
                        }
                    }

                    int pointsAdded = expandRight->at(i-1).pointsInCenter + expandLeft->at(j).pointsInCenter + bottomPointsAdded + expandRight->at(i-1).upperSubrange->size() + expandLeft->at(j).upperSubrange->size();
                    int pointsRemoved = upperLevel.pointsRemovedAtFurthest;
                    double currentSimilarity = (double) (pointsAdded - pointsRemoved) / (double) (initialPoints + pointsAdded);

                    if (currentSimilarity > similarity) {
                        currentSimilarity = similarity;
                        fairMinimum[0] = expandLeft->at(j).expandedPoint.coordinates[0];
                        fairMinimum[1] = minY;
                        fairMaximum[0] = expandRight->at(i-1).expandedPoint.coordinates[0];
                        fairMaximum[1] = maxY;
                    }
                }                

                //expansions with 1 shrink on left

                if (j > 0 && lowerLevel.finished) {
                    double minY = lowerLevel.furthest.coordinates[1];
                    int bottomPointsAdded = lowerLevel.intermediates->size() + 1 + expandRight->at(i).lowerSubrange->size() + expandLeft->at(j-1).lowerSubrange->size();

                    double maxY;
                    vector<Point> *pointsInRange = new vector<Point>();
                    copy_if(points->begin(), points->end(), back_inserter(pointsInRange), [expandLeft, expandRight, i, j, upperLevel, minY](Point pnt){
                        return (pnt.coordinates[0] >= expandLeft->at(j-1).expandedPoint.coordinates[0] && pnt.coordinates[1] >= minY && pnt.coordinates[0] <= expandRight->at(i).expandedPoint.coordinates[0] && pnt.coordinates[1] <upperLevel.furthest.coordinates[1]);
                    });
                    maxY = pointsInRange->at(0).coordinates[1];
                    for (Point pnt : *pointsInRange) {
                        if (pnt.coordinates[1] >= maxY) {
                            maxY = pnt.coordinates[1];
                        }
                    }

                    int pointsAdded = expandRight->at(i).pointsInCenter + expandLeft->at(j-1).pointsInCenter + bottomPointsAdded + expandRight->at(i).upperSubrange->size() + expandLeft->at(j-1).upperSubrange->size();
                    int pointsRemoved = upperLevel.pointsRemovedAtFurthest;
                    double currentSimilarity = (double) (pointsAdded - pointsRemoved) / (double) (initialPoints + pointsAdded);

                    if (currentSimilarity > similarity) {
                        currentSimilarity = similarity;
                        fairMinimum[0] = expandLeft->at(j-1).expandedPoint.coordinates[0];
                        fairMinimum[1] = minY;
                        fairMaximum[0] = expandRight->at(i).expandedPoint.coordinates[0];
                        fairMaximum[1] = maxY;
                    }
                }
                i--;
                j++;
            }

            int i = 0;
            while (i < expandRight->size() - 1 && expandRight->at(i).improvement + upperLevel.improvementAtFurthest + lowerLevel.improvementUnderClosest -1 < threshold) {
                i++;
            }
            int j = 0;
            while (j < shrinkRight->size() && shrinkRight->at(j).improvement + expandRight->at(i).improvement + upperLevel.improvementAtFurthest + lowerLevel.improvementUnderClosest - 1 < threshold) {
                j++;
            }
           
            while (i >= 0 && j < shrinkRight->size()) {
                //expansions with 0 shrinks

                double minY;
                int bottomPointsAdded = 0;
                if (expandRight->at(i).lowerSubrange->empty() && shrinkRight->at(j).pointsRemovedFromBottom == 0) {
                    minY = lowerLevel.closest.coordinates[1];
                } else {
                    sort(lowerLevel.intermediates->begin(), lowerLevel.intermediates->end(), pointSorter(1, false));
                    vector<Point> *intermediates = new vector<Point>();
                    copy_if(lowerLevel.intermediates->begin(), lowerLevel.intermediates->end(), back_inserter(intermediates), [shrinkRight, j](Point pnt) {
                        return pnt.coordinates[0] >= shrinkRight->at(j).shrunkPoint.coordinates[0];
                    });
                    sort(expandRight->at(i).lowerSubrange->begin(), expandRight->at(i).lowerSubrange->end(), pointSorter(1, false));
                    findY(minY, bottomPointsAdded, expandRight->at(i).lowerSubrange, intermediates, expandRight->at(i).improvementAboveBottom + shrinkRight->at(j).improvementAboveBottom, expandRight->at(i).improvement + shrinkRight->at(j).improvement, blueInRange < redInRange, true);
                }


                double maxY;
                vector<Point> *pointsInRange = new vector<Point>();
                copy_if(points->begin(), points->end(), back_inserter(pointsInRange), [expandRight, shrinkRight, i, j, upperLevel, minY](Point pnt){
                    return (pnt.coordinates[0] >= shrinkRight->at(j).shrunkPoint.coordinates[0] && pnt.coordinates[1] >= minY && pnt.coordinates[0] <= expandRight->at(i).expandedPoint.coordinates[0] && pnt.coordinates[1] < upperLevel.furthest.coordinates[1]);
                });
                maxY = pointsInRange->at(0).coordinates[1];
                for (Point pnt : *pointsInRange) {
                    if (pnt.coordinates[1] >= maxY) {
                        maxY = pnt.coordinates[1];
                    }
                }
                
                int pointsAdded = lowerLevel.pointsUnderClosest - shrinkRight->at(j).pointsRemovedFromCenter + shrinkRight->at(j).pointsRemovedFromStartingRange  + expandRight->at(i).pointsInCenter + bottomPointsAdded + expandRight->at(i).upperSubrange->size();
                int pointsRemoved = upperLevel.pointsRemovedAtFurthest + shrinkRight->at(j).pointsRemovedFromStartingRange;
                double currentSimilarity = (double) (pointsAdded - pointsRemoved) / (double) (initialPoints + pointsAdded) ;

                if (currentSimilarity > similarity) {
                    currentSimilarity = similarity;
                    fairMinimum[0] = shrinkRight->at(j).shrunkPoint.coordinates[0];
                    fairMinimum[1] = minY;
                    fairMaximum[0] = expandRight->at(i).expandedPoint.coordinates[0];
                    fairMaximum[1] = maxY;
                }

                //expansions with 1 shrink on right

                if (i > 0 && lowerLevel.finished) {
                    double minY = lowerLevel.furthest.coordinates[1];
                    int bottomPointsAdded = lowerLevel.intermediates->size() + 1 - shrinkRight->at(j).pointsRemovedFromBottom + expandRight->at(i-1).lowerSubrange->size();
                    double maxY;
                    vector<Point> *pointsInRange = new vector<Point>();
                    copy_if(points->begin(), points->end(), back_inserter(pointsInRange), [expandRight, shrinkRight, i, j, upperLevel, minY](Point pnt){
                        return (pnt.coordinates[0] >= shrinkRight->at(j).shrunkPoint.coordinates[0] && pnt.coordinates[1] >= minY && pnt.coordinates[0] <= expandRight->at(i-1).expandedPoint.coordinates[0] && pnt.coordinates[1] < upperLevel.furthest.coordinates[1]);
                    });
                    maxY = pointsInRange->at(0).coordinates[1];
                    for (Point pnt : *pointsInRange) {
                        if (pnt.coordinates[1] >= maxY) {
                            maxY = pnt.coordinates[1];
                        }
                    }
                    
                    int pointsAdded = lowerLevel.pointsUnderClosest - shrinkRight->at(j).pointsRemovedFromCenter + shrinkRight->at(j).pointsRemovedFromStartingRange  + expandRight->at(i-1).pointsInCenter + bottomPointsAdded + expandRight->at(i-1).upperSubrange->size();
                    int pointsRemoved = upperLevel.pointsRemovedAtFurthest + shrinkRight->at(j).pointsRemovedFromStartingRange;
                    double currentSimilarity = (double) (pointsAdded - pointsRemoved) / (double) (initialPoints + pointsAdded) ;

                    if (currentSimilarity > similarity) {
                        currentSimilarity = similarity;
                        fairMinimum[0] = shrinkRight->at(j).shrunkPoint.coordinates[0];
                        fairMinimum[1] = minY;
                        fairMaximum[0] = expandRight->at(i-1).expandedPoint.coordinates[0];
                        fairMaximum[1] = maxY;
                    }
                }                

                //expansions with 1 shrink on left

                if (j > 0 && lowerLevel.finished) {
                    double minY = lowerLevel.furthest.coordinates[1];
                    int bottomPointsAdded = lowerLevel.intermediates->size() + 1 + expandRight->at(i).upperSubrange->size() - shrinkRight->at(j-1).pointsRemovedFromBottom;
                    
                    double maxY;
                    vector<Point> *pointsInRange = new vector<Point>();
                    copy_if(points->begin(), points->end(), back_inserter(pointsInRange), [expandRight, shrinkRight, i, j, upperLevel, minY](Point pnt){
                        return (pnt.coordinates[0] >= shrinkRight->at(j-1).shrunkPoint.coordinates[0] && pnt.coordinates[1] >= minY && pnt.coordinates[0] <= expandRight->at(i).expandedPoint.coordinates[0] && pnt.coordinates[1] < upperLevel.furthest.coordinates[1]);
                    });
                    maxY = pointsInRange->at(0).coordinates[1];
                    for (Point pnt : *pointsInRange) {
                        if (pnt.coordinates[1] >= maxY) {
                            maxY = pnt.coordinates[1];
                        }
                    }
                    
                    int pointsAdded = lowerLevel.pointsUnderClosest - shrinkRight->at(j-1).pointsRemovedFromCenter + shrinkRight->at(j-1).pointsRemovedFromStartingRange  + expandRight->at(i).pointsInCenter + bottomPointsAdded + expandRight->at(i).upperSubrange->size();
                    int pointsRemoved = upperLevel.pointsRemovedAtFurthest + shrinkRight->at(j-1).pointsRemovedFromStartingRange;
                    double currentSimilarity = (double) (pointsAdded - pointsRemoved) / (double) (initialPoints + pointsAdded) ;

                    if (currentSimilarity > similarity) {
                        currentSimilarity = similarity;
                        fairMinimum[0] = shrinkRight->at(j-1).shrunkPoint.coordinates[0];
                        fairMinimum[1] = minY;
                        fairMaximum[0] = expandRight->at(i).expandedPoint.coordinates[0];
                        fairMaximum[1] = maxY;
                    }
                }
                i--;
                j++;
            }

            int i = 0;
            while (i < expandLeft->size() - 1 && expandLeft->at(i).improvement + upperLevel.improvementAtFurthest + lowerLevel.improvementUnderClosest -1 < threshold) {
                i++;
            }
            int j = 0;
            while (j < shrinkLeft->size() && shrinkLeft->at(j).improvement + expandRight->at(i).improvement + upperLevel.improvementAtFurthest + lowerLevel.improvementUnderClosest - 1 < threshold) {
                j++;
            }
           
            while (i >= 0 && j < shrinkLeft->size()) {
                //expansions with 0 shrinks

                double minY;
                int bottomPointsAdded = 0;
                if (expandLeft->at(i).lowerSubrange->empty() && shrinkLeft->at(j).pointsRemovedFromBottom == 0) {
                    minY = lowerLevel.closest.coordinates[1];
                } else {
                    sort(lowerLevel.intermediates->begin(), lowerLevel.intermediates->end(), pointSorter(1, false));
                    vector<Point> *intermediates = new vector<Point>();
                    copy_if(lowerLevel.intermediates->begin(), lowerLevel.intermediates->end(), back_inserter(intermediates), [shrinkLeft, j](Point pnt) {
                        return pnt.coordinates[0] <= shrinkLeft->at(j).shrunkPoint.coordinates[0];
                    });
                    sort(expandLeft->at(i).lowerSubrange->begin(), expandLeft->at(i).lowerSubrange->end(), pointSorter(1, false));
                    findY(minY, bottomPointsAdded, expandLeft->at(i).lowerSubrange, intermediates, expandLeft->at(i).improvementAboveBottom + shrinkLeft->at(j).improvementAboveBottom, expandLeft->at(i).improvement + shrinkLeft->at(j).improvement, blueInRange < redInRange, true);
                }

                double maxY;
                vector<Point> *pointsInRange = new vector<Point>();
                copy_if(points->begin(), points->end(), back_inserter(pointsInRange), [expandLeft, shrinkLeft, i, j, upperLevel, minY](Point pnt){
                    return (pnt.coordinates[0] >= expandLeft->at(i).expandedPoint.coordinates[0] && pnt.coordinates[1] >= minY && pnt.coordinates[0] <= shrinkLeft->at(j).shrunkPoint.coordinates[0] && pnt.coordinates[1] < upperLevel.furthest.coordinates[1]);
                });
                maxY = pointsInRange->at(0).coordinates[1];
                for (Point pnt : *pointsInRange) {
                    if (pnt.coordinates[1] >= maxY) {
                        maxY = pnt.coordinates[1];
                    }
                }
                
                int pointsAdded = lowerLevel.pointsUnderClosest - shrinkLeft->at(j).pointsRemovedFromCenter + shrinkLeft->at(j).pointsRemovedFromStartingRange  + expandLeft->at(i).pointsInCenter + bottomPointsAdded + expandLeft->at(i).upperSubrange->size();
                int pointsRemoved = upperLevel.pointsRemovedAtFurthest + shrinkLeft->at(j).pointsRemovedFromStartingRange;
                double currentSimilarity = (double) (pointsAdded - pointsRemoved) / (double) (initialPoints + pointsAdded) ;

                if (currentSimilarity > similarity) {
                    currentSimilarity = similarity;
                    fairMinimum[0] = expandLeft->at(i).expandedPoint.coordinates[0];
                    fairMinimum[1] = minY;
                    fairMaximum[0] = shrinkLeft->at(j).shrunkPoint.coordinates[0];
                    fairMaximum[1] = maxY;
                }

                //expansions with 1 shrink on right

                if (i > 0 && lowerLevel.finished) {
                    double minY = lowerLevel.furthest.coordinates[1];
                    int bottomPointsAdded = lowerLevel.intermediates->size() + 1 + expandLeft->at(i-1).lowerSubrange->size() - shrinkLeft->at(j).pointsRemovedFromBottom;

                    double maxY;
                    vector<Point> *pointsInRange = new vector<Point>();
                    copy_if(points->begin(), points->end(), back_inserter(pointsInRange), [expandLeft, shrinkLeft, i, j, upperLevel, minY](Point pnt){
                        return (pnt.coordinates[0] >= expandLeft->at(i-1).expandedPoint.coordinates[0] && pnt.coordinates[1] >= minY && pnt.coordinates[0] <= shrinkLeft->at(j).shrunkPoint.coordinates[0] && pnt.coordinates[1] < upperLevel.furthest.coordinates[1]);
                    });
                    maxY = pointsInRange->at(0).coordinates[1];
                    for (Point pnt : *pointsInRange) {
                        if (pnt.coordinates[1] >= maxY) {
                            maxY = pnt.coordinates[1];
                        }
                    }
                    
                    int pointsAdded = lowerLevel.pointsUnderClosest - shrinkLeft->at(j).pointsRemovedFromCenter + shrinkLeft->at(j).pointsRemovedFromStartingRange  + expandLeft->at(i-1).pointsInCenter + bottomPointsAdded + expandLeft->at(i-1).upperSubrange->size();
                    int pointsRemoved = upperLevel.pointsRemovedAtFurthest + shrinkLeft->at(j).pointsRemovedFromStartingRange;
                    double currentSimilarity = (double) (pointsAdded - pointsRemoved) / (double) (initialPoints + pointsAdded) ;

                    if (currentSimilarity > similarity) {
                        currentSimilarity = similarity;
                        fairMinimum[0] = expandLeft->at(i-1).expandedPoint.coordinates[0];
                        fairMinimum[1] = minY;
                        fairMaximum[0] = shrinkLeft->at(j).shrunkPoint.coordinates[0];
                        fairMaximum[1] = maxY;
                    }
                }                

                //expansions with 1 shrink on left

                if (j > 0 && lowerLevel.finished) {
                    double minY = lowerLevel.furthest.coordinates[1];
                    int bottomPointsAdded = lowerLevel.intermediates->size() + 1 + expandLeft->at(i).upperSubrange->size() - shrinkLeft->at(j-1).pointsRemovedFromBottom;
                    
                    double maxY;
                    vector<Point> *pointsInRange = new vector<Point>();
                    copy_if(points->begin(), points->end(), back_inserter(pointsInRange), [expandLeft, shrinkLeft, i, j, upperLevel, minY](Point pnt){
                        return (pnt.coordinates[0] >= expandLeft->at(i).expandedPoint.coordinates[0]  && pnt.coordinates[1] >= minY && pnt.coordinates[0] <= shrinkLeft->at(j-1).shrunkPoint.coordinates[0]&& pnt.coordinates[1] < upperLevel.furthest.coordinates[1]);
                    });
                    maxY = pointsInRange->at(0).coordinates[1];
                    for (Point pnt : *pointsInRange) {
                        if (pnt.coordinates[1] >= maxY) {
                            maxY = pnt.coordinates[1];
                        }
                    }
                    
                    int pointsAdded = lowerLevel.pointsUnderClosest - shrinkLeft->at(j-1).pointsRemovedFromCenter + shrinkLeft->at(j-1).pointsRemovedFromStartingRange  + expandLeft->at(i).pointsInCenter + bottomPointsAdded + expandLeft->at(i).upperSubrange->size();
                    int pointsRemoved = upperLevel.pointsRemovedAtFurthest + shrinkLeft->at(j-1).pointsRemovedFromStartingRange;
                    double currentSimilarity = (double) (pointsAdded - pointsRemoved) / (double) (initialPoints + pointsAdded) ;

                    if (currentSimilarity > similarity) {
                        currentSimilarity = similarity;
                        fairMinimum[0] = expandLeft->at(i).expandedPoint.coordinates[0];
                        fairMinimum[1] = minY;
                        fairMaximum[0] = shrinkLeft->at(j-1).shrunkPoint.coordinates[0];
                        fairMaximum[1] = maxY;
                    }                }
                i--;
                j++;
            }
            int i = 0;
            while (i < shrinkLeft->size() - 1 && shrinkLeft->at(i).improvement + upperLevel.improvementAtFurthest + lowerLevel.improvementUnderClosest -1 < threshold) {
                i++;
            }
            int j = 0;
            while (j < shrinkRight->size() && shrinkRight->at(j).improvement + expandLeft->at(i).improvement + upperLevel.improvementAtFurthest + lowerLevel.improvementUnderClosest - 1 < threshold) {
                j++;
            }
           
            while (i >= 0 && j < shrinkRight->size()) {
                //expansions with 0 shrinks

                double minY;
                int bottomPointsAdded = 0;
                if (shrinkLeft->at(i).pointsRemovedFromBottom == 0 && shrinkRight->at(j).pointsRemovedFromBottom == 0) {
                    minY = lowerLevel.closest.coordinates[1];
                } else {
                    sort(lowerLevel.intermediates->begin(), lowerLevel.intermediates->end(), pointSorter(1, false));
                    vector<Point> *intermediates = new vector<Point>();
                    copy_if(lowerLevel.intermediates->begin(), lowerLevel.intermediates->end(), back_inserter(intermediates), [shrinkRight, j, shrinkLeft, i](Point pnt) {
                        return pnt.coordinates[0] >= shrinkRight->at(j).shrunkPoint.coordinates[0] && pnt.coordinates[1] <= shrinkLeft->at(i).shrunkPoint.coordinates[0];
                    });
                    findY(minY, bottomPointsAdded, new vector<Point>(), intermediates, shrinkLeft->at(i).improvementAboveBottom + shrinkRight->at(j).improvementAboveBottom, shrinkLeft->at(i).improvement + shrinkRight->at(j).improvement, blueInRange < redInRange, true);
                }

                double maxY;
                vector<Point> *pointsInRange = new vector<Point>();
                copy_if(points->begin(), points->end(), back_inserter(pointsInRange), [shrinkLeft, shrinkRight, i, j, upperLevel, minY](Point pnt){
                    return (pnt.coordinates[0] >= shrinkRight->at(j).shrunkPoint.coordinates[0] && pnt.coordinates[1] >= minY && pnt.coordinates[0] <= shrinkLeft->at(i).shrunkPoint.coordinates[0] && pnt.coordinates[1] < upperLevel.furthest.coordinates[1]);
                });
                maxY = pointsInRange->at(0).coordinates[1];
                for (Point pnt : *pointsInRange) {
                    if (pnt.coordinates[1] >= maxY) {
                        maxY = pnt.coordinates[1];
                    }
                }
                
                int pointsAdded = lowerLevel.pointsUnderClosest - shrinkRight->at(j).pointsRemovedFromCenter + shrinkRight->at(j).pointsRemovedFromStartingRange - shrinkLeft->at(i).pointsRemovedFromCenter + shrinkLeft->at(i).pointsRemovedFromStartingRange + bottomPointsAdded;
                int pointsRemoved = upperLevel.pointsRemovedAtFurthest + shrinkRight->at(j).pointsRemovedFromStartingRange + shrinkLeft->at(i).pointsRemovedFromStartingRange;
                double currentSimilarity = (double) (pointsAdded - pointsRemoved) / (double) (initialPoints + pointsAdded) ;

                if (currentSimilarity > similarity) {
                    currentSimilarity = similarity;
                    fairMinimum[0] = shrinkRight->at(j).shrunkPoint.coordinates[0];
                    fairMinimum[1] = minY;
                    fairMaximum[0] = shrinkLeft->at(i).shrunkPoint.coordinates[0];
                    fairMaximum[1] = maxY;
                }

                //expansions with 1 shrink on right

                if (i > 0 && lowerLevel.finished) {
                    double minY = lowerLevel.furthest.coordinates[1];
                    int bottomPointsAdded = lowerLevel.intermediates->size() + 1 - shrinkRight->at(j).pointsRemovedFromBottom - shrinkLeft->at(i-1).pointsRemovedFromBottom;

                    double maxY;
                    vector<Point> *pointsInRange = new vector<Point>();
                    copy_if(points->begin(), points->end(), back_inserter(pointsInRange), [shrinkLeft, shrinkRight, i, j, upperLevel, minY](Point pnt){
                        return (pnt.coordinates[0] >= shrinkRight->at(j).shrunkPoint.coordinates[0] && pnt.coordinates[1] >= minY && pnt.coordinates[0] <= shrinkLeft->at(i-1).shrunkPoint.coordinates[0] && pnt.coordinates[1] < upperLevel.furthest.coordinates[1]);
                    });
                    maxY = pointsInRange->at(0).coordinates[1];
                    for (Point pnt : *pointsInRange) {
                        if (pnt.coordinates[1] >= maxY) {
                            maxY = pnt.coordinates[1];
                        }
                    }
                    
                    int pointsAdded = lowerLevel.pointsUnderClosest - shrinkRight->at(j).pointsRemovedFromCenter + shrinkRight->at(j).pointsRemovedFromStartingRange  - shrinkLeft->at(i-1).pointsRemovedFromCenter + shrinkLeft->at(i-1).pointsRemovedFromStartingRange + bottomPointsAdded;
                    int pointsRemoved = upperLevel.pointsRemovedAtFurthest + shrinkRight->at(j).pointsRemovedFromStartingRange + shrinkLeft->at(i-1).pointsRemovedFromStartingRange;
                    double currentSimilarity = (double) (pointsAdded - pointsRemoved) / (double) (initialPoints + pointsAdded) ;

                    if (currentSimilarity > similarity) {
                        currentSimilarity = similarity;
                        fairMinimum[0] = shrinkRight->at(j).shrunkPoint.coordinates[0];
                        fairMinimum[1] = minY;
                        fairMaximum[0] = shrinkLeft->at(i-1).shrunkPoint.coordinates[0];
                        fairMaximum[1] = maxY;
                    }
                }                

                //expansions with 1 shrink on left

                if (j > 0 && lowerLevel.finished) {
                    double minY = lowerLevel.furthest.coordinates[1];
                    int bottomPointsAdded = lowerLevel.intermediates->size() + 1 - shrinkRight->at(j-1).pointsRemovedFromBottom - shrinkLeft->at(i).pointsRemovedFromBottom;

                    double maxY;
                    vector<Point> *pointsInRange = new vector<Point>();
                    copy_if(points->begin(), points->end(), back_inserter(pointsInRange), [shrinkLeft, shrinkRight, i, j, upperLevel, minY](Point pnt){
                        return (pnt.coordinates[0] >= shrinkRight->at(j-1).shrunkPoint.coordinates[0] && pnt.coordinates[1] >= minY && pnt.coordinates[0] <= shrinkLeft->at(i).shrunkPoint.coordinates[0] && pnt.coordinates[1] < upperLevel.furthest.coordinates[1]);
                    });
                    maxY = pointsInRange->at(0).coordinates[1];
                    for (Point pnt : *pointsInRange) {
                        if (pnt.coordinates[1] >= maxY) {
                            maxY = pnt.coordinates[1];
                        }
                    }
                    
                    int pointsAdded = lowerLevel.pointsUnderClosest - shrinkRight->at(j-1).pointsRemovedFromCenter + shrinkRight->at(j-1).pointsRemovedFromStartingRange  - shrinkLeft->at(i).pointsRemovedFromCenter + shrinkLeft->at(i).pointsRemovedFromStartingRange + bottomPointsAdded;
                    int pointsRemoved = upperLevel.pointsRemovedAtFurthest + shrinkRight->at(j-1).pointsRemovedFromStartingRange + shrinkLeft->at(i).pointsRemovedFromStartingRange;
                    double currentSimilarity = (double) (pointsAdded - pointsRemoved) / (double) (initialPoints + pointsAdded) ;

                    if (currentSimilarity > similarity) {
                        currentSimilarity = similarity;
                        fairMinimum[0] = shrinkRight->at(j-1).shrunkPoint.coordinates[0];
                        fairMinimum[1] = minY;
                        fairMaximum[0] = shrinkLeft->at(i).shrunkPoint.coordinates[0];
                        fairMaximum[1] = maxY;
                    }
                }                
                i--;
                j++;
            }
        }
    }
}