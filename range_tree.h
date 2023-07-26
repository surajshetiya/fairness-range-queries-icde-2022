#include <string>
#include <algorithm>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stack>
#include <utility>
#include <math.h>
#include <unordered_set>
#include <limits>

using namespace std;

struct Point {
    bool isBlue;
    bool inRange;
    double *coordinates;
    int dimensions;

    Point() {
        coordinates = nullptr;
    }

    Point(const Point &other) {
        isBlue = other.isBlue;
        inRange = other.inRange;
        coordinates = new double[other.dimensions];
        for (int i = 0; i < other.dimensions; i++) {
            coordinates[i] = other.coordinates[i];
        }
        dimensions = other.dimensions;
    }

    Point& operator=(const Point &other) {
        isBlue = other.isBlue;
        inRange = other.inRange;
        dimensions = other.dimensions;
        if (coordinates != nullptr) {
            delete [] coordinates;
        }
        coordinates = new double[other.dimensions];
        for (int i = 0; i < other.dimensions; i++) {
            coordinates[i] = other.coordinates[i];
        }
        return *this;
    }

    ~Point() {
        if (coordinates != nullptr) {
            delete [] coordinates;
        }
        coordinates = nullptr;
    }
};

struct avlNode {
    avlNode *left;
    avlNode *right;
    int height;
    Point point;
};

struct rangeNode {
    rangeNode *left;
    rangeNode *right;
    rangeNode *parent;
    rangeNode *nextDimension;
    int dimenison;
    int depth;
    int height;
    int size;
    Point point;
};

struct range {
    double *start;
    double *end;
    bool inclusiveXMin;
    bool inclusiveXMax;
    bool inclusiveYMin;
    bool inclusiveYMax;

    range() {
        start = nullptr;
        end = nullptr;
    }

    range(double *s, double *e, bool ixmin, bool ixmax, bool iymin, bool iymax) {
        //start = s;
        //end = e;
        start = new double[2];
        end = new double[2];
        start[0] = s[0];
        end[0] = e[0];
        start[1] = s[1];
        end[1] = e[1];
        inclusiveXMin = ixmin;
        inclusiveXMax = ixmax;
        inclusiveYMin = iymin;
        inclusiveYMax = iymax; 
    }

    range& operator=(const range& other) {
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
        start = new double[2];
        end = new double[2];
        start[0] = other.start[0];
        end[0] = other.end[0];
        start[1] = other.start[1];
        end[1] = other.end[1];
        return *this;
    }

    range(const range& other) {
        inclusiveXMin = other.inclusiveXMin;
        inclusiveXMax = other.inclusiveXMax;
        inclusiveYMin = other.inclusiveYMin;
        inclusiveYMax = other.inclusiveYMax;
        start = new double[2];
        end = new double[2];
        start[0] = other.start[0];
        end[0] = other.end[0];
        start[1] = other.start[1];
        end[1] = other.end[1];
    }

    ~range() {
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

class AVLTree {
    avlNode *root;
    int dimension;
    
    private: 
    void outputTree() {
        vector<Point> *traversal = preOrderTraversal();
        for (Point point : *traversal) {
            cout << "{" << point.coordinates[0] << ", " << point.coordinates[1] << "}, ";
        }
        cout << endl;
    }

    private:
        void computeHeight(avlNode *node) {
            node->height = max(getHeight(node->left), getHeight(node->right)) + 1;
        }

        int getHeight(avlNode *node) {
            return node == nullptr ? -1 : node->height;
        }

        avlNode *rightRotate(avlNode *node) {
            avlNode *A = node->left;
            avlNode *B = node->left->right;

            node->left = B;
            A->right = node;

            computeHeight(node);
            computeHeight(A);

            return A;
        }

        avlNode *leftRotate(avlNode *node) {
            avlNode *A = node->right;
            avlNode *B = node->right->left;

            node->right = B;
            A->left = node;

            computeHeight(node);
            computeHeight(A);

            return A;
        }

        avlNode *insertRecursive(Point point, avlNode *current) {
            if (current == nullptr) {
                avlNode *node = new avlNode;
                node->left = nullptr;
                node->right = nullptr;
                node->height = 0;
                node->point = point;
                return node;
            }
            if (point.coordinates[this->dimension] <= current->point.coordinates[this->dimension]) {
                current->left = insertRecursive(point, current->left);
            } else {
                current->right = insertRecursive(point, current->right);
            }

            computeHeight(current);

            if (getHeight(current->left) > getHeight(current->right) + 1 && point.coordinates[this->dimension] <= current->left->point.coordinates[this->dimension]) {
                return rightRotate(current);
            } else if (getHeight(current->left) > getHeight(current->right) + 1) {
                current->left = leftRotate(current->left);
                return rightRotate(current); 
            } else if (getHeight(current->right) > getHeight(current->left) + 1 && point.coordinates[this->dimension] > current->right->point.coordinates[this->dimension]) {
                return leftRotate(current);
            } else if (getHeight(current->right) > getHeight(current->left) + 1) {
                current->right = rightRotate(current->right);
                return leftRotate(current);
            }
            return current;
        }

        void preOrderTraversalRecursive(vector<Point> *result, avlNode *current) {
            if (current == nullptr) {
                return;
            }
            (*result).push_back(current->point);
            preOrderTraversalRecursive(result, current->left);
            preOrderTraversalRecursive(result, current->right);
        }
        
        void deleteNode(avlNode *node) {
            if (node == nullptr) {
                return;
            }
            deleteNode(node->left);
            deleteNode(node->right);
            delete node;
        }

    public:
        ~AVLTree() {
            deleteNode(this->root);
        }

        AVLTree(int dimension) {
            this->root = nullptr;
            this->dimension = dimension;
        };

        void insert(Point point) {
            this->root = insertRecursive(point, this->root);
        };

        vector<Point> *preOrderTraversal() {
            vector<Point> *result = new vector<Point>();
            result->reserve(pow(getHeight(this->root), 2));
            preOrderTraversalRecursive(result, this->root);
            return result;   
        }
};

class RangeTree {
    rangeNode *root;
    rangeNode *invertedRoot;
    int dimensions;

    private:
        void computeHeight(rangeNode *node) {
            node->height = max(getHeight(node->left), getHeight(node->right)) + 1;
        }

        int getHeight(rangeNode *node) {
            return node == nullptr ? -1 : node->height;
        }

        void computeSize(rangeNode *node) {
            node->size = getSize(node->left) + getSize(node->right) + 1;
        }

        int getSize(rangeNode *node) {
            return node == nullptr ? 0 : node->size;
        }
        vector<Point> *readPoints(string filename, double *minimum, double *maximum, int dimensions) {
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
                point.coordinates = new double[dimensions];
                for (int i = 0; i < dimensions; i++) {
                    ss >> point.coordinates[i];
                }
                int isBlue;
                //1 is blue?
                ss >> isBlue;
                point.isBlue = isBlue == 1;
                point.inRange = true;
                point.dimensions = dimensions;
                for (int i = 0; i < dimensions; i++) {
                    if (point.coordinates[i] < minimum[i] || point.coordinates[i] > maximum[i]) {
                        point.inRange = false;
                    }    
                }
                (*points).push_back(point);
                getline(file, line, '\n');
            }
            return points;  
        }

        rangeNode *insertPoint(rangeNode* current, Point point, int dimension, int depth) {
            if (current == nullptr) {
                rangeNode *node = new rangeNode;
                node->left = nullptr;
                node->right = nullptr;
                node->nextDimension = nullptr; 
                node->parent = nullptr;
                node->dimenison = dimension;
                node->point = point;
                node->depth = depth;
                node->height = 0;
                node->size = 1;
                return node;
            }
            if (point.coordinates[dimension] <= current->point.coordinates[dimension]) {
                current->left = insertPoint(current->left, point, dimension, depth+1);
                current->left->parent = current;
            } else {
                current->right = insertPoint(current->right, point, dimension, depth+1);
                current->right->parent = current;
            }
            computeHeight(current);
            computeSize(current);
            return current;
        }

        rangeNode *findPoint(rangeNode *current, double target, int dimension) {
            if (target < current->point.coordinates[dimension] && current->left == nullptr) {
                return current;
            } else if (target < current->point.coordinates[dimension]) {
                return findPoint(current->left, target, dimension);
            } else if (target > current->point.coordinates[dimension] && current->right == nullptr) {
                return current;
            } else if (target > current->point.coordinates[dimension]) {
                return findPoint(current->right, target, dimension);
            } else {
                return current;
            }
        }

        void buildNextDimensionTree(rangeNode *current, vector<Point> *descendants, int currentDimension) {
            if (current == nullptr) {
                return;
            }
            vector<Point> *leftDescendants = new vector<Point>();
            leftDescendants->reserve(getSize(current->left));
            vector<Point> *rightDescendants = new vector<Point>();
            rightDescendants->reserve(getSize(current->right));
            buildNextDimensionTree(current->left, leftDescendants, currentDimension);
            buildNextDimensionTree(current->right,rightDescendants, currentDimension);
 
            vector<Point> allDescendants;
            allDescendants.reserve((*leftDescendants).size() + (*rightDescendants).size() + 1);
            allDescendants.insert(allDescendants.end(), (*leftDescendants).begin(), (*leftDescendants).end());
            allDescendants.push_back(current->point);
            allDescendants.insert(allDescendants.end(), (*rightDescendants).begin(), (*rightDescendants).end());

            AVLTree tree(1);

            for (Point point : allDescendants) {
                tree.insert(point);
            }
            vector<Point> *nextDimensionTree = tree.preOrderTraversal();
            for (Point point : *nextDimensionTree) {
                current->nextDimension = insertPoint(current->nextDimension, point, 1, 0);
            }
            for (Point point : allDescendants) {
                (*descendants).push_back(point);
            }
            delete nextDimensionTree;
            delete leftDescendants;
            delete rightDescendants;
        }

        vector<Point> *findPointsInYSubtree(rangeNode *current, double minY, double maxY, bool inclusiveYMin, bool inclusiveYMax) {
            vector<Point> *points = new vector<Point>();
            // points->reserve(getSize(current) + 1);

            if (current == nullptr) {
                return points;
            }

            rangeNode *startY = findPoint(current, minY, 1);
            while (startY != nullptr && (inclusiveYMin ? startY->point.coordinates[1] < minY : startY->point.coordinates[1] <= minY)) {
                startY = next(startY);
            }
            rangeNode *endY = findPoint(current, maxY, 1);
            while (endY != nullptr && (inclusiveYMax ? endY->point.coordinates[1] > maxY : endY->point.coordinates[1] >= maxY)) {
                endY = prev(endY);
            }

            if (startY == nullptr || endY == nullptr || startY->point.coordinates[1] > endY->point.coordinates[1]) {
                return points;
            }

            while (startY != endY) {
               points->push_back(startY->point);
                startY = next(startY);
            }

            points->push_back(endY->point);

            return points;
        }

        rangeNode *next(rangeNode *current) {
            if (current->right != nullptr) {
                rangeNode *iter = current->right;
                while (iter->left != nullptr) {
                    iter = iter->left;
                }
                return iter;
            } else if (current->parent != nullptr) {
                rangeNode *iter = current->parent;
                rangeNode *prev = current;
                while (iter != nullptr && iter->right == prev) {
                    prev = iter;
                    iter = iter->parent;
                }
                return iter;
            }
            return nullptr;
        }


        rangeNode *prev(rangeNode *current) {
            if (current->left != nullptr) {
                rangeNode *iter = current->left;
                while (iter->right != nullptr) {
                    iter = iter->right;
                }
                return iter;            } else if (current->parent != nullptr) {
                rangeNode *iter = current->parent;
                rangeNode *prev = current;
                while (iter != nullptr && iter->left == prev) {
                    prev = iter;
                    iter = iter->parent;
                }
                return iter;
            }
            return nullptr;
        }

        rangeNode *findFork(rangeNode *first, rangeNode *second) {
            while (first->depth > second->depth) {
                first = first->parent;
            }
            while (second->depth > first->depth) {
                second = second->parent;
            }
            while (first != nullptr && second != nullptr) {
                if (first == second) {
                    return first;
                }
                first = first->parent;
                second = second->parent;
            }
            throw new runtime_error("Unable to find fork");
        }

        //the included trees on the left side of the fork
        vector<rangeNode*> *leftSideTrees(rangeNode *start, rangeNode *fork) {
            vector<rangeNode*> *leftTrees = new vector<rangeNode*>();
            leftTrees->reserve(getHeight(fork) + 1);
            while (start != fork) {
                if (start->right != nullptr) {
                    leftTrees->push_back(start->right->nextDimension);
                }
                //skip the ones that are less than start.
                while (start->parent->right == start) {
                    start = start->parent;
                }
                start = start->parent;
            }
            return leftTrees;
        }
        
        //the included individual nodes on the left side of the fork
        vector<rangeNode*> *leftSideNodes(rangeNode *start, rangeNode *fork) {
            vector<rangeNode*> *leftNodes = new vector<rangeNode*>();
            leftNodes->reserve(getHeight(fork) + 1);
            while (start != fork) {
                leftNodes->push_back(start);
                //skip the ones that are less than start.
                while (start->parent->right == start) {
                    start = start->parent;
                }
                start = start->parent;
            }
            return leftNodes;
        }

        //the included trees on the right side of the fork        
        vector<rangeNode*> *rightSideTrees(rangeNode *end, rangeNode *fork) {
            vector<rangeNode*> *rightTrees = new vector<rangeNode*>();
            rightTrees->reserve(getHeight(fork) + 1);
            while (end != fork) {
                if (end->left != nullptr) {
                    rightTrees->push_back(end->left->nextDimension);
                }
                //skip the ones that are greater than end.
                while (end->parent->left == end) {
                    end = end->parent;
                }
                end = end->parent;
            }
            return rightTrees;
        }
        
        //the included individual nodes on the right side of the fork
        vector<rangeNode*> *rightSideNodes(rangeNode *end, rangeNode *fork) {
            vector<rangeNode*> *rightNodes = new vector<rangeNode*>();
            rightNodes->reserve(getHeight(fork) + 1);
            while (end != fork) {
                rightNodes->push_back(end);
                //skip the ones that are greater than end.
                while (end->parent->left == end) {
                    end = end->parent;
                }
                end = end->parent;
            }
            return rightNodes;
        }

        vector<Point> *inversePoints(vector<Point> *original) {
            vector<Point> *result = new vector<Point>();
            result->reserve(original->size() + 1);
            for (Point p: *original) {
                Point point;
                point.dimensions = p.dimensions;
                point.coordinates = new double[p.dimensions];
                for (int i = 2; i < p.dimensions; i++) {
                    point.coordinates[i] = p.coordinates[i];
                }
                point.coordinates[0] = p.coordinates[1];
                point.coordinates[1] = p.coordinates[0];
                point.isBlue = p.isBlue;
                point.inRange = p.inRange;
                result->push_back(point);
            }
            return result;
        }

        Point * findNearestInTree(rangeNode* root, double *start, double *end, int dimension) {
            bool min = start[dimension] < end[dimension];
            double startC0 = start[1-dimension] < end[1-dimension] ? start[1-dimension] : end[1-dimension];
            double endC0 = start[1-dimension] < end[1-dimension] ? end[1-dimension] : start[1-dimension];
            double startC1 = min ? start[dimension] : end[dimension];
            double endC1 = min ? end[dimension] : start[dimension];

            rangeNode *startD1 = findPoint(root, startC0, 0);

            while (startD1 != nullptr && startD1->point.coordinates[0] < startC0 ) {
                startD1 = next(startD1);
            }

            rangeNode *endD1 = findPoint(root, endC0, 0);
            while (endD1 != nullptr && endD1->point.coordinates[0] > endC0) {
                endD1 = prev(endD1);
            }

            if (startD1 == nullptr || endD1 == nullptr || startD1->point.coordinates[0] > endD1->point.coordinates[0]) {
                return nullptr;
            }

            rangeNode *fork = findFork(startD1, endD1);

            vector<Point>candidates;
            candidates.reserve(getSize(fork));

            vector<rangeNode*> *lsn = leftSideNodes(startD1, fork);
            for (rangeNode *node: *lsn) {
                if (min ? (node->point.coordinates[1] > startC1 && node->point.coordinates[1] <= endC1) : 
                    (node->point.coordinates[1] >= startC1 && node->point.coordinates[1] < endC1)) {
                    candidates.push_back(node->point);
                }
            }

            vector<rangeNode*> *rsn = rightSideNodes(endD1, fork);
            for (rangeNode *node: *rsn) {
                if (min ? (node->point.coordinates[1] > startC1 && node->point.coordinates[1] <= endC1) : 
                    (node->point.coordinates[1] >= startC1 && node->point.coordinates[1] < endC1)) {
                    candidates.push_back(node->point);
                }
            }

            if (min ? (fork->point.coordinates[1] > startC1 && fork->point.coordinates[1] <= endC1) : 
                (fork->point.coordinates[1] >= startC1 && fork->point.coordinates[1] < endC1)) {
                candidates.push_back(fork->point);            

            }

            vector<rangeNode*> *lst = leftSideTrees(startD1, fork);
            for (rangeNode *node : *lst) {
                rangeNode *current = findPoint(node, min ? startC1 : endC1, 1);
                if (min) {
                    while (current != nullptr && current->point.coordinates[1] <= startC1) {
                        current = next(current);
                    }
                } else {
                    while (current != nullptr && current->point.coordinates[1] >= endC1) {
                        current = prev(current);
                    }
                }
                if (current == nullptr || (min && current->point.coordinates[1] > endC1) || (!min && current->point.coordinates[1] < startC1)) {
                    continue;
                }
                candidates.push_back(current->point);
            }

            vector<rangeNode*> *rst = rightSideTrees(endD1, fork);
            for (rangeNode *node : *rst) {
                rangeNode *current = findPoint(node, min ? startC1 : endC1, 1);
                if (min) {
                    while (current != nullptr && current->point.coordinates[1] <= startC1) {
                        current = next(current);
                    }
                } else {
                    while (current != nullptr && current->point.coordinates[1] >= endC1) {
                        current = prev(current);
                    }
                }
                if (current == nullptr || (min && current->point.coordinates[1] > endC1) || (!min && current->point.coordinates[1] < startC1)) {
                    continue;
                }
                candidates.push_back(current->point);
            }
            delete rsn;
            delete lsn;
            delete rst;
            delete lst;

            if (candidates.empty()) {
                return nullptr;
            }
            Point max = candidates[0];
            for (Point p : candidates) {
                if (min ? p.coordinates[1] < max.coordinates[1] : p.coordinates[1] > max.coordinates[1]) {
                    max = p;
                }
            }
            Point *result = new Point();
            result->coordinates = new double[this->dimensions];
            result->coordinates[0] = dimension == 1 ? max.coordinates[0] : max.coordinates[1];
            result->coordinates[1] = dimension == 1 ? max.coordinates[1] : max.coordinates[0];
            for (int i = 2; i < this->dimensions; i++) {
                result->coordinates[i] = max.coordinates[i];
            }
            result->isBlue = max.isBlue;
            result->inRange = max.inRange;

            return result;
        }

        void deleteNode(rangeNode* current) {
            if (current == nullptr) {
                return;
            }
            deleteNode(current->left);
            current->left = nullptr;
            deleteNode(current->right);
            current->right = nullptr;
            deleteNode(current->nextDimension);
            current->nextDimension = nullptr;
            delete current;
        }

    public:
        RangeTree(string filename, double *minimum, double *maximum) : RangeTree(filename, minimum, maximum, 2){}

        RangeTree(string filename, double *minimum, double *maximum, int dimensions) {
            this->root = nullptr;
            this->invertedRoot = nullptr;
            this->dimensions = dimensions;

            vector<Point> *points = readPoints(filename, minimum, maximum, dimensions);
            AVLTree tree(0);
            for (Point point : *points) {
                tree.insert(point);
            }
            vector<Point> *zeroDimensionTree = tree.preOrderTraversal();
            for (Point point : *zeroDimensionTree) {
                this->root = insertPoint(this->root, point, 0, 0);
            }
            vector<Point> *descendants = new vector<Point>();
            descendants->reserve(getSize(this->root));
            buildNextDimensionTree(this->root, descendants, 0);

            vector<Point> *invertedPoints = inversePoints(points);
            AVLTree invertedTree(0);
            for (Point point : *invertedPoints) {
                invertedTree.insert(point);
            }
            vector<Point> *zeroDimensionInvertedTree = invertedTree.preOrderTraversal();
            for (Point point : *zeroDimensionInvertedTree) {
                this->invertedRoot = insertPoint(this->invertedRoot, point, 0, 0);
            }
            vector<Point> *invertedDescendants = new vector<Point>();
            invertedDescendants->reserve(getSize(this->invertedRoot));
            buildNextDimensionTree(this->invertedRoot, invertedDescendants, 0);

            delete zeroDimensionTree;
            delete descendants;
            delete points;
            delete zeroDimensionInvertedTree;
            delete invertedDescendants;
            delete invertedPoints;
        }

        ~RangeTree() {
            deleteNode(this->root);
            this->root = nullptr;
            deleteNode(this->invertedRoot);
            this->invertedRoot = nullptr;
        }

        vector<Point> *findPointsInRange(range r) {
            double startXDimension = r.start[0] < r.end[0] ? r.start[0] : r.end[0];
            double endXDimension = r.start[0] > r.end[0] ? r.start[0] : r.end[0];
            double startYDimension = r.start[1] < r.end[1] ? r.start[1] : r.end[1];
            double endYDimension = r.start[1] > r.end[1] ? r.start[1] : r.end[1];

            rangeNode *startX = findPoint(root, startXDimension, 0);

            while (startX != nullptr && (r.inclusiveXMin ? startX->point.coordinates[0] < startXDimension : startX->point.coordinates[0] <= startXDimension)) {
                startX = next(startX);
            }
            rangeNode *endX = findPoint(root, endXDimension, 0);
            while (endX != nullptr && (r.inclusiveXMax ? endX->point.coordinates[0] > endXDimension : endX->point.coordinates[0] >= endXDimension)) {
                endX = prev(endX);
            }

            if (startX == nullptr || endX == nullptr || startX->point.coordinates[0] > endX->point.coordinates[0]) {
                return new vector<Point>();
            }

            rangeNode *fork = findFork(startX, endX);

            vector<Point> *pointsInRange = new vector<Point>;
            // pointsInRange->reserve(getSize(fork));

            vector<rangeNode*> *lsn = leftSideNodes(startX, fork);
            for (rangeNode *node : *lsn) {
                if ((r.inclusiveYMin ? node->point.coordinates[1] >= startYDimension : node->point.coordinates[1] > startYDimension) &&
                        (r.inclusiveYMax ? node->point.coordinates[1] <= endYDimension : node->point.coordinates[1] < endYDimension)) {
                    pointsInRange->push_back(node->point);
                }
            }
            delete lsn;

            vector<rangeNode*> *rsn = rightSideNodes(endX, fork);
            for (rangeNode *node : *rsn) {
                if ((r.inclusiveYMin ? node->point.coordinates[1] >= startYDimension : node->point.coordinates[1] > startYDimension) &&
                        (r.inclusiveYMax ? node->point.coordinates[1] <= endYDimension : node->point.coordinates[1] < endYDimension)) {
                    pointsInRange->push_back(node->point);
                }
            }
            delete rsn;

            if ((r.inclusiveYMin ? fork->point.coordinates[1] >= startYDimension : fork->point.coordinates[1] > startYDimension) &&
                    (r.inclusiveYMax ? fork->point.coordinates[1] <= endYDimension : fork->point.coordinates[1] < endYDimension)) {
                pointsInRange->push_back(fork->point);
            }

            vector<rangeNode*> *lst = leftSideTrees(startX, fork);
            for (rangeNode *node : *lst) {
                vector<Point> *ySubTree = findPointsInYSubtree(node, startYDimension, endYDimension, r.inclusiveYMin, r.inclusiveYMax);
                for (Point p : *ySubTree) {
                    pointsInRange->push_back(p);
                }
                delete ySubTree;
            }
            delete lst;

            vector<rangeNode*> *rst = rightSideTrees(endX, fork);
            for (rangeNode *node : *rst) {
                vector<Point> *ySubTree = findPointsInYSubtree(node, startYDimension, endYDimension, r.inclusiveYMin, r.inclusiveYMax);
                for (Point p : *ySubTree) {
                    pointsInRange->push_back(p);
                }
                delete ySubTree;
            }
            delete rst;
            return pointsInRange;
        }

        vector<Point> *findSkyline(double *start, double *end) {
            bool xMin = start[0] < end[0];
            bool yMin = start[1] < end[1];
            stack<range> ranges;
            range startRange;
            startRange.start = start;
            startRange.end = end;
            startRange.inclusiveXMin = true;
            startRange.inclusiveYMin = true;
            startRange.inclusiveXMax = true;
            startRange.inclusiveYMax = true;
            vector<Point> *skyline = new vector<Point>();
            ranges.push(startRange);
            while (!ranges.empty()) {
                range r = ranges.top();
                ranges.pop();
                vector<Point> *points = findPointsInRange(r);
                if ((*points).size() == 0) {
                    delete points;
                    continue;
                }
                Point minY;
                double pos1 = yMin ? min(r.start[0], r.end[0]) : max(r.start[0], r.end[0]);
                double pos2 = yMin ? min(r.start[1], r.end[1]) : max(r.start[1], r.end[1]);
                if ((*points)[0].coordinates[0] == start[0] && (*points)[0].coordinates[1] == start[1]) {
                    if ((*points).size() == 1) {
                        delete points;
                        continue;
                    }
                    minY = (*points)[1]; 
                } else {
                    minY = (*points)[0];
                }
                for (int i = 1; i < (*points).size(); i++) {
                    if ((*points)[i].coordinates[0] == start[0] && (*points)[i].coordinates[1] == start[1]) {
                        continue;
                    }
                    if (yMin) {
                        if (((*points)[i]).coordinates[1] < minY.coordinates[1]) {
                            minY = (*points)[i];
                        }
                    } else {
                        if (((*points)[i]).coordinates[1] > minY.coordinates[1]) {
                            minY = (*points)[i];
                        }
                    }
                }
                (*skyline).push_back(minY);
                double *newEnd = new double[2];
                newEnd[0] = xMin ? min(r.start[0], r.end[0]) : max(r.start[0], r.end[0]);
                newEnd[1] = yMin ? max(r.start[1], r.end[1]) : min(r.start[1], r.end[1]);
                range newRange;
                newRange.start = new double[2];
                newRange.start[0] = minY.coordinates[0];
                newRange.start[1] = minY.coordinates[1];
                newRange.end = newEnd;
                //these are kind of hacked, need to be recalculated for D dimensions
                
                newRange.inclusiveXMin = xMin ? true : false;
                newRange.inclusiveXMax = xMin ? false : true;
                newRange.inclusiveYMin = yMin ? false : true;
                newRange.inclusiveYMax = yMin ? true : false;
                ranges.push(newRange);
                delete points;
        }
        return skyline;
    }

    Point *findNearestEnhanced(double *start, double *end, int dimension) {
            range r = {start, end, false, false, false, false};
            vector<Point> *pointsInRange = findPointsInRange(r);
            
            bool *increasing = new bool[this->dimensions];
            for (int i = 2; i < this->dimensions; i++) {
                increasing[i] = start[i] < end[i];
            }

            double smallestDifference = numeric_limits<double>::max();
            Point bestPoint;
            bool pointsFound = false;
            
            for (Point p : *pointsInRange) {
                bool inRange = true;
                for (int i = 2; i < this->dimensions; i++) {
                    if (increasing[i] ? (p.coordinates[i] <= start[i] || p.coordinates[i] >= end[i]) : (p.coordinates[i] >= start[i] || p.coordinates[i] <= end[i])) {
                        inRange = false;
                    }
                }
                if (inRange && (abs(p.coordinates[dimension] - start[dimension]) < smallestDifference)) {
                    smallestDifference = abs(p.coordinates[dimension] - start[dimension]);
                    bestPoint = p;
                    pointsFound = true;
                }
            }

            if (!pointsFound) {
                return nullptr;
            }

            Point *result = new Point();
            result->isBlue = bestPoint.isBlue;
            result->inRange = bestPoint.inRange;
            result->dimensions = bestPoint.dimensions;
            result->coordinates = new double[bestPoint.dimensions];
            for (int i = 0; i < bestPoint.dimensions; i++) {
                result->coordinates[i] = bestPoint.coordinates[i];
            }
            return result;
    }

        vector<Point> *findSkylineEnhanced(double *start, double *end) {
            bool *dMin = new bool[this->dimensions];
            for (int i = 0; i < dimensions; i++) {
                dMin[i] = start[i] < end[i];
            }

            range r = {start, end, true, true, true, true};
            vector<Point> *pointsInRange = findPointsInRange(r);

            for (int i = 2; i < this->dimensions; i++) {
                vector<Point> *pointsInRangeNext = new vector<Point>();
                for (Point p : *pointsInRange) {
                    if (dMin[i] ? (p.coordinates[i] > start[i] && p.coordinates[i] < end[i]) : (p.coordinates[i] < start[i] && p.coordinates[i] > end[i])) {
                        pointsInRangeNext->push_back(p);
                    }
                }
                delete pointsInRange;
                pointsInRange = pointsInRangeNext;
            }

            vector<Point> *skylinePoints = new vector<Point>();

            for (Point p : *pointsInRange){
                bool addToSkyline = true;
                for (Point o : *pointsInRange) {
                    bool pointsEqual = true;
                    for (int i = 0; i < this->dimensions; i++) {
                        if (p.coordinates[i] != o.coordinates[i]) {
                            pointsEqual = false;
                        }
                    }
                    if (pointsEqual) {
                        continue;
                    }
                    bool dominated = true;
                    for (int i = 0; i < this->dimensions; i++) {
                        if (dMin[i] ? (p.coordinates[i] < o.coordinates[i]) : (p.coordinates[i] > o.coordinates[i])) {
                            dominated = false;
                        }
                    }
                    if (dominated) {
                        addToSkyline = false;
                    }
                }
                if (addToSkyline) {
                    skylinePoints->push_back(p);
                }
            }

            return skylinePoints;
        }


        Point *findNearest(double *start, double *end, int dimension) {
            if (dimension == 0) {
                return findNearestInTree(this->invertedRoot, start, end, dimension);
            } else {
                return findNearestInTree(this->root, start, end, dimension);
            }
        }
};
