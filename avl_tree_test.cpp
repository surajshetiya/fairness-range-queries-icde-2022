#include "range_tree.h"
#include <assert.h>
#include <iostream>

using namespace std;

bool compare2dArray(double a[], double b[]) {
    return a[0] == b[0] && a[1] == b[1];
}

int leftRotate() {
    Point point1;
    point1.coordinates = new double[1];
    point1.coordinates[0] = 1.0;
    Point point2;
    point2.coordinates = new double[1];
    point2.coordinates[0] = 2.0;
    Point point3;
    point3.coordinates = new double[1];
    point3.coordinates[0] = 3.0;
    
    AVLTree tree(0);
    tree.insert(point1);
    tree.insert(point2);
    tree.insert(point3);

    vector<Point> *traversal = tree.preOrderTraversal();

    if (!compare2dArray((*traversal)[0].coordinates, point2.coordinates)) {
        return false;
    } else if (!compare2dArray((*traversal)[1].coordinates, point1.coordinates)) {
        return false;
    } else if (!compare2dArray((*traversal)[2].coordinates, point3.coordinates)) {
        return false;
    } else {
        return true;
    }
}

int twoLeftRotates() {
    Point point1;
    point1.coordinates = new double[1];
    point1.coordinates[0] = 1.0;
    Point point2;
    point2.coordinates = new double[1];
    point2.coordinates[0] = 2.0;
    Point point3;
    point3.coordinates = new double[1];
    point3.coordinates[0] = 3.0;
    Point point4;
    point4.coordinates = new double[1];
    point4.coordinates[0] = 4.0;
    Point point5;
    point5.coordinates = new double[1];
    point5.coordinates[0] = 5.0;
    
    AVLTree tree(0);
    tree.insert(point1);
    tree.insert(point2);
    tree.insert(point3);
    tree.insert(point4);
    tree.insert(point5);


    vector<Point> *traversal = tree.preOrderTraversal();

    if (!compare2dArray((*traversal)[0].coordinates, point2.coordinates)) {
        return false;
    } else if (!compare2dArray((*traversal)[1].coordinates, point1.coordinates)) {
        return false;
    } else if (!compare2dArray((*traversal)[2].coordinates, point4.coordinates)) {
        return false;
    } else if (!compare2dArray((*traversal)[3].coordinates, point3.coordinates)) {
        return false;
    } else if (!compare2dArray((*traversal)[4].coordinates, point5.coordinates)) {
        return false;
    } else {
        return true;
    }
}


int rightRotate() {
    Point point1;
    point1.coordinates = new double[1];
    point1.coordinates[0] = 3.0;
    Point point2;
    point2.coordinates = new double[1];
    point2.coordinates[0] = 2.0;
    Point point3;
    point3.coordinates = new double[1];
    point3.coordinates[0] = 1.0;
    
    AVLTree tree(0);
    tree.insert(point1);
    tree.insert(point2);
    tree.insert(point3);

    vector<Point> *traversal = tree.preOrderTraversal();

    if (!compare2dArray((*traversal)[0].coordinates, point2.coordinates)) {
        return false;
    } else if (!compare2dArray((*traversal)[1].coordinates, point3.coordinates)) {
        return false;
    } else if (!compare2dArray((*traversal)[2].coordinates, point1.coordinates)) {
        return false;
    } else {
        return true;
    }
}


int leftRightRotate() {
    Point point1;
    point1.coordinates = new double[1];
    point1.coordinates[0] = 1.0;
    Point point2;
    point2.coordinates = new double[1];
    point2.coordinates[0] = 3.0;
    Point point3;
    point3.coordinates = new double[1];
    point3.coordinates[0] = 2.0;
    
    AVLTree tree(0);
    tree.insert(point1);
    tree.insert(point2);
    tree.insert(point3);

    vector<Point> *traversal = tree.preOrderTraversal();

    if (!compare2dArray((*traversal)[0].coordinates, point3.coordinates)) {
        return false;
    } else if (!compare2dArray((*traversal)[1].coordinates, point1.coordinates)) {
        return false;
    } else if (!compare2dArray((*traversal)[2].coordinates, point2.coordinates)) {
        return false;
    } else {
        return true;
    }
}

int rightLeftRotate() {
    Point point1;
    point1.coordinates = new double[1];
    point1.coordinates[0] = 3.0;
    Point point2;
    point2.coordinates = new double[1];
    point2.coordinates[0] = 1.0;
    Point point3;
    point3.coordinates = new double[1];
    point3.coordinates[0] = 2.0;
    
    AVLTree tree(0);
    tree.insert(point1);
    tree.insert(point2);
    tree.insert(point3);

    vector<Point> *traversal = tree.preOrderTraversal();

    if (!compare2dArray((*traversal)[0].coordinates, point3.coordinates)) {
        return false;
    } else if (!compare2dArray((*traversal)[1].coordinates, point2.coordinates)) {
        return false;
    } else if (!compare2dArray((*traversal)[2].coordinates, point1.coordinates)) {
        return false;
    } else {
        return true;
    }
}

int main() {
    int testsFailed = 0;

    if (!leftRotate()) {
        cout << "Left rotate failed!" << endl;
        testsFailed++;
    }
    if (!twoLeftRotates()) {
        cout << "Two left rotates failed!" << endl;
        testsFailed++;
    }
    if (!rightRotate()) {
        cout << "Right rotate failed!" << endl;
        testsFailed++;
    }
    if (!leftRightRotate()) {
        cout << "Left right rotate failed!" << endl;
        testsFailed++;
    }
    if (!rightLeftRotate()) {
        cout << "Right left rotate failed!" << endl;
        testsFailed++;
    }

    if (testsFailed == 0) {
        cout << "All tests passed!" << endl;
    }
}