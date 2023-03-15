#include <iostream>
#include <cmath>
#include "vector.h"

using namespace std;

int main() {
    double v[] = {1, 1, 1};
    if (norm(v, 3) - sqrt(3.0) < pow(10, -12)) {
        cout << "Test norm(): passed" << endl;
    } else {
        cout << "Test norm(): failed" << endl;
    }
    double v1[] = {1, 1, 2};
    double v2[] = {1, 2, 3};
    if(dot(v1,v2,3,3)-(v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2])<pow(10,-12)){
        cout << "Test dot(): passed" << endl;
    }
    else{
        cout << "Test norm(): failed" << endl;
    }
    vector<double> c= cross(v1,v2,3,3);
    if(c[0]-(v1[0]*v2[0]-v1[2]*v2[1])<pow(10,12) && c[1]-(v1[2]*v2[0]-v1[0]*v2[2])<pow(10,12) & c[2]-(v1[0]*v2[1]-v1[1]*v2[0])<pow(10,12) ){
        cout << "Test cross(): passed" << endl;
    }
    else{
        cout << "Test norm(): failed" << endl;
    }
}