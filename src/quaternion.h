#pragma once

#include <vector>

using namespace std;

class Quaternion
{
public:

    Quaternion(float aw, float ai, float aj, float ak);
    Quaternion(float phi, float psi, float theta);
    Quaternion(float angle, vector<float> axis);
    
    Quaternion(Quaternion& q);
    Quaternion(){};
    Quaternion(const Quaternion& q);


    ~Quaternion(){};

    static void generateConeRotations(vector<Quaternion>& rotations,float coneAngle, float coneSampling, float inplaneAngle, float inplaneSampling);
    static void generateConeRotations(vector<Quaternion>& rotations, vector<float>& phiAngles, float coneAngle, float coneSampling, float inplaneAngle, float inplaneSampling);

    static void slerp(vector<Quaternion>& output, Quaternion& q1, Quaternion& q2, float stepSize, float interpolationLimit = 1.0f);

    static Quaternion mult(Quaternion& q1,Quaternion& q2);
    vector<float> rotate(Quaternion& q,vector<float>& point);
    vector<float> pointRotate(vector<float>& point);
    Quaternion inverseQuaternion();
    void getRotationMatrix(vector<float>& matrix);

    void normalize();
    static float jointAngle(Quaternion& q1,Quaternion& q2);
    static Quaternion add(Quaternion& q1, Quaternion& q2);
    static Quaternion mult(Quaternion& q1, float value);
    vector<float> getEulerAngles();
    void getAxisAngle(float& angle, vector<float>& axis);

    static void getAngularDistance(Quaternion q1, Quaternion q2, float& angleInDegrees, float& distance);

    bool isEqual(Quaternion& q);

    float w;
    float i;
    float j;
    float k;

private:

    void operator*(float value);
    Quaternion operator+(Quaternion& q);

};
