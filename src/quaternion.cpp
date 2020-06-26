#include "quaternion.h"
#include "defines.h"
#include "vectorRoutines.h"

using namespace vr;

Quaternion::Quaternion(float angle, vector<float> axis)
{
    vr::normalize(axis);
    float halfAngle = angle / 2.0f;

    w = cos(halfAngle);
    i = axis[0] * sin(halfAngle);
    j = axis[1] * sin(halfAngle);
    k = axis[2] * sin(halfAngle);

    normalize();
}

Quaternion::Quaternion(const Quaternion& q)
{
    w=q.w;
    i=q.i;
    j=q.j;
    k=q.k;
}

Quaternion::Quaternion(float aw, float ai, float aj, float ak)
{
    w = aw;
    i = ai;
    j = aj;
    k = ak;
}

Quaternion::Quaternion(Quaternion& q)
{
    w=q.w;
    i=q.i;
    j=q.j;
    k=q.k;
}

Quaternion::Quaternion(float phi, float psi, float theta)
{
    double rphi = phi*M_PI/180;
    double rpsi = psi*M_PI/180;
    double rtheta = theta*M_PI/180;

    w=(float)(cos((rpsi+rphi)/2)*cos(rtheta/2));
    i=(float)(cos((rpsi-rphi)/2)*sin(rtheta/2));
    j=(float)(sin((rpsi-rphi)/2)*sin(rtheta/2));
    k=(float)(sin((rpsi+rphi)/2)*cos(rtheta/2));

    if(fabs(w)<QEPS)
        w=0.0f;

    if(fabs(i)<QEPS)
        i=0.0f;

    if(fabs(j)<QEPS)
        j=0.0f;

    if(fabs(k)<QEPS)
        k=0.0f;

}

vector<float> Quaternion::getEulerAngles()
{
    float psi;
    float theta;
    float phi;

    if(fabs(i)<QEPS)
        i=0;

    if(fabs(j)<QEPS)
            j=0;

    if(fabs(k)<QEPS)
            k=0;

    if(i==0.0f && j==0.0f)
    {
        if(k==0.0f)
        {
            phi=0.0f;
            psi=0.0f;
            theta=0.0f;
        }
        else if(w<0.0f)
        {
            psi=-2.0f*asin(k);
            phi=0.0f;
            theta=0.0f;
        }
        else
        {
            psi=2*asin(k);
            theta=0.0f;
            phi=0.0f;
        }
    }
    else if(w==0.0f & k==0.0f)
    {
         if(j==-1.0f || j==1.0f)
         {
            psi=0.0f;
            theta=acos(-i*i-j*j);
            phi=-2*asin(j);
         }
        else if(i==0.0f || j==0.0f)
        {
            phi=0.0f;
            psi=0.0f;
            theta=acos(-i*i-j*j);
        }
        else if(i<0.0f & j<0.0f)
        {
            psi=0.0f;
            theta=acos(-i*i-j*j);
            phi=2.0f*asin(j);
        }
        else if(i<0.0f)
        {
            psi=0.0f;
            theta=acos(-i*i-j*j);
            phi=-2.0f*asin(-j);
        }
        else
        {
            psi=0.0f;
            theta=acos(-i*i-j*j);
            phi=-2.0f*asin(j);
        }
    }
    else
    {
        psi = atan2((i*k+j*w),-(j*k-i*w));
        theta = acos(-i*i-j*j+k*k+w*w);
        phi = atan2((i*k-j*w),(j*k+i*w));
    }

    psi = psi*180.0f/M_PI;
    theta = theta*180.0f/M_PI;
    phi = phi*180.0f/M_PI;

    // TODO: this is a nasty hack - needs to be figured out why the nan is sometimes returned
    if(isnan(theta))
        theta=180.0f;

    return vector<float>{phi,psi,theta};
}

void Quaternion::getAxisAngle(float& angle, vector<float>& axis)
{
    axis.resize(3);

    if (w == 1.0f)
    {
        angle = 0.0f;
        axis = { 0.0f, 0.0f, 0.0f };
    }
    else
    {
        angle = 2.0f * acos(w);
        axis[0] = i / sqrt(1.0f - w*w);
        axis[1] = j / sqrt(1.0f - w*w);
        axis[2] = k / sqrt(1.0f - w*w);
    }
}

float Quaternion::jointAngle(Quaternion& q1,Quaternion& q2)
{
    float angle = q1.w*q2.w + q1.i*q2.i +q1.j*q2.j +q1.k*q2.k;

    angle = acos(angle);

    return angle;
}

void Quaternion::operator*(float value)
{
    w*=value;
    i*=value;
    j*=value;
    k*=value;
}

Quaternion Quaternion::operator+(Quaternion& q)
{
    Quaternion res;

    res.w*=q.w + w;
    res.i*=q.i + i ;
    res.j*=q.j + j;
    res.k*=q.k + k;

    return res;
}

Quaternion Quaternion::mult(Quaternion& q1, float value)
{
    Quaternion q(q1);

    q.w*=value;
    q.i*=value;
    q.j*=value;
    q.k*=value;

    return q;
}

Quaternion Quaternion::add(Quaternion& q1, Quaternion& q2)
{
    Quaternion q;

    q.w = q1.w + q2.w;
    q.i = q1.i + q2.i;
    q.j = q1.j + q2.j;
    q.k = q1.k + q2.k;

    return q;
}

Quaternion Quaternion::mult(Quaternion& q1,Quaternion& q2)
{
    Quaternion q;

    q.w = q1.w*q2.w - q1.i*q2.i - q1.j*q2.j - q1.k*q2.k;
    q.i = q1.w*q2.i + q1.i*q2.w - q1.j*q2.k + q1.k*q2.j;
    q.j = q1.w*q2.j + q1.i*q2.k + q1.j*q2.w - q1.k*q2.i;
    q.k = q1.w*q2.k - q1.i*q2.j + q1.j*q2.i + q1.k*q2.w;

    q.normalize();

    return q;
}

vector<float> Quaternion::rotate(Quaternion& q,vector<float>& point)
{

    float coeff1 = 2*(q.i*point[0] + q.j*point[1] + q.k*point[2]);
    float coeff2 = q.w*q.w - (q.i*q.i + q.j*q.j + q.k*q.k);

    vector<float> result(3);

    result[0] = coeff1*q.i + coeff2*point[0] + 2*q.w*(q.j*point[2] - q.k*point[1]);
    result[1] = coeff1*q.j + coeff2*point[1] + 2*q.w*(q.k*point[0] - q.i*point[2]);
    result[2] = coeff1*q.k + coeff2*point[2] + 2*q.w*(q.i*point[1] - q.j*point[0]);

    return result;
}

vector<float> Quaternion::pointRotate(vector<float>& point)
{
    return rotate(*this,point);
}

void Quaternion::normalize()
{
    float size = sqrt(w*w + i*i + j*j + k*k);

    w/=size;
    i/=size;
    j/=size;
    k/=size;
}

bool Quaternion::isEqual(Quaternion& q)
{
    if(fabs(i-q.i)<QEPS && fabs(j-q.j)<QEPS && fabs(k-q.k)<QEPS &&fabs(w-q.w)<QEPS)
        return true;
    else
        return false;

}

void Quaternion::getAngularDistance(Quaternion q1, Quaternion q2, float& angleInDegrees, float& distance)
{
    q1.normalize();
    q2.normalize();

    float sum = q1.i*q2.i+q1.j*q2.j+q1.k*q2.k+q1.w*q2.w;

    angleInDegrees = acos(2.0f*sum*sum - 1.0f)*180.f/M_PI;

    if(isnan(angleInDegrees))
    {
        angleInDegrees = 0.0f;
    }

    distance = 1.0f - sum*sum;
}

void Quaternion::getRotationMatrix(vector<float>& matrix)
{
    matrix.resize(9);

    matrix[0] = 1.0f - 2.0f*(j*j + k*k);
    matrix[1] = 2.0f*(i*j - w*k);
    matrix[2] = 2.0f*(i*k + w*j);

    matrix[3] = 2.0f*(i*j + w*k);
    matrix[4] = 1.0f - 2.0f*(i*i + k*k);
    matrix[5] = 2.0f*(j*k - w*i);

    matrix[6] = 2.0f*(i*k - w*j);
    matrix[7] = 2.0f*(j*k + w*i);
    matrix[8] = 1.0f - 2.0f*(i*i + j*j);


    /*matrix[0] = 1.0f - 2.0f*(j*j + k*k);
    matrix[3] = 2.0f*(i*j - w*k);
    matrix[6] = 2.0f*(i*k + w*j);

    matrix[1] = 2.0f*(i*j + w*k);
    matrix[4] = 1.0f - 2.0f*(i*i + k*k);
    matrix[7] = 2.0f*(j*k - w*i);

    matrix[2] = 2.0f*(i*k - w*j);
    matrix[5] = 2.0f*(j*k + w*i);
    matrix[8] = 1.0f - 2.0f*(i*i + j*j);*/
}

void Quaternion::slerp(vector<Quaternion>& output, Quaternion& q1, Quaternion& q2, float stepSize,float interpolationLimit)
{

    float angler=jointAngle(q1,q2);

    if(fabs(angler/(2*M_PI))<EPS)
        return;
    double s;
    for( s = 0; s<interpolationLimit-EPS; s += stepSize)
    {
        float coeff1 = sin((1-s)*angler)/sin(angler);
        float coeff2 = sin(s*angler)/sin(angler);

        Quaternion nq1(q1);
        Quaternion nq2(q2);
        nq1*coeff1;
        nq2*coeff2;
        Quaternion res(add(nq1,nq2));
        output.push_back(res);
    }
}

Quaternion Quaternion::inverseQuaternion()
{
    return Quaternion(w,-i,-j,-k);
}

void Quaternion::generateConeRotations(vector<Quaternion>& rotations,float coneAngle, float coneSampling, float inplaneAngle, float inplaneSampling)
{
    rotations.push_back(Quaternion(1.0f,0.0f,0.0f,0.0f));
    //cout << rotations.back().w << " " << rotations.back().i << " " << rotations.back().j << " " << rotations.back().k << endl;

    // Prepare inplane rotations
    vector<Quaternion> phiRotations;

    for(float phi = -inplaneAngle/2.0f; phi <= inplaneAngle/2.0f; phi+=inplaneSampling)
    {
        phiRotations.push_back(Quaternion(phi,0.0f,0.0f));
        if(abs(phi)>EPS)
        {
            rotations.push_back(Quaternion::mult(phiRotations.back(),rotations[0]));
            //cout << rotations.back().w << " " << rotations.back().i << " " << rotations.back().j << " " << rotations.back().k << endl;
        }
    }
    // Prepare cone rotations
    float coneOpening = coneAngle*M_PI/ (2.0f*180.f);

    Quaternion cq1(coneOpening, vector<float>{1.0f, 0.0f, 0.0f});
    Quaternion cq2(0.0f, vector<float>{1.0f, 0.0f, 0.0f});
   
    float verticalSamplingStep = 1.0f / floor(coneAngle/2.0f / coneSampling);

    vector<Quaternion> verticalSamplingPoints;
    Quaternion::slerp(verticalSamplingPoints,cq1,cq2,verticalSamplingStep);

    vector<float> radii;
    for (size_t v = 0; v < verticalSamplingPoints.size(); v++)
    {
        float axisAngle;
        vector<float> axis;
        verticalSamplingPoints[v].getAxisAngle(axisAngle, axis);
        radii.push_back(sin(axisAngle));
    }

    float arcLength =  2.0f*M_PI*coneSampling / 360.0f;

    for(size_t v=0; v<verticalSamplingPoints.size(); v++)
    {
        Quaternion hq1(0.0f, vector<float>{0.0f, 0.0f, 1.0f});
        Quaternion hq2(M_PI/2.0f, vector<float>{0.0f, 0.0f, 1.0f});

        double numberOfSamples = 2.0*M_PI*fabs((double)radii[v]) / arcLength;
        float horizontalSamplingStep = 4.0f / round(numberOfSamples);
        
        vector<Quaternion> horizontalSamplingPoints;

        Quaternion::slerp(horizontalSamplingPoints,hq1,hq2,horizontalSamplingStep,4.0);

        for(size_t h =0; h < horizontalSamplingPoints.size(); h++)
        {
            Quaternion coneQuat=Quaternion::mult(verticalSamplingPoints[v],horizontalSamplingPoints[h]);
            for(size_t p=0; p< phiRotations.size();p++)
            {
                rotations.push_back(Quaternion::mult(phiRotations[p],coneQuat));
                //rotations.push_back(Quaternion::mult(coneQuat,phiRotations[p]));
                //cout << rotations.back().w << " " << rotations.back().i << " " << rotations.back().j << " " << rotations.back().k << endl;
            }
        }
    }
}

void Quaternion::generateConeRotations(vector<Quaternion>& rotations, vector<float>& phiAngles, float coneAngle, float coneSampling, float inplaneAngle, float inplaneSampling)
{
    rotations.push_back(Quaternion(1.0f,0.0f,0.0f,0.0f));

    for(float phi = -inplaneAngle/2.0f; phi <= inplaneAngle/2.0f; phi+=inplaneSampling)
    {
        phiAngles.push_back(phi);
    }

    // Prepare cone rotations
    float coneOpening = coneAngle*M_PI/ (2.0f*180.f);

    Quaternion cq1(coneOpening, vector<float>{1.0f, 0.0f, 0.0f});
    Quaternion cq2(0.0f, vector<float>{1.0f, 0.0f, 0.0f});

    float verticalSamplingStep = 1.0f / floor(coneAngle/2.0f / coneSampling);

    vector<Quaternion> verticalSamplingPoints;
    Quaternion::slerp(verticalSamplingPoints,cq1,cq2,verticalSamplingStep);

    vector<float> radii;
    for (size_t v = 0; v < verticalSamplingPoints.size(); v++)
    {
        float axisAngle;
        vector<float> axis;
        verticalSamplingPoints[v].getAxisAngle(axisAngle, axis);
        radii.push_back(sin(axisAngle));
    }

    float arcLength =  2.0f*M_PI*coneSampling / 360.0f;

    for(size_t v=0; v<verticalSamplingPoints.size(); v++)
    {
        Quaternion hq1(0.0f, vector<float>{0.0f, 0.0f, 1.0f});
        Quaternion hq2(M_PI/2.0f, vector<float>{0.0f, 0.0f, 1.0f});

        double numberOfSamples = 2.0*M_PI*fabs((double)radii[v]) / arcLength;
        float horizontalSamplingStep = 4.0f / round(numberOfSamples);

        vector<Quaternion> horizontalSamplingPoints;

        Quaternion::slerp(horizontalSamplingPoints,hq1,hq2,horizontalSamplingStep,4.0);

        for(size_t h =0; h < horizontalSamplingPoints.size(); h++)
        {
            rotations.push_back(Quaternion::mult(verticalSamplingPoints[v],horizontalSamplingPoints[h]));

        }
    }
}
