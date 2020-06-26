#pragma once
#include <vector>
#include <complex>
#include <iostream>
using namespace std;

namespace vr
{
    //template<typename typeName>
    vector<float> operator*(vector<float>& left, vector<float>& right);
    vector<complex<float>> operator*(vector<complex<float>>& left, vector<complex<float>>& right);

    //template<typename typeName>
    vector<float> operator+(vector<float>& vec, float value);

   // template<typename typeName>
    void operator*(vector<float>& vec, float value);

    vector<float> operator-(vector<float>& vec);

    //template<typename typeName>
    vector<float> operator+(vector<float>& left, vector<float>& right);

    //template<typename typeName>
    vector<float> operator-(vector<float>& left, vector<float>& right);

   // template<typename typeName>
    //vector<float> operator/(vector<float>& left, vector<float>& right);

    vector<float> vec_div(vector<float>& left, vector<float>& right, float eps);
    void vec_div(vector<float>& vec, float value);
    void vec_div(vector<complex<double>>& vec1, vector<float>& vec2, float eps);

   // template<typename typeName>
   // vector<float> operator/(vector<float>& left, float divisor);

   // template<typename typeName>
    float vec_sum(vector<float>& vec);

    //template<typename typeName>
    //void vec_mult(vector<float>& vec, typeName value);

   // template<typename typeName>
    void vec_mult(vector<float>& vec1, vector<float>& vec2);

    void vec_conj_mult(vector<complex<float>>& vec1, vector<complex<float>>& vec2, vector<complex<float>>& res);

    void vec_conj_mult(vector<complex<float>>& vec1, vector<complex<float>>& vec2);

    void vec_mult(vector<complex<float>>& vec1, vector<complex<float>>& vec2);

    void vec_mult(std::vector<std::complex<double>>& vec1, std::vector<float>& vec2, std::vector<float>& vec3);

    //template<typename typeName>
    void vec_mult(vector<float>& vec, float value);
    void vec_mult(vector<complex<float>>& vec, float value);
    void vec_mult(vector<complex<double>>& vec, double value);
    void vec_mult(vector<complex<float>>& vec, float value);
    void vec_mult(vector<complex<double>>& left, vector<float>& right);

    float vec_mean(vector<float>& data);
    float vec_variance(vector<float>& data, float mean);
    void  vec_normalize(vector<float>& data);

    float dot_product(vector<float>& vec1, vector<float>& vec2);
    
    void normalize(vector<float>& vec);
    //template<typename typeName>
   // vector<complex<float>> operator*(vector<complex<float>>& left, vector<float>& right);

    //template<typename typeName>
    vector<complex<double>> operator*(vector<complex<double>>& left, vector<double>& right);

    //vector<complex<double>> operator*(vector<complex<double>>& left, vector<float>& right);



   // template vector<float> operator+<float>(vector<float>& vec, float value);
   // template vector<float> operator*<float>(vector<float >& left, vector<float>& right);
   // template vector<float> operator*<float>(vector<float>& left, vector<float>& right);
   // template vector<float> operator+<float>(vector<float>& vec, float value);
   // template vector<float> operator+<float>(vector<float>& left, vector<float>& right);
   // template vector<float> operator-<float>(vector<float>& left, vector<float>& right);
   // template vector<float> operator/<float>(vector<float>& left, vector<float>& right);
   // template vector<float> operator/<float>(vector<float>& left, float divisor);
   // template float vec_sum<float>(vector<float>& vec);
   // template void vec_mult<float>(vector<float>& vec, float value);
   // template void vec_mult<double>(vector<complex<double>>& vec, double value);
    //template void vec_mult<float>(vector<float>& vec1, vector<float>& vec2);
    //template void vec_mult<float>(vector<complex<float>>& vec, float value);
    //template vector<complex<float>> operator*<float>(vector<complex<float>>& left, vector<float>& right);
    //template vector<complex<double>> operator*<float>(vector<complex<double>>& left, vector<float>& right);
}
