#include <complex>
#include <iostream>
#include "vectorRoutines.h"

using namespace std;

    void vr::normalize(vector<float>& vec)
    {
        float norm = sqrt(dot_product(vec, vec));
        vec_div(vec, norm);
    }

    float vr::dot_product(vector<float>& vec1, vector<float>& vec2)
    {
        float result=0;

        for(size_t i=0; i< vec1.size(); i++)
        {
            result =  result + vec1[i]*vec2[i];
        }

        return result;

    }

    //template<typename typeName>
    vector<float> vr::operator*(vector<float>& left, vector<float>& right)
    {
        // check sizes of vectors first
        if (left.size() != right.size())
        {
            cout << endl << "ERROR: Vectors do not have the same size!" << endl;
            exit(EXIT_FAILURE);
        }

        vector<float> multResult;
        multResult.resize(left.size());
        for (size_t i = 0; i < multResult.size(); i++)
        {
            multResult[i] = left[i] * right[i];
        }

        return multResult;
    }

    vector<complex<float>> vr::operator*(vector<complex<float>>& left, vector<complex<float>>& right)
    {
        // check sizes of vectors first
        if (left.size() != right.size())
        {
            cout << endl << "ERROR: Vectors do not have the same size!" << endl;
            exit(EXIT_FAILURE);
        }

        vector<complex<float>> multResult;
        multResult.resize(left.size());
        for (size_t i = 0; i < multResult.size(); i++)
        {
            multResult[i] = left[i] * right[i];
        }

        return multResult;
    }

   // template<typename float>
    vector<float> vr::operator+(vector<float>& vec, float value)
    {
        vector<float> addResult;
        addResult.resize(vec.size());

        for (size_t i = 0; i < vec.size(); i++)
        {
            addResult[i] = vec[i] + value;
        }

        return addResult;
    }

    //template<typename typeName>
    void vr::operator*(vector<float>& vec, float value)
    {
        for (size_t i = 0; i < vec.size(); i++)
        {
            vec[i] = vec[i] * value;
        }
    }

    vector<float> vr::operator-(vector<float>& vec)
    {
        vector<float> result(vec.size());
        
        for (size_t i = 0; i < vec.size(); i++)
            result[i] = -vec[i];

        return result;
    }

    //template<typename typeName>
    vector<float> vr::operator+(vector<float>& left, vector<float>& right)
    {
        // check sizes of vectors first
        if (left.size() != right.size())
        {
            cout << endl << "ERROR: Vectors do not have the same size!" << endl;
            exit(EXIT_FAILURE);
        }

        vector<float> addResult;
        addResult.resize(left.size());
        for (size_t i = 0; i < addResult.size(); i++)
        {
            addResult[i] = left[i] + right[i];
        }

        return addResult;
    }

    //template<typename typeName>
    vector<float> vr::operator-(vector<float>& left, vector<float>& right)
    {
        // check sizes of vectors first
        if (left.size() != right.size())
        {
            cout << endl << "ERROR: Vectors do not have the same size!" << endl;
            exit(EXIT_FAILURE);
        }

        vector<float> addResult;
        addResult.resize(left.size());
        for (size_t i = 0; i < addResult.size(); i++)
        {
            addResult[i] = left[i] - right[i];
        }

        return addResult;
    }

    //template<typename typeName>
    /*vector<float> vr::operator/(vector<float>& left, vector<float>& right)
    {
        // check sizes of vectors first
        if (left.size() != right.size())
        {
            cout << endl << "ERROR: Vectors do not have the same size!" << endl;
            exit(EXIT_FAILURE);
        }

        vector<float> divResult;
        divResult.resize(left.size());
        for (size_t i = 0; i < divResult.size(); i++)
        {
            if (right[i] != 0.0f)
                divResult[i] = left[i] / right[i];
            else
            {
                cout << endl << "ERROR: Vector division by zero!" << endl;
                exit(EXIT_FAILURE);
            }
        }

        return divResult;
    }*/

    vector<float> vr::vec_div(vector<float>& left, vector<float>& right, float eps)
    {
        // check sizes of vectors first
       if (left.size() != right.size())
       {
           cout << endl << "ERROR: Vectors do not have the same size!" << endl;
           exit(EXIT_FAILURE);
       }

       vector<float> divResult;
       divResult.resize(left.size());
       for (size_t i = 0; i < divResult.size(); i++)
       {
           if ( fabs(right[i]) < eps )
               divResult[i] = 0.0f;
           else
           {
               divResult[i] = left[i] / right[i];
           }
       }

       return divResult;
    }

    void vr::vec_div(vector<float>& vec, float value)
    {
        if (value == 0)
        {
            cout << endl << "ERROR: Vector division by zero!" << endl;
            exit(EXIT_FAILURE);
        }

        for (size_t i = 0; i < vec.size(); i++)
        {
            vec[i] = vec[i] / value;
        }
    }


    //template<typename typeName>
   /* vector<float> vr::operator/(vector<float>& left, float divisor)
    {
        if (divisor == 0)
        {
            cout << endl << "ERROR: Vector division by zero!" << endl;
            exit(EXIT_FAILURE);
        }

        vector<float> divResult;
        divResult.resize(left.size());

        for (size_t i = 0; i < divResult.size(); i++)
        {
            divResult[i] = left[i] / divisor;
        }

        return divResult;
    }*/

   // template<typename typeName>
    float vr::vec_sum(vector<float>& vec)
    {
        float sum = 0;
        for (size_t i = 0; i < vec.size(); i++)
        {
            sum += vec[i];
        }

        return sum;
    }

    //template<typename typeName>
    void vr::vec_mult(vector<float>& vec, float value)
    {
        for (size_t i = 0; i < vec.size(); i++)
        {
            vec[i] = vec[i] * value;
        }
    }

    //template<typename typeName>
    void vr::vec_mult(vector<float>& vec1, vector<float>& vec2)
    {
        for (size_t i = 0; i < vec1.size(); i++)
        {
            vec1[i] = vec1[i] * vec2[i];
        }
    }

    void vr::vec_conj_mult(vector<complex<float>>& vec1, vector<complex<float>>& vec2)
    {
        for (size_t i = 0; i < vec1.size(); i++)
        {
            vec1[i] = conj(vec1[i]) *vec2[i];
        }
    }

    void vr::vec_conj_mult(vector<complex<float>>& vec1, vector<complex<float>>& vec2, vector<complex<float>>& res)
    {
        if(res.size()==0)
        {
            res.resize(vec1.size());
        }

        for (size_t i = 0; i < vec1.size(); i++)
        {
            res[i] = conj(vec1[i]) *vec2[i];
        }
    }

    void vr::vec_mult(vector<complex<float>>& vec1, vector<complex<float>>& vec2)
    {
        for (size_t i = 0; i < vec1.size(); i++)
        {
            vec1[i] = vec1[i] *vec2[i];
        }
    }

   // template<typename typeName>
    void vr::vec_mult(vector<complex<float>>& vec, float value)
    {
        for (size_t i = 0; i < vec.size(); i++)
        {
            vec[i] = vec[i] * value;
        }
    }

 	void vr::vec_mult(vector<complex<double>>& vec, double value)
    {
        for (size_t i = 0; i < vec.size(); i++)
        {
            vec[i] = vec[i] * value;
        }
    }

    void vr::vec_mult(vector<complex<double>>& vec1, vector<float>& vec2, vector<float>& vec3)
    {
        // check sizes of vectors first
        if (vec1.size() != vec2.size() || vec1.size() != vec3.size())
        {
            cout << endl << "ERROR: Vectors do not have the same size!" << endl;
            exit(EXIT_FAILURE);
        }

        for (size_t i = 0; i < vec1.size(); i++)
        {
            vec1[i] = complex<double>(vec1[i].real() * vec2[i] *vec3[i], vec1[i].imag()* vec2[i] *vec3[i]);
        }
    }

    void vr::vec_div(vector<complex<double>>& vec1, vector<float>& vec2, float eps)
    {
        // check sizes of vectors first
          if (vec1.size() != vec2.size())
          {
              cout << endl << "ERROR: Vectors do not have the same size!" << endl;
              exit(EXIT_FAILURE);
          }

          for (size_t i = 0; i < vec1.size(); i++)
          {
              if ( fabs(vec2[i]) < eps )
                  vec1[i] = complex<double>(0.0,0.0);
              else
              {
                  vec1[i] =  complex<double>(vec1[i].real()/ (double)vec2[i], vec1[i].imag()/(double)vec2[i]);
              }
          }
    }

    //template<typename typeName>
   /* vector<complex<float>> vr::operator*(vector<complex<float>>& left, vector<float>& right)
    {
        // check sizes of vectors first
        if (left.size() != right.size())
        {
            cout << endl << "ERROR: Vectors do not have the same size!" << endl;
            exit(EXIT_FAILURE);
        }

        vector<complex<float>> multResult;
        multResult.resize(left.size());

        for (size_t i = 0; i < multResult.size(); i++)
        {
            multResult[i] = complex<float>(left[i].real() * right[i], left[i].imag() * right[i]);
        }

        return multResult;
    }*/

    //template<typename typeName>
    vector<complex<double>> vr::operator*(vector<complex<double>>& left, vector<double>& right)
    {
        // check sizes of vectors first
        if (left.size() != right.size())
        {
            cout << endl << "ERROR: Vectors do not have the same size!" << endl;
            exit(EXIT_FAILURE);
        }

        vector<complex<double>> multResult;
        multResult.resize(left.size());

        for (size_t i = 0; i < multResult.size(); i++)
        {
            multResult[i] = complex<float>(left[i].real() * right[i], left[i].imag() * right[i]);
        }

        return multResult;
    }

    void vr::vec_mult(vector<complex<double>>& left, vector<float>& right)
    {
        // check sizes of vectors first
        if (left.size() != right.size())
        {
            cout << endl << "ERROR: Vectors do not have the same size!" << endl;
            exit(EXIT_FAILURE);
        }

        for (size_t i = 0; i < left.size(); i++)
        {
            left[i] = complex<float>(left[i].real() * (double)right[i], left[i].imag() * (double)right[i]);
        }

    }

    float vr::vec_mean(vector<float>& data)
    {
        float vectorSum = vr::vec_sum(data);
        float vectorMean = vectorSum / (float)data.size();

        return vectorMean;
    }

    float vr::vec_variance(vector<float>& data, float mean)
    {
        float variance=0.0f;

        for(size_t i=0;i<data.size();i++)
        {
            float diff = data[i] - mean;
            variance = variance + diff*diff;
        }

        variance = variance/ (float)data.size();

        return variance;

    }

    void  vr::vec_normalize(vector<float>& data)
    {
        float mean = vec_mean(data);
        float stdDev = sqrt(vec_variance(data, mean));
        
        for (size_t i = 0; i<data.size(); i++)
        {
            data[i] = (data[i] - mean) / stdDev;
        }
    }

    /*vector<complex<double>> vr::operator*(vector<complex<double>>& left, vector<float>& right)
    {
        // check sizes of vectors first
        if (left.size() != right.size())
        {
            cout << endl << "ERROR: Vectors do not have the same size!" << endl;
            exit(EXIT_FAILURE);
        }

        vector<complex<double>> multResult;
        multResult.resize(left.size());

        for (size_t i = 0; i < multResult.size(); i++)
        {
            multResult[i] = complex<float>(left[i].real() * (double)right[i], left[i].imag() * (double)right[i]);
        }

        return multResult;
    }*/

