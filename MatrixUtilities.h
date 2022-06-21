//
// Created by burak on 22.03.2022.
//



#ifndef HW_1_MATRIXUTILITIES_H
#define HW_1_MATRIXUTILITIES_H

#include <cstdlib>
#include <iostream>
#include <chrono>
#include <vector>
#include <random>
#include <cmath>

#define MAXRANDOMVALUE 10.0
#define MINRANDOMVALUE 1.0

namespace matrixutilies {

    bool isPrime(int n)
{
    if (n <= 1)
        return false;
    for (int i = 2; i < n; i++)
        if (n % i == 0)
            return false;
 
    return true;
}

void FindDivisionCoeff(int nprocs, int x[3]){
    int divider=2;
    x[0]=1;
    x[1]=1;
    x[2]=1;

    int ind=0;
    while(1 <nprocs)
    {
        while(nprocs%divider==0){
            nprocs/=divider;

            x[ind]*=divider;

            ind++;
        }
        divider++;
        while (!isPrime(divider))
        {
            divider++;
        }

        

    }



}


    template<typename A>
    void CreateVector(A *&vector,int size){
        vector= new A[size];
        std::uniform_real_distribution<double>  unif(MINRANDOMVALUE,MAXRANDOMVALUE);

        std::default_random_engine re(std::chrono::system_clock::now().time_since_epoch().count());

        for (int i = 0; i < size; i++)
        {
            vector[i]=static_cast<A>(unif(re));
        }

    }

    template<typename A>
    A* CreateVector(int size){
        A* vector= new A[size];
        std::uniform_real_distribution<double>  unif(MINRANDOMVALUE,MAXRANDOMVALUE);    
        std::default_random_engine re(std::chrono::system_clock::now().time_since_epoch().count());

        for (int i = 0; i < size; i++)
        {
            vector[i]=static_cast<A>(unif(re));
            

        }
        return vector;

    }
     template<typename A>
    A* CreateLocalVector(const int size,const int seed){
        A* vector= new A[size];
        std::uniform_real_distribution<double>  unif(MINRANDOMVALUE,MAXRANDOMVALUE);
        //std::default_random_engine re;
        std::default_random_engine re(seed);

        for (int i = 0; i < size; i++)
        {
            vector[i]=static_cast<A>(unif(re));
            //vector[i]=static_cast<A>(i);


        }
        return vector;

    }

     template<typename A>
    A* CreateLocalVector(const int size,const int seed,int padding){
        A* vector= new A[size];
        std::uniform_real_distribution<double>  unif(MINRANDOMVALUE,MAXRANDOMVALUE);
        //std::default_random_engine re;
        std::default_random_engine re(seed);

        for (int i = 0; i < size-padding; i++)
        {
            vector[i]=static_cast<A>(unif(re));
            //vector[i]=static_cast<A>(i);


        }
        for (size_t i = size-padding; i < size; i++)
        {
            vector[i]=0;
        
        }
        
        return vector;

    }




    template<typename T>
    T* CreateMatrix1Dynamic(const int size) {
        T* matrix;
        matrix = new T [size*size];

        return matrix;
    }

    template<typename T>
    T* CreateMatrix1Dynamic(const int m,const int n) {
        T* matrix;
        matrix = new T [n*m];

        return matrix;
    }

    template<typename T>
    T** CreateRandomMatrix2DDynamic(const int size) {
        T** matrix;
        matrix = new T* [size];
        for (int i = 0; i < size; i++)
        {
            matrix[i] = new T[size];
        }
        std::uniform_real_distribution<double>  unif(MINRANDOMVALUE,MAXRANDOMVALUE);
        //std::default_random_engine re;
        std::default_random_engine re(std::chrono::system_clock::now().time_since_epoch().count());

        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                matrix[i][j]=static_cast<T>(unif(re));
               // matrix[i][j]=static_cast<T>(j);

            }
        }

        return matrix;
    }


   template<typename T>
    T** CreateRandomMatrix2DDynamic(const int m,const int n) {
        T** matrix;
        matrix = new T* [m];
        for (int i = 0; i < m; i++)
        {
            matrix[i] = new T[n];
        }

        std::uniform_real_distribution<double>  unif(MINRANDOMVALUE,MAXRANDOMVALUE);
        std::default_random_engine re(std::chrono::system_clock::now().time_since_epoch().count());
        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < n; ++j) {
                matrix[i][j]=static_cast<T>(unif(re));
            }
        }
        return matrix;
    }

   template<typename T>
    T** CreateRandomMatrix2DDynamicPadding(const int size,const int padding) {
        T** matrix;
        matrix = new T* [size];
        for (int i = 0; i < size; i++)
        {
            matrix[i] = new T[size];
        }
        std::uniform_real_distribution<double>  unif(MINRANDOMVALUE,MAXRANDOMVALUE);
        //std::default_random_engine re;
        std::default_random_engine re(std::chrono::system_clock::now().time_since_epoch().count());

        for (int i = 0; i < size-padding; ++i) {
            for (int j = 0; j < size; ++j) {
                matrix[i][j]=static_cast<T>(unif(re));
               // matrix[i][j]=static_cast<T>(j);

            }
        }
        for (int i = size-padding; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
               // matrix[i][j]=static_cast<T>(unif(re));
                matrix[i][j]=0;

            }
        }
        return matrix;

    }



    template<typename A>
    void PrintMatrix( A **matrix,int size,int rank){
        std::cout<< "\n--------------Rank: "<<rank<<"---------------\n";

    
        for (int i = 0; i < size; i++)
        {
            for (int j = 0; j < size; j++)
            {

                std::cout<<matrix[i][j]<<"\t";

            }
            std::cout<<"\n";

        }
        std::cout<< "\n------------------------------------"<<std::endl;


    }

    template<typename A>
    void PrintVector( const A *vector,int size,int rank){

        std::cout<< "\n--------------Rank: "<<rank<<"---------------\n";

        for (int j = 0; j < size; j++)
        {

            std::cout<<vector[j]<<"\t";

        }
        std::cout<< "\n-------------------------------------"<<std::endl;
    }

    template<typename A>
    A* MultipleMatrixVector(A *vector, A** matrix,int size){
        A* result = new A[size];
        for (size_t i = 0; i < size; i++)
        {
            result[i]=0;
            for (int j = 0; j < size; j++)
            {
                result[i]+=matrix[i][j]*vector[j];
            }

        }

        return result;
    }

        template<typename A>
    void MultipleMatrixVector(A *vector_in, A** matrix,A* result,int size){
      
        for (size_t i = 0; i < size; i++)
        {
            result[i]=0;  
            int count=0;         
            for (int j = 0; j < size; j++)
            {
                
                result[i]+=matrix[i][j]*vector_in[j];
            }

        }
    }
    template<typename A>
    void SpecialMulitpleMatrix(A *vector_in,A* result,int size){

        int plate_length =std::cbrt(size) ;
       // plate_length =20;

   /* for (int i = 0; i < size; i++)
    {

        result[i]=0;
        result[i]+= 6*vector_in[i];


        if (i%plate_length != 0) //left
        {
            result[i]+= -1*vector_in[i-1];
        }

        if (i%plate_length != plate_length-1) //right
        {

            result[i]+= -1*vector_in[i+1];
        }

        if (std::floor((i%(plate_length*plate_length))/plate_length)!=0 ) //back
        {
            result[i]+= -1*vector_in[i-plate_length];
        }
        
        if (std::floor((i%(plate_length*plate_length))/plate_length)!=plate_length-1 ) //front
        {
            
             result[i]+= -1*vector_in[i+plate_length];
        }

        if (std::floor(i/(plate_length*plate_length)) != 0) //Top
        {
             result[i]+= -1*vector_in[i-plate_length*plate_length];
        }


       if (std::floor(i/(plate_length*plate_length)) != plate_length-1) //bottom
        {
             result[i]+= -1*vector_in[i+plate_length*plate_length];

        }
    }*/

    
     
    
        for (int z = 1; z < plate_length-1; z++)
        {
            int Z_= z*plate_length*plate_length;
            for (int y = 1; y < plate_length-1; y++)
            {   int Y_= y*plate_length;
                for (int x = 1; x < plate_length-1; x++)
                {                                        
                        result[(x+Z_+Y_)] = 6* vector_in[(x+Z_+Y_)] - vector_in[(x+1+Z_+Y_)] 
                        - vector_in[(x-1+Z_+Y_)] - vector_in[(x+Z_+Y_+plate_length)] - vector_in[(x+Z_+Y_-plate_length)]
                        - vector_in[(x+Z_+Y_+plate_length*plate_length)] - vector_in[(x+Z_+Y_-plate_length*plate_length)] ;
                }
                
            }
            
        }

   {
       for (int z = 0, i=0; z < plate_length; z++)
       {
           int x=0;
           for (int y = 0; y < plate_length; y++)
           {
                i=z*plate_length*plate_length+y*plate_length+x;

                result[i]=0;
                result[i]+= 6*vector_in[i];

                result[i]+= -1*vector_in[i+1];

                if (std::floor((i%(plate_length*plate_length))/plate_length)!=0 ) //back
                {
                    result[i]+= -1*vector_in[i-plate_length];
                }
                
                if (std::floor((i%(plate_length*plate_length))/plate_length)!=plate_length-1 ) //front
                {
                    
                    result[i]+= -1*vector_in[i+plate_length];
                }

                if (std::floor(i/(plate_length*plate_length)) != 0) //Top
                {
                    result[i]+= -1*vector_in[i-plate_length*plate_length];
                }


                if (std::floor(i/(plate_length*plate_length)) != plate_length-1) //bottom
                {
                    result[i]+= -1*vector_in[i+plate_length*plate_length];

                }

           }
           
       }    
    
    }    

        for (int z = 0, i=0; z < plate_length; z++)
       {
           int x=plate_length-1;
           for (int y = 0; y < plate_length; y++)
           {    
                i=z*plate_length*plate_length+y*plate_length+x;

                result[i]=0;
                result[i]+= 6*vector_in[i];

                    result[i]+= -1*vector_in[i-1];

                if (std::floor((i%(plate_length*plate_length))/plate_length)!=0 ) //back
                {
                    result[i]+= -1*vector_in[i-plate_length];
                }
                
                if (std::floor((i%(plate_length*plate_length))/plate_length)!=plate_length-1 ) //front
                {
                    
                    result[i]+= -1*vector_in[i+plate_length];
                }

                if (std::floor(i/(plate_length*plate_length)) != 0) //Top
                {
                    result[i]+= -1*vector_in[i-plate_length*plate_length];
                }


            if (std::floor(i/(plate_length*plate_length)) != plate_length-1) //bottom
                {
                    result[i]+= -1*vector_in[i+plate_length*plate_length];

                }

           }
           
       }    

        for (int y = 0, i=0; y < plate_length; y++)
       {
           int z=0;
           for (int x = 0; x < plate_length; x++)
           {
                i=z*plate_length*plate_length+y*plate_length+x;

                result[i]=0;
                result[i]+= 6*vector_in[i];

                if (i%plate_length != 0) //left
                {
                    result[i]+= -1*vector_in[i-1];
                }

                if (i%plate_length != plate_length-1) //right
                {

                    result[i]+= -1*vector_in[i+1];
                }

                if (std::floor((i%(plate_length*plate_length))/plate_length)!=0 ) //back
                {
                    result[i]+= -1*vector_in[i-plate_length];
                }
                
                if (std::floor((i%(plate_length*plate_length))/plate_length)!=plate_length-1 ) //front
                {
                    
                    result[i]+= -1*vector_in[i+plate_length];
                }


                    result[i]+= -1*vector_in[i+plate_length*plate_length];


           }
           
       } 



    for (int y = 0, i=0; y < plate_length; y++)
       {
           int z=plate_length-1;
           for (int x = 0; x < plate_length; x++)
           {
                i=z*plate_length*plate_length+y*plate_length+x;

                result[i]=0;
                result[i]+= 6*vector_in[i];

                if (i%plate_length != 0) //left
                {
                    result[i]+= -1*vector_in[i-1];
                }

                if (i%plate_length != plate_length-1) //right
                {

                    result[i]+= -1*vector_in[i+1];
                }

                if (std::floor((i%(plate_length*plate_length))/plate_length)!=0 ) //back
                {
                    result[i]+= -1*vector_in[i-plate_length];
                }
                
                if (std::floor((i%(plate_length*plate_length))/plate_length)!=plate_length-1 ) //front
                {
                    
                    result[i]+= -1*vector_in[i+plate_length];
                }


                    result[i]+= -1*vector_in[i-plate_length*plate_length];


           }
           
       }  


        for (int z = 0, i=0; z < plate_length; z++)
       {
           int y=0;
           for (int x = 0; x < plate_length; x++)
           {
                i=z*plate_length*plate_length+y*plate_length+x;

                result[i]=0;
                result[i]+= 6*vector_in[i];

                if (i%plate_length != 0) //left
                {
                    result[i]+= -1*vector_in[i-1];
                }

                if (i%plate_length != plate_length-1) //right
                {

                    result[i]+= -1*vector_in[i+1];
                }


                    
                    result[i]+= -1*vector_in[i+plate_length];


            if (std::floor(i/(plate_length*plate_length)) != 0) //Top
                {
                    result[i]+= -1*vector_in[i-plate_length*plate_length];
                }


            if (std::floor(i/(plate_length*plate_length)) != plate_length-1) //bottom
                {
                    result[i]+= -1*vector_in[i+plate_length*plate_length];

                }


           }
           
       } 

        for (int z = 0, i=0; z < plate_length; z++)
       {
           int y=plate_length-1;
           for (int x = 0; x < plate_length; x++)
           {
                i=z*plate_length*plate_length+y*plate_length+x;

                result[i]=0;
                result[i]+= 6*vector_in[i];

                if (i%plate_length != 0) //left
                {
                    result[i]+= -1*vector_in[i-1];
                }

                if (i%plate_length != plate_length-1) //right
                {

                    result[i]+= -1*vector_in[i+1];
                }


                    
                    result[i]+= -1*vector_in[i-plate_length];


            if (std::floor(i/(plate_length*plate_length)) != 0) //Top
                {
                    result[i]+= -1*vector_in[i-plate_length*plate_length];
                }


            if (std::floor(i/(plate_length*plate_length)) != plate_length-1) //bottom
                {
                    result[i]+= -1*vector_in[i+plate_length*plate_length];

                }


           }
           
       } 
    
    }  
    

   /* 
template<typename A>
void RedBlackMatrixMultiplication(A *X,A* result,int size){
    int plate_length =std::cbrt(size) ;
    int middle= size-std::floor(size/2)+1+ (size%2==1 ? +1 : 0);  
   // int ind= middle > i? i*2: (i-middle)*2+1 ;
    for (int i = 0,ind =0; i < middle; i++, ind+=2) //Red
    {
        result[i]=0;
        result[i]+= 6*X[i];        
        if (ind%plate_length != 0) //left
            result[i]+= -1*X[middle+((i-1)-1)/2];
        if (ind%plate_length != plate_length-1) //right
            result[i]+= -1*X[middle+((i+1)-1)/2];
        if (std::floor((ind%(plate_length*plate_length))/plate_length)!=0 ) //back
            result[i]+= -1*X[middle + ((i-plate_length)-1)/2];
        if (std::floor((ind%(plate_length*plate_length))/plate_length)!=plate_length-1 ) //front            
             result[i]+= -1*X[middle + ((i+plate_length)-1)/2];
        if (std::floor(ind/(plate_length*plate_length)) != 0) //Top        {
             result[i]+= -1*X[middle + ((i-plate_length*plate_length)-1)/2];
        if (std::floor(ind/(plate_length*plate_length)) != plate_length-1) //bottom
             result[i]+= -1*X[middle + ((i+plate_length*plate_length)-1)/2];

    }

    for (int i = middle,ind =1; i < size; i++, ind+=2) //Black
    {
        result[i]=0;
        result[i]+= 6*X[i];        
        if (ind%plate_length != 0) //left
            result[i]+= -1*X[(i-1)/2];
        if (ind%plate_length != plate_length-1) //right
            result[i]+= -1*X[(i+1)/2];
        if (std::floor((ind%(plate_length*plate_length))/plate_length)!=0 ) //back
            result[i]+= -1*X[ (i-plate_length)/2];
        if (std::floor((ind%(plate_length*plate_length))/plate_length)!=plate_length-1 ) //front            
             result[i]+= -1*X[ (i+plate_length)/2];
        if (std::floor(ind/(plate_length*plate_length)) != 0) //Top        {
             result[i]+= -1*X[ (i-plate_length*plate_length)/2];
        if (std::floor(ind/(plate_length*plate_length)) != plate_length-1) //bottom
             result[i]+= -1*X[(i+plate_length*plate_length)/2];

    }   
    

}

*/
    template<typename A>
    A MultipleVectorVector(A *vector1, A* vector2,int size){
        A result =0.0;
        for (int i = 0; i < size; i++)
        {        
            result+=vector1[i]*vector2[i];
        }

        return result;


    
    }
    template<typename A>
    void MultipleScalerVector(A scalerValue, A* vector,int size){
        for (int i = 0; i < size; i++)
        {
            vector[i]*=scalerValue;
        }       
        
    }

    template<typename A>
    A* MultipleScalerVector(A scalerValue, A* vector,int size,bool inplace){
         A* result= new A[size];
        for (int i = 0; i < size; i++)
        {
            result[i]=scalerValue*vector[i];
        }  
        return  result; 
        
    }

    template<typename A>
    void SumVectorVector(A* vectorout, A* vectorin,int size){
        for (int i = 0; i < size; i++)
        {
            vectorout[i]+=vectorin[i];
        }       
        
    }
    template<typename A>
    void SubtractVectorVector(A* vectorout, A* vectorin,int size){
        
        for (int i = 0; i < size; i++)
        {
            vectorout[i]-=vectorin[i];
        } 
           

    }



    template<typename A>
    double Norm(A *vector,int size,bool squared=true){
        double result{};
        for (size_t i = 0; i < size; i++)
        {
            result+=vector[i]*vector[i];

        }
        
        if (squared)
        {
            return sqrt(result);

        }
        else{
            return result;
        }
        
        

    }

    template<typename A>
    void DeleteMatrix(A **matrix,int size){

    if (matrix==nullptr)
    {
        return;
    }
    
        for (int i = 0; i < size; i++)
        {
            delete[] matrix[i];
        }
        delete[] matrix;

    }

    template<typename A>
    void DeleteVector(A *vector){

        if (vector != nullptr)
        {
            delete[] vector;
        }
        
        

    }

};


#endif //HW_1_MATRIXUTILITIES_H
