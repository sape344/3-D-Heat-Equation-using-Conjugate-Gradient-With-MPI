#include "MatrixUtilities.h"

struct HeatConditions
{
    double top{};
    double bottom{};
    double left{};
    double right{};
    double front{};
    double back{};
    double initial_heat{};

};

template<typename T>
T* ConjugateGradient(T** A,T*B,int size,double tol=1e-12);


template<typename T>
void CreateHeatMap(T**& A, T*&b,HeatConditions& cond,int plate_length);


int main(int argc, char const *argv[])
{

    double** A= nullptr;
    double* b= nullptr;
    int plate_length= 25;

    HeatConditions cond;
    cond.top=100.0;
    cond.left=100.0;
    cond.bottom=100.0;
    cond.right=100.0;
    cond.front=100.0;
    cond.back=100.0;
    cond.initial_heat=100.0;

    CreateHeatMap(A, b,cond, plate_length);
    auto t_start = std::chrono::high_resolution_clock::now();
    auto x = ConjugateGradient(A,b, plate_length*plate_length*plate_length);
    auto t_end = std::chrono::high_resolution_clock::now();
    double sum=0.0;
    std::cout<<"Wall clock : "<<std::chrono::duration<double, std::milli>(t_end-t_start).count()<<"\n";
    matrixutilies::DeleteMatrix(A,plate_length*plate_length*plate_length);
    matrixutilies::DeleteVector(b);
    matrixutilies::DeleteVector(x);


    







    int* arr=new int[9];

    for(int i=0;i<9;i++){
        arr[i]=i-4;
        std::cout<<i-4<<"\t";
    }
    std::cout<<"\n";
    std::cout<<"Result: "<< matrixutilies::Norm(arr,9)<<"\n";
    return 0;
}


template<typename T>
T* ConjugateGradient(T** A,T*B,int size,double tol)
{


    T* x =matrixutilies::CreateVector<T>(size); //initial guess

    auto residual = matrixutilies::MultipleMatrixVector(x,A,size); // local
    matrixutilies::SubtractVectorVector(residual,B,size); // local
    auto direction = matrixutilies::MultipleScalerVector(-1.0,residual,size,true); // local
    auto residual_norm= matrixutilies::Norm(residual,size); // MPI reduce to get global norm

    int num_iter=0;
    //Iteration is started 
    auto A_direction=matrixutilies::MultipleMatrixVector(direction,A,size); // local
    auto sum_resual= matrixutilies::MultipleVectorVector(residual,residual,size); // MPI_reduce to get global sum_resual
    auto alpha =sum_resual / matrixutilies::MultipleVectorVector(direction,A_direction,size); // MPI_Reduce to get global alpha
    for (int i = 0; i < size; i++)
    {
        x[i]+= alpha*direction[i]; //local
        residual[i]+= alpha*A_direction[i];
    }

    auto beta = matrixutilies::MultipleVectorVector(residual,residual,size)/sum_resual; // MPI_Reduce to get global beta

    for (int i = 0; i < size; i++)
    {
        direction[i] = -residual[i] +beta*direction[i]; //local
    }
    
    num_iter+=1;
    residual_norm= matrixutilies::Norm(residual,size); // MPI_Reduce to get global residual_norm

    while(residual_norm > tol)
    {
      


        //auto t_start = std::chrono::high_resolution_clock::now();
        //matrixutilies::MultipleMatrixVector(direction,A,A_direction,size); //en çok süreyi burası harcıyor
        matrixutilies::SpecialMulitpleMatrix(direction,A_direction,size); //en çok süreyi burası harcıyor //local

        //auto t_end = std::chrono::high_resolution_clock::now();        
        //std::cout<<"inside Wall clock : "<<std::chrono::duration<double, std::milli>(t_end-t_start).count()<<"\n";
        
        sum_resual= matrixutilies::MultipleVectorVector(residual,residual,size); // MPI_Reduce to get global sum_resual
        alpha =sum_resual / matrixutilies::MultipleVectorVector(direction,A_direction,size); // MPI_Reduce to get global alpha



        for (int i = 0; i < size; i++)
        {
            x[i]+= alpha*direction[i]; //local
            residual[i]+= alpha*A_direction[i]; //local
        }


        beta = matrixutilies::MultipleVectorVector(residual,residual,size)/sum_resual;
        for (int i = 0; i < size; i++)
        {
            direction[i] = -residual[i] +beta*direction[i];
        }
    
         num_iter+=1;
        residual_norm= matrixutilies::Norm(residual,size);

        std::cout<<"Residual norm is "<<residual_norm<<std::endl;
        std::cout<<"Step size "<<num_iter<<"\n";

    }
    for (int i = 0; i < size; i++)
    {
        std::cout<<"x["<<i<<"]="<<x[i]<<"\n";
    }
    

    std::cout<<"Step size "<<num_iter<<"\n";
    matrixutilies::DeleteVector(direction);
    matrixutilies::DeleteVector(A_direction);
    matrixutilies::DeleteVector(residual);



    return x;

}


template<typename T>
void CreateHeatMap(T**& A, T*&b,HeatConditions& cond,int plate_length){  
    
    const double coeff= 1;  
    int size= plate_length*plate_length*plate_length;
    A = new T* [size];
    for (int i = 0; i < size; i++)
    {
        A[i] = new T[size];
    }

    b = new T [size];

    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            A[i][j]=0.0;
        }
        
    }


    

    for (int i = 0; i < size; i++)
    {
        b[i]=0.0;
        A[i][i]=1-(1-6*coeff);

        if (i%plate_length != 0) //left
        {
            A[i][i-1]=-coeff;
        }
        else{
            b[i]+=cond.left;
        }

        if (i%plate_length != plate_length-1) //right
        {
            A[i][i+1]=-coeff;
        }
        else{
            b[i]+=cond.right;
        }

        if (std::floor((i%(plate_length*plate_length))/plate_length)!=0 ) //back
        {
            A[i][i-plate_length]=-coeff;
        }
        else{
            b[i]+=cond.back;
        }
        
        if (std::floor((i%(plate_length*plate_length))/plate_length)!=plate_length-1 ) //front
        {
            A[i][i+plate_length]=-coeff;
        }
        else{
            b[i]+=cond.front;
        }

        if (std::floor(i/(plate_length*plate_length)) != 0) //Top
        {
            A[i][i-plate_length*plate_length]=-coeff;
        }
        else
        {
             b[i]+=cond.top;
        }

        if (std::floor(i/(plate_length*plate_length)) != plate_length-1) //bottom
        {
            A[i][i+plate_length*plate_length]=-coeff;
        }
        else
        {
             b[i]+=cond.bottom;
        }
        

    }
    



}