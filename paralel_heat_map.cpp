#include "mpi.h"
#include "MatrixUtilities.h"


#define LEFT 0
#define RIGHT 1
#define UP 2
#define DOWN 3
#define FRONT 4
#define BACK 5

#define SIZE 0
#define DISP 1

#define X 0
#define Y 1
#define Z 2

#define SEND 0
#define RECV 1

typedef double DefType;
  
    
template<typename A>
void RedBlackMatrixMultiplication(A *x,A* result, A* comm_buffer,const int comm_info[6][2] ,int local_sizes[3]);
void UpdateGhostPoints(int neighbors[6], int comm_buffer_info[6][2],DefType* comm_buffer ,const MPI_Datatype& XY_TYPE,const MPI_Datatype& XZ_TYPE,const MPI_Datatype& YZ_TYPE, const MPI_Comm& cart_comm);
void UpdateGhostPoints(int neighbors[6],int comm_buffer_info[6][2] ,DefType** comm_buffer ,const MPI_Datatype& XY_TYPE,const MPI_Datatype& XZ_TYPE,const MPI_Datatype& YZ_TYPE, const MPI_Comm& cart_comm);
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






int index(const int& x, const int& y, const int& z, const int local_sizes[3]);
int index(const int& index, const int& size );
void FillLocalB(DefType* local_b, int local_sizes[3], int neighbors[6] ,HeatConditions& cond);
template<typename A>
void UpdateCommBuffer(const A *Xk,A* comm_buffer,const int comm_info[6][2] ,int local_sizes[3]);


int main(int argc, char *argv[])
{

    if (argc!=2)
    {
        std::cout<<"Please Enter a matrix size\n";
        return -1;
    }
    MPI_Init(&argc, &argv);
    int rank, nprocs;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    const int N= std::atoi(argv[1]); //matrix size is NxNxN

    int cart_dims[3];
    matrixutilies::FindDivisionCoeff(nprocs,cart_dims);

    //MPI Cartesian topology
    int periods[3]={0,0,0} ,coords[3],reorder=1,neighbors[6];
    MPI_Comm cart_comm;
    MPI_Cart_create(MPI_COMM_WORLD, 3, cart_dims, periods, reorder, &cart_comm);    
    MPI_Cart_coords(cart_comm, rank, 3, coords);

    MPI_Cart_shift(cart_comm, 0, 1, &neighbors[LEFT], &neighbors[RIGHT]); //left and right
    MPI_Cart_shift(cart_comm, 1, 1, &neighbors[UP], &neighbors[DOWN]); //up and down
    MPI_Cart_shift(cart_comm, 2, 1, &neighbors[FRONT], &neighbors[BACK]); //front and back

    int local_sizes[3];
    local_sizes[X] = N/cart_dims[X];
    local_sizes[Y] = N/cart_dims[Y];
    local_sizes[Z] = N/cart_dims[Z];
    int local_sizes_total=local_sizes[X]*local_sizes[Y]*local_sizes[Z];

    // Communication Types
    MPI_Datatype XY_TYPE,XZ_TYPE,YZ_TYPE;
	MPI_Type_contiguous(local_sizes[X] *local_sizes[Y], MPI_DOUBLE, &XY_TYPE); // Front and back
	MPI_Type_commit(&XY_TYPE);
    MPI_Type_contiguous(local_sizes[X] *local_sizes[Z], MPI_DOUBLE, &XZ_TYPE); //Up and Down
    MPI_Type_commit(&XZ_TYPE);
    MPI_Type_contiguous(local_sizes[Y] *local_sizes[Z], MPI_DOUBLE, &YZ_TYPE);  // Left and Right
    MPI_Type_commit(&YZ_TYPE);

    // Allocate memory for local matrices
    DefType** comm_buffer= new DefType* [2]; 
    DefType* local_b= new DefType[local_sizes_total];
    DefType* local_x= new DefType[local_sizes_total];
    std::fill(local_x, local_x+local_sizes_total, 0.0);

    

    //Data size and displacement
    int comm_buffer_info[6][2];  
    comm_buffer_info[LEFT][DISP]= 0;
    comm_buffer_info[LEFT][SIZE]= neighbors[LEFT] != -2 ? local_sizes[Y]*local_sizes[Z]: 0;
    comm_buffer_info[RIGHT][DISP]= comm_buffer_info[LEFT][DISP]+ comm_buffer_info[LEFT][SIZE];    
    comm_buffer_info[RIGHT][SIZE]= neighbors[RIGHT] != -2 ? local_sizes[Y]*local_sizes[Z]: 0;
    comm_buffer_info[UP][DISP]= comm_buffer_info[RIGHT][DISP]+ comm_buffer_info[RIGHT][SIZE];
    comm_buffer_info[UP][SIZE]= neighbors[UP] != -2 ? local_sizes[X]*local_sizes[Z]: 0;
    comm_buffer_info[DOWN][DISP]= comm_buffer_info[UP][DISP]+ comm_buffer_info[UP][SIZE];
    comm_buffer_info[DOWN][SIZE]= neighbors[DOWN] != -2 ? local_sizes[X]*local_sizes[Z]: 0;
    comm_buffer_info[FRONT][DISP]= comm_buffer_info[DOWN][DISP]+ comm_buffer_info[DOWN][SIZE];
    comm_buffer_info[FRONT][SIZE]= neighbors[FRONT] != -2 ? local_sizes[X]*local_sizes[Y]: 0;
    comm_buffer_info[BACK][DISP]= comm_buffer_info[FRONT][DISP]+ comm_buffer_info[FRONT][SIZE];  
    comm_buffer_info[BACK][SIZE]= neighbors[BACK] != -2 ? local_sizes[X]*local_sizes[Y]: 0;

    
    //Allocate memory for communication buffer
    comm_buffer[SEND] = new DefType[comm_buffer_info[BACK][SIZE]+ comm_buffer_info[BACK][DISP]]; 
    comm_buffer[RECV] = new DefType[comm_buffer_info[BACK][SIZE]+ comm_buffer_info[BACK][DISP]]; 

    std::fill(comm_buffer[SEND], comm_buffer[SEND]+comm_buffer_info[BACK][SIZE]+ comm_buffer_info[BACK][DISP], 0.0);
    std::fill(comm_buffer[RECV], comm_buffer[RECV]+comm_buffer_info[BACK][SIZE]+ comm_buffer_info[BACK][DISP], 0.0);

    HeatConditions cond{100.0, 0.0, 0.0, 0.0,100.0, 100.0, 0.0}; //Initial conditions

    bool cube_inside = comm_buffer_info[BACK][SIZE]+ comm_buffer_info[BACK][DISP] == (2*local_sizes[X]*local_sizes[Y]) + 2*local_sizes[X]*local_sizes[Z] + 2*local_sizes[Y]*local_sizes[Z] ? true : false;
    //Heat Map Craeating
    if (cube_inside)
    {
       //the process inside cube
        std::fill(local_b, local_b+local_sizes_total, cond.initial_heat);
    }
    else
    {
       //process on the surface of the cube
        FillLocalB( local_b, local_sizes,  neighbors , cond);

        
    }



   
    DefType *residual = new DefType[local_sizes_total];
    DefType* A_direction = new DefType[local_sizes_total];


    const double tolerence = 1e-7;
    const double tolerence_m4 = 1e-4;
    auto t_start = std::chrono::high_resolution_clock::now();

    //Conjugate Gradient Method Started 
    
    UpdateCommBuffer(local_x,comm_buffer[SEND],comm_buffer_info , local_sizes);
    UpdateGhostPoints( neighbors, comm_buffer_info,comm_buffer , XY_TYPE,  XZ_TYPE,  YZ_TYPE,  cart_comm);
    RedBlackMatrixMultiplication(local_x, residual,comm_buffer[RECV],comm_buffer_info , local_sizes);


    matrixutilies::SubtractVectorVector(residual, local_b, local_sizes_total);
    DefType* direction = matrixutilies::MultipleScalerVector<DefType>(-1.0,residual,local_sizes_total,true);
    double residual_norm = matrixutilies::Norm(residual, local_sizes_total,false);
    MPI_Allreduce(&residual_norm, &residual_norm, 1, MPI_DOUBLE, MPI_SUM,  MPI_COMM_WORLD);
    int num_iterations = 0;

    
    UpdateCommBuffer(direction,comm_buffer[SEND],comm_buffer_info , local_sizes);
    UpdateGhostPoints( neighbors, comm_buffer_info,comm_buffer , XY_TYPE,  XZ_TYPE,  YZ_TYPE,  cart_comm);
    RedBlackMatrixMultiplication(direction, A_direction,comm_buffer[RECV],comm_buffer_info , local_sizes);



    double sum_resual= matrixutilies::MultipleVectorVector(residual,residual,local_sizes_total);
    MPI_Allreduce(&sum_resual, &sum_resual, 1, MPI_DOUBLE, MPI_SUM,  MPI_COMM_WORLD);
    double dir_adir_mulp= matrixutilies::MultipleVectorVector(direction,A_direction,local_sizes_total);
    MPI_Allreduce(&dir_adir_mulp, &dir_adir_mulp, 1, MPI_DOUBLE, MPI_SUM,  MPI_COMM_WORLD);

    double alpha = sum_resual/dir_adir_mulp;
    
    for (int i = 0; i < local_sizes_total; i++)
    {
        local_x[i]+= alpha*direction[i]; //local
        residual[i]+= alpha*A_direction[i];
    }


    double beta = matrixutilies::MultipleVectorVector(residual,residual,local_sizes_total);
    MPI_Allreduce(&beta, &beta, 1, MPI_DOUBLE, MPI_SUM,  cart_comm);
    beta = beta/sum_resual;

    for (int i = 0; i < local_sizes_total; i++)
        direction[i] = -residual[i] +beta*direction[i]; 
    num_iterations++;
    residual_norm = matrixutilies::Norm(residual, local_sizes_total,false);
    MPI_Allreduce(&residual_norm, &residual_norm, 1, MPI_DOUBLE, MPI_SUM,  cart_comm);
    residual_norm=std::sqrt(residual_norm);

    while (residual_norm > tolerence)
    {
        
        UpdateCommBuffer(direction,comm_buffer[SEND],comm_buffer_info , local_sizes);
        UpdateGhostPoints( neighbors, comm_buffer_info,comm_buffer , XY_TYPE,  XZ_TYPE,  YZ_TYPE,  cart_comm);
        RedBlackMatrixMultiplication(direction, A_direction,comm_buffer[RECV],comm_buffer_info , local_sizes);

        sum_resual= matrixutilies::MultipleVectorVector(residual,residual,local_sizes_total);
        MPI_Allreduce(&sum_resual, &sum_resual, 1, MPI_DOUBLE, MPI_SUM,  cart_comm);
        dir_adir_mulp= matrixutilies::MultipleVectorVector(direction,A_direction,local_sizes_total);
        MPI_Allreduce(&dir_adir_mulp, &dir_adir_mulp, 1, MPI_DOUBLE, MPI_SUM,  cart_comm);
        alpha = sum_resual/dir_adir_mulp;
        for (int i = 0; i < local_sizes_total; i++)
        {
            local_x[i]+= alpha*direction[i]; //local
            residual[i]+= alpha*A_direction[i];
        }
        beta = matrixutilies::MultipleVectorVector(residual,residual,local_sizes_total);
        MPI_Allreduce(&beta, &beta, 1, MPI_DOUBLE, MPI_SUM,  cart_comm);
        beta = beta/sum_resual;
        for (int i = 0; i < local_sizes_total; i++)
        {
            direction[i] = -residual[i] +beta*direction[i]; 
        }
        num_iterations++;
        residual_norm = matrixutilies::Norm(residual, local_sizes_total,false);
        MPI_Allreduce(&residual_norm, &residual_norm, 1, MPI_DOUBLE, MPI_SUM,  cart_comm);
        residual_norm=std::sqrt(residual_norm);



    }
    //Conjugate Gradient Method Ended


    if (rank ==0 )
    {
        auto t_end = std::chrono::high_resolution_clock::now();

        std::cout<<"Cube size ,"<<N<<", Core size ,"<<nprocs<<", Wall clock : ,"<<std::chrono::duration<double, std::milli>(t_end-t_start).count()<<", ms"<<"\n";

    }
    






    delete[] A_direction;
    delete[] residual;
    delete[] comm_buffer[RECV];
    delete[] comm_buffer[SEND];
    delete[] comm_buffer;
    delete[] local_b;
    delete[] local_x;


    MPI_Type_free(&XY_TYPE);
    MPI_Type_free(&XZ_TYPE);
    MPI_Type_free(&YZ_TYPE);
    MPI_Comm_free(&cart_comm);


    MPI_Finalize();

    return 0;
}



void UpdateGhostPoints(int neighbors[6],int comm_buffer_info[6][2] ,DefType* comm_buffer ,const MPI_Datatype& XY_TYPE,const MPI_Datatype& XZ_TYPE,const MPI_Datatype& YZ_TYPE, const MPI_Comm& cart_comm)
{

    if ( neighbors[LEFT] != -2)
    {
    //    std::cout<<"I am "<<rank<<" and my left neighbor is "<<neighbors[LEFT]<<std::endl;
        MPI_Sendrecv_replace(comm_buffer+comm_buffer_info[LEFT][DISP], 1, YZ_TYPE, neighbors[LEFT], LEFT, neighbors[LEFT], RIGHT, cart_comm,  MPI_STATUS_IGNORE);

    }

    if ( neighbors[RIGHT] != -2)
    {
    //  std::cout<<"I am "<<rank<<" and my RIGHT neighbor is "<<neighbors[RIGHT]<<std::endl;
        MPI_Sendrecv_replace(comm_buffer+comm_buffer_info[RIGHT][DISP], 1, YZ_TYPE, neighbors[RIGHT], RIGHT, neighbors[RIGHT], LEFT, cart_comm,  MPI_STATUS_IGNORE);
    }


    if ( neighbors[UP] != -2)
    {
    //   std::cout<<"I am "<<rank<<" and my UP neighbor is "<<neighbors[UP]<<std::endl;
        MPI_Sendrecv_replace(comm_buffer+comm_buffer_info[UP][DISP], 1, XZ_TYPE, neighbors[UP], UP, neighbors[UP], DOWN, cart_comm, MPI_STATUS_IGNORE);
    }

    if ( neighbors[DOWN] != -2){
    //   std::cout<<"I am "<<rank<<" and my DOWN neighbor is "<<neighbors[DOWN]<<std::endl;
        MPI_Sendrecv_replace(comm_buffer+comm_buffer_info[DOWN][DISP], 1, XZ_TYPE, neighbors[DOWN], DOWN, neighbors[DOWN], UP, cart_comm,  MPI_STATUS_IGNORE);
    }

        
    if ( neighbors[FRONT] != -2){
    //    std::cout<<"I am "<<rank<<" and my FRONT neighbor is "<<neighbors[FRONT]<<std::endl;
        MPI_Sendrecv_replace(comm_buffer+comm_buffer_info[FRONT][DISP], 1, XY_TYPE, neighbors[FRONT], FRONT, neighbors[FRONT], BACK, cart_comm,  MPI_STATUS_IGNORE);
    }

    if ( neighbors[BACK] != -2){
    //  std::cout<<"I am "<<rank<<" and my BACK neighbor is "<<neighbors[BACK]<<std::endl;
        MPI_Sendrecv_replace(comm_buffer+comm_buffer_info[BACK][DISP], 1, XY_TYPE, neighbors[BACK], BACK, neighbors[BACK], FRONT, cart_comm,  MPI_STATUS_IGNORE);
    }
}

void UpdateGhostPoints(int neighbors[6],int comm_buffer_info[6][2] ,DefType** comm_buffer ,const MPI_Datatype& XY_TYPE,const MPI_Datatype& XZ_TYPE,const MPI_Datatype& YZ_TYPE, const MPI_Comm& cart_comm)
{

    MPI_Request requests[6];
    MPI_Status status;
    if ( neighbors[LEFT] != -2){
        MPI_Isend(comm_buffer[SEND]+comm_buffer_info[LEFT][DISP],1,YZ_TYPE,neighbors[LEFT],LEFT,cart_comm,&requests[0]);

    }   


    if ( neighbors[RIGHT] != -2)
    {
        MPI_Isend(comm_buffer[SEND]+comm_buffer_info[RIGHT][DISP],1,YZ_TYPE,neighbors[RIGHT],RIGHT,cart_comm,&requests[1]);
    }


    if ( neighbors[UP] != -2)
    {
        MPI_Isend(comm_buffer[SEND]+comm_buffer_info[UP][DISP],1,XZ_TYPE,neighbors[UP],UP,cart_comm,&requests[2]);
       
    }

    if ( neighbors[DOWN] != -2){
        MPI_Isend(comm_buffer[SEND]+comm_buffer_info[DOWN][DISP],1,XZ_TYPE,neighbors[DOWN],DOWN,cart_comm,&requests[3]);
    }

        
    if ( neighbors[FRONT] != -2){

        MPI_Isend(comm_buffer[SEND]+comm_buffer_info[FRONT][DISP],1,XY_TYPE,neighbors[FRONT],FRONT,cart_comm,&requests[4]);

   }

    if ( neighbors[BACK] != -2){

        MPI_Isend(comm_buffer[SEND]+comm_buffer_info[BACK][DISP],1,XY_TYPE,neighbors[BACK],BACK,cart_comm,&requests[5]);

    }

    if ( neighbors[LEFT] != -2){
        MPI_Recv(comm_buffer[RECV]+comm_buffer_info[LEFT][DISP],1,YZ_TYPE,neighbors[LEFT],RIGHT,cart_comm,&status);

    }   


    if ( neighbors[RIGHT] != -2)
    {
        MPI_Recv(comm_buffer[RECV]+comm_buffer_info[RIGHT][DISP],1,YZ_TYPE,neighbors[RIGHT],LEFT,cart_comm,&status);
    }


    if ( neighbors[UP] != -2)
    {
        MPI_Recv(comm_buffer[RECV]+comm_buffer_info[UP][DISP],1,XZ_TYPE,neighbors[UP],DOWN,cart_comm,&status);
       
    }

    if ( neighbors[DOWN] != -2){
        MPI_Recv(comm_buffer[RECV]+comm_buffer_info[DOWN][DISP],1,XZ_TYPE,neighbors[DOWN],UP,cart_comm,&status);
    }

        
    if ( neighbors[FRONT] != -2){

        MPI_Recv(comm_buffer[RECV]+comm_buffer_info[FRONT][DISP],1,XY_TYPE,neighbors[FRONT],BACK,cart_comm,&status);

   }

    if ( neighbors[BACK] != -2){

        MPI_Recv(comm_buffer[RECV]+comm_buffer_info[BACK][DISP],1,XY_TYPE,neighbors[BACK],FRONT,cart_comm,&status);

    }

   
    



  /*  MPI_Request requests[12] ;
    if ( neighbors[LEFT] != -2){
        MPI_Isend(comm_buffer[SEND]+comm_buffer_info[LEFT][DISP],1,YZ_TYPE,neighbors[LEFT],LEFT,cart_comm,&requests[0]);
        MPI_Irecv(comm_buffer[RECV]+comm_buffer_info[LEFT][DISP],1,YZ_TYPE,neighbors[LEFT],RIGHT,cart_comm,&requests[1]);

    }   


    if ( neighbors[RIGHT] != -2)
    {
        MPI_Isend(comm_buffer[SEND]+comm_buffer_info[RIGHT][DISP],1,YZ_TYPE,neighbors[RIGHT],RIGHT,cart_comm,&requests[2]);
        MPI_Irecv(comm_buffer[RECV]+comm_buffer_info[RIGHT][DISP],1,YZ_TYPE,neighbors[RIGHT],LEFT,cart_comm,&requests[3]);
    }


    if ( neighbors[UP] != -2)
    {
        MPI_Isend(comm_buffer[SEND]+comm_buffer_info[UP][DISP],1,XZ_TYPE,neighbors[UP],UP,cart_comm,&requests[4]);
        MPI_Irecv(comm_buffer[RECV]+comm_buffer_info[UP][DISP],1,XZ_TYPE,neighbors[UP],DOWN,cart_comm,&requests[5]);
       
    }

    if ( neighbors[DOWN] != -2){
        MPI_Isend(comm_buffer[SEND]+comm_buffer_info[DOWN][DISP],1,XZ_TYPE,neighbors[DOWN],DOWN,cart_comm,&requests[6]);
        MPI_Irecv(comm_buffer[RECV]+comm_buffer_info[DOWN][DISP],1,XZ_TYPE,neighbors[DOWN],UP,cart_comm,&requests[7]);
    }

        
    if ( neighbors[FRONT] != -2){

        MPI_Isend(comm_buffer[SEND]+comm_buffer_info[FRONT][DISP],1,XY_TYPE,neighbors[FRONT],FRONT,cart_comm,&requests[8]);
        MPI_Irecv(comm_buffer[RECV]+comm_buffer_info[FRONT][DISP],1,XY_TYPE,neighbors[FRONT],BACK,cart_comm,&requests[9]);

   }

    if ( neighbors[BACK] != -2){

        MPI_Isend(comm_buffer[SEND]+comm_buffer_info[BACK][DISP],1,XY_TYPE,neighbors[BACK],BACK,cart_comm,&requests[10]);
        MPI_Irecv(comm_buffer[RECV]+comm_buffer_info[BACK][DISP],1,XY_TYPE,neighbors[BACK],FRONT,cart_comm,&requests[11]);

    }

    for (int i = 0; i < 12; i+=1)
    {
        MPI_Wait(requests[i], MPI_STATUS_IGNORE);
    }*/
    


}


/*
int index(const int& x, const int& y, const int& z, const int local_sizes[3] ){
    return x + local_sizes[X] * (y + local_sizes[Y] * z);
}
*/

int index(const int& index, const int& middle ){
    //int middle= size-std::floor(size/2)+1+ (size%2==1 ? +1 : 0); 
    return index;
   /* int result;
    if((index) %2==0 ){
        result= index/2;
    }else{
        result= ((index-1)/2)+middle-1;

    }

    return result;*/
}

int index(const int& x, const int& y, const int& z, const int local_sizes[3] ){
    int size = local_sizes[0]*local_sizes[1]*local_sizes[2];
    int middle= size-std::floor(size/2)+1+ (size%2==1 ? +1 : 0); 
    int index_1= x + local_sizes[X] * (y + local_sizes[Y] * z);

    return index(index_1,middle);

}


int BlackIndex(const int& index, const int& middle ){

    return ((index-1)/2)+middle-1;
    
}

int RedIndex(const int& index ){

    return (index)/2;
    
}



void FillLocalB(DefType* local_b, int local_sizes[3], int neighbors[6] ,HeatConditions& cond){ 
    std::fill(local_b, local_b+(local_sizes[0]*local_sizes[1]*local_sizes[2]), cond.initial_heat);

    if (neighbors[LEFT]==-2)
    {
        for (int j = 0; j < local_sizes[Z]; j++)
        {
            for (int k = 0; k < local_sizes[Y]; k++)
            {
                local_b[index(0,k,j,local_sizes)]+=cond.left;
            }
        }
    }
    if (neighbors[RIGHT]==-2)
    {
        for (int j = 0; j < local_sizes[Z]; j++)
        {
            for (int k = 0; k < local_sizes[Y]; k++)
            {
                local_b[index(local_sizes[X]-1,k,j,local_sizes)]+=cond.right;
            }
        }
    }
    if (neighbors[UP]==-2)
    {
        for (int j = 0; j < local_sizes[Z]; j++)
        {
            for (int k = 0; k < local_sizes[X]; k++)
            {
                local_b[index(k,local_sizes[Y]-1,j,local_sizes)]+=cond.top;
            }
        }
    }
    if (neighbors[DOWN]==-2)
    {
        for (int j = 0; j < local_sizes[Z]; j++)
        {
            for (int k = 0; k < local_sizes[X]; k++)
            {
                local_b[index(k,0,j,local_sizes)]+=cond.bottom;
            }
        }
    }
    if (neighbors[FRONT]==-2)
    {
        for (int j = 0; j < local_sizes[Y]; j++)
        {
            for (int k = 0; k < local_sizes[X]; k++)
            {
                local_b[index(k,j,0,local_sizes)]+=cond.front;
            }
        }
    }
    if (neighbors[BACK]==-2)
    {
        for (int j = 0; j < local_sizes[X]; j++)
        {
            for (int k = 0; k < local_sizes[Y]; k++)
            {
                local_b[index(k,j,local_sizes[Z]-1,local_sizes)]+=cond.back;
            }
        }
    }              

}




template<typename A>
void RedBlackMatrixMultiplication(A *Xk,A* result, A* comm_buffer,const int comm_info[6][2] ,int local_sizes[3]){
    
    int size = local_sizes[0]*local_sizes[1]*local_sizes[2];
    int middle= size-std::floor(size/2)+1+ (size%2==1 ? +1 : 0);  
    for ( int rb = 0; rb <= 1; rb+=1) //Red Black
    {
        for (int z = 1; z < local_sizes[Z]-1; z++)
        {
            int Z_= z*local_sizes[X]*local_sizes[Y];
            for (int y = 1; y < local_sizes[Y]-1; y++)
            {   int Y_= y*local_sizes[X];
                for (int x = 1; x < local_sizes[X]-1; x++)
                {
                    if((x+y+z)%2 == rb)
                    {                        
                        result[index(x+Z_+Y_,middle)] = 6* Xk[index(x+Z_+Y_,middle)] - Xk[index(x+1+Z_+Y_,middle)] 
                        - Xk[index(x-1+Z_+Y_,middle)] - Xk[index(x+Z_+Y_+local_sizes[X],middle)] - Xk[index(x+Z_+Y_-local_sizes[X],middle)]
                        - Xk[index(x+Z_+Y_+local_sizes[X]*local_sizes[Y],middle)] - Xk[index(x+Z_+Y_-local_sizes[X]*local_sizes[Y],middle)] ;
                    }
                }
                
            }
            
        }
        
    }

    //LEFT 
    {
        int x=0;
        int ind=0;
        for ( int rb = 0; rb <=1; rb+=1) //Red Black
        {
            for (int z = 0; z < local_sizes[Z]; z++) 
                for (int y = 0; y < local_sizes[Y]; y++)
                {   
                    if((x+y+z) %2== rb)
                    {   
                        ind=z*local_sizes[X]*local_sizes[Y]+y*local_sizes[X]+x;                 
                        result[index(ind,middle)] = 6* Xk[index(ind,middle)];
                        result[index(ind,middle)] -= Xk[index(ind+1,middle)];
                        if (y != 0)
                            result[index(ind,middle)] -= Xk[index(ind-local_sizes[X],middle)];
                        if (y != local_sizes[Y]-1)
                            result[index(ind,middle)] -= Xk[index(ind+local_sizes[X],middle)];
                        if (z != 0)
                            result[index(ind,middle)] -= Xk[index(ind-local_sizes[X]*local_sizes[Y],middle)];
                        if (z != local_sizes[Z]-1)
                            result[index(ind,middle)] -= Xk[index(ind+local_sizes[X]*local_sizes[Y],middle)];
                        if (comm_info[LEFT][SIZE])
                        {
                            result[index(ind,middle)] -=comm_buffer[comm_info[LEFT][DISP]+y+z*local_sizes[Y]];
                          //  comm_buffer[comm_info[LEFT][DISP]+y+z*local_sizes[Y]] = result[index(ind,middle)] ;
                        }    
                }
            }
            
        }
    }

    

    //RIGHT 
        {
            int x=local_sizes[X]-1;
            int ind=0;
            for ( int rb = 0; rb <=1; rb+=1) //Red Black
            {
                for (int z = 0; z < local_sizes[Z]; z++) 
                    for (int y = 0; y < local_sizes[Y]; y++)
                    {   
                        if((x+y+z) %2== rb)
                        {   
                            ind=z*local_sizes[X]*local_sizes[Y]+y*local_sizes[X]+x;                
                            result[index(ind,middle)] = 6* Xk[index(ind,middle)];
                            result[index(ind,middle)] -= Xk[index(ind-1,middle)];
                            if (y != 0)
                                result[index(ind,middle)] -= Xk[index(ind-local_sizes[X],middle)];
                            if (y != local_sizes[Y]-1)
                                result[index(ind,middle)] -= Xk[index(ind+local_sizes[X],middle)];
                            if (z != 0)
                                result[index(ind,middle)] -= Xk[index(ind-local_sizes[X]*local_sizes[Y],middle)];
                            if (z != local_sizes[Z]-1)
                                result[index(ind,middle)] -= Xk[index(ind+local_sizes[X]*local_sizes[Y],middle)];
                            if (comm_info[RIGHT][SIZE])
                            {
                                result[index(ind,middle)] -=comm_buffer[comm_info[RIGHT][DISP]+y+z*local_sizes[Y]];
                             //   comm_buffer[comm_info[RIGHT][DISP]+y+z*local_sizes[Y]] = result[index(ind,middle)] ;
                        }    
                    }
                }
                
            }
        }

    
    //UP 
    {
        int y=local_sizes[Y]-1;
        int ind=0;
        for ( int rb = 0; rb <=1; rb+=1) //Red Black
        {
            for (int z = 0; z < local_sizes[Z]; z++) 
                for (int x = 0; x < local_sizes[X]; x++)
                {   
                    if((x+y+z) %2== rb)
                    {   
                        ind=z*local_sizes[X]*local_sizes[Y]+y*local_sizes[X]+x;                
                        result[index(ind,middle)] = 6* Xk[index(ind,middle)];
                        if (x !=0)
                            result[index(ind,middle)] -= Xk[index(ind-1,middle)];
                        if (x != local_sizes[X]-1)
                            result[index(ind,middle)] -= Xk[index(ind+1,middle)];
                        result[index(ind,middle)] -= Xk[index(ind-local_sizes[X],middle)];
                        if (z != 0)
                            result[index(ind,middle)] -= Xk[index(ind-local_sizes[X]*local_sizes[Y],middle)];
                        if (z != local_sizes[Z]-1)
                            result[index(ind,middle)] -= Xk[index(ind+local_sizes[X]*local_sizes[Y],middle)];
                        if (comm_info[UP][SIZE])
                        {
                            result[index(ind,middle)] -=comm_buffer[comm_info[UP][DISP]+x+z*local_sizes[X]];
                         //   comm_buffer[comm_info[UP][DISP]+x+z*local_sizes[X]] = result[index(ind,middle)] ;
                        }    
                    }
                }
                
            }
        }

    

    //DOWN
    {
        int y=0;
        int ind=0;
        for ( int rb = 0; rb <=1; rb+=1) //Red Black
        {
            for (int z = 0; z < local_sizes[Z]; z++) 
                for (int x = 0; x < local_sizes[X]; x++)
                {   
                    if((x+y+z) %2== rb)
                    {   
                        ind=z*local_sizes[X]*local_sizes[Y]+y*local_sizes[X]+x;                
                        result[index(ind,middle)] = 6* Xk[index(ind,middle)];
                        if (x !=0)
                            result[index(ind,middle)] -= Xk[index(ind-1,middle)];
                        if (x != local_sizes[X]-1)
                            result[index(ind,middle)] -= Xk[index(ind+1,middle)];
                        result[index(ind,middle)] -= Xk[index(ind+local_sizes[X],middle)];
                        if (z != 0)
                            result[index(ind,middle)] -= Xk[index(ind-local_sizes[X]*local_sizes[Y],middle)];
                        if (z != local_sizes[Z]-1)
                            result[index(ind,middle)] -= Xk[index(ind+local_sizes[X]*local_sizes[Y],middle)];
                        if (comm_info[DOWN][SIZE])
                        {
                            result[index(ind,middle)] -=comm_buffer[comm_info[DOWN][DISP]+x+z*local_sizes[X]];
                         //   comm_buffer[comm_info[DOWN][DISP]+x+z*local_sizes[X]] = result[index(ind,middle)] ;
                        }    
                    }
                }
                
            }
        }

    

    //FRONT
    {
        int z=0;
        int ind=0;
        for ( int rb = 0; rb <=1; rb+=1) //Red Black
        {
            for (int y = 0; y < local_sizes[Y]; y++) 
                for (int x = 0; x < local_sizes[X]; x++)
                {   
                    if((x+y+z) %2== rb)
                    {   
                        ind=z*local_sizes[X]*local_sizes[Y]+y*local_sizes[X]+x;                
                        result[index(ind,middle)] = 6* Xk[index(ind,middle)];
                        if (x !=0)
                            result[index(ind,middle)] -= Xk[index(ind-1,middle)];
                        if (x != local_sizes[X]-1)
                            result[index(ind,middle)] -= Xk[index(ind+1,middle)];
                        if (y != 0)
                            result[index(ind,middle)] -= Xk[index(ind-local_sizes[X],middle)];
                        if(y != local_sizes[Y]-1)
                            result[index(ind,middle)] -= Xk[index(ind+local_sizes[X],middle)];

                        result[index(ind,middle)] -= Xk[index(ind+local_sizes[X]*local_sizes[Y],middle)];
                        if (comm_info[FRONT][SIZE])
                        {
                            result[index(ind,middle)] -=comm_buffer[comm_info[FRONT][DISP]+x+y*local_sizes[X]];
                          //  comm_buffer[comm_info[FRONT][DISP]+x+y*local_sizes[X]] = result[index(ind,middle)] ;
                        }    
                    }
                }
                
            }
        }

    

   
        //BACK
    {
        int z=local_sizes[Z]-1;
        int ind=0;
        for ( int rb = 0; rb <=1; rb+=1) //Red Black
        {
            for (int y = 0; y < local_sizes[Y]; y++) 
                for (int x = 0; x < local_sizes[X]; x++)
                {   
                    if((x+y+z) %2== rb)
                    {   
                        ind=z*local_sizes[X]*local_sizes[Y]+y*local_sizes[X]+x;                
                        result[index(ind,middle)] = 6* Xk[index(ind,middle)];
                        if (x !=0)
                            result[index(ind,middle)] -= Xk[index(ind-1,middle)];
                        if (x != local_sizes[X]-1)
                            result[index(ind,middle)] -= Xk[index(ind+1,middle)];
                        if (y != 0)
                            result[index(ind,middle)] -= Xk[index(ind-local_sizes[X],middle)];
                        if(y != local_sizes[Y]-1)
                            result[index(ind,middle)] -= Xk[index(ind+local_sizes[X],middle)];

                        result[index(ind,middle)] -= Xk[index(ind-local_sizes[X]*local_sizes[Y],middle)];
                        if (comm_info[BACK][SIZE])
                        {
                            result[index(ind,middle)] -=comm_buffer[comm_info[BACK][DISP]+x+y*local_sizes[X]];
                           // comm_buffer[comm_info[BACK][DISP]+x+y*local_sizes[X]] = result[index(ind,middle)] ;
                        }    
                    }
                }
                
            }
        }

    




    

}



template<typename A>
void UpdateCommBuffer(const A *Xk,A* comm_buffer,const int comm_info[6][2] ,int local_sizes[3]){
    int size = local_sizes[0]*local_sizes[1]*local_sizes[2];
    int middle= size-std::floor(size/2)+1+ (size%2==1 ? +1 : 0);  


    if (comm_info[LEFT][SIZE])   //LEFT 
    {
        int x=0;
        int ind=0;

        for (int z = 0; z < local_sizes[Z]; z++){
            for (int y = 0; y < local_sizes[Y]; y++)
            {                      
                comm_buffer[comm_info[LEFT][DISP]+y+z*local_sizes[Y]] = Xk[index(ind,middle)] ;
            }
        }            
        
    }


    if (comm_info[RIGHT][SIZE])   //RIGHT 
    {
        int x=local_sizes[X]-1;
        int ind=0;

        for (int z = 0; z < local_sizes[Z]; z++){ 
            for (int y = 0; y < local_sizes[Y]; y++)
            {                      
                comm_buffer[comm_info[RIGHT][DISP]+y+z*local_sizes[Y]] = Xk[index(ind,middle)] ;
            }
        } 
    }           
        
    


    //UP 
    if (comm_info[UP][SIZE]) 
    {
        int y=local_sizes[Y]-1;
        int ind=0;

        for (int z = 0; z < local_sizes[Z]; z++){
            for (int x = 0; x < local_sizes[X]; x++)
            {   
                ind=z*local_sizes[X]*local_sizes[Y]+y*local_sizes[X]+x;
                comm_buffer[comm_info[UP][DISP]+x+z*local_sizes[X]] = Xk[index(ind,middle)];                        
            }
        } 

    }


    //DOWN 
    if (comm_info[DOWN][SIZE]) 
    {
        int y=0;
        int ind=0;

        for (int z = 0; z < local_sizes[Z]; z++){
            for (int x = 0; x < local_sizes[X]; x++)
            {   
                ind=z*local_sizes[X]*local_sizes[Y]+y*local_sizes[X]+x;
                comm_buffer[comm_info[DOWN][DISP]+x+z*local_sizes[X]] = Xk[index(ind,middle)];                        
            }
        } 

    }

        //FRONT
    if (comm_info[FRONT][SIZE]) 
    {
        int z=0;
        int ind=0;
        for (int y = 0; y < local_sizes[Y]; y++){ 
            for (int x = 0; x < local_sizes[X]; x++)
            {   
             
                ind=z*local_sizes[X]*local_sizes[Y]+y*local_sizes[X]+x;
                comm_buffer[comm_info[FRONT][DISP]+x+y*local_sizes[X]] = Xk[index(ind,middle)] ;                
    
            }
            
                
        }
    }


    //BACK
    if (comm_info[BACK][SIZE]) 
    {
        int z=local_sizes[Z]-1;
        int ind=0;
        for (int y = 0; y < local_sizes[Y]; y++){ 
            for (int x = 0; x < local_sizes[X]; x++)
            {   
                ind=z*local_sizes[X]*local_sizes[Y]+y*local_sizes[X]+x;
                comm_buffer[comm_info[BACK][DISP]+x+y*local_sizes[X]] = Xk[index(ind,middle)] ;                
            }
   
        }
    }

  
              

}




