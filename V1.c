#include "exalgorithms.h"
#include "mpi.h"
#include <float.h>
#include <math.h>
#include "readfile.h"



/** 
 * The function that calculates the knnresult struct, given the:
 *  1)double* X        -->    The corpus points array
 *  2)double* Y        -->    The query points array
 *  3)int rank         -->    The rank of the process that called the function
 *  4)ind d            -->    The number of dimensions of the space
 *  5)int k            -->    The number of nearest neighbours that I want to calculate
 *  6)int comm_size    -->    The number of processes called
 *  7)int source_rank  -->    The rank of the process that sent the X array
 *  8)int corpus_total -->    The size of the starting array that we divided into blocks 
 * 
 * This function is similar to the kNN_block function that was implemented in V0.c
 * The main difference is that the rows n and m, stored in the arrays X and Y respectively,
 * are calculated inside the function depending on the size of the communicator,
 * the rank of the process that called the function and the rank of the process that sent the Y array .-
 */
knnresult kNNV1(double * X, double * Y, int rank, int d, int k, int comm_size, int source_rank, int corpus_total){

    //Calculating the rows of the subprocesses (i.e. the number of points stored in them)
    int subprocess_rows = corpus_total/comm_size;

    //Calculating the rows of the root process
    int root_process_rows = corpus_total - (comm_size-1)*subprocess_rows;
    
    //The number of rows in Y
    int m;
    if(rank==0){
        m = root_process_rows;
    }
    else{
        m = subprocess_rows;
    }
    
    //The number of rows in X
    int n;
    if(source_rank==0){
        n = root_process_rows;
    }
    else{
        n = subprocess_rows;
    }


    int Y_parts = 10;


    //Initializing the struct to be returned
    knnresult result;
    result.ndist = malloc(m*k*sizeof(double));
    result.nidx = malloc(m*k*sizeof(int));
    result.k = k;
    result.m = m;


    if(result.ndist==NULL){
        printf("Error: Couldn't allocate memory.\n");
        exit(-1);
    }
    if(result.nidx==NULL){
        printf("Error: Couldn't allocate memory.\n");
        exit(-1);
    }

    //Helpful variable, that will later assist on filling the result.ndist, result.nidx vectors (with the kNN of each query).
    int final_vectors_indices = 0;

    //Pointer that assists on blocking the big Y vector into smaller Y_small ones.
    double* Y_small = Y;

    //In the following loop, we block the Y array and fill the result.ndist and result.nidx arrays with the appropriate values.
    for(int i=0; i<Y_parts; i++){



        //FIRST STEP    -> Creating the blocks. We assign m/Y_parts query points to each iteration, except of the last which takes all the remaining ones
        //This is done for the cases when the m is not an integral multiple of the Y_parts
        int Y_small_rows = m/Y_parts;

        //Runs in the last iteration
        if(i==Y_parts -1){
            Y_small_rows = m - (Y_parts -1)*(m/Y_parts);
        }


        //SECOND STEP   -> Creating the D_small and indices_small arrays, by only working with the values of the Y_small array.
        double* D_small = create_D(X, Y_small, Y_small_rows , n, d);
        int* indices_small = malloc(Y_small_rows*n*sizeof(int));

        if(indices_small==NULL){
            printf("Error: Couldn't allocate memory.\n");
            exit(-1);
        }

        //The index changes due to the fact that the X is blocked
        int index_change = 0;

        if(source_rank!=0){
            index_change = root_process_rows + subprocess_rows*(source_rank-1);
        }
 

        for(int j=0; j<Y_small_rows; j++){
            for(int q=0; q<n; q++){
                indices_small[n*j +q] = q + index_change;
            }
        }


        //THIRD STEP    -> Using quickselect and quicksort to find the kNN for every point of the Y block we are currently working with.

        //Pointers that will allow me to use quick select on the correct parts of the D_small/indices_small arrays.
        double* D_pointer = D_small;
        int* indices_pointer = indices_small;


        for(int j=0; j<Y_small_rows; j++){

            //Find the kth neighbour after that the first k elements in D_pointer are the kNN
            quick_select(D_pointer, indices_pointer, n, k);

            //Sort the neighbours
            quick_sort(D_pointer, indices_pointer, k);

            //Store the neighbours in the struct to be returned
            for (int q=0; q<k; q++){

                result.ndist[final_vectors_indices] = D_pointer[q];
                result.nidx[final_vectors_indices] = indices_pointer[q];
                final_vectors_indices++;
        
            }
            //Moving the pointers along so that they point to the next row.
            D_pointer = D_pointer+n;
            indices_pointer = indices_pointer+n;
        }

        //Freeing the memory that was allocated in create_D function
        free(D_small);
        free(indices_small);

        //Changing the Y_small pointer so that it points in the next block of the Y matrix.
        Y_small = Y_small + Y_small_rows* d;

    }

    return result;

}

/**
 * Function that finds the rank of the process that initially sent array Y and then calls the kNNV1 function.
 * We need to find the rank of the source in order to specify the rows in Y (calculated in kNNV1) .-
 */
knnresult compute_kNN(double * X_local, double * Y, int rank, int d, int k, int comm_size, int corpus_total, int iteration_num){

    int source_rank;
    if(rank>=iteration_num){
        source_rank = rank-iteration_num;
    }
    else{
        source_rank = comm_size -iteration_num+rank;
    }

    knnresult return_struct = kNNV1(Y, X_local, rank, d, k, comm_size, source_rank, corpus_total);
    
    return return_struct;
}


/** 
 * Function that splits the corpus array into p blocks, where p is the number of processes called.
 * 
 * Each MPI process has:
 *  1)an X array (X_local) that is a block of the starting corpus array
 *  2)a Y array from which it calculates the distances from X and keeps the k smallest
 *  3)and a Z array that stores the points received from another process
 * 
 * The output of the function is a knnresult struct that contains the information of the kNN
 * for the points in X_local
 */
knnresult distrAllkNN(double* X, int n, int d, int k){

    /*Initialising all the variables that the MPI processes need.*/
    int comm_size,rank;
    int subprocess_rows;                //The amount of rows the subprocesses will have
    int root_process_rows;              //The amount of rows the root process will have
    int current_rows;                   //The amount of rows the current process has

    double* X_local;                    //Array with the local points X(corpus)
    double* Y;                          //Array with the points Y(query) that X_local is being compared to
    double* Z;                          //Array with the incoming points
    
    knnresult result_local;             //The struct containing the results of the process (It changes every time we get a new Z)
    knnresult result_temp;              //The struct containing the result of the kNN search between the X_local and the Y array.

    MPI_Request send_request;
    MPI_Request recv_request;
    MPI_Status send_status;
    MPI_Status recv_status;
    /*Finished initialising*/


    /*Getting info about the processes*/
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    subprocess_rows = n/comm_size;                            //We assign a standard n/comm_size rows to each subprocess
    root_process_rows = n - (comm_size-1)*(n/comm_size);      //We assign the rest to the root process.    
    if(rank==0){
        current_rows = root_process_rows;
    }
    else{
        current_rows = subprocess_rows;
    }
    /*Finished getting info and allocating the memory*/


    /*Using Scatterv to send the X blocks to each MPI process, freeing the memory we don't need and allocating the Y, Z*/
    //Each process has a unique X_local array, which contains a block from the original X array
    X_local = malloc(current_rows*d*sizeof(double));

    if(X_local==NULL){
        printf("Error: Couldn't allocate memory.\n");
        exit(-1);
    }

    //Variables used in MPI_Scatterv
    int* send_counts = malloc(comm_size*sizeof(int));
    int* displacements = malloc(comm_size*sizeof(int));

    if(send_counts==NULL){
        printf("Error: Couldn't allocate memory.\n");
        exit(-1);
    }
    if(displacements==NULL){
        printf("Error: Couldn't allocate memory.\n");
        exit(-1);
    }

    send_counts[0] = root_process_rows*d;
    displacements[0]=0;
    for(int i=1; i<comm_size; i++){
        send_counts[i] = subprocess_rows*d;
        displacements[i] = (root_process_rows +(i-1)*subprocess_rows)*d;
    }

    MPI_Scatterv(X, send_counts, displacements, MPI_DOUBLE, X_local, root_process_rows*d, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    free(send_counts);
    free(displacements);


    Y = malloc(root_process_rows*d*sizeof(double));

    if(Y==NULL){
        printf("Error: Couldn't allocate memory.\n");
        exit(-1);
    }

    //In the beginning the Y array is the same as X_local
    for(int i=0; i<root_process_rows*d; i++){
        Y[i] = X_local[i];
    }

    if(rank==0){
        free(X);
    }

    Z = malloc(root_process_rows*d*sizeof(double));

    if(Z==NULL){
        printf("Error: Couldn't allocate memory.\n");
        exit(-1);
    }
    /*Finished sending, freeing and allocating*/


    /** Below we make the message passing comm_size -1 times, and we use the merge_structs so that the result_local has the kNN results,
    *   which are refreshed every time we get a a new block **/
    for(int i=0; i<comm_size -1; i++){
        
        //The information moves in a circle, meaning a process sends points to the one with the next rank and receives points from the previous
        //To complete the circle the last process (with rank==comm_size) sends the points to the root process (with rank==0)
        if(rank==comm_size -1){
            MPI_Isend(Y, root_process_rows*d, MPI_DOUBLE, 0, 4, MPI_COMM_WORLD, &send_request);
            MPI_Irecv(Z, root_process_rows*d, MPI_DOUBLE, comm_size-2, 4, MPI_COMM_WORLD, &recv_request);
        }
        else if(rank==0){
            MPI_Isend(Y, root_process_rows*d, MPI_DOUBLE, 1, 4, MPI_COMM_WORLD, &send_request);
            MPI_Irecv(Z, root_process_rows*d, MPI_DOUBLE, comm_size-1, 4, MPI_COMM_WORLD, &recv_request);
        }

        else{
            MPI_Isend(Y, root_process_rows*d, MPI_DOUBLE, rank+1, 4, MPI_COMM_WORLD, &send_request);
            MPI_Irecv(Z, root_process_rows*d, MPI_DOUBLE, rank-1, 4, MPI_COMM_WORLD, &recv_request);
        }   


        if(i==0){
            result_local = compute_kNN(X_local, Y, rank, d, k, comm_size, n, i);
        }
        else{
            result_temp = compute_kNN(X_local, Y, rank, d, k, comm_size, n, i);
            merge_structs(&result_local, result_temp);
        }

        MPI_Wait(&send_request, &send_status);
        MPI_Wait(&recv_request, &recv_status); 
        
        pointer_swap(&Y, &Z);
    }

    result_temp = compute_kNN(X_local, Y, rank, d, k, comm_size, n, comm_size -1);
    merge_structs(&result_local, result_temp);
    /*Finished the message passing and the kNN computations.*/

    free(X_local);
    free(Y);
    free(Z);

    return result_local;
}



int main(int argc, char* argv[]){

    //Initialising the environment
    MPI_Init(&argc, &argv);
    int comm_size;
    int rank;
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    //End of initialization

    int n = 0;
    int d = 0;
    int k = 0;

    //If flag == 0 then the user gave a file as an input
    int flag = 0;

    //Used to read the input files
    int numbered = 0;

    double* X;


    if(argc < 3){
        if(rank==0){
            printf("\nNot enough arguments!\nUsage: mpirun -np [No. processes] %s [Path to file] [No. kNN]\n", argv[0]);
            printf("\t*****OR*****");
            printf("\nUsage: mpirun -np [No. processes] %s [No. query points] [No. dimensions] [No. kNN]\n\n", argv[0]);
        }
        exit(-1);
    }


    if (argc==3){
        //Call by reference in order to change their values in the function
        find_sizes(argv[1], &n, &d, &numbered);

        k = atoi(argv[2]);
    }
    else if (argc==4){
        n = atoi(argv[1]);
        d = atoi(argv[2]); 
        k = atoi(argv[3]);

        //The user wants to use a random array
        flag = 1;
    }
    else{
        printf("Too many arguments\n");
        exit(-1);
    }
    


    //Check if n and k are valid
    if(n<k){
        if(rank==0){
            printf("\nCannot find kNN because k exceeds the number of query points.\n\n");
        }
        exit(-1);
    }

    if(rank==0){
        if(!flag){
            printf("\nCurrently in %s running %s with n = %d, d = %d, k = %d and p = %d.\n", argv[0], argv[1], n, d, k, comm_size);
        }
        else{
            printf("\nCurrently in %s running a random array with n = %d, d = %d, k = %d and p = %d.\n", argv[0], n, d, k, comm_size);
        }
    }


    //X is allocated only in the process with root==0
    if(rank==0){            
        X = malloc(n*d*sizeof(double));

        if(X==NULL){
            printf("Error: Couldn't allocate memory.\n");
            exit(-1);
        }


        if(flag==0){
            X = getX(argv[1], n, d, numbered); 
            printf("\nArray in file successfully read\n");  
        }
        else{
            srand(time(NULL));

            for(int i=0; i<n*d; i++){
                X[i] = (double)rand()/100000000;
            }
            printf("\nRandom array complete\n");
        }
    }

    struct timespec begin, end;

    if (rank==0){

        // Starting the clock
        clock_gettime(CLOCK_MONOTONIC, &begin);

    }
    
    knnresult a = distrAllkNN(X, n, d, k);

    if(rank==0){

        // Stopping the clock
        clock_gettime(CLOCK_MONOTONIC, &end);
        long seconds = end.tv_sec - begin.tv_sec;
        long nanoseconds = end.tv_nsec - begin.tv_nsec;
        double elapsed = seconds + nanoseconds*1e-9;

        printf("\nTime elapsed: %.5f seconds.\n\n", elapsed);
    }


    free(a.ndist);
    free(a.nidx);

    if(rank==0){
        if(!flag){
            printf("Done with %s that run %s with n = %d, d = %d, k = %d and p = %d.\n\n", argv[0], argv[1], n, d, k, comm_size);
        }
        else{
            printf("Done with %s that run a random array with n = %d, d = %d, k = %d and p = %d.\n\n", argv[0], n, d, k, comm_size);
        }
    }

    MPI_Finalize();
    return 0;

}