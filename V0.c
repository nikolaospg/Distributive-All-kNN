#include "exalgorithms.h"
#include "readfile.h"


/** 
 * The function that calculates the knnresult struct, given the:
 *  1)double* X -->     The corpus points array
 *  2)double* Y -->     The query points array
 *  3)int n     -->     The number of corpus points
 *  4)int m     -->     The number of query points
 *  5)ind d     -->     The number of dimensions of the space
 *  6)int k     -->     The number of nearest neighbours that I want to calculate
 * 
 * The way it is done is:
 *  First, we compute the initial D matrix, by calling the create_D function. We also create an m*n element vector, called indices, 
 *  that contains the indices of the respective elements in the D matrix.
 *  Then, by using the quick_select and quick_sort function, which are written in exalgorithms.h, we search for the kth nearest neighbour.
 *  The quick_select finds the kth smallest element in the array and puts it in its correct position.
 *  After that, since all the elements before k are smaller than the element there, these are our kNN.
 *  However, before storing the elements we sort them with quicksort.
 *  We create the struct to be returned and store each neighbour in the knnresult.ndist and the indices in knnresult.nidx. .-
 */
knnresult kNN(double * X, double * Y, int n, int m, int d, int k){

    //Creating the D using our function create_D
    double* D = create_D(X, Y, m ,n, d);

    //Initializing the indices vectors
    int* indices = malloc(m*n*sizeof(int));

    //Check if allocation failed
    if(indices==NULL){
        printf("Error: Couldn't allocate memory.");
        exit(-1);
    }

    for(int i=0; i<m; i++){
        for(int j=0; j<n; j++){
            indices[n*i +j] = j;
        }
    }


    //Creating the mxk arrays, so that they contain only the knn and the respective indices.
    //The D_pointer and indices_pointer are pointers we use to call our select function.
    double* D_pointer = D;
    int* indices_pointer = indices;

    //Initializing the struct to be returned
    knnresult result;
    result.k = k;
    result.m = m;
    result.ndist = malloc(m*k*sizeof(double));
    result.nidx = malloc(m*k*sizeof(double));

    if(result.ndist==NULL){
        printf("Error: Couldn't allocate memory.");
        exit(-1);
    }    
    if(result.nidx==NULL){
        printf("Error: Couldn't allocate memory.");
        exit(-1);
    }


    //Iterating for every row
    for(int i=0; i<m; i++){       

        //Find the kth neighbour after that the first k elements in D_pointer are the kNN
        quick_select(D_pointer, indices_pointer, n, k);

        //Sort the neighbours
        quick_sort(D_pointer, indices_pointer, k);

        //Store the neighbours in the struct to be returned
        for(int j=0; j<k; j++){
            result.ndist[j + k*i] = D_pointer[j];
            result.nidx[j + k*i] = indices_pointer[j];
        }

        //Moving the help pointers one row ahead to continue with the next iteration
        D_pointer = D_pointer+n;
        indices_pointer = indices_pointer+n;
    }


    //Freeing the old and larger Arrays and returning the struct.
    free(D);
    free(indices);
    
    return result;
}


/** 
 * The kNN_block function acts as an improvement to the previous kNN function. The difference is that
 * with this one we can actually block the Y array and not have to create large intermidiate arrays,
 * The two functions are similar, and the outputs/inputs are exactly the same.
 * For a better understanding of what is going on, we recommend you to study the comments of the kNN function .-
 */
knnresult kNN_block(double * X, double * Y, int n, int m, int d, int k){

    //The number of blocks to divide Y
    int Y_parts = 10;    

    //Initializing the struct to be returned
    knnresult result;
    result.ndist = malloc(m*k*sizeof(double));
    result.nidx = malloc(m*k*sizeof(int));
    result.k = k;
    result.m = m;

    if(result.ndist==NULL){
        printf("Error: Couldn't allocate memory.");
        exit(-1);
    }
    if(result.nidx==NULL){
        printf("Error: Couldn't allocate memory.");
        exit(-1);
    }

    //Helpful variable, that will later assist on filling the result.ndist, result.nidx vectors (with the kNN of each query).
    int final_vectors_indices = 0;

    //Pointer that assists on blocking the big Y vector into smaller Y_small ones.
    double* Y_small = Y;

    //In the following loop, we block the Y array and fill the result.ndist and result.nidx arrays with the appropriate values.
    for(int i=0; i<Y_parts; i++){


        //FIRST STEP    -> Creating the blocks. We assign m/Y_parts query points to each iteration, except of the last which takes all the remaining ones
        //This is done for the cases when m is not an integral multiple of Y_parts
        int Y_small_rows = m/Y_parts;

        //Adjusts the size of Y_small for the last iteration
        if(i==Y_parts-1){
            Y_small_rows = m-(Y_parts -1)*(m/Y_parts);
        }


        //SECOND STEP   -> Creating the D_small and indices_small arrays, by only working with the values of the Y_small array.
        double* D_small = create_D(X, Y_small, Y_small_rows, n, d);
        int* indices_small = malloc(Y_small_rows*n*sizeof(int));

        if(indices_small==NULL){
            printf("Error: Couldn't allocate memory.");
            exit(-1);
        }

        for(int j=0; j<Y_small_rows; j++){
            for(int q=0; q<n; q++){
              indices_small[n*j +q] = q;
            }
        }


        //THIRD STEP    -> Using quickselect and quicksort to find the kNN for every point of the Y block we are currently working with.


        //Pointers that will allow me to use quickselect on the correct parts of the D_small/indices_small arrays.
        double* D_pointer = D_small;                  
        int* indices_pointer = indices_small;
        
        
        for(int j=0; j<Y_small_rows; j++){

            //Find the kth neighbour after that the first k elements in D_pointer are the kNN
            quick_select(D_pointer, indices_pointer, n, k);

            //Sort the neighbours
            quick_sort(D_pointer, indices_pointer, k);

            //Store the neighbours in the struct to be returned
            for(int q=0; q<k; q++){
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



int main(int argc, char* argv[]){
  
    int n = 0;
    int d = 0;
    int k = 0;

    //If flag == 0 then the user gave a file as an input
    int flag = 0;

    //Used to read the input files
    int numbered = 0;

    double* X;


    if(argc < 3){
        printf("\nNot enough arguments!\nUsage: %s [Path to file] [No. kNN]\n", argv[0]);
        printf("\t*****OR*****");
        printf("\nUsage: %s [No. query points] [No. dimensions] [No. kNN]\n\n", argv[0]);
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
        printf("\nCannot find kNN because k exceeds the number of query points.\n\n");
        exit(-1);
    }


    if(!flag){
        printf("\nCurrently in %s running %s with n = %d, d = %d and k = %d.\n", argv[0], argv[1], n, d, k);
    }
    else{
        printf("\nCurrently in %s running a random array with n = %d, d = %d and k = %d.\n", argv[0], n, d, k);
    }


    //X is allocated only in the process with root==0           
    X = malloc(n*d*sizeof(double));

    if(X==NULL){
        printf("Error: Couldn't allocate memory for X.\n");
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


    struct timespec begin, end;

    // Starting the clock
    clock_gettime(CLOCK_MONOTONIC, &begin);

    knnresult b = kNN_block(X, X, n, n, d, k);

    // Stopping the clock
    clock_gettime(CLOCK_MONOTONIC, &end);
    long seconds = end.tv_sec - begin.tv_sec;
    long nanoseconds = end.tv_nsec - begin.tv_nsec;
    double elapsed = seconds + nanoseconds*1e-9;

    printf("\nTime elapsed: %.5f seconds.\n\n", elapsed);

    free(X);
    free(b.nidx);
    free(b.ndist);

    if(!flag){
        printf("Done with %s that run %s with n = %d, d = %d and k = %d.\n\n", argv[0], argv[1], n, d, k);
    }
    else{
        printf("Done with %s that run a random array with n = %d, d = %d and k = %d.\n\n", argv[0], n, d, k);
    }

    
    return 0;
}