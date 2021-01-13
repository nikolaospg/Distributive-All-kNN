/** 
 * This file contains many useful functions that we used to implement the programs.
 * 
 * Contains:
 *  i) General purpose functions:
 *      3) void pointer_swap(double** a, double** b)
 *      4) double point_distance( double* X, double* Y, int size)
 * 
 *  ii) Functions that calculate the array D (given by the exercise):
 *      5) double* Array_reduce(double* A, int rows, int cols)
 *      6) double* create_D(double* X, double* Y, int m, int n, int d)
 * 
 *  iii) Quick select algorithms:
 *      7) void vector_element_swap(void *a, int index1, int index2, int type_flag)
 *      8) int vector_partition(double* A, int* indices, int size)
 *      9) double quick_select(double *A, int* indices, int size, int k)
 *      10) void quick_sort(double *A, int* indices, int size)
 *      11) void merge_structs(knnresult* A, knnresult B)
 *      12) int vector_partitionV2(double* A, int* indices_global, int* indices_local, int size)
 *      13) double quick_selectV2(double *A, int* indices_global, int* indices_local, int size, int k) 
 */

#include <stdio.h>
#include <stdlib.h>
#include <cblas.h>
#include <time.h>
#include <string.h>


//Struct given by the exercise
typedef struct knnresult{
    int    * nidx;    //!< Indices (0-based) of nearest neighbors [m-by-k]
    double * ndist;   //!< Distance of nearest neighbors          [m-by-k]
    int      m;       //!< Number of query points                 [scalar]
    int      k;       //!< Number of nearest neighbors            [scalar]
} knnresult;



/*****************************************************************************/
/*                        General purpose functions                          */
/*****************************************************************************/

/** 
 * This simple function changes the values of 2 pointers.
 * The parameteres are pointers to pointers, because we wanted to pass the pointer variables by reference
 */
void pointer_swap(double** a, double** b){
    double* temp=*a;
    *a=*b;
    *b=temp;
}


/*Simple function, which helps find the square of a distance between two vectors, X and Y, with a certain size.*/
double point_distance( double* X, double* Y, int size){         
   
	double count = 0;
	double diff;
	for(int i=0; i<size; i++){
		diff = X[i] - Y[i];
		count += diff*diff;
	}
	return count;
}



/*****************************************************************************/
/*                  Functions that calculate the array D                     */
/*****************************************************************************/

/** 
 * The function below calculates the sum(X.^2,2) in terms of the matlab expression, and returns a vector containing the values of the Matrix to Vector reduction.
 * 
 * Inputs: double* A --> The array that I want to make the calculation with. (In our case X or Y)
 *         int rows  --> The number of rows of the A array (number of queries m in case of Y, number of corpus points n in case of X)
 *         int cols  --> The number of columns of the A array (the dimension of the metric space)
 * 
 * Output: The mx1 array Y_reduce 
 * 
 *      **WARNING : The returned vector represents a column in the case of Y, or a mx1 array, and a row in the case of X, or an 1xd Array !** .-
 */
double* Array_reduce(double* A, int rows, int cols){

    double* ret_column = calloc(rows,sizeof(double));

    if(ret_column==NULL){
        printf("Error: Couldn't allocate memory.\n");
        exit(-1);
    }

    int row_index;      //Index helping me later to access the appropriate elements

    //Calculating the sum(A(i,j)^2 for each j element of the row and for every row i, and filling the ret_column vector 
    for(int i=0; i<rows; i++){
        row_index=cols*i;

        for(int j=0; j<cols; j++){
            ret_column[i]+= A[row_index+j]*A[row_index+j];     
        }
    }

    return ret_column;
}


/** 
 * This function creates the mxn distance matrix.
 * 
 * Inputs: double* X --> The X array (with the total number of corpus points)
 *         double* Y --> The Y array (the one with the queries)
 *         int m --> the number of queries
 *         int n --> the number of elements of X
 *         int d --> the dimension of the metric space
 * 
 * Outputs: The D array (returned with a m*nx1 vector)
 * 
 * The way it is done is:
 * 
 * First, we compute the X-reduced and Y_reduced vectors, which are the third and first terms of the MATLAB expression we had been given.
 * Then, with a double loop, we compute the "sum" of the first and the third terms of the expression. It is named as D (the initial values of D)
 * Then, we find the actual distance array D, by using the cblas_dgemm function. With the proper arguments, we calculate the expression:
 * C := alpha*op(A)*op(B) + beta*C where  alpha<-(-2), op(A)<-Y, op(B)<- X.', beta=1 and C<-D.
 * This is exactly the MATLAB expression implemented with the cblas function.
 * This way we get the form of D (mxn). Later, with another function, we get the KNN for every row .- 
 */
double* create_D(double* X, double* Y, int m, int n, int d){

    double* D = calloc(m*n, sizeof(double));
    if(D==NULL){
        printf("Error: Couldn't allocate memory.\n");
        exit(-1);
    }
    double* Y_reduced = Array_reduce(Y, m, d);  //Calculating the reductions, which are the first and third term of the expression
    double* X_reduced = Array_reduce(X, n, d);
    int index;                                  //Index helping me to change the appropriate elements later
  
    //Calculating the "sum" of the X_reduced,Y_reduce.
    //We are careful because the Y_reduce vector represents a column, while the X_reduced vector represents a row 
    for(int i=0; i<m; i++){
        index=n*i;
        for(int j=0; j<n; j++){
            D[index +j] =X_reduced[j] + Y_reduced[i];     
        }
    }

    free(Y_reduced);
    free(X_reduced);

    //Calculating the parameters and calling the cblas_dgemm to actually calculate the D matrix.
    int lda=d;
    int ldb=d;
    int ldc=n;

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, m, n ,d ,-2, Y, lda, X, ldb, 1, D, ldc);

    return D;
}



/*****************************************************************************/
/*                          Quick select algorithms                          */
/*****************************************************************************/

/** 
 * This function swaps two elements of a vector.
 * Inputs: void* a         --> pointer to the elements of the vector
 *         int index1      --> the index of the first element
 *         int index2      --> the index of the second element
 *         int type_flag   --> A flag so that we can specify the type of the elements.
 *     If type_flag==0 then the elements are integers, else the elements are doubles.  
 * 
 * The pointer is a void pointer so that we can actually work with both types of elements .-
 */
void vector_element_swap(void *a, int index1, int index2, int type_flag){

    if(type_flag == 0){
        int* b=a;
        int temp;              
        temp=b[index1];
        b[index1]=b[index2];
        b[index2]=temp;
    }

    else{
        double* b=a;
        double temp;               
        temp=b[index1];
        b[index1]=b[index2];
        b[index2]=temp;
    }

}



/** 
 * The two functions below serve as a method to apply the quickselect algorithm.
 * When we swap the elements of the array, we also swap the indices of the elements in the starting array.
 * To achieve this, we demand an array with the indices as a parameter for the functions.
 * This way we change both of these vectors simultaneously.
 * The first function is used to partition the array (with the first element as a pivot),
 * and the second one is a non recursive quickselect implementation .-
 */
 

/*  This function returns the index of the pivot in the sorted array. */
int vector_partition(double* A, int* indices, int size){

    //Pivot is the first value of the vector A
    double pivot=A[0];
    int i=0;
    int j=size;

    if(size == 1){
        return 0;
    }

    while(i<j){
    
        do{
            i++;
        }while(A[i]<=pivot && i<size-1);

        do{
            j--;
        }while(A[j]>=pivot && j>0);

        if(i<j){
            vector_element_swap(A, i, j, 1);
            vector_element_swap(indices, i, j, 0);
        }
    }
    
    //Place the pivot in the correct position in the sorted vector
    vector_element_swap(A, j, 0, 1);
    vector_element_swap(indices, j, 0, 0);
   
    return j;

}



/** 
 * This function applies the quick_select method for a specific k on the A array.
 * After this function is called the kth element is in its correct position and the elements on the left are smaller than it.
 * It returns the value of the kth element in array A.
 */
double quick_select(double *A, int* indices, int size, int k){

    int ret = vector_partition(A, indices, size);

    //While loop controling the outcome of the algorithm.
    while(ret != (k-1) && size>1 ){

        //Partitioning on the left of the pivot
        if(k-1<ret){
            size = ret;
            ret = vector_partition(A, indices, size);
        }

        //Partitioning on the right of the pivot (changing the parameters in an appropriate way so that the correct sub-vector is partitioned)
        else if(k-1>ret){
            k = k - ret -1;
            size = size - ret -1;
            indices = indices + ret +1;
            A = A + ret +1;
            ret = vector_partition(A, indices, size);
        }
    }

    return A[ret];

}


/**
 * Takes an array A with a specific size and sorts it.
 * The indices vector changes accordingly.
 */
void quick_sort(double *A, int* indices, int size){

    if(size>1){

        int ret = vector_partition(A, indices, size);

        double* left = A;
        int* left_indices = indices;
        quick_sort(left, left_indices, ret+1);

        double* right = A+ret+1;
        int* right_indices = indices +ret+1;
        quick_sort(right, right_indices, size-ret-1);
    }
    
}


/**
 * This is the function that allowed us to implement the V1. It helps us to keep finding the kNN to the blocks of X every process was assigned
 * everytime we get a new block from the other processes.
 * 
 * Inputs: knnresult* A --> The result struct we have, which will be used until the end when the process will have computed the distances from all the X points.
 *         knnresult B  --> The result struct from the knn computation between the local X block and the block we just acquired from another process.
 * 
 * One could say that this function "merges" the two structures and alters the first one (the arrays it points to, to be specific),
 * so that it is the one who holds the actual kNN up until that point.
 * This is the reason why the first one should be passed be reference, so that it can change.
 * By calling the function each time we get a new block, we manage to have the final results in the A struct .-  
 */
void merge_structs(knnresult* A, knnresult B){
    
    //Creating the variables to help us implement the loops beneath.
    int k = A->k;
    int m = A->m;
    int count = 0;
    int a_index = 0;
    int b_index = 0;
    double* ret_values = malloc(k*m*sizeof(double));    //These are the new arrays that the A struct will point to
    int* ret_indices = malloc(k*m*sizeof(int));         //and they are of equal size to the previous arrays (A->ndist and A->nidx)

    if(ret_values==NULL || ret_indices==NULL){
        printf("Error: Couldn't allocate memory.\n");
        exit(-1);
    }
    //Finished creating the variables


    //Condition that checks if there is a mismatch of rows in A and B
    //If B has less rows than A the last m-B.m rows of A must remain unchanged
    if (B.m < m){

        //Iterating every row of A that will not change
        for (int i=B.m; i<m; i++){

            //Iterating every element of the i row
            for(int j=0; j<k; j++){
                
                //Copying the last rows in the new arrays
                ret_values[k*i + j] = A->ndist[k*i + j];
                ret_indices[k*i+j] = A->nidx[k*i+ j];
            }
            
        }

        //Adjust m
        m = B.m;
    }
    //Else if B has more rows than A then no adjustment is needed


    //For every (remaining) row of A
    for(int i=0; i<m; i++){

        //We make a first test, for lower complexity reasons.
        //Here, we compare the last element of A distances and the first one of the B distances.
        //If the A dist element is smaller than the B dist, then all of the elements of the A dist array are
        //smaller than the B counterparts, so we just fill the row with A values.
        if(A->ndist[k*i+ k-1] < B.ndist[k*i]){
            for(int j=0; j<k; j++){
                ret_values[k*i+j] = A->ndist[k*i+ j];
                ret_indices[k*i+j] = A->nidx[k*i+ j];
            }
            continue;
        }
        //End of test

        //Filling the rows of the new arrays. With every comparison, we get the smaller element and put it to the ret_values arrays.
        //If they are the same, we put both of the values.
        //We always check if the count is smaller than k to fill elements into the ret_values arrays
        else{
            count = 0;
            a_index = 0;
            b_index = 0;
            while(count < k){

                //For the case the A element is larger
                if(A->ndist[a_index +k*i] > B.ndist[b_index +k*i]){          
                    ret_values[count +k*i] = B.ndist[b_index + k*i];
                    ret_indices[count+k*i] = B.nidx[b_index+ k*i];
                    b_index++;
                    count++;
                }

                //For the case the B element is largerd
                else if(A->ndist[a_index+ k*i] < B.ndist[b_index+ k*i]){
                    ret_values[count+ k*i] = A->ndist[a_index+ k*i];
                    ret_indices[count+ k*i] = A->nidx[a_index+ k*i];
                    a_index++;
                    count++;
                }

                //For the case the elements are equal
                else{
                    ret_values[count+ k*i] = B.ndist[b_index+ k*i];
                    ret_indices[count+ k*i] = B.nidx[b_index+ k*i];
                    count++;
                    b_index++;
                    if(count < k){
                        ret_values[count+ k*i] = A->ndist[a_index+ k*i];
                        ret_indices[count+ k*i] = A->nidx[a_index+ k*i];
                        count++;
                        a_index++;
                    }
                }
            }
        }
        //Finished filling the row.
       
    }
    
    //Freeing the old arrays and changing the A pointers, so that they point to the new arrays.
    free(B.ndist);
    free(B.nidx);
    free(A->nidx);
    free(A->ndist);
    A->ndist = ret_values;
    A->nidx = ret_indices;
    //End of changing. The A result struct is now refreshed
}



/**
 * The two functions below are small changes in the partition and quick_select functions that we used to implement V2.
 * Here, we swap the elements of another array as well, which is called indices_local.
 * The reason why we do so, is because we want to insert elements in the tree, without having to use searching functions to insert them in a proper way.
 * Βy passing the indices_global array, we can easily see which elements are on the left side of the pivot and which ones are on the right, helping us make
 * the space partitions for the tree creation .-
 */
int vector_partitionV2(double* A, int* indices_global, int* indices_local, int size){

    //Pivot is the first value of the vector A
    double pivot=A[0];
    int i=0;
    int j=size;

    if(size == 1){
        return 0;
    }

    while(i<j){
    
        do{
            i++;
        }while(A[i]<=pivot && i<size-1);

        do{
            j--;
        }while(A[j]>=pivot && j>0);

        if(i<j){
            vector_element_swap(A, i, j, 1);
            vector_element_swap(indices_global, i, j, 0);
            vector_element_swap(indices_local, i, j, 0);
        }
    }
    
    //Place the pivot in the correct position in the sorted vector
    vector_element_swap(A,j, 0, 1);
    vector_element_swap(indices_global, j, 0, 0);
    vector_element_swap(indices_local, j, 0, 0);
   
    return j;

}


double quick_selectV2(double *A, int* indices_global, int* indices_local, int size, int k){

    int ret = vector_partitionV2(A, indices_global, indices_local, size);

    //While loop controling the outcome of the algorithm.
    while(ret != (k-1) && size>1 ){

        //Partitioning on the left of the pivot
        if(k-1<ret){
            size = ret;
            ret = vector_partitionV2(A, indices_global, indices_local, size);
        }

        //Partitioning on the right of the pivot (changing the parameters in an appropriate way so that the correct sub-vector is partitioned)
        else if(k-1>ret){
            k = k - ret -1;
            size = size - ret -1;
            indices_global = indices_global + ret +1;
            indices_local = indices_local + ret +1;
            A = A + ret +1;
            ret = vector_partitionV2(A, indices_global, indices_local, size);
        }
    }


    return A[ret];

}