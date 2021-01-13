#include "exalgorithms.h"
#include "mpi.h"
#include <float.h>
#include <math.h>
#include "readfile.h"

//This global variable is used to easily change the number of elements in each bucket when calling the tree_create recursively. Check the tree_create function.
int bucket_index = 0;


/**
 * The fuction the tree_fun calls to create the tree. Take a look at tree_fun first.
 * It is called recursively.
 * 
 * Inputs:
 *  1)double* X_sub             --> The metric space to make the tree of
 *  2)int* indices_global       --> Indices of the points (that we take from X)
 *  3)int x_rows                --> The amount of rows the X array has
 *  4)int d                     --> Number of dimensions      
 *  5)int root_key              --> The key-number that the root of the subtree has
 *  6)double* bucket_pointer    --> Pointer to the tree array (used to "fill" the buckets)
 *  7)int level_count           --> The level in the tree that the root node has
 *  8)int final_levels          --> Total levels of the tree (taking into account that the buckets reduce it)
 */
void tree_create(double* X_sub, int* indices_global, int x_rows, int d, int root_key, double* bucket_pointer, int level_count, int final_levels){
   
    //If clause which is useful to stop the recursive calls.
    if(level_count<final_levels){

        //Calculating useful variables for the tree creation.
        int vp_elements_num = pow(2, final_levels-1) -1;
        int bucket_number = vp_elements_num +1;

        double* vp_elements_pointer = bucket_pointer + bucket_number;                          //These pointers are used to access the various elements/parts of the array.
        double* leaf_elements_pointer = vp_elements_pointer + vp_elements_num*(3+d);
        //Finished calculating the variables.
        

        //Creating the root point of this sub tree into the array
        double* root_pointer = vp_elements_pointer + (3+d)*root_key;

        root_pointer[0] = indices_global[0];        //As we said, the first number is the index of this point
        root_pointer[1] = root_key;                 //The second one is the key
    
        //Choose the first point of X_sub as a vantage point
        for(int i=0; i<d; i++){
            root_pointer[3+i] = X_sub[i];           //Copying the coordinates of the point into the tree
        }
        //Finished creating for the moment. I will later have to put the median distance in root_pointer[2] as well.
  

        //Creating the D array and the ones to be used in the next call
        double* D = create_D(X_sub+d, X_sub, 1, x_rows-1, d);

        int* indices_local = malloc((x_rows-1)* sizeof(int));         //These indices will help me create the X_left and X_right arrays without the need to use searching methods.
        if(indices_local==NULL){
            printf("Error: Couldn't allocate memory.\n");
            exit(-1);
        }
        for(int i=0; i<x_rows-1; i++){
            indices_local[i] = i+1;
        }

        //Finding and storing the median
        root_pointer[2] = quick_selectV2(D, indices_global+1, indices_local, x_rows-1, (x_rows-1)/2+1);

        //Calculating the X_left, X_right and every other useful array, which will be used in the recursive call.
        int left_rows = (x_rows-1)/2;
	    int right_rows = x_rows-1 - left_rows;

	    double* X_left = malloc(left_rows*d*sizeof(double));
	    double* X_right = malloc(right_rows*d*sizeof(double));

        int * indices_left = malloc(left_rows*sizeof(int));
        int * indices_right = malloc(right_rows*sizeof(int));

        if(X_left==NULL || X_right==NULL || indices_left==NULL || indices_right==NULL){
            printf("Error: Couldn't allocate memory.\n");
            exit(-1);
        }

        int index_change;

        for(int i=0; i<left_rows; i++){
            indices_left[i] = indices_global[i+1];
            index_change = indices_local[i] *d;

            for(int j=0; j<d; j++){
                X_left[j +i*d] = X_sub[index_change+j];
            }
        }

        for(int i=0; i<right_rows; i++){
            indices_right[i] = indices_global[i+ left_rows +1];
		    index_change = indices_local[i + left_rows]*d;

		    for(int j=0; j<d; j++){
			    X_right[i*d +j] = X_sub[index_change+j ];
		    }
	    }
        //Finished creating the D array and every other parameter to call the function again.

        
        
        //Working in the case of the points that spawn buckets, with B>=1 points.
        if(level_count==final_levels-1){
            
            //Filling the two buckets, that are "children" of this node which is the root of the subtree
            //The use of the global variable could have been avoided, but serves as an easy and sharp equivalent to passing more function arguments.
            if(bucket_index==0){
                bucket_pointer[0] = left_rows;            //Filling the elements in the first part of the array, which give info about the points in each bucket.
                bucket_pointer[1] = bucket_pointer[0] +right_rows;
                
            }
            //The special case of index==0 was taken, because we do not have a special zero value element in the beginning of the array (which may be generally used in similar situations)
            else{
                bucket_pointer[bucket_index] = bucket_pointer[bucket_index-1] +left_rows;
                bucket_pointer[bucket_index+1] = bucket_pointer[bucket_index] +right_rows;
                int leaf_index_change = bucket_pointer[bucket_index-1];
                leaf_elements_pointer = leaf_elements_pointer+leaf_index_change*(2+d);
            }

            //Filling the leaf elements by carefully changing the leaf_elements_pointer, using the bucket number elements.
            for(int i=0; i<left_rows; i++){
                leaf_elements_pointer[0] = indices_left[i];
                leaf_elements_pointer[1] = root_key;

                for(int j=0; j<d; j++){
                    leaf_elements_pointer[2+j] = X_left[i*d+j];
                }

                leaf_elements_pointer += 2+d;
            }

            for(int i=0; i<right_rows; i++){
                leaf_elements_pointer[0] = indices_right[i];
                leaf_elements_pointer[1] = root_key;
       
                for(int j=0;j<d; j++){
                    leaf_elements_pointer[2+j] = X_right[i*d+j];
                }
              
                leaf_elements_pointer += 2+d;
            }
           
            bucket_index = bucket_index+2;   

            //Freeing whatever is not useful
            if(root_key!=0){
                free(X_sub);
            }
            free(X_left);
            free(X_right);
            free(D);
            free(indices_global);
            free(indices_left);
            free(indices_right);
            free(indices_local);
        }
        //Finished with the buckets/ leaves


        else{
            //Freeing what I don't need, and recursively call the function
            if(root_key!=0){
                free(X_sub);
            }
            free(indices_global);
            free(D);
            free(indices_local);

            level_count++;
            
            tree_create(X_left, indices_left, left_rows, d, root_key*2+1, bucket_pointer, level_count, final_levels);
            tree_create(X_right, indices_right, right_rows, d, root_key*2+2, bucket_pointer, level_count, final_levels);
        }
    }
}


/**
 * This is the function the user can use to create the tree. The way we decided to create the tree was to use a SERIAL REPRESENTATION of this data structure. 
 * The reason for this choice is that with this way we can EASILY pass the whole tree AT ONCE, using the MPI. This implementation was able to allow us to write the algorithms for
 * the tree creation and search with a reasonable time complexity and a space complexity not that much larger than the one in the previous versions(V0, V1).
 * The space complexity could maybe have been a little bit smaller, but the code would have been more grumpy and the complexity improvement was quite insignificant.
 * 
 * The basic idea behind the tree creation is the idea you proposed to us about using buckets of elements with B>1 elements. We did use buckets of elements, but they do not hold a specific 
 * number of elements in each tree creation. The criterion by which we choose the amount of elements is so that the final vp tree will have the "form of a balanced binary tree". So, we create
 * buckets in the last level of the tree with the appropriate size, so that our tree HAS A BALANCED FORM. Of course it also helped in the sense that we generally had to make less distance computations
 * to reach the bottom of the tree
 * 
 * The tree is represented in an array of doubles. The array has three different parts.
 * 
 * In the first part, the number of elements is equal to the number of buckets the tree has. In each element we have stored the amount of points each bucket has, but in a cumulative fashion
 * (i.e. in the first element we have stored the amount of points up until the first bucket, in the kth element the amount of points up until the kth bucket and so on).
 * Helpful so that we can easily access the bucket elements.
 * 
 * In the second part, we store the vantage points. Along with the d dimensions, we also store 3 more numbers associated with the point.
 * The first one is the index of this point. We store this, because the knnstruct needs the indices of the elements and not just the distances.
 * The second one is the key this point has in the tree structure. Useful to easily access the other family "nodes"
 * The third one is the median distance this vantage point has when the space partition takes place. It is later also used in the search algorithm, so we store it now so that we dont have to compute it again
 * The rest d elements are the coordinates of the vp.
 * 
 * In the third part, we store the leaf points, which are all in buckets (balanced form).Along with the d dimensions, we also store 2 more numbers associated with the point
 * The first one is the index of the point (useful for the knn struct).
 * The second one is key THAT THE FATHER of the point has. A number useful when performing the search.
 * The rest d elements are the coordinates of the leaf point. The leaf points are accessed using the elements in the first part of the array.
 * 
 * 
 * The tree_fun just calculates some parameters that are used and then calls the tree_create function which uses recursion to create the tree.
 */
double* tree_fun(double* X, int* indices, int d, int n){

    //Calculating the number of levels the tree will have.
    int levels = 0;       //Τhe levels the simple vp tree would have had
    int final_levels;
    int result = n;

    do{
        result = result/2;
        levels++;
    }while(result>0);

    if(levels<1){
        printf("Something went wrong with the number of levels... Exiting!\n");
        exit(-1);
    }
    else if(levels<3){
        final_levels = levels;
    }
    else if(levels<10){
        final_levels = levels-1;
    }
    else{
        final_levels = levels-2;
    }
    //Finished calculating the number of levels of the tree.


    //Calculating the number of elements, the size of the tree array and creating it
    int vp_elements_num = pow(2, final_levels-1) -1;
    int leaf_elements_num = n - vp_elements_num;
    int bucket_number = vp_elements_num +1;
    int total_size = bucket_number + vp_elements_num*(3+d) + leaf_elements_num*(2+d);

    double* tree = calloc(total_size, sizeof(double));
    double* bucket_pointer = tree;
    //Finished calculating and creating

    tree_create(X, indices, n, d, 0, bucket_pointer, 1, final_levels);

    return tree;
}


/**
 * Struct used as a return value in the tree_search function below
 */
typedef struct treeSearch{
    int bucket_id;
    double* dist;
} treeSearch;


/** 
 * This function takes a specific tree, a specific query point and a specific root and takes the path until the bottom of the tree, saving the distances and the indices in the meantime.
 * It stores the values in the result_struct which is passed as an argument.
 * If the root passed as an argument is the first one (id==0), then it just creates the indices and distances arrays of the struct, else it "merges" the distances calculated by this
 * specific search with the ones that already are in the result_struct.
 * 
 * The return value of the function is a struct that contains 2 variables:
 *  1)int bucket_id --> This is the key of the bucket in which we land. Or to state it in a better sense, the key that the bucket would have had if we stored the bucket keys. 
 *                      Used to find the key of the parent (useful later).
 *  2)double* dist  --> This pointer holds all the distances of the query point from the vantage points it encountered.
 *                      The way the distances are stored is that the first element of the array is the distance first vantage point and so on.
 *                      Later when we go up the tree we access the array from the last element to the first in order to acces the correct distance (more in whole_search function).
 * 
 * This method is later used in the whole search recursive method .-
 */
treeSearch tree_search(double* tree, double* query, int root_id, int d, int levels_num, int vp_elements_num, double* bucket_pointer, knnresult* result_struct){
    
    /*Calculating and initialising many variables that will be later used. Their use is what their names imply.*/
    double* vp_elements_pointer = bucket_pointer + vp_elements_num + 1;
    double* leaf_elements_pointer = vp_elements_pointer + vp_elements_num*(3+d);

    //Arrays to be used to form the knnresult struct
    double* distances = malloc((levels_num-1)*sizeof(double));      //Initialy allocating as a levels_num size, will later reallocate it for as much is needed (when I learn the amount of elements in the bucket).
    int* indices = malloc((levels_num-1)*sizeof(int));


    if(distances==NULL || indices==NULL){
        printf("Error: Couldn't allocate memory.\n");
        exit(-1);
    }
    
    double* vp_pointer_temp;
    double* leaf_pointer_temp;
    double distance;                    //This variable will temporarily hold the distance between two points (the one will be the query, and the other some point of the tree)
    
    int bucket_points;
    int leaf_index;
    int bucket_index;
    int k = result_struct->k;
    int current_node_id = root_id;

    //Initialising the struct to be returned
    treeSearch return_struct;
    return_struct.dist = malloc((levels_num-1)*sizeof(double));
    if(return_struct.dist==NULL){
        printf("Error: Couldn't allocate memory.\n");
        exit(-1);
    }


    int count = 0;                      //Variable counting the amount of elements I calculated the distances from
    /*Finished calculating and initialising.*/


    //Moving down the tree and stopping when a bucket is found
    while(current_node_id < vp_elements_num){
       
        vp_pointer_temp = vp_elements_pointer + (3+d)*current_node_id;      //Now pointing at the root "sub-array"
        distance = point_distance(query, vp_pointer_temp+3, d );            //+3 Because the first three elements are not coordinates.
        distances[count] = distance;                                        //Placing the results                       
        indices[count] = vp_pointer_temp[0];

        return_struct.dist[count] = distance;

        count++;
         

        //Following the path of the left child
        if(distance < (int)*(vp_pointer_temp+2) ){
            current_node_id = 2*current_node_id +1;               //left child
        }
        //Following the path of the right child
        else{
            current_node_id = 2*current_node_id +2;               //right child
        }
    }
    //Finished working with the vantage points


    return_struct.bucket_id = current_node_id;


    //Working with the right bucket, in the case where I have to search there
    if(current_node_id%2==0){
            
        bucket_index = current_node_id - vp_elements_num -1;              //The index i will use to access the correct bucket
        leaf_index = bucket_pointer[bucket_index];
        bucket_points = bucket_pointer[bucket_index+1] -leaf_index;

        distances = realloc(distances, (levels_num-1+bucket_points)*sizeof(double));          //Reallocating so that the arrays have the proper size.
        indices = realloc(indices, (levels_num-1+bucket_points)*sizeof(int));

        leaf_pointer_temp = leaf_elements_pointer + leaf_index*(2+d);
       
        for(int i=0; i<bucket_points; i++){
            distance = point_distance(query, leaf_pointer_temp +2, d);
            distances[count] = distance;
            indices[count] = leaf_pointer_temp[0];
            count++;
            leaf_pointer_temp += 2+d;
        }
    }
    //Finished working with the right bucket.


    //Working with the left bucket, in the case where I have to search there
    else if(current_node_id%2==1){
        
        bucket_index = current_node_id - vp_elements_num;              //The index I will use to access the correct bucket
     
        if(bucket_index==0){
            leaf_index = 0;
        }
        else{
            leaf_index = bucket_pointer[bucket_index-1];
        }
      
        bucket_points = bucket_pointer[bucket_index] - leaf_index;

        distances = realloc(distances, (levels_num-1+bucket_points)*sizeof(double));      //Reallocating so that the arrays have the proper size.
        indices = realloc(indices, (levels_num-1+bucket_points)*sizeof(int));


        leaf_pointer_temp = leaf_elements_pointer+leaf_index*(2+d);

        for(int i=0; i<bucket_points; i++){
            distance = point_distance(query, leaf_pointer_temp +2, d );
            distances[count] = distance;
            indices[count] = leaf_pointer_temp[0];
            count++;
            leaf_pointer_temp += 2+d;
        }
    }
    //Finished working with the left bucket.

    
    //Creating the kNN arrays (distances and indices)
    double* struct_distances = malloc(k*sizeof(double));
    int* struct_indices = malloc(k*sizeof(int));
    if(struct_distances==NULL || struct_indices==NULL){
        printf("Error: Couldn't allocate memory.\n");
        exit(-1);
    }


    //This is in case I still have not gathered k distances, so I have to fill up the rest of the elements with DBL_MAX.
    if(k>=count+1){
        quick_select(distances, indices, count, count);

        quick_sort(distances, indices, count);

        for(int i=0; i<count; i++){
			struct_distances[i] = distances[i];
			struct_indices[i] = indices[i];
		}

		for(int i=count; i<k; i++){
            struct_distances[i] = DBL_MAX;
		    struct_indices[i] = -1;
	    }
    }

    else{
        quick_select(distances, indices, count, k);

        quick_sort(distances, indices, k);

        for(int i=0; i<k; i++){
			struct_distances[i] = distances[i];
			struct_indices[i] = indices[i];
		}
	}
    //Finished creating the knn arrays


    //If root_id==0 I just store the arrays in the struct
    if(root_id==0){
        result_struct->nidx = struct_indices;
        result_struct->ndist = struct_distances;
    }
    
    //Else I "merge" the arrays with the ones I had in the old struct.
    else{
        knnresult local;
        local.k = k;
        local.m = result_struct->m;
        local.nidx = struct_indices;
        local.ndist = struct_distances;

        merge_structs(result_struct, local);
    }


    free(distances);
    free(indices);


    return return_struct;
}

/**
 * This function is called by the kNN to alter the struct until it reaches its final correct form.
 * It is a recursive function, and it calls the tree_search function everytime there is an intersection.
 *  The way we judge if there is an intersection is by comparing three distances:
 *  1)The median distance of the vantage point
 *  2)The distance between the vantage point and the query (already calculated in tree_search)
 *  3)The maximum distance stored in the result_struct
 * The way we compare them can be seen below. Special care was taken, because we have to compare THE DISTANCES, and NOT THE SQUARES OF THE DISTANCES, which are generally stored in the arrays.-
 */
void whole_search(double* tree, double* query, int root_id, int d, int levels_num, int vp_elements_num, double* bucket_pointer, knnresult* result_struct, int root_level){
   
    //Calculating and initialising useful variables for the search
    double* vp_elements_pointer = bucket_pointer + vp_elements_num+1;
    double* leaf_elements_pointer = vp_elements_pointer + vp_elements_num*(3+d);
    
    int current_root_id;
    
    double distance;
    double* current_point;

    treeSearch search;
    // Finished initialising
   
   //If clause used to control the recursion
    if(root_level==levels_num){
        //Updating the kNN according to the points contained in the bucket
        search = tree_search(tree, query, root_id, d ,levels_num, vp_elements_num, bucket_pointer, result_struct);

        free(search.dist);
    }


    else if(root_level<levels_num && root_level>0){
        //Calling the tree_search on the root of the subtree.
        search = tree_search(tree, query, root_id, d , levels_num, vp_elements_num, bucket_pointer, result_struct);

        current_root_id = search.bucket_id;
       
        //In the for iteration we begin to make comparisons with the nodes, until we reach the top of the tree. If there is an intersection, then we call the same whole_search on the appropriate child
        for(int i=0; i<levels_num - root_level; i++){
            
            //Access the array from the last element to the first to take the appropriate distance
            distance = search.dist[levels_num - root_level -1 - i];
          
            //Currently working with a right child
            if(current_root_id%2==0){
                current_root_id = current_root_id/2-1;          //Getting the father Node.
                current_point = vp_elements_pointer + current_root_id*(3+d);


                //The condition that tells us when there is an intersection is
                //Distmax + median_distance > dist(query, vp)
                //We get those values and create the if clause.
                if(sqrt(result_struct->ndist[result_struct->k -1]) +sqrt(current_point[2]) > sqrt(distance) ){
                    whole_search(tree, query, 2*current_root_id +1, d, levels_num, vp_elements_num, bucket_pointer, result_struct, levels_num -i);  
                }
            }


            //Currenty working with a left child
            else{
                current_root_id = current_root_id/2;            //Getting the father Node.
                current_point = vp_elements_pointer + current_root_id*(3+d);

                //The condition that tells us when there is an intersection is
                //Distmax + median_distance - dist(query, vp)> 0
                //We get those values and create the if clause.
                if(sqrt(result_struct->ndist[result_struct->k -1]) +sqrt(current_point[2]) -sqrt(distance) > 0){
                    whole_search(tree, query, 2*current_root_id +2, d, levels_num, vp_elements_num, bucket_pointer,result_struct, levels_num -i);
                }
            }
        }

        free(search.dist);

    }
}


/**
 * The knnV2 calls the whole_search to alter a struct that it already has created , and returns it.
 * It calculates the kNN for one specific point at a time and then creates the kNN for the whole Y.
 * */
knnresult knnV2(double* tree, double* Y, int d, int n, int k, int m){

    /*Calculating useful stuff, that will help us to make the search*/
    int levels = 0;       //the levels the simple vp tree would have had
    int final_levels;
    int result = n;
    do{
        result = result/2;
        levels++;
    }while(result>0);

    if(levels<1){
        printf("Something went wrong with the number of levels... Exiting!\n");
        exit(-1);
    }
    else if(levels<3){
        final_levels = levels;
    }
    else if(levels<10){
        final_levels = levels-1;
    }
    else{
        final_levels = levels-2;
    }

    int vp_elements_num = pow(2, final_levels-1) -1;
    double* bucket_pointer = tree;
    /*Finished calculating */

    knnresult ret_struct;
    ret_struct.k = k;
    ret_struct.m = m;
    double* ndist = malloc(k*m*sizeof(double));
    int* nidx = malloc(k*m*sizeof(int));
    if(ndist==NULL || nidx==NULL){
        printf("Error: Couldn't allocate memory.\n");
        exit(-1);
    }

    //Keeps the results of a single query point
    knnresult* small;

    //Iterate for every row(point) of Y
    for(int i=0; i<m; i++){
        small = malloc(sizeof(knnresult));
        if(small==NULL){
            printf("Error: Couldn't allocate memory.\n");
            exit(-1);
        }
        small->k = k;
        small->m = 1;

        whole_search(tree, Y+i*d, 0, d, final_levels, vp_elements_num, bucket_pointer, small, 1);

        //Passing the kNN to the array to be returned
        for(int j=0; j<k; j++){
            ndist[j+i*k] = small->ndist[j];
            nidx[j+i*k] = small->nidx[j];
        }

        free(small->ndist);
        free(small->nidx);
        free(small);
    }

    ret_struct.nidx = nidx;
    ret_struct.ndist = ndist;
   
    return ret_struct;
}


knnresult distrAllkNN2(double* X, int n, int d, int k){

    /*Initialising all the variables that the MPI processes need.*/
    int comm_size,rank;
    int subprocess_rows;                         //The amount of rows the subprocesses will have
    int root_process_rows;                       //The amount of rows the root process will have
    int current_rows;                            //The amount of rows the current process has
    double* X_local;        
    double* Y_tree;                              //The Y and Y arrays will hold trees that will be passed among the processes
    double* Z_tree;
    int root_tree_size;
    int sub_tree_size;
    knnresult result_local;                 //The struct containing the results of the process (It changes every time we get a new Z)
    knnresult result_temp;                  //The struct containing the result of the kNN search between the X_local and the Y array.

    MPI_Request send_request;
    MPI_Request recv_request;
    MPI_Status send_status;
    MPI_Status recv_status;
    /*Finished initialising*/


    /*Getting info about the processes*/
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    subprocess_rows=n/comm_size;                                //We assign a standard n/comm_size rows to each subprocess
    root_process_rows = n - (comm_size-1)*(n/comm_size);        //We assign the rest to the root process.    
    if(rank==0){
        current_rows = root_process_rows;
    }
    else{
        current_rows = subprocess_rows;
    }
    /*Finished getting info and allocating the memory*/



    /*Calculating the size of the trees, by having in mind that the process will have root_process_rows (it is the maximum)*/
    //First we calculate the size that the root process tree will have
    int levels = 0;       
    int final_levels;
    int result = root_process_rows;
    do{
        result = result/2;
        levels++;
    }while(result>0);
    
    if(levels<1){
        printf("Something went wrong with the number of levels... Exiting!\n");
        exit(-1);
    }
    else if(levels<3){
        final_levels = levels;
    }
    else if(levels<10){
        final_levels = levels-1;
    }
    else{
        final_levels = levels-2;
    }

    int vp_elements_num = pow(2, final_levels-1) -1;
    int leaf_elements_num = root_process_rows - vp_elements_num;
    int bucket_number = vp_elements_num +1;
    root_tree_size = bucket_number + vp_elements_num*(3+d) + leaf_elements_num*(2+d);

    //Now we calculate the size for the sub processes as well
    result = subprocess_rows;
    levels = 0;

    do{
        result=result/2;
        levels++;
    }while(result>0);
    
    if(levels<1){
        printf("Something went wrong with the number of levels... Exiting!\n");
        exit(-1);
    }
    else if(levels<3){
        final_levels = levels;
    }
    else if(levels<10){
        final_levels = levels-1;
    }
    else{
        final_levels = levels-2;
    }

    vp_elements_num = pow(2, final_levels-1) -1;
    leaf_elements_num = subprocess_rows- vp_elements_num;
    bucket_number = vp_elements_num +1;
    sub_tree_size = bucket_number + vp_elements_num*(3+d) + leaf_elements_num*(2+ d);
    /*Finished calculating the size of the trees*/


    /*Using Scatterv to send the X blocks to each MPI process and freeing the memory I don't need*/
    X_local = malloc(current_rows*d*sizeof(double));               //Allocating X_local memory
    int* send_counts = malloc(comm_size*sizeof(int));
    int* displacements = malloc(comm_size*sizeof(int));

    if(X_local==NULL || send_counts==NULL || displacements==NULL){
        printf("Error: Couldn't allocate memory.\n");
        exit(-1);
    }

    send_counts[0] = root_process_rows*d;
    displacements[0] = 0;
    for(int i=1; i<comm_size; i++){
        send_counts[i] = subprocess_rows*d;
        displacements[i] = (root_process_rows +(i-1)*subprocess_rows)*d;
    }

    MPI_Scatterv(X, send_counts, displacements, MPI_DOUBLE, X_local, root_process_rows*d, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    if(rank==0){
        free(X);
    }
    /*Finished sending, and freeing */



    /*Creating the trees in each process  */
    int* indices = malloc(current_rows*sizeof(int));
    if(indices==NULL){
        printf("Error: Couldn't allocate memory.\n");
        exit(-1);
    }

    if(rank==0){
        for(int i=0;i<current_rows; i++){
            indices[i] = i;
        }
    }
    else{
        for(int i=0;i<current_rows; i++){
            indices[i] = i+ root_process_rows +subprocess_rows*(rank-1) ;
        }
    }


    Y_tree = tree_fun(X_local, indices, d, current_rows);
    Z_tree = malloc(root_tree_size*sizeof(double));

    if(Z_tree==NULL){
        printf("Error: Couldn't allocate memory.\n");
        exit(-1);
    }

    //Reallocating the Y_trees on the sub processes so that they take as much memory as the root process tree.
    if(rank !=0){
        Y_tree = realloc(Y_tree, root_tree_size*sizeof(double));
    }

   
    result_local = knnV2(Y_tree, X_local, d, current_rows, k+1, current_rows);
    /*Finished creating the trees and got the first knn struct*/


    
    /*Below we make the message passing comm_size -1 times, and we use the merge_structs so that the result_local has the kNN results, which are refreshed every time we get a a new block*/
    for(int i=0; i<comm_size -1; i++){
        if(rank==comm_size -1){
            MPI_Isend(Y_tree, root_tree_size, MPI_DOUBLE, 0, 4, MPI_COMM_WORLD, &send_request);
            MPI_Irecv(Z_tree, root_tree_size, MPI_DOUBLE, comm_size-2, 4, MPI_COMM_WORLD, &recv_request);
        }
        else if(rank==0){
            MPI_Isend(Y_tree, root_tree_size, MPI_DOUBLE, 1, 4, MPI_COMM_WORLD, &send_request);
            MPI_Irecv(Z_tree, root_tree_size, MPI_DOUBLE, comm_size-1, 4, MPI_COMM_WORLD, &recv_request);
        }

        else{
            MPI_Isend(Y_tree, root_tree_size, MPI_DOUBLE, rank+1, 4, MPI_COMM_WORLD, &send_request);
            MPI_Irecv(Z_tree, root_tree_size, MPI_DOUBLE, rank-1, 4, MPI_COMM_WORLD, &recv_request);
        }   

        if (i!=0){
            if(rank ==i){
                result_temp = knnV2(Y_tree, X_local, d, root_process_rows, k+1, current_rows);
                merge_structs(&result_local, result_temp);
            }
            else{
                result_temp = knnV2(Y_tree, X_local, d, subprocess_rows, k+1, current_rows);
                merge_structs(&result_local, result_temp);
            }
        }

        MPI_Wait(&send_request, &send_status);
        MPI_Wait(&recv_request, &recv_status); 
        pointer_swap(&Y_tree, &Z_tree);
    }

    if(rank==comm_size-1){
        result_temp = knnV2(Y_tree, X_local, d, root_process_rows, k+1, current_rows);
        merge_structs(&result_local, result_temp);
    }
    else{
        result_temp = knnV2(Y_tree, X_local, d, subprocess_rows, k+1, current_rows);
        merge_structs(&result_local, result_temp);
    }
    /*Finished the message passing and the kNN computations.*/

    
    free(X_local);
    free(Y_tree);
    free(Z_tree);
    
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
    
    knnresult a = distrAllkNN2(X, n, d, k);


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