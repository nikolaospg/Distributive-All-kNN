#include <stdlib.h>
#include <stdio.h>
#include <string.h>

//Reads the name of the input file and sets the variables n, d and numbered
void find_sizes(char* filename, int* n, int* d, int* numbered){
    
    if (strstr(filename, "ColorHistogram.asc") != NULL){
        *n = 68040;
        *d = 32;
        *numbered = 1;
    }
    else if (strstr(filename, "ColorMoments.asc") != NULL){
        *n = 68040;
        *d = 9;
        *numbered = 1;
    }
    else if (strstr(filename, "CoocTexture.asc") != NULL){
        *n = 68040;
        *d = 16;
        *numbered = 1;
    }
    else if (strstr(filename, "LayoutHistogram.asc") != NULL){
        *n = 68040;
        *d = 32;
        *numbered = 1;
    }
    else if (strstr(filename, "features.csv") != NULL){
        *n = 106574;
        *d = 518;
        *numbered = 4;
    }
    else if (strstr(filename, "MiniBooNE_PID.txt") != NULL){
        *n = 130064;
        *d = 50;
        *numbered = 3;
    }
    else if (strstr(filename, "BBC.txt") != NULL){
        *n = 17720;
        *d = 17;
        *numbered = 2;
    }
    else if (strstr(filename, "CNN.txt") != NULL){
        *n = 22545;
        *d = 17;
        *numbered = 2;
    }
    else if (strstr(filename, "CNNIBN.txt") != NULL){
        *n = 33117;
        *d = 17;
        *numbered = 2;
    }
    else if (strstr(filename, "NDTV.txt") != NULL){
        *n = 17051;
        *d = 17;
        *numbered = 2;
    }
    else if (strstr(filename, "TIMESNOW.txt") != NULL){
        *n = 39252;
        *d = 17;
        *numbered = 2;
    }
    else{
        printf("No file found.\n");
        exit(1);
    }
}


double* getX(char* filename, int n, int d, int numbered){

    double* X = NULL;

    X = (double*)malloc(n*d*sizeof(double));

    if (X==NULL){
        printf("Couldn't allocate memory for X\n");
        exit(-1);
    }

    FILE *matFile = fopen(filename, "r");

	if (matFile == NULL){
		printf("Couldn't open file\n");
        exit(1);
	}

	double num,temp;

	int i=0,j=0;
	
    if(numbered==1){
		while (fscanf(matFile,"%lf", &num) != EOF){
			if(i%(d+1)!=0){
    			X[j]=num;
				j++;
			}
			i++;
		}
	}
    else if(numbered==2){
		for(i=0; i<n; i++){
			fscanf(matFile,"%lf",&num);
			for(j=0; j<d; j++){
	    		if(fscanf(matFile," %lf:%lf",&temp,&num)==EOF) break;
		    		X[i*d+j]=num;
        	}
			fscanf(matFile,"%*[^\n]\n");
		}
	}
    else if(numbered==3){
        //Skip first line
		fscanf(matFile,"%*[^\n]\n");
        
		for(i=0; i<n; i++){
			for(j=0; j<d; j++){
				if(fscanf(matFile,"%lf",&num)==EOF) break;
				X[i*d+j]=num;
       		}
			fscanf(matFile,"%*[^\n]\n");
		}
    }
    
    else{
        //Skip first four lines
		for(int skip=0; skip<4; skip++){
			fscanf(matFile,"%*[^\n]\n");
		}

		for(i=0; i<n; i++){
			fscanf(matFile,"%lf",&num);
			for(j=0; j<d; j++){
				if(fscanf(matFile,",%lf",&num)==EOF) break;
				X[i*d+j]=num;
       		}
			fscanf(matFile,"%*[^\n]\n");
		}
	}

	fclose(matFile);

    return X;
}

//Thanks to Chris Pavlidis for the information