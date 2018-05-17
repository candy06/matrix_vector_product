/*
********************************************************************************
*                            Polytech'Nice Sophia
*
*
* Filename : main.c
* Programmer : Loïc ROSE
* Description : algorithm to do a matrix vector product using OpenMPI
********************************************************************************
*/

/*
********************************************************************************
*                               INCLUDE FILES
********************************************************************************
*/

# include <stdio.h>
# include <stdlib.h>
# include <mpi.h>
# include <omp.h>

/*
********************************************************************************
*                               CONSTANTS & MACRO
********************************************************************************
*/

# define ROOT_PROCESS 0

/*
********************************************************************************
*                               GLOBAL VARIABLES
********************************************************************************
*/

int N;
int currentProcessID, P, nextNode, prevNode, lastNode;
int * matrix = NULL;
int * vector = NULL;
int * submatrix = NULL;
int * computedResult = NULL;
int * result = NULL;

/*
********************************************************************************
*                               GetNumberOfLines()
*
* Description : method that, given a file path, open the file and count the
                number of lines in it
* Arguments   : the filename (char *)
* Returns     : number of lines in the file
********************************************************************************
*/

int GetNumberOfLines(const char * filename) {
  FILE * matrixFile = fopen(filename, "r");
  int lines = 0;
  while (!feof(matrixFile)) {
    char c = fgetc(matrixFile);
    if (c == '\n') {
      lines++;
    }
  }
  fclose(matrixFile);
  return lines;
}

/*
********************************************************************************
*                               FillMatrixOrVector()
*
* Description : method that, given a file path, open the file and fill the
                given array of int with the integers read in it
* Arguments   : the filename (char *) and the matrix or vector (int *)
* Returns     : void
********************************************************************************
*/

void FillMatrixOrVector(const char * filename, int * array) {
  FILE * file = fopen(filename, "r");
  int i = 0;
  while (!feof(file)) {
    fscanf(file, "%d", &array[i++]);
  }
  fclose(file);
}

/*
********************************************************************************
*                               Broadcast()
*
* Description : method that, given an array and its size, will broadcast it
                to all the nodes of the ring
* Arguments   : the data we want to broadcast (int *) and its size (int)
* Returns     : void
********************************************************************************
*/

void Broadcast(int * data, int dataSize) {
  if (currentProcessID == ROOT_PROCESS) {
    MPI_Send(data, dataSize, MPI_INT, nextNode, 0, MPI_COMM_WORLD);
  } else if (currentProcessID == lastNode) {
    MPI_Recv(data, dataSize, MPI_INT, prevNode, 0,
      MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  } else {
    MPI_Recv(data, dataSize, MPI_INT, prevNode, 0,
      MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Send(data, dataSize, MPI_INT, nextNode, 0, MPI_COMM_WORLD);
  }
}

/*
********************************************************************************
*                               Scatter()
*
* Description : methat that, given a source data will send a different block of
                size 'sizeDst' to all the nodes of the ring (i.e each node will
                have a portion of the source data) and store it in the array
                dataDst
* Arguments   : the source and destination arrays (int *) and their sizes (int)
* Returns     : void
********************************************************************************
*/

void Scatter(int * dataSrc, int sizeSrc, int * dataDst, int sizeDst) {
  if (currentProcessID == ROOT_PROCESS) {
    for (int i = sizeDst ; i < sizeSrc ; i += sizeDst) {
      MPI_Send(&dataSrc[i], sizeDst, MPI_INT, nextNode, 0, MPI_COMM_WORLD);
    }
    submatrix = &dataSrc[currentProcessID];
  } else {
    MPI_Recv(dataDst, sizeDst, MPI_INT, prevNode, 0,
      MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    for (int i = 0 ; i < (ROOT_PROCESS - 1 - currentProcessID + P) % P ; i++) {
      MPI_Send(dataDst, sizeDst, MPI_INT, nextNode, 0, MPI_COMM_WORLD);
      MPI_Recv(dataDst, sizeDst, MPI_INT, prevNode, 0,
        MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
  }
}

/*
********************************************************************************
*                               Gather()
*
* Description : send all the computed results (product of the submatrix and
                vector of each processor) to the root process
* Arguments   : the computed result array (int *) and its size
* Returns     : void
********************************************************************************
*/

void Gather(int * computedResult, int computedResultSize) {
  if (currentProcessID == ROOT_PROCESS) {
    // Copy the first computed elements of P0 into the beginning of result
    for (int k = 0 ; k < computedResultSize ; k++) {
      result[k] = computedResult[k];
    }
    for (int i = 0 ; i < (P - 1) % P ; i++) {
      MPI_Recv(&result[computedResultSize+i*computedResultSize],
        computedResultSize, MPI_INT, prevNode, 0,
        MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
  } else {
    MPI_Send(computedResult, computedResultSize, MPI_INT,
      nextNode, 0, MPI_COMM_WORLD);
    for (int i = 0 ; i < (currentProcessID - 1) % P ; i++) {
      MPI_Recv(computedResult, computedResultSize, MPI_INT, prevNode, 0,
        MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Send(computedResult, computedResultSize, MPI_INT,
        nextNode, 0, MPI_COMM_WORLD);
    }
  }
}

/*
********************************************************************************
*                               Compute()
*
* Description : compute the prduct (sub)matrix vector and store the result in
                the array computedResult
* Arguments   : the submatrix, vector and computedResult arrays (int *) and
                their sizes
* Returns     : void
* Note        : THIS METHOD IS PARALLEILIZED WITH OMP
********************************************************************************
*/

void Compute(int * submatrix, int submatrixSize, int * vector, int vectorSize,
                int * computedResult, int computedResultSize) {
  int tmp = 0;
  int i, j;
  # pragma omp parallel for private(tmp) private(j)
  for (i = 0 ; i < computedResultSize ; i++) {
    for (j = 0 ; j < vectorSize ; j++) {
      if (i == 0) {
        tmp += submatrix[j] * vector[j];
      } else {
        tmp += submatrix[j+vectorSize*i] * vector[j];
      }
    }
    computedResult[i] = tmp;
    tmp = 0;
  }
}

/*
********************************************************************************
*                               DisplayFinalResult()
*
* Description : method that display the array given in parameter in the correct
                format (for the exercise)
* Arguments   : the array (int *) and its size (int)
* Returns     : void
********************************************************************************
*/

void DisplayFinalResult(int * result, int sizeResult) {
  for (int i = 0 ; i < sizeResult ; i++) {
    printf("%d\n", result[i]);
  }
}

/*
********************************************************************************
*                               GetNumberOfLines()
*
* Description : the classic main method of a C program
* Arguments   : the first argument is the filepath of the file that has the
                matrix and the second one represents the filepath of the file
                that has the vector
* Returns     : 0
********************************************************************************
*/

int main(int argc, char *argv[]) {


  // Get the file paths of the files for the matrix and the vector
  char * matrixFilePath = argv[1];
  char * vectorFilePath = argv[2];

  // Initialize the context for MPI
  MPI_Init(NULL, NULL);

  // Set the number of processes to the variable P
  MPI_Comm_size(MPI_COMM_WORLD, &P);
  // Set the current process number to the variable currentProcessID
  MPI_Comm_rank(MPI_COMM_WORLD, &currentProcessID);

  // Index to know the next/previous node of the current one and the last node
  nextNode = (currentProcessID + 1) % P;
  prevNode = (currentProcessID - 1 + P) % P;
  lastNode = (ROOT_PROCESS + P - 1) % P;

  // Read the vector file (smaller) to get N in the ROOT process
  if (currentProcessID == ROOT_PROCESS) {
    N = GetNumberOfLines(vectorFilePath);
  }

  // Broadcast N to every processes
  Broadcast(&N, 1);

  // Allocate memory for the result array
  result = malloc(N*sizeof(int));

  // Allocate memory for the matrix and the vector
  matrix = malloc(N*N*sizeof(int));
  vector = malloc(N*sizeof(int));

  // Allocate memory for the submatrix of each process
  submatrix = malloc((N/P)*N*sizeof(int));

  // Allocate memory for the temporary results computed by each process
  computedResult = malloc((N/P)*sizeof(int));

  // Read files and fill matrix/vector arrays in the ROOT process
  if (currentProcessID == ROOT_PROCESS) {
    FillMatrixOrVector(matrixFilePath, matrix);
    FillMatrixOrVector(vectorFilePath, vector);
  }

  // Broadcast the vector to every nodes
  Broadcast(vector, N);

  // Scatter the matrix
  Scatter(matrix, N*N, submatrix, N/P*N);

  // Do the matrix/vector multiplication in every nodes
  Compute(submatrix, (N/P)*N, vector, N, computedResult, N/P);

  // Gather all the sub results in the root process
  Gather(computedResult, N/P);

  // Display the final result in the root process
  if (currentProcessID == ROOT_PROCESS) {
    DisplayFinalResult(result, N);
  }

  MPI_Finalize();

  return 0;
}
