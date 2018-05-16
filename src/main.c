# include <stdio.h>
# include <stdlib.h>
# include <mpi.h>

# define ROOT_PROCESS 0

int N;
int currentProcessID, P, nextNode, prevNode, lastNode;
int * matrix = NULL;
int * vector = NULL;
int * submatrix = NULL;
int * computedResult = NULL;
int * result = NULL;

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

void FillMatrixOrVector(const char * filename, int * array) {
  FILE * file = fopen(filename, "r");
  int i = 0;
  while (!feof(file)) {
    fscanf(file, "%d", &array[i++]);
  }
  fclose(file);
}

void Broadcast(int * data, int dataSize) {
  if (currentProcessID == ROOT_PROCESS) {
    MPI_Send(data, dataSize, MPI_INT, nextNode, 0, MPI_COMM_WORLD);
  } else if (currentProcessID == lastNode) {
    MPI_Recv(data, dataSize, MPI_INT, prevNode, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  } else {
    MPI_Recv(data, dataSize, MPI_INT, prevNode, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Send(data, dataSize, MPI_INT, nextNode, 0, MPI_COMM_WORLD);
  }
}

void Scatter(int * dataSrc, int sizeSrc, int * dataDst, int sizeDst) {
  if (currentProcessID == ROOT_PROCESS) {
    for (int i = sizeDst ; i < sizeSrc ; i += sizeDst) {
      MPI_Send(&dataSrc[i], sizeDst, MPI_INT, nextNode, 0, MPI_COMM_WORLD);
    }
    submatrix = &dataSrc[currentProcessID];
  } else {
    MPI_Recv(dataDst, sizeDst, MPI_INT, prevNode, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    for (int i = 0 ; i < (ROOT_PROCESS - 1 - currentProcessID + P) % P ; i++) {
      MPI_Send(dataDst, sizeDst, MPI_INT, nextNode, 0, MPI_COMM_WORLD);
      MPI_Recv(dataDst, sizeDst, MPI_INT, prevNode, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
  }
}


void Gather(int * computedResult, int computedResultSize) {
  if (currentProcessID == ROOT_PROCESS) {
    result = computedResult;
    for (int i = 0 ; i < (P - 1) % P ; i++) {
      //printf("%d\n", computedResultSize+i*computedResultSize);
      MPI_Recv(&result[computedResultSize+i*computedResultSize], computedResultSize, MPI_INT, prevNode, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
  } else {
    MPI_Send(computedResult, computedResultSize, MPI_INT, nextNode, 0, MPI_COMM_WORLD);
    for (int i = 0 ; i < (currentProcessID - 1) % P ; i++) {
      MPI_Recv(computedResult, computedResultSize, MPI_INT, prevNode, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Send(computedResult, computedResultSize, MPI_INT, nextNode, 0, MPI_COMM_WORLD);
    }
  }
}


void Compute(int * submatrix, int submatrixSize, int * vector, int vectorSize, int * computedResult, int computedResultSize) {
  int index = 0;
  int tmp = 0;
  for (int i = 0 ; i < computedResultSize ; i++) {
    for (int j = 0 ; j < vectorSize ; j++) {
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

void DisplayFinalResult(int * result, int sizeResult) {
  for (int i = 0 ; i < sizeResult ; i++) {
    printf("%d\n", result[i]);
  }
}

int main(int argc, char *argv[]) {


  char * matrixFilePath = argv[1];
  char * vectorFilePath = argv[2];

  MPI_Init(NULL, NULL);

  MPI_Comm_size(MPI_COMM_WORLD, &P);
  MPI_Comm_rank(MPI_COMM_WORLD, &currentProcessID);
  nextNode = (currentProcessID + 1) % P;
  prevNode = (currentProcessID - 1 + P) % P;
  lastNode = (ROOT_PROCESS + P - 1) % P;

  // Read the vector file (smaller) to get N in the ROOT process
  if (currentProcessID == ROOT_PROCESS) {
    N = GetNumberOfLines(vectorFilePath);
  }

  // Broadcast N to every processes
  Broadcast(&N, 1);
  result = malloc(N*sizeof(int));

  // Allocate memory for the matrix and the vector
  matrix = malloc(N*N*sizeof(int));
  vector = malloc(N*sizeof(int));

  // Allocate enough memory for the submatrix of each process
  submatrix = malloc((N/P)*N*sizeof(int));

  // Allocate enough memory for the temporary results computed by each process
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

  Compute(submatrix, (N/P)*N, vector, N, computedResult, N/P);

  /*
  for (int i = 0 ; i < N/P ; i++) {
     if (i == 0) {
       printf("From process ID %d : %d ", currentProcessID, computedResult[i]);
     } else if (i == N/P - 1) {
       printf("%d\n", computedResult[i]);
     } else {
       printf("%d ", computedResult[i]);
     }
  }*/

  //MPI_Gather(computedResult, N/P, MPI_INT, result, N/P, MPI_INT, ROOT_PROCESS, MPI_COMM_WORLD);
  Gather(computedResult, N/P);


  if (currentProcessID == ROOT_PROCESS) {
    DisplayFinalResult(result, N);
  }

  MPI_Finalize();

  return 0;
}
