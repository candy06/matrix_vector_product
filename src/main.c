# include <stdio.h>
# include <stdlib.h>
# include <mpi.h>

# define ROOT_PROCESS 0

int N;
int currentProcessID, P, nextNode, prevNode, lastNode;
int * matrix = NULL;
int * vector = NULL;

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

int main(int argc, char *argv[]) {


  char * matrixFilePath = argv[1];
  char * vectorFilePath = argv[2];

  MPI_Init(NULL, NULL);

  MPI_Comm_size(MPI_COMM_WORLD, &P);
  MPI_Comm_rank(MPI_COMM_WORLD, &currentProcessID);
  nextNode = (currentProcessID + 1) % P;
  prevNode = (currentProcessID - 1 + P) % P;
  lastNode = (ROOT_PROCESS + P - 1) % P;

  if (currentProcessID == ROOT_PROCESS) {
    N = GetNumberOfLines(vectorFilePath);
  }

  Broadcast(&N, 1);

  matrix = malloc(N*N*sizeof(int));
  vector = malloc(N*sizeof(int));

  if (currentProcessID == ROOT_PROCESS) {
    FillMatrixOrVector(matrixFilePath, matrix);
    FillMatrixOrVector(vectorFilePath, vector);
  }

  Broadcast(vector, N*N);

  //printf("ProcessID %d : N = %d\n", currentProcessID, N);
  //Broadcast(vector, N);


  for (int i = 0 ; i < N ; i++) {
    (i == N - 1) ? printf("%d \n", vector[i]) : printf("%d ", vector[i]);
  }



  MPI_Finalize();

  return 0;
}
