# include <stdio.h>
# include <stdlib.h>
# include <mpi.h>

# define ROOT_PROCESS 0

int N;
int currentProcessID, numberOfProcesses;
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

int main(int argc, char *argv[]) {


  char * matrixFilePath = argv[1];
  char * vectorFilePath = argv[2];

  MPI_Init(NULL, NULL);

  MPI_Comm_size(MPI_COMM_WORLD, &numberOfProcesses);
  MPI_Comm_rank(MPI_COMM_WORLD, &currentProcessID);

  if (currentProcessID == ROOT_PROCESS) {
    N = GetNumberOfLines(vectorFilePath);
    matrix = malloc(N*N*sizeof(int));
    vector = malloc(N*sizeof(int));
    FillMatrixOrVector(matrixFilePath, matrix);
    FillMatrixOrVector(vectorFilePath, vector);
    for (int i = 0 ; i < N*N ; i++) {
      (i == N*N - 1) ? printf("%d \n", matrix[i]) : printf("%d ", matrix[i]);
    }
    for (int i = 0 ; i < N ; i++) {
      (i == N - 1) ? printf("%d \n", vector[i]) : printf("%d ", vector[i]);
    }
  }


  MPI_Finalize();

  return 0;
}
