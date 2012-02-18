#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

int main(int argc, char **argv) {

  double *data = NULL;
  int N = 10, i;

  data = (double*)malloc(N * sizeof(double));
  if (data == NULL)
    return 1;

#pragma omp parallel shared(data) private(i)
  {
    int num = omp_get_thread_num();

#pragma omp for
    for (i = 0; i < N; ++i) {
      data[i] = num;
    }
  }
  /* end of the parallel block */

  for (i = 0; i < N; ++i) {
    printf("data[%d] = %f\n", i, data[i]);
  }

  free(data);


  return 0;
}
