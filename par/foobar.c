#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv) {

  double *data = NULL;
  int N = 16, i, num = 0;

  data = (double*)malloc(N * sizeof(double));
  if (data == NULL)
    return 1;

  /* beginning of the parallel block */
#pragma omp parallel shared(data) private(i)
  {
    int my_num;

#pragma omp critical
    my_num = num++;

#pragma omp for
    for (i = 0; i < N; ++i) {
      data[i] = my_num;
    }
  } /* end of the parallel block */

  for (i = 0; i < N; ++i)
    printf("data[%d] = %f\n", i, data[i]);

  free(data);

  return 0;
}
