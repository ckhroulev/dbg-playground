#include <math.h>
#include <stdio.h>

int foobar(double* data, int length, double* output) {
  int i;

#pragma omp parallel default(shared) private(i)
  {

#pragma omp for
    for (i = 0; i < length; ++i) {
      output[i] = sqrt(data[i]);
    }

  }

  return 0;
}
