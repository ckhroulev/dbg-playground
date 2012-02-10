int foo(double *data, int len) {
  int j;
  for (j = 0; j < len; ++j)
    data[j] += 1;

  return 0;
}
