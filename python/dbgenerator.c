#include <Python.h>
#include <arrayobject.h>

static PyObject* dbg_upslope_area(PyObject *self, PyObject *args);

/*
 * Functions in this module.
 */
static PyMethodDef DrainageBasinMethods[] = {
  {"upslope_area",
   dbg_upslope_area,
   METH_VARARGS,
   "compute logit"},
  {NULL,
   NULL,
   0,
   NULL}
};

/* Upslope area computation */
static PyObject* dbg_upslope_area(PyObject *self, PyObject *args)
{
  PyObject *arg1=NULL, *arg2=NULL, *out=NULL;
  PyObject *dem=NULL, *mask=NULL, *output=NULL;

  int ndims;
  npy_intp *dims;
  double *dem_data, *mask_data;

  if (!PyArg_ParseTuple(args, "OOO!", &arg1, &arg2,
                        &PyArray_Type, &out))
    return NULL;

  dem = PyArray_FROM_OTF(arg1, NPY_DOUBLE, NPY_IN_ARRAY);
  if (dem == NULL) return NULL;

  mask = PyArray_FROM_OTF(arg2, NPY_DOUBLE, NPY_IN_ARRAY);
  if (mask == NULL) goto fail;

  output = PyArray_FROM_OTF(out, NPY_DOUBLE, NPY_INOUT_ARRAY);
  if (output == NULL) goto fail;

  /* code that makes use of arguments */
  ndims = PyArray_NDIM(dem);
  if (ndims != PyArray_NDIM(mask))
    goto fail;

  if (ndims != PyArray_NDIM(output))
    goto fail;

  dims = PyArray_DIMS(dem);
  dem_data = (double*)PyArray_DATA(dem);
  mask_data = (double*)PyArray_DATA(mask);

  /* You will probably need at least
     nd = PyArray_NDIM(<..>)    -- number of dimensions
     dims = PyArray_DIMS(<..>)  -- npy_intp array of length nd
     showing length in each dim.
     dptr = (double *)PyArray_DATA(<..>) -- pointer to data.

     If an error occurs goto fail.
  */

  Py_DECREF(dem);
  Py_DECREF(mask);
  Py_DECREF(output);
  Py_INCREF(Py_None);
  return Py_None;

 fail:
  Py_XDECREF(dem);
  Py_XDECREF(mask);
  PyArray_XDECREF_ERR(output);
  return NULL;
}

/* This initiates the module using the above definitions. */
#if PY_VERSION_HEX >= 0x03000000
static struct PyModuleDef moduledef = {
  PyModuleDef_HEAD_INIT,
  "spam",
  NULL,
  -1,
  SpamMethods,
  NULL,
  NULL,
  NULL,
  NULL
};

PyObject *PyInit_spam(void)
{
  PyObject *m;
  m = PyModule_Create(&moduledef);
  if (!m) {
    return NULL;
  }
  return m;
}
#else
PyMODINIT_FUNC initdbgenerator(void)
{
  PyObject *m;

  m = Py_InitModule("dbgenerator", DrainageBasinMethods);
  if (m == NULL) {
    return;
  }

  import_array();               /* initialize numpy */
}
#endif
