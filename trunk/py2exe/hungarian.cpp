/*
	hungarian
	Harold Cooper <hbc@mit.edu>
	hungarian.cpp 2004-09-02
*/

#include "Python.h"
#include "numpy/arrayobject.h"
#include "asp.h"
#include <stdio.h>

static PyObject *
hungarian(PyObject *self, PyObject *args)
//hungarian(costs)
{
  //PyObject *ocosts;
  PyArrayObject *costs;
  npy_intp n;
  long *rowsol;
  long *colsol;
  cost *buf,**ccosts;

  if (!PyArg_ParseTuple(args, "O!", &PyArray_Type,&costs))
    return NULL;

  if (costs->nd!=2)
    {
      PyErr_SetString(PyExc_ValueError,"hungarian() requires a 2-D matrix");
      goto error;
    }
  n = costs->dimensions[0];
  if (costs->dimensions[1]!=n)
    {
      PyErr_SetString(PyExc_ValueError,"hungarian() requires a square matrix");
      goto error;
    }


  //get inputted matrix as a 1-D C array:
  buf = (cost*)costs->data;
  //buf = (cost *)NA_OFFSETDATA(costs);

  //copy inputted matrix into a 2-dimensional C array:
  ccosts = (cost **)malloc(sizeof(cost *)*n);
  if(!ccosts)
    {
      PyErr_NoMemory();
      goto error;
    }
  for(int i=0;i<n;i++)
    ccosts[i] = &buf[i*n];

  rowsol = (long *)malloc(sizeof(long)*n);
  colsol = (long *)malloc(sizeof(long)*n);
  if(!(rowsol&&colsol))
    {
      PyErr_NoMemory();
      goto error;
    }

  //run hungarian!:
  asp(n,ccosts,rowsol,colsol);

  //NA_InputArray() incremented costs, but now we're done with it, so let it get GC'ed:
  //Py_XDECREF(costs);

  return Py_BuildValue("(OO)",PyArray_SimpleNewFromData(1,&n,PyArray_LONG,(char*)rowsol),PyArray_SimpleNewFromData(1,&n,PyArray_LONG,(char*)colsol));

 error:
  //Py_XDECREF(costs);
  return NULL;
}

static PyMethodDef HungarianMethods[] = {
    {"hungarian",  hungarian, METH_VARARGS,
     "Solves the linear assignment problem using the Hungarian\nalgorithm.\n\nhungarian() takes a single argument - a square cost matrix - and returns a\ntuple of the form\n(row_assigns,col_assigns)."},
    {NULL, NULL, 0, NULL}        /* Sentinel (terminates structure) */
};

PyMODINIT_FUNC
inithungarian(void)
{
    (void) Py_InitModule("hungarian", HungarianMethods);
	import_array();
}

int
main(int argc, char *argv[])
{
    /* Pass argv[0] to the Python interpreter */
    Py_SetProgramName(argv[0]);

    /* Initialize the Python interpreter.  Required. */
    Py_Initialize();

    /* Add a static module */
    inithungarian();

	return 0;
}
