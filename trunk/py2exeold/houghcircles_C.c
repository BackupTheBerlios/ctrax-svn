#include "Python.h"
#include "numpy/arrayobject.h"
#include <stdlib.h>
#include <math.h>

#define SQUARED(XX) ((XX)*(XX))

static int binround(double x, double * edges, int nbins);
static void houghcircles_main(double * acc, double * x, double * y, double * w,
		       double* bina, double* binb, double * binr,
		       int npts, int na, int nb,
		       int nr);

static PyObject* houghcircles_C(PyObject *self,PyObject *args){

  int npts, na, nb, nr;  
  PyArrayObject * xarray, *yarray, *warray, *binedgesaarray,
    *bincentersbarray, *bincentersrarray, *accarray;
  double * x, *y, *w, *binedgesa, *bincentersb, *bincentersr, *acc;
  npy_intp dims[3];
  
  if (!PyArg_ParseTuple(args, "O!O!O!O!O!O!", 
			&PyArray_Type,&xarray,
			&PyArray_Type,&yarray,
			&PyArray_Type,&warray,
			&PyArray_Type,&binedgesaarray,
			&PyArray_Type,&bincentersbarray,
			&PyArray_Type,&bincentersrarray))
    return NULL;

  npts = xarray->dimensions[0];
  na = binedgesaarray->dimensions[0]-1;
  nb = bincentersbarray->dimensions[0];
  nr = bincentersrarray->dimensions[0];

  dims[0] = nr; dims[1] = nb; dims[2] = na;
  accarray = (PyArrayObject *) PyArray_SimpleNew(3,&dims[0],NPY_DOUBLE);
  //accarray=(PyArrayObject *) PyArray_FromDims(3,dims,NPY_DOUBLE);

  /* get inputs as 1-D C arrays*/ 
  x = (double*)xarray->data;
  y = (double*)yarray->data;
  w = (double*)warray->data;
  binedgesa = (double*)binedgesaarray->data;
  bincentersb = (double*)bincentersbarray->data;
  bincentersr = (double*)bincentersrarray->data;
  acc = (double*)PyArray_DATA(accarray);

  houghcircles_main(acc,x,y,w,binedgesa,bincentersb,bincentersr,
		    npts,na,nb,nr);

  return PyArray_Return(accarray);
}

static void houghcircles_main(double * acc, double * x, double * y, double * w,
			      double* bina, double* binb, double * binr,
			      int npts, int na, int nb,
			      int nr){
  
  int ib, ir, i;
  double rsquared, xpt, ypt;
  int pt, iacc;
  double d;
  int soln1, soln2;
  double *binrsquared;

  for(i = 0; i < na*nb*nr; i++) acc[i] = 0;

  binrsquared = (double*)malloc(sizeof(double)*nr);
  for(ir = 0; ir < nr; ir++){
    binrsquared[ir] = binr[ir]*binr[ir];
  }

  for(pt = 0; pt < npts; pt++){
    iacc = 0;
    xpt = x[pt];
    ypt = y[pt];
    for(ir = 0; ir < nr; ir++){
      rsquared = binrsquared[ir];
      for(ib = 0; ib < nb; ib++,iacc+=na){

	/* solve for a such that (x - a)^2 + (y - b)^2 = r^2
	   x^2 - 2ax + a^2 + (y - b)^2 - r^2 = 0
	   a^2 - 2ax + (x^2 + (y - b)^2 - r^2) = 0
	   a = x +/- sqrt(x^2 - (x^2 + (y - b)^2 - r^2))
	   = x +/- sqrt( r^2 - (y - b)^2 )
	   in order for there to be a solution for this b, we 
	   require r >= abs(y - b)*/

	d = SQUARED(ypt - binb[ib]);
	if( d > rsquared){
	  /* no solution */
	  continue;
	}
	d = sqrt(rsquared - d);
	soln1 = binround(xpt+d,bina,na);
	soln2 = binround(xpt-d,bina,na);
	
	/* one solution */
	if(soln1 == soln2){
	  /* out of bounds? */
	  if(soln1 < 0){
	    continue;
	  }
	  acc[iacc+soln1]+=w[pt];
	  continue;
	}

	/* two solutions */
	if(soln1 >= 0){
	  acc[iacc+soln1]+=w[pt];
	}
	if(soln2 >= 0){
	  acc[iacc+soln2]+=w[pt];
	}

      }
    }
    
  }
  
  free(binrsquared);
  binrsquared = 0;

  return;

}

static int binround(double x, double * edges, int nbins){

  int bin;

  /* check for out of bounds */
  if( (x < edges[0]) || (x >= edges[nbins]) ){
    return(-1);
  }

  /* for now, just do a linear search */
  for(bin = 1; bin <= nbins; bin++){
    if( x < edges[bin] ){
      return(bin-1);
    }
  } 
  return(-1);

}

static PyMethodDef Houghcircles_CMethods[] = {
    {"houghcircles_C",  houghcircles_C, METH_VARARGS},
    {NULL, NULL}        /* Sentinel (terminates structure) */
};

PyMODINIT_FUNC
inithoughcircles_C(void)
{
    (void) Py_InitModule("houghcircles_C", Houghcircles_CMethods);
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
    inithoughcircles_C();
	
    return 0;
}
