#include <Python.h>
#include "geom2d.h"

static PyObject* geom2d_add(PyObject *self, PyObject *args)
{
    Point p;
    Point q;
    if (!PyArg_ParseTuple(args, "ffff", &(p.x), &(p.y), &(q.x), &(q.y)))
        return NULL;
    geom2d_point_add(&p, q);
    return Py_BuildValue("ff", p.x, p.y);
}

//-------------------------------------------------------------------
// Method table
//-------------------------------------------------------------------
static PyMethodDef Geom2dMethods[] = {
    {"add", geom2d_add, METH_VARARGS, "Add two points."},
    {NULL, NULL, 0, NULL} // sentinel
};

//-------------------------------------------------------------------
// Initialization function
//-------------------------------------------------------------------
PyMODINIT_FUNC initgeom2d(void)
{
    (void) Py_InitModule("geom2d", Geom2dMethods);
}

