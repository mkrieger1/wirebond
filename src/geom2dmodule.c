#include <Python.h>
#include "geom2d.h"


//-------------------------------------------------------------------
// Method definitions
//-------------------------------------------------------------------
static PyObject* geom2d_add(PyObject *self, PyObject *args)
{
    Point p;
    Point q;
    if (!PyArg_ParseTuple(args, "ffff", &(p.x), &(p.y), &(q.x), &(q.y)))
        return NULL;
    geom2d_point_add(&p, q);
    return Py_BuildValue("ff", p.x, p.y);
}

static PyObject* geom2d_sub(PyObject *self, PyObject *args)
{
    Point p;
    Point q;
    if (!PyArg_ParseTuple(args, "ffff", &(p.x), &(p.y), &(q.x), &(q.y)))
        return NULL;
    geom2d_point_sub(&p, q);
    return Py_BuildValue("ff", p.x, p.y);
}

static PyObject* geom2d_mul(PyObject *self, PyObject *args)
{
    Point p;
    float a;
    if (!PyArg_ParseTuple(args, "fff", &(p.x), &(p.y), &a))
        return NULL;
    geom2d_point_mul(&p, a);
    return Py_BuildValue("ff", p.x, p.y);
}

static PyObject* geom2d_div(PyObject *self, PyObject *args)
{
    Point p;
    float a;
    if (!PyArg_ParseTuple(args, "fff", &(p.x), &(p.y), &a))
        return NULL;
    if (geom2d_point_div(&p, a)) {
        PyErr_SetNone(PyExc_ZeroDivisionError);
        return NULL;
    }
    return Py_BuildValue("ff", p.x, p.y);
}

//-------------------------------------------------------------------
// Method table
//-------------------------------------------------------------------
static PyMethodDef Geom2dMethods[] = {
    {"add", geom2d_add, METH_VARARGS, "Add two points."},
    {"sub", geom2d_sub, METH_VARARGS, "Subtract a point from another."},
    {"mul", geom2d_mul, METH_VARARGS, "Multiply a point with a scalar."},
    {"div", geom2d_div, METH_VARARGS, "Divide a point by a scalar."},
    {NULL, NULL, 0, NULL} // sentinel
};

//-------------------------------------------------------------------
// Initialization function
//-------------------------------------------------------------------
PyMODINIT_FUNC initgeom2d(void)
{
    (void) Py_InitModule("geom2d", Geom2dMethods);
}

