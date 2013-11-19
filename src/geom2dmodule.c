#include <Python.h>
#include "geom2d.h"

static PyObject* geom2d_add(PyObject *self, PyObject *args)
{
    Point *p;
    Point q;
    if (PyArg_ParseTuple(args, "ffff", &(p->x), &(p->y), &(q.x), &(q.y))) {
        geom2d_point_add(p, q);
    }
    return Py_BuildValue("ff", &(p->x), &(p->y));
}

static PyMethodDef Geom2dMethods[] = {
    {"add", geom2d_add, METH_VARARGS, "Add two points."},
    {NULL, NULL, 0, NULL} // sentinel
};
