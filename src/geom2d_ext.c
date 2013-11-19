#include <Python.h>
#include "geom2d.h"


//===================================================================
// Method definitions
//===================================================================

//--------------------------------------------------------------------
// Point
//--------------------------------------------------------------------

// functions that manipulate a Point
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

static PyObject* geom2d_neg(PyObject *self, PyObject *args)
{
    Point p;
    if (!PyArg_ParseTuple(args, "ff", &(p.x), &(p.y)))
        return NULL;
    geom2d_point_neg(&p);
    return Py_BuildValue("ff", p.x, p.y);
}

static PyObject* geom2d_rotate(PyObject *self, PyObject *args)
{
    Point p;
    float angle;
    if (!PyArg_ParseTuple(args, "fff", &(p.x), &(p.y), &angle))
        return NULL;
    geom2d_point_rotate(&p, angle);
    return Py_BuildValue("ff", p.x, p.y);
}

static PyObject* geom2d_normalize(PyObject *self, PyObject *args)
{
    Point p;
    if (!PyArg_ParseTuple(args, "ff", &(p.x), &(p.y)))
        return NULL;
    if (geom2d_point_normalize(&p)) {
        PyErr_SetNone(PyExc_ZeroDivisionError);
        return NULL;
    }
    return Py_BuildValue("ff", p.x, p.y);
}

// functions that calculate a value from a Point
static PyObject* geom2d_quadrant(PyObject *self, PyObject *args)
{
    Point p;
    int quadrant;
    if (!PyArg_ParseTuple(args, "ff", &(p.x), &(p.y)))
        return NULL;
    quadrant = geom2d_point_quadrant(p);
    return Py_BuildValue("i", quadrant);
}

static PyObject* geom2d_abs(PyObject *self, PyObject *args)
{
    Point p;
    float value;
    if (!PyArg_ParseTuple(args, "ff", &(p.x), &(p.y)))
        return NULL;
    value = geom2d_point_abs(p);
    return Py_BuildValue("f", value);
}

static PyObject* geom2d_angle(PyObject *self, PyObject *args)
{
    Point p;
    float angle;
    if (!PyArg_ParseTuple(args, "ff", &(p.x), &(p.y)))
        return NULL;
    angle = geom2d_point_angle(p);
    return Py_BuildValue("f", angle);
}

static PyObject* geom2d_dot(PyObject *self, PyObject *args)
{
    Point p;
    Point q;
    float dot;
    if (!PyArg_ParseTuple(args, "ffff", &(p.x), &(p.y), &(q.x), &(q.y)))
        return NULL;
    dot = geom2d_point_dot(p, q);
    return Py_BuildValue("f", dot);
}

//--------------------------------------------------------------------
// Line
//--------------------------------------------------------------------

static PyObject* geom2d_dist_point_line_(PyObject *self, PyObject *args)
{
    Point p;
    Line L;
    if (!PyArg_ParseTuple(args, "ffffff", &(p.x), &(p.y),
                          &(L.start.x), &(L.start.y), &(L.end.x), &(L.end.y)))
        return NULL;
    Point dist;
    float pos;
    if (geom2d_dist_point_line(p, L, &dist, &pos)) {
        PyErr_SetString(PyExc_ZeroDivisionError, "No valid line.");
        return NULL;
    }
    return Py_BuildValue("fff", dist.x, dist.y, pos);
}

//===================================================================
// Method table
//===================================================================
static PyMethodDef Geom2dMethods[] = {
    {"add",       geom2d_add, METH_VARARGS, "Add two points."},
    {"sub",       geom2d_sub, METH_VARARGS, "Subtract a point from another."},
    {"mul",       geom2d_mul, METH_VARARGS, "Multiply a point with a scalar."},
    {"div",       geom2d_div, METH_VARARGS, "Divide a point by a scalar."},
    {"neg",       geom2d_neg, METH_VARARGS, "Negate coordinates."},
    {"rotate",    geom2d_rotate, METH_VARARGS, "Rotate by an angle."},
    {"normalize", geom2d_normalize, METH_VARARGS, "Divide by length."},
    {"quadrant",  geom2d_quadrant, METH_VARARGS, "Calculate quadrant."},
    {"abs",       geom2d_abs, METH_VARARGS, "Calculate length."},
    {"angle",     geom2d_angle, METH_VARARGS, "Calculate angle."},
    {"dot",       geom2d_dot, METH_VARARGS, "Calculate dot product."},
    {"dist_point_line", geom2d_dist_point_line_, METH_VARARGS, "Calculate point-line distance."},
    {NULL, NULL, 0, NULL} // sentinel
};

//===================================================================
// Initialization function
//===================================================================
PyMODINIT_FUNC initgeom2d(void)
{
    (void) Py_InitModule("geom2d", Geom2dMethods);
}

