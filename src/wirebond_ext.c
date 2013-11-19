#include <Python.h>
#include "geom2d.h"

//===================================================================
// Method definitions
//===================================================================


//-------------------------------------------------------------------
// repulsion point <-> point
//-------------------------------------------------------------------
static PyObject* repulsion_point_point(PyObject *self, PyObject *args)
{
    Point p;
    Point q;
    float min;
    if (!PyArg_ParseTuple(args, "fffff", &(p.x), &(p.y), &(q.x), &(q.y), &min))
        return NULL;

    Point dist = p;
    geom2d_point_sub(&dist, q);

    Point force;
    float violation = min - geom2d_point_abs(dist);
    if (violation > 0) {
        force = dist;
        if (geom2d_point_normalize(&force)) {
            PyErr_SetNone(PyExc_ZeroDivisionError);
            return NULL;
        }
        geom2d_point_mul(&force, violation);
    }
    return Py_BuildValue("fff", violation, force.x, force.y);
}


//-------------------------------------------------------------------
// repulsion point <-> line
//-------------------------------------------------------------------
static PyObject* repulsion_point_line(PyObject *self, PyObject *args)
{
    Point p;
    Line L;
    float min;
    if (!PyArg_ParseTuple(args, "fffffff", &(p.x), &(p.y),
                                           &(L.start.x), &(L.start.y),
                                           &(L.end.x), &(L.end.y), &min))
        return NULL;

    Point force;
    float violation = 0.0;
    Point dist;
    float pos;
    if (geom2d_dist_point_line(p, L, &dist, &pos)) {
        PyErr_SetNone(PyExc_ZeroDivisionError);
        return NULL;
    }
    if (pos > 0 && pos < 1) {
        violation = min - geom2d_point_abs(dist);
        if (violation > 0) {
            force = dist;
            if (geom2d_point_normalize(&force)) {
                PyErr_SetNone(PyExc_ZeroDivisionError);
                return NULL;
            }
            geom2d_point_mul(&force, violation);
        } // violation > 0
    } // pos > 0 && pos < 1
    return Py_BuildValue("fff", violation, force.x, force.y);
}


//===================================================================
// Method table
//===================================================================
static PyMethodDef WirebondExtMethods[] = {
    {"repulsion_point_point", repulsion_point_point, METH_VARARGS,
     "Helper function for repulsion_pboard."},
    {"repulsion_point_line", repulsion_point_line, METH_VARARGS,
     "Helper function for repulsion_pboard_wire and repulsion_pchip_wire."},
    {NULL, NULL, 0, NULL} // sentinel
};


//===================================================================
// Initialization function
//===================================================================
PyMODINIT_FUNC initwirebond_ext(void)
{
    (void) Py_InitModule("wirebond_ext", WirebondExtMethods);
}

