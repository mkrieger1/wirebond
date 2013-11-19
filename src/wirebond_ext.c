#include <Python.h>
#include "geom2d.h"

//===================================================================
// Method definitions
//===================================================================

static PyObject* repulsion_pboard(PyObject *self, PyObject *args)
{
    Point p; // bond.pboard
    Point q; // otherbond.pboard
    float min; // min_dist_pboard
    if (!PyArg_ParseTuple(args, "fffff", &(p.x), &(p.y), &(q.x), &(q.y), &min))
        return NULL;

    Point dist = p;
    geom2d_point_sub(&dist, q);

    Point force;
    float violation = min - geom2d_point_abs(dist);
    if (violation > 0) {
        force = dist;
        geom2d_point_normalize(&force);
        geom2d_point_mul(&force, violation);
    }
    return Py_BuildValue("fff", violation, force.x, force.y);
}


//===================================================================
// Method table
//===================================================================
static PyMethodDef WirebondExtMethods[] = {
    {"repulsion_pboard", repulsion_pboard, METH_VARARGS,
                         "Helper function for repulsion_pboard."},
    {NULL, NULL, 0, NULL} // sentinel
};

//===================================================================
// Initialization function
//===================================================================
PyMODINIT_FUNC initwirebond_ext(void)
{
    (void) Py_InitModule("wirebond_ext", WirebondExtMethods);
}

