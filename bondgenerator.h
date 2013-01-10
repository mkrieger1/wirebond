#ifndef __LIB_BONDGENERATOR_H
#define __LIB_BONDGENERATOR_H

#include <math.h>


//--------------------------------------------------------------------
// point in 2-D euclidian space
//--------------------------------------------------------------------
typedef struct {
    float x;
    float y;
} Point;

static void add(Point *p, Point q);
static void sub(Point *p, Point q);
static void mul(Point *p, float a);
static void div(Point *p, float a);
static void neg(Point *p);

static float abs(Point p);
static float angle(Point p);
static float dot(Point p, Point q);

static void rotate(Point *p, float phi);
static Point rotated(Point p, float phi);

static void normalize(Point *p);
static Point normalized(Point p);


// forces are also represented by 2-D points
typedef Point Force;

// 'vectors' are also represented by 2-D points
typedef Point Vector;


//--------------------------------------------------------------------
// straight line between two Points
//--------------------------------------------------------------------
typedef struct {
    Point start;
    Point end;
} Line;

static Point footpoint(Point p, Line L);

static Point intersect_line_x(Line L, float x);
static Point intersect_line_y(Line L, float y);

static Point intersect_line_circle(Line L, float radius, int quadrant);


//--------------------------------------------------------------------
// representation of a bond wire as two connected Points
//--------------------------------------------------------------------
typedef struct {
    Point pchip;  // position on the chip (usually fixed)
    Point pboard; // position on the PCB (movable)

    float length; // distance pchip->pboard in µm
    float angle;  // polar angle of the line pchip->pboard

    Point force;  // total force acting on pboard
} Bond;

static void calc_pboard(Bond *b);
static void calc_length_angle(Bond *b);

static void set_pboard(Bond *b, Point pboard);
static void set_length(Bond *b, float length);
static void set_angle(Bond *b, float angle);

static void stretch(Bond *b, float length);
static void rotate(Bond *b, float angle);
static void move(Bond *b, Vector v);

static void add_force(Bond *b, Force f);
static void apply_force(Bond *b);

static void snap_rectangle(Bond *b, float x, float y, float radius);


//--------------------------------------------------------------------
// pair of bond wires
//--------------------------------------------------------------------
typedef struct {
    Bond* bonds[2];
    float min_dist[3]; // {pboard-pboard, pboard-wire, pchip-wire}
} BondPair;

static Vector dist_pboard_pboard(const BondPair *P, int order);
static Vector dist_pboard_wire(const BondPair *P, int order);
static Vector dist_pchip_wire(const BondPair *P, int order);

static void repulsion_pboard_pboard(BondPair *P, float damp);
static void repulsion_pboard_wire(BondPair *P, float damp);
static void repulsion_pchip_wire(BondPair *P, float damp);


#endif

