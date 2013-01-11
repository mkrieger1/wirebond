#ifndef __LIB_BONDGENERATOR_H
#define __LIB_BONDGENERATOR_H

#include <math.h>


//====================================================================
// point in 2-D euclidian space
//====================================================================
typedef struct {
    float x;
    float y;
} Point;

void add(Point *p, Point q);
void sub(Point *p, Point q);
void mul(Point *p, float a);
void div(Point *p, float a, int *exists);
void neg(Point *p);

float abs_pt(Point p);
float angle(Point p);
float dot(Point p, Point q);

void rotate_pt(Point *p, float angle);
Point rotated_pt(Point p, float angle);

void normalize(Point *p, int *exists);
Point normalized(Point p, int *exists);


// forces are also represented by 2-D points
typedef Point Force;

// 'vectors' are also represented by 2-D points
typedef Point Vector;


//====================================================================
// straight line between two Points
//====================================================================
typedef struct {
    Point start;
    Point end;
} Line;

Vector to_vector(Line L);

Point footpoint(Point p, Line L, int *exists);

Point intersect_line_x(Line L, float x, int *exists);
Point intersect_line_y(Line L, float y, int *exists);

Point intersect_line_circle(Line L, float radius, int quadrant,
                            int *exists);


//====================================================================
// representation of a bond wire as two connected Points
//====================================================================
typedef struct {
    Point pchip;  // position on the chip (usually fixed)
    Point pboard; // position on the PCB (movable)

    float length; // distance pchip->pboard in µm
    float angle;  // polar angle of the line pchip->pboard

    Point force;  // total force acting on pboard
} Bond;

void calc_pboard(Bond *b);
void calc_length_angle(Bond *b);

void set_pboard(Bond *b, Point pboard);
void set_length(Bond *b, float length);
void set_angle(Bond *b, float angle);

void stretch(Bond *b, float length);
void rotate_bd(Bond *b, float angle);
void move(Bond *b, Vector v);

void add_force(Bond *b, Force f);
void apply_force(Bond *b);

void snap_rectangle(Bond *b, float x, float y, float radius);


//====================================================================
// pair of bond wires
//====================================================================
typedef struct {
    Bond* bonds[2];
    float min_dist[3]; // {pboard-pboard, pboard-wire, pchip-wire}
} BondPair;

Vector dist_pboard_pboard(const BondPair *P, int order);
Vector dist_pboard_wire(const BondPair *P, int order);
Vector dist_pchip_wire(const BondPair *P, int order);

void repulsion_pboard_pboard(BondPair *P, float damp);
void repulsion_pboard_wire(BondPair *P, float damp);
void repulsion_pchip_wire(BondPair *P, float damp);


#endif

