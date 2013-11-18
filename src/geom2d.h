#ifndef __GEOM2D_H
#define __GEOM2D_H

#include <math.h>


//====================================================================
// Point/vector in 2-D euclidian space
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

int quadrant(Point p);
float abs_pt(Point p);
float angle(Point p);
float dot(Point p, Point q);

void rotate_pt(Point *p, float angle);
Point rotated_pt(Point p, float angle);

void normalize(Point *p, int *exists);
Point normalized(Point p, int *exists);


//====================================================================
// Line through two points
//====================================================================
typedef struct {
    Point start;
    Point end;
} Line;

Point to_vector(Line L);

Point footpoint(Point p, Line L, int *exists);

Point intersect_line_x(Line L, float x, int *exists);
Point intersect_line_y(Line L, float y, int *exists);

Point intersect_line_circle(Line L, Point center, float radius,
                            int quadrant, int *exists);


#endif

