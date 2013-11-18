#ifndef __GEOM2D_H
#define __GEOM2D_H

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

int quadrant(Point p);
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

Point intersect_line_circle(Line L, Point center, float radius,
                            int quadrant, int *exists);


#endif

