#ifndef __GEOM2D_H
#define __GEOM2D_H


//====================================================================
// Point/vector in 2-D euclidian space
//====================================================================
typedef struct {
    float x;
    float y;
} Point;

// functions that manipulate a Point and return 0 on success
int geom2d_point_add(Point *p, Point q);
int geom2d_point_sub(Point *p, Point q);
int geom2d_point_mul(Point *p, float a);
int geom2d_point_div(Point *p, float a);
int geom2d_point_neg(Point *p);

int geom2d_point_rotate(Point *p, float angle);
int geom2d_point_normalize(Point *p);

// functions that take one or more Points and return a calculated value
int geom2d_point_quadrant(Point p);
float geom2d_point_abs(Point p);
float geom2d_point_angle(Point p);
float geom2d_point_dot(Point p, Point q);


//====================================================================
// Line through two points
//====================================================================
typedef struct {
    Point start;
    Point end;
} Line;

int geom2d_dist_point_line(Point p, Line L, Point *dist, float *pos);

/*
Point intersect_line_x(Line L, float x, int *exists);
Point intersect_line_y(Line L, float y, int *exists);

Point intersect_line_circle(Line L, Point center, float radius,
                            int quadrant, int *exists);
*/

#endif

