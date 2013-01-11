#include "bondgenerator.h"


//====================================================================
// point in 2-D euclidian space
//====================================================================

// add q to p (in-place)
void add(Point *p, Point q)
{
    p->x += q.x; p->y += q.y;
}

// subtract q from p (in-place)
void sub(Point *p, Point q)
{
    p->x -= q.x; p->y -= q.y;
}

// multiply p by a (in-place)
void mul(Point *p, float a)
{
    p->x *= a; p->y *= a;
}

// divide p by a (in-place)
void div(Point *p, float a, int *exists)
{
    if (a != 0.0)
    {
        *exists = 0;
    }
    else
    {
        *exists = 1;
        p->x /= a; p->y /= a;
    }
}

// negate p (in-place)
void neg(Point *p)
{
    mul(p, -1.0);
}

// calculate L2 norm
float abs_pt(Point p)
{
    return sqrt(p.x*p.x + p.y*p.y);
}

// calculate polar angle
float angle(Point p)
{
    return atan2(p.y, p.x);
}

// calculate dot-product of p and q
float dot(Point p, Point q)
{
    return p.x*q.x + p.y*q.y;
}

// rotate p (in-place)
void rotate_pt(Point *p, float angle)
{
    float x = cos(angle)*p->x - sin(angle)*p->y;
    float y = sin(angle)*p->x + cos(angle)*p->y;
    p->x = x;
    p->y = y;
}

// return rotated copy
Point rotated_pt(Point p, float angle)
{
    rotate_pt(&p, angle);
    return p;
}

// normalize p (in-place)
void normalize(Point *p, int *exists)
{
    div(p, abs_pt(*p), exists);
}

// return normalized copy
Point normalized(Point p, int *exists)
{
    normalize(&p, exists);
    return p;
}

//====================================================================
// straight line between two Points
//====================================================================

// return vector start->end
Vector to_vector(Line L)
{
    sub(&(L.end), L.start);
    return L.end;
}

// return footpoint of p on L
Point footpoint(Point p, Line L, int *exists)
{
    sub(&p, L.start);
    Vector dir = to_vector(L);
    normalize(&dir, exists);
    float t = dot(p, dir);
    mul(&dir, t);
    add(&(L.start), dir);
    return L.start;
}

// return point on L with given x
Point intersect_line_x(Line L, float x, int *exists)
{
    if (L.start.x == L.end.x) // parallel
    {
        *exists = 0;
        return L.end; // we have to return something
    }
    else
    {
        *exists = 1;
        Point p = L.start;
        Vector dir = to_vector(L);
        float a = (x-p.x)/dir.x;
        mul(&dir, a);
        add(&p, dir);
        return p;
    }
}

// return point on L with given y
Point intersect_line_y(Line L, float y, int *exists)
{
    if (L.start.y == L.end.y) // parallel
    {
        *exists = 0;
        return L.end; // we have to return something
    }
    else
    {
        *exists = 1;
        Point p = L.start;
        Vector dir = to_vector(L);
        float a = (y-p.y)/dir.y;
        mul(&dir, a);
        add(&p, dir);
        return p;
    }
}

// return intersection point of L with circle around origin
// source: http://mathworld.wolfram.com/Circle-LineIntersection.html
Point intersect_line_circle(Line L, float radius, int quadrant,
                            int *exists)
{
    Point p1 = L.start;
    Point p2 = L.end;
    float length = abs_pt(to_vector(L));
    float D = p1.x*p2.y - p1.y*p2.x;
    float DD = (radius*length)*(radius*length) - D*D;
    if (DD < 0)
    {
        *exists = 0;
        return L.end; // we have to return something
    }  
    else
    {
        Point q1, q2;

        float dx = p2.x-p1.x;
        float dy = p2.y-p1.y;
        int sgn = (dy < 0) ? -1 : 1;
        float sDD = sqrt(DD);
        float lsq = length*length;

        q1.x = ( D*dy + sgn*dx*sDD) / lsq;
        q1.y = (-D*dx + abs(dx)*sDD)/ lsq;
        q2.x = ( D*dy - sgn*dx*sDD) / lsq;
        q2.y = (-D*dx - abs(dx)*sDD)/ lsq;

        // quadrant matching: 0 = NE, 1 = NW, 2 = SW, 3 = SE
        int qsgnx = (quadrant == 0 || quadrant == 3);
        int qsgny = (quadrant == 0 || quadrant == 1);
        int sgnq1x = (q1.x >= 0);
        int sgnq1y = (q1.y >= 0);
        int sgnq2x = (q2.x >= 0);
        int sgnq2y = (q2.y >= 0);

        if ((qsgnx == sgnq1x) && (qsgny == sgnq1y))
        {
            *exists = 1;
            return q1;
        }
        else if ((qsgnx == sgnq2x) && (qsgny == sgnq2y))
        {
            *exists = 1;
            return q2;
        }
        else
        {
            *exists = 0;
            return L.end; // we have to return something
        }
    }
}








