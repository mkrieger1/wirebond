#include "bondgenerator.h"


//--------------------------------------------------------------------
// point in 2-D euclidian space
//--------------------------------------------------------------------

void add(Point *p, Point q)
{
    p->x += q.x;
    p->y += q.y;
}

void sub(Point *p, Point q)
{
    p->x -= q.x;
    p->y -= q.y;
}

void mul(Point *p, float a)
{
    p->x *= a;
    p->y *= a;
}

void div(Point *p, float a)
{
    if (a != 0.0)
    {
        p->x /= a;
        p->y /= a;
    }
}

void neg(Point *p)
{
    mul(p, -1.0);
}

float abs_pt(Point p)
{
    return sqrt(p.x*p.x + p.y*p.y);
}

float angle(Point p)
{
    return atan2(p.y, p.x);
}

float dot(Point p, Point q)
{
    return p.x*q.x + p.y*q.y;
}

void rotate_pt(Point *p, float angle)
{
    float x = cos(angle)*p->x - sin(angle)*p->y;
    float y = sin(angle)*p->x + cos(angle)*p->y;
    p->x = x;
    p->y = y;
}

Point rotated_pt(Point p, float angle)
{
    rotate_pt(&p, angle);
    return p;
}

void normalize_pt(Point *p)
{
    div(p, abs_pt(*p));
}

Point normalized_pt(Point p)
{
    normalize_pt(&p);
    return p;
}


//--------------------------------------------------------------------
// straight line between two Points
//--------------------------------------------------------------------
Vector to_vector(Line L)
{
    sub(&(L.end), L.start);
    return L.end;
}

Point footpoint(Point p, Line L)
{
    sub(&p, L.start);
    Vector dir = to_vector(L);
    normalize_pt(&dir);
    float t = dot(p, dir);
    mul(&dir, t);
    add(&(L.start), dir);
    return L.start;
}

