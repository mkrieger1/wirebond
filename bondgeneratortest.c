#include "bondgenerator.h"
#include <stdio.h>

void print_pt(Point p)
{
    printf("x = %.2f  y = %.2f\n", p.x, p.y);
}


int main()
{
    int exists = 0;

    printf("Testing bondgenerator C library.\n");

    Point p;
    p.x = 3.0;
    p.y = 2.0;

    Point q;
    q.x = -1.0;
    q.y =  4.0;

    printf("Point p: ");
    print_pt(p);

    printf("Point q: ");
    print_pt(q);

    printf("add q: ");
    add(&p, q);
    print_pt(p);

    printf("sub q: ");
    sub(&p, q);
    print_pt(p);

    printf("mul 3: ");
    mul(&p, 3.0);
    print_pt(p);

    printf("div 2: ");
    div(&p, 2.0, &exists);
    print_pt(p);

    printf("neg: ");
    neg(&p);
    print_pt(p);

    printf("abs: %.2f\n", abs_pt(p));
    printf("angle: %.2f\n", angle(p)*180.0/M_PI);
    printf("dot p q: %.2f\n", dot(p, q));

    printf("rotate 90°: ");
    rotate_pt(&p, M_PI/2);
    print_pt(p);

    printf("q2 = rotated q\n");
    Point q2 = rotated_pt(q, M_PI/2);
    print_pt(q2);
    print_pt(q);

    printf("normalize p: ");
    normalize(&p, &exists);
    print_pt(p);

    printf("r = normalized q\n");
    Point r = normalized(q, &exists);
    print_pt(r);
    print_pt(q);

    printf("footpoint:\n");
    Point a, b;
    a.x = 3.0; a.y = 4.0;
    b.x = 5.0; b.y = 6.0;
    Line L;
    L.start = a;
    L.end = b;
    Point c;
    c.x = 4.5;
    c.y = 4.5;
    Point F = footpoint(c, L, &exists);
    print_pt(F);

    return 0;
}

