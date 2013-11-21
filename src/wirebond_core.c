#include "wirebond_core.h"

//====================================================================
// representation of a bond wire as two connected Points
//====================================================================
void calc_pboard(Bond *b)
{
    Vector wire = {1.0, 0.0};
    rotate_pt(&wire, b->angle);
    mul(&wire, b->length);
    b->pboard = b->pchip;
    add(&(b->pboard), wire);
}

void calc_length_angle(Bond *b)
{
    Vector wire = b->pboard;
    sub(&wire, b->pchip);
    b->length = abs_pt(wire);
    b->angle = angle(wire);
}

void set_pboard(Bond *b, Point pboard)
{
    b->pboard = pboard;
    calc_length_angle(b);
}

void set_length(Bond *b, float length)
{
    b->length = length;
    calc_pboard(b);
}

void set_angle(Bond *b, float angle)
{
    b->angle = angle;
    calc_pboard(b);
}

void stretch(Bond *b, float length)
{
    set_length(b, b->length + length);
}

void rotate_bd(Bond *b, float angle)
{
    set_angle(b, b->angle + angle);
}

void move(Bond *b, Vector v)
{
    Point pboard_new = b->pboard;
    add(&pboard_new, v);
    set_pboard(b, pboard_new);
}

void add_force(Bond *b, Force f)
{
    add(&(b->force), f);
}

void apply_force(Bond *b)
{
    Vector wire = b->pboard;
    sub(&wire, b->pchip);
    add(&wire, b->force);
    set_angle(b, angle(wire));
    b->force = (Force) {0.0, 0.0};
}

void snap_rectangle(Bond *b, float x, float y, float radius)
{
    int exists_x, exists_y, exists_c;
    Vector bond = b->pboard;
    sub(&bond, b->pchip);
    int q = quadrant(bond);
    int sgnx = ((q == 0 || q == 1)) ? 1 : -1;
    int sgny = ((q == 0 || q == 3)) ? 1 : -1;
    Line wire = (Line) {b->pboard, b->pchip};
    pboard_x = intersect_line_x(wire, sgnx*x, &exists_x);
    pboard_y = intersect_line_y(wire, sgny*y, &exists_y);
    Point center = (Point) {sgnx*x, sgny*y};
    pboard_c = intersect_line_circle(wire, center, radius, q, &exists_c);
    /// ....
}


//====================================================================
// pair of bond wires
//====================================================================
Vector dist_pboard_pboard(const BondPair *P, int order)
{
    // TODO
}

Vector dist_pboard_wire(const BondPair *P, int order)
{
    // TODO
}

Vector dist_pchip_wire(const BondPair *P, int order)
{
    // TODO
}


void repulsion_pboard_pboard(BondPair *P, float damp)
{
    // TODO
}

void repulsion_pboard_wire(BondPair *P, float damp)
{
    // TODO
}

void repulsion_pchip_wire(BondPair *P, float damp)
{
    // TODO
}

