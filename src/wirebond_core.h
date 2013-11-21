#ifndef __WIREBOND_CORE_H
#define __WIREBOND_CORE_H

#include <geom2d.h>


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
void move(Bond *b, Point v);

void add_force(Bond *b, Point f);
void apply_force(Bond *b);

void snap_rectangle(Bond *b, float x, float y, float radius);


//====================================================================
// pair of bond wires
//====================================================================
typedef struct {
    Bond* bonds[2];
    float min_dist[3]; // {pboard-pboard, pboard-wire, pchip-wire}
} BondPair;

Point dist_pboard_pboard(const BondPair *P, int order);
Point dist_pboard_wire(const BondPair *P, int order);
Point dist_pchip_wire(const BondPair *P, int order);

void repulsion_pboard_pboard(BondPair *P, float damp);
void repulsion_pboard_wire(BondPair *P, float damp);
void repulsion_pchip_wire(BondPair *P, float damp);


#endif

