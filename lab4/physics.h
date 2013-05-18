#ifndef _physics_h
#define _physics_h

#include "coordinate.h"

#define STEP_SIZE 1.0 /* the step size use in the integration */

int feuler(pcord_t *a, double time_step);

double wall_collide(pcord_t *p, cord_t wall);


double collide(pcord_t *p1, pcord_t *p2);

void interact(pcord_t *p1, pcord_t *p2, double t);


#endif
