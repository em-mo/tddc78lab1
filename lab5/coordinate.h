#ifndef _coordinate_h
#define _coordinate_h

struct cord {
    double x0 ;
    double x1 ;
    double y0 ;
    double y1 ;
} ;

struct part_cord {
    double x ;
    double y ;
    double vx ;
    double vy ;
} ;

typedef struct cord cord_t ;
typedef struct part_cord pcord_t ;

#endif
