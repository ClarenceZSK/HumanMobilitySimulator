#ifndef COMMON_H
#define COMMON_H

#define PI 3.1415926
#define INVALID_POINT (Point(IMPOSSIBLE_VALUE, IMPOSSIBLE_VALUE))
#define BIG_DISTANCE (double)0xffffffff
#define IMPOSSIBLE_VALUE double(0xffffff)
#define BIG_VALUE (double)0xffffffff

#define myabs(a) ((a) > 0 ? (a) : (-(a)))
#define CEILING(a) ((int)(a) == (a) ? (a) : ((a)+1))    //only for positive value!
#define EQUAL_DOUBLE_ERROR double(0.0001)
#define DOUBLE_EQUAL(a, b) (myabs((a)-(b)) < EQUAL_DOUBLE_ERROR)
#define DOUBLE_EQUAL_PT(a, b) (myabs((a.x)-(b.x)) < EQUAL_DOUBLE_ERROR && myabs((a.y)-(b.y)) < EQUAL_DOUBLE_ERROR)
#define TOPOFILE_NAME "FieldTopology"

#define mymax(a,b)    (((a) > (b)) ? (a) : (b))
#define mymin(a,b)    (((a) < (b)) ? (a) : (b))

#define PS_ONLINE 0
#define PS_TOLEFT 1
#define PS_TORIGHT 2

extern double COMM_RANGE;
extern int NUM_SENSORS;
extern int QUASI_UDG;
extern double QUASI_UDG_ALPHA;
extern double GRID_SIZE;

extern int TIME_UNIT;
extern int TIME_UNIT_COUNTER;
extern int UPDATE_TIME_UNIT;        //update human move for every UPDATE_TIME_UNIT time units

//extern double TIME_COUNTER;
extern double RHO_0;
#endif