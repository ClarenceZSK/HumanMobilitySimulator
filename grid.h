//
//  grid.h
//  HumanMobilitySimulator
//
//  Created by Shengkai Zhang on 7/3/14.
//  Copyright (c) 2014 Shengkai Zhang. All rights reserved.
//

#ifndef HumanMobilitySimulator_grid_h
#define HumanMobilitySimulator_grid_h

#include <iostream>
#include <set>
#include <map>
#include <vector>
#include <algorithm>
#include "math.h"
#include "geometry.h"
#include "common.h"
#include "random.h"
#include "maplocation.h"
#include "field.h"

class Field;

class Grid
{
public:
	Grid();
	~Grid();
    //Attributes
	int index;
    double size;
    Point location;
    int rho0;    //parameter rho_0 in the model
    bool hasHuman;
    int humanNum;
    int areaID;
    bool broadcastCovered;
    set<Grid*> gridNeighbors;
    Field *field;
    
    
    //Methods
    void init(int index_, Point location_, int rho_);
    void searchingNeighbors();
	
	inline bool operator == (const Grid &right) const
	{
		if(index == right.index)
			return true;
		else
			return false;
	}
};

#endif
