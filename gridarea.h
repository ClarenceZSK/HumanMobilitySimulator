//
//  gridarea.h
//  HumanMobilitySimulator
//
//  Created by Shengkai Zhang on 12/3/14.
//  Copyright (c) 2014 Shengkai Zhang. All rights reserved.
//

#ifndef HumanMobilitySimulator_gridarea_h
#define HumanMobilitySimulator_gridarea_h

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
#include "grid.h"

class Field;
class Grid;

class GridArea
{
public:
	GridArea();
	~GridArea();
    //Attributes
	int index;
    Field *field;
    vector<Grid*> gridArea;
    double humanDensity;
    double averageDistanceToCenter;
    
    //mobility topology analysis
    int gridAreaVisitedTimes;                           //gridAreaVisitedTimes determins the gridarea size
    map<GridArea*, int> mobilityTopologyNeighbors;          //neighbor grid area pointer, jump times
    map<GridArea*, double> mobilityTransportProb;
    
    
    //Methods
	
	inline bool operator == (const GridArea &right) const
	{
		if(index == right.index)
			return true;
		else
			return false;
	}
};

#endif
