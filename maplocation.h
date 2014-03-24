//
//  maplocation.h
//  HumanMobilitySimulator
//
//  Created by Shengkai Zhang on 1/3/14.
//  Copyright (c) 2014 Shengkai Zhang. All rights reserved.
//

#ifndef HumanMobilitySimulator_maplocation_h
#define HumanMobilitySimulator_maplocation_h

#include <iostream>
#include <set>
#include <map>
#include <vector>
#include <algorithm>
#include "math.h"
#include "geometry.h"
#include "common.h"
#include "random.h"

class Maplocation
{
public:
    Maplocation();
	Maplocation(Point p, int f);
	~Maplocation();
    //Attributes
	Point position;
    int visitFrequency;
    
    inline bool operator == (const Maplocation &right) const
	{
		if(position == right.position)
			return true;
		else
			return false;
	}
    
    inline bool operator != (const Maplocation &right) const
	{
		if(position != right.position)
			return true;
		else
			return false;
	}
    
};

#endif
