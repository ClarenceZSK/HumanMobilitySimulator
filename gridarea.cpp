//
//  gridarea.cpp
//  HumanMobilitySimulator
//
//  Created by Shengkai Zhang on 12/3/14.
//  Copyright (c) 2014 Shengkai Zhang. All rights reserved.
//

#include "gridarea.h"

GridArea::GridArea()
{
	index = 0;
    field = NULL;
    gridAreaVisitedTimes = 0;
    humanDensity = 0;
    averageDistanceToCenter = 0;
}

GridArea::~GridArea()
{
}

