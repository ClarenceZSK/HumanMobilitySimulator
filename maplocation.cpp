//
//  maplocation.cpp
//  HumanMobilitySimulator
//
//  Created by Shengkai Zhang on 1/3/14.
//  Copyright (c) 2014 Shengkai Zhang. All rights reserved.
//

#include "maplocation.h"

Maplocation::Maplocation()
{
    position.x = 0;
    position.y = 0;
	visitFrequency = 1;
}

Maplocation::Maplocation(Point p, int f)
{
    position = p;
	visitFrequency = f;
}

Maplocation::~Maplocation()
{
}