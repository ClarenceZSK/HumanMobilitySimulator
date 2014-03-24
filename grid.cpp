//
//  grid.cpp
//  HumanMobilitySimulator
//
//  Created by Shengkai Zhang on 7/3/14.
//  Copyright (c) 2014 Shengkai Zhang. All rights reserved.
//

#include "grid.h"

Grid::Grid()
{
	index = 0;
    field = NULL;
    size = 1;
    rho0 = 0;
    hasHuman = false;
    humanNum = 0;
    areaID = -1;
    broadcastCovered = false;
	gridNeighbors.clear();
	
}

Grid::~Grid()
{
}

void Grid::init(int index_, Point location_, int rho0_)
{
	index = index_;
	location = location_;
    rho0 = rho0_;
}

void Grid::searchingNeighbors()
{
    //only adjacent grids are neighbors
    int maxX = field->maxX;
    int id = index;
    if (id+1 < maxX) {
        gridNeighbors.insert(&field->gridPool[id+1]);
        if (id+1+maxX < field->nGrid)
            gridNeighbors.insert(&field->gridPool[id+1+maxX]);
        if (id+1-maxX >= 0)
            gridNeighbors.insert(&field->gridPool[id+1-maxX]);
    }
    if (id-1 >= 0) {
        gridNeighbors.insert(&field->gridPool[id-1]);
        if (id-1+maxX < field->nGrid)
            gridNeighbors.insert(&field->gridPool[id-1+maxX]);
        if (id-1-maxX >= 0)
            gridNeighbors.insert(&field->gridPool[id-1-maxX]);
    }
    if (id+maxX < field->nGrid) {
        gridNeighbors.insert(&field->gridPool[id+maxX]);
    }
    if (id-maxX >= 0) {
        gridNeighbors.insert(&field->gridPool[id+1]);
    }
}
