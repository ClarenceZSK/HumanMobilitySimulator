//
//  human.h
//  HumanMobilitySimulator
//
//  Created by Shengkai Zhang on 28/2/14.
//  Copyright (c) 2014 Shengkai Zhang. All rights reserved.
//

#ifndef HumanMobilitySimulator_human_h
#define HumanMobilitySimulator_human_h

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
class Maplocation;
class Grid;

class Human
{
public:
	Human();
	~Human();
    //Attributes
	int index;
    vector<pair<double, int> > waitingTimeSequence;   //time, visitedLocations.size()
	Maplocation currentMapPosition;
    vector<Maplocation> visitedLocations;           //this vector becomes redundant, please ignore all references
    multimap<int, Maplocation> userVisitedRecord;   //visited times, maplocation
    multimap<double, Point, greater<double> > normalizedFrequencyLocations;
    int gridNum;
    //move control
    bool readyToMove;
    bool waitNow;
    
    double powerLawT;   //waiting time
    double timeCounter;
    double timePassed;
    
    UniformRandVar uRandR;          //for powerLawR
    UniformRandVar uRandT;          //for powerLawT
    UniformRandVar uRandM;          //for moveP
    UniformRandVar uRandV;          //for revisitP
    UniformRandVar uRandA;          //for angle
    
    //mobility topology analysis
    vector<pair<Grid*, Grid*> > jumpGridRecord;     //record jump steps
    
    Field *field;
    Grid *grid;
    
    //for output St and fk
    int startIndex;
    vector<double> visitFrequency;
    
    //for output rg
    vector<pair<double, double> > rg_t; //time passed, rg value
    int maxRGT;                         //for normalizing
    void computeRGT();
    
    //Methods
    void init(int index_, Point position_);
    
    //void moveModel1(double randr, double randt, double angle, double moveP, double revisitP);
    //void moveModel2(double randr, double randt, double angle, double moveP, double revisitP);
    //void moveModel2SingleTrack(double randr, double randt, double angle, double moveP, double revisitP);
    
    void moveModel1();
    void moveModel2();
    void moveModel2SingleTrack();
    
    void normalizeLocations();
    int inGrid(double gridSize);       //return the grid number in which the human is
    void updateGridInfo();
    
    /*
    bool mapLocationSort (Maplocation ml1,Maplocation ml2) {
        return (ml1.visitFrequency > ml2.visitFrequency);
    }
    */
	
	inline bool operator == (const Human &right) const
	{
		if(index == right.index)
			return true;
		else
			return false;
	}
};

#endif
