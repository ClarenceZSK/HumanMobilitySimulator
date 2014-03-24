#ifndef SENSORFIELD_H
#define SENSORFIELD_H

#pragma warning(disable: 4786)

#include <iostream>
#include <set>
#include <map>
#include <vector>
#include <algorithm>
#include "math.h"
#include "geometry.h"
#include "common.h"
#include "random.h"
#include "assert.h"
#include "human.h"
#include "grid.h"
#include "gridarea.h"

using namespace std;

class Human;
class GridArea;

class Field
{
public:
	Field(int maxX_, int minX_, int maxY_, int minY_, 
		int nHuman_);
	~Field();
    
    //Attibutes
    int maxX;
    int minX;
    int maxY;
    int minY;
    int nHuman;
    int nGrid;
    int nGridAreas;
    int moves;
    Human *humanPool;
    Grid *gridPool;
    GridArea *gridAreaPool;
    
    multiset<double, greater<double> > deltaR;
    
    //Methods
	void InitDeployment(double gridSize);
    int MergeGirdAreas();               //return total number of areas
    void EraseVectorElement(vector<Human*> &dl_vector, Human *s);
    void ComputeAverageDistance();
    void ComputeHumanDensity();
    
    void DumpData(char *filename1, char *filename2);
    void DumpCCDF(int x, char *filename1, char *filename2);
    void DumpGrid(char *filename);
    void DumpMobilityTopologyAnalysis(char *filename1, char *filename2, char *filename3, char *filename4, char *filename5);
    void DumpDistanceMap(char *filename1);
    
    void DumpRG(char *filename1, char *filename2);
    
    void DumpDeltaR(char *filename);
	
};


#endif