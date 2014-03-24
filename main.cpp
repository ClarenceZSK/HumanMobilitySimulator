//
//  main.cpp
//  HumanMobilitySimulator
//
//  Created by Shengkai Zhang on 27/2/14.
//  Copyright (c) 2014 Shengkai Zhang. All rights reserved.
//

#include <stdio.h>
#include "common.h"
#include "random.h"
#include "geometry.h"
#include "field.h"
#include "human.h"

using namespace std;

Field *globalField = NULL;
double GRID_SIZE = 1;
double RHO_0 = 1;
//double TIME_COUNTER = 0;
//int TIME_UNIT = 10000;
int TIME_UNIT = 360000;      //100 hours = 360000 seconds
int TIME_UNIT_COUNTER = 0;
int UPDATE_TIME_UNIT = 60;


int main(int argc, char *argv [])
{
    bool singleHumanTrack = true;   //tracking single human's mobility pattern, outputing St, fk
	
    //UniformRandVar uniRand(0, 1);
    //UniformRandVar uniRandAngle(0, 360);
    //double randr = 0;
    //double randt = 0;
    //double angle = 0;
    //double moveP = 0;
    //double revisitP = 0;
    int moveSteps = 1000;
   
    globalField = new Field(1000, 0, 1000, 0, 1000); //maxX(maxX_), minX(minX_), maxY(maxY_), minY(minY_), nHuman(nHuman_)
    globalField->moves = moveSteps;
    globalField->InitDeployment(GRID_SIZE);  //init virtual locations, the parameter is the grid size
    
    bool needDrawing = false;
    
    if (singleHumanTrack) {
        //int humanIndex = 0;
        for (int i = 0; i < globalField->nHuman; ++i) {
            for (int s = 0; s < moveSteps; ++s) {
                Human *h = &globalField->humanPool[i];
                h->moveModel2SingleTrack();                //(randr, randt, angle, moveP, revisitP);
                globalField->humanPool[i].normalizeLocations();
            }
            cout << "Processing human " << i << endl;
        }
        
        globalField->DumpData("st.txt", "fk.txt");
        
        for (int i = 0; i < globalField->nHuman; ++i) {
            Human *h = &globalField->humanPool[i];
            h->startIndex = 1;
            h->computeRGT();
        }
        
        globalField->DumpRG("rgt.txt", "Prg.txt");
    }
    
    else if (!needDrawing) {
        for (TIME_UNIT_COUNTER = 0; TIME_UNIT_COUNTER < TIME_UNIT; TIME_UNIT_COUNTER += UPDATE_TIME_UNIT) {
            for (int i = 0; i < globalField->nHuman; ++i) {
                //randr = uniRand.value();
                //randt = uniRand.value();
                //angle = uniRandAngle.value();
                //moveP = uniRand.value();
                //revisitP = uniRand.value();
				Human *h = &globalField->humanPool[i];
                //double randr = h->uRandR.value();
                //double randt = h->uRandT.value();
                //double angle = h->uRandA.value();
                //double moveP = h->uRandM.value();
                //double revisitP = h->uRandV.value();
                h->moveModel2();                        //(randr, randt, angle, moveP, revisitP);
				
                
            }
            cout << "Time " << TIME_UNIT_COUNTER << endl;
			         
            if (TIME_UNIT_COUNTER >= 30000) {
                globalField->DumpGrid("grid30000.txt");
                break;
            }
        }
    
        //globalField->DumpDeltaR("deltaR.txt");
        int totalAreas = globalField->MergeGirdAreas();
        cout << "Total areas " << totalAreas << endl;
        globalField->ComputeHumanDensity();
        globalField->ComputeAverageDistance();
        
        
        globalField->DumpCCDF(0, "areaCCDF0.txt", "densityCCDF.txt");   //1 indicates that pick humanNum > 1 grids
        globalField->DumpCCDF(1, "areaCCDF1.txt", "densityCCDF.txt");
        globalField->DumpCCDF(2, "areaCCDF2.txt", "densityCCDF.txt");
        globalField->DumpCCDF(5, "areaCCDF5.txt", "densityCCDF.txt");
        globalField->DumpCCDF(10, "areaCCDF10.txt", "densityCCDF.txt");
        
        globalField->DumpMobilityTopologyAnalysis("gridAreaSizeCCDF.txt", "gridAreaDegreeCCDF.txt", "transportProbCCDF.txt", "size.txt", "link.txt");
        globalField->DumpDistanceMap("distanceMap.txt");
    }
	
	return 1;
}
