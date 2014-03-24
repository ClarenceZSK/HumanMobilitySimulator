#include "field.h"


Field::Field(int maxX_, int minX_, int maxY_, int minY_, 
			 int nHuman_):
maxX(maxX_), minX(minX_), maxY(maxY_), minY(minY_), 
nHuman(nHuman_)

{
    moves = 0;
    nGrid = 0;
    nGridAreas = 0;
	humanPool = new Human[nHuman];
	gridAreaPool = NULL;
	gridPool = NULL;
    deltaR.clear();
}

Field::~Field()
{
	delete [] humanPool;
    delete [] gridPool;
    delete [] gridAreaPool;
}


void Field::EraseVectorElement(vector<Human*> &dl_vector, Human *s)
{
	for(vector<Human*>::iterator dvpos = dl_vector.begin(); dvpos != dl_vector.end(); ++dvpos)
	{
		Human *d = *dvpos;
		if(d == s)
		{
			dl_vector.erase(dvpos);
			break;
		}
	}
}


void Field::InitDeployment(double gridSize)
{
    UniformRandVar uniRand(2500, 7500);
    
    //define grids
    nGrid = maxX/gridSize;
    nGrid = nGrid*nGrid;
    gridPool = new Grid[nGrid];
    for (int i = 0; i < nGrid; ++i) {
        int j = i*gridSize/maxX;
        if (j == 1) {
            bool stop = 1;
        }
        double offset = gridSize/2.0;
        int x = i%maxX;
        Point p(x*gridSize+offset, j*gridSize+offset);
        gridPool[i].init(i, p, RHO_0);      //int index_, Point location_, int rho0_
        gridPool[i].size = gridSize;
        gridPool[i].field = this;
    }
    
	//deploy human
    for (int i = 0; i < nHuman; ++i) {
       
        Point p(maxX/2.0, maxY/2.0);
        humanPool[i].init(i, p);
        humanPool[i].field = this;
    }
    
    //update the grid state that whether it has human
    int gR = humanPool[0].inGrid(gridSize);
	if (gR < 0 || gR >= 1000000)
	{
		cout << "Grid number is wrong" << endl;
		exit(2);
	}
    gridPool[gR].hasHuman = true;
    
    //init gridNum for humanPool
    for (int i = 0; i < nHuman; ++i) {
        humanPool[i].gridNum = gR;
        gridPool[gR].humanNum++;
        humanPool[i].grid = &gridPool[gR];
    }
    
    //get neighbors for each grid
    for (int i = 0; i < nGrid; ++i) {
        gridPool[i].searchingNeighbors();
    }
}

int Field::MergeGirdAreas()
{
    //int counter = 0;
	if (gridAreaPool != NULL) {
		delete [] gridAreaPool;
		gridAreaPool = NULL;
	}
	for(int i = 0; i < nGrid; ++i)
	{
		gridPool[i].broadcastCovered = 0;
		gridPool[i].areaID = -1;
	}

    typedef set<Grid*>::iterator SCITER;
    int areaID = 0;
    for (int i = 0; i < nGrid; ++i)
    {
        if(gridPool[i].areaID < 0 && !gridPool[i].broadcastCovered && gridPool[i].hasHuman)       //not be assigned gridarea ID, start to merge
        {
            set<Grid*> *frontier = new set<Grid*>;
            gridPool[i].broadcastCovered = 1;
            gridPool[i].areaID = areaID;
            frontier->insert(&gridPool[i]);
            while (!frontier->empty())
            {
                set<Grid*> *newFrontier = new set<Grid*>;
                for(set<Grid*>::iterator iter = frontier->begin(); iter != frontier->end(); ++iter)
                {
                    Grid *g = (*iter);
                    for(SCITER gnpos = g->gridNeighbors.begin(); gnpos != g->gridNeighbors.end(); ++gnpos)
                    {
                        Grid *n = *gnpos;
                        if (n != &gridPool[i] && !n->broadcastCovered && n->hasHuman)
                        {
                            n->areaID = areaID;
                            newFrontier->insert(n);
                        }
                        n->broadcastCovered = true;
                    }
                }
                delete frontier;
                frontier = newFrontier;
            }
            delete frontier;
            ++areaID;
        }
    }
    for(int i = 0; i < nGrid; ++i)
	{
		gridPool[i].broadcastCovered = 0;
	}
    
    //Now, we know how many grid areas in total with areaID
    //Initiate grid areas
    nGridAreas = areaID;
    gridAreaPool = new GridArea[nGridAreas];
    for (int j = 0; j < nGrid; ++j) {
        int aID = gridPool[j].areaID;
        if (aID >= 0) {
            gridAreaPool[aID].gridArea.push_back(&gridPool[j]);
            gridAreaPool[aID].index = aID;
            gridAreaPool[aID].field = this;
        }
    }
    
    return nGridAreas;
}

void Field::ComputeAverageDistance()
{
    Point center(maxX/2.0, maxY/2.0);
    double totalDistance = 0;
    int totalNumGrid = 0;
    for (int i = 0; i < nGridAreas; ++i) {
        GridArea *ga = &gridAreaPool[i];
        totalNumGrid = ga->gridArea.size();
        totalDistance = 0;
        for (int j = 0; j < totalNumGrid; ++j) {
            totalDistance += ga->gridArea[j]->location.distance(center);
        }
        ga->averageDistanceToCenter = totalDistance/totalNumGrid;
    }
}

void Field::ComputeHumanDensity()
{
    for (int i = 0; i < nGridAreas; ++i)
    {
        GridArea *ga = &gridAreaPool[i];
        int totalHuman = 0;
        int totalGridArea = ga->gridArea.size();
        for(int j = 0; j < totalGridArea; ++j)
        {
            totalHuman += ga->gridArea[j]->humanNum;
        }
        ga->humanDensity = (double)totalHuman/(totalGridArea*(gridPool[0].size*gridPool[0].size));
    }
}



void Field::DumpData(char *filename1, char *filename2)
{
	if (!filename1 || !filename2)
		return;

	FILE *fp;
	fp = fopen(filename1, "w");
	
    if(!fp) {
		fprintf(stderr, "open file %s to write failed!\n", filename1);
		return;
	}
    
	double time_ = 0;
    double visitedPlaces = 0;
    int updateTimeUnit = 1;
    for (TIME_UNIT_COUNTER = 0; TIME_UNIT_COUNTER < TIME_UNIT; TIME_UNIT_COUNTER += updateTimeUnit) {
        
        visitedPlaces = 0;
        for (int i = 0; i < nHuman; ++i) {
            
            Human *h = &humanPool[i];
            for (int j = h->startIndex; j < h->waitingTimeSequence.size(); ++j) {
                
                time_ = h->waitingTimeSequence[j].first;
            
                if (TIME_UNIT_COUNTER < time_ ) {
                    visitedPlaces += h->waitingTimeSequence[j].second;
                    break;
                }
                else {
                    ++h->startIndex;
                    if (h->startIndex >= h->waitingTimeSequence.size()) {
                        h->startIndex = h->waitingTimeSequence.size()-1;
                        visitedPlaces += h->waitingTimeSequence[h->startIndex].second;
                    }
                }
            }
        }
        visitedPlaces = visitedPlaces/nHuman;
        fprintf(fp, "%d %lf\n", TIME_UNIT_COUNTER, visitedPlaces); //time point, visited places
        //if (TIME_UNIT_COUNTER/updateTimeUnit >= 10) {
        //    updateTimeUnit = updateTimeUnit*10;
        //}
    }
	fclose(fp);
    cout << "File St.txt is outputed." << endl;
    
    
    fp = fopen(filename2, "w");
    if(!fp) {
		fprintf(stderr, "open file %s to write failed!\n", filename2);
		return;
	}
    
    //double visitFrequency = 0;
    //Point location;
    multimap<double, Point, greater<double> >::iterator nflpos;
    int minCount = 0xffff;
    //int count = 0;
    for (int i = 0; i < nHuman; ++i) {
        //count = 0;
        Human *h = &humanPool[i];
        for (nflpos = h->normalizedFrequencyLocations.begin(); nflpos != h->normalizedFrequencyLocations.end(); ++nflpos) {
            //++count;
            h->visitFrequency.push_back(nflpos->first);
            //location = nflpos->second;
            //fprintf(fp, "%lf %lf %lf\n", visitFrequency, location.x, location.y); //visit frequency, location.x, location.y
            //fprintf(fp, "%d %lf\n", count, visitFrequency); //id, visit frequency,
        }
        if (minCount > h->normalizedFrequencyLocations.size()) {
            minCount = h->normalizedFrequencyLocations.size();
        }
    }
    minCount = 1000;
    double totalFrequency = 0;
    for (int j = 0; j < minCount; ++j) {
        totalFrequency = 0;
        for (int i = 0; i < nHuman; ++i) {
            
            Human *h = &humanPool[i];
            totalFrequency += h->visitFrequency[j];
            
        }
        totalFrequency = totalFrequency/nHuman;
        fprintf(fp, "%d %lf\n", j, totalFrequency); //id, visit frequency,
    }
	fclose(fp);
    cout << "File fk.txt is outputed." << endl;
    
}

void Field::DumpRG(char *filename1, char *filename2)        //rgt.txt, Prg.txt
{
	if (!filename1 || !filename2)
		return;
    
	FILE *fp;
	fp = fopen(filename1, "w");         //rgt.txt
	
    if(!fp) {
		fprintf(stderr, "open file %s to write failed!\n", filename1);
		return;
	}
    
	double time_ = 0;
    map<double, double>::iterator rg_pos;
    double averRG = 0;
    int updateTimeUnit = 1;
    for (TIME_UNIT_COUNTER = 0; TIME_UNIT_COUNTER < TIME_UNIT; TIME_UNIT_COUNTER += updateTimeUnit) {
        
        averRG = 0;
        for (int i = 0; i < nHuman; ++i) {
            
            Human *h = &humanPool[i];
            for (int j = h->startIndex; j < h->rg_t.size(); ++j) {
                
                time_ = h->rg_t[j].first;
                
                if (TIME_UNIT_COUNTER < time_ ) {
                    averRG += h->rg_t[j].second/h->maxRGT;
                    break;
                }
                else {
                    ++h->startIndex;
                    if (h->startIndex >= h->rg_t.size()) {
                        h->startIndex = h->rg_t.size()-1;
                        averRG += h->rg_t[h->startIndex].second/h->maxRGT;
                    }
                }
            }
        }
        averRG = averRG/nHuman;
        fprintf(fp, "%d %lf\n", TIME_UNIT_COUNTER, averRG); //time point, rg_t
        //if (TIME_UNIT_COUNTER/updateTimeUnit >= 10) {
        //    updateTimeUnit = updateTimeUnit*10;
        //}
    }
	fclose(fp);
    cout << "File rgt.txt is outputed." << endl;
    
    
    fp = fopen(filename2, "w");         //output Prg.txt
    if(!fp) {
		fprintf(stderr, "open file %s to write failed!\n", filename2);
		return;
	}
    
    vector<double> forOutput;
    
    for (int i = 0; i < nHuman; ++i)
    {
        Human *h = &humanPool[i];
        for (int j = 0; j < h->rg_t.size(); ++j) {
            forOutput.push_back(h->rg_t[j].second);
        }
    }
    sort(forOutput.begin(), forOutput.end());
    
    //binding forOutput, binding 1, 2, 4, 8
    //multiset<double>::iterator mapgapos = forOutput.begin();
    int division = 10;
    for (int i = 0; i < forOutput.size(); ++i) {
        double value  = forOutput[i];
        double tag = value/division;
        if (tag < 1) {
            
            if (tag >= 0 && tag < 0.15)
                forOutput[i] = 1*division/10;
            else if (tag >= 0.15 && tag < 0.3)
                forOutput[i] = 2*division/10;
            else if (tag >= 0.3 && tag < 0.6)
                forOutput[i] = 4*division/10;
            else if (tag >= 0.6 && tag < 0.9)
                forOutput[i] = 8*division/10;
            else if (tag >= 0.9 && tag < 1)
                forOutput[i] = 10*division/10;
            
        }
        else {
            --i;
            division = division*10;
        }
    }
    
    int counter = 0;
    //mapgapos = forOutput.begin();
    int i = 0;
    multimap<double, double> PDF;              //area, portion
    while (i < forOutput.size())
    {
        double rg1 = forOutput[i];
        ++counter;
        ++i;
        if (i == forOutput.size()) {
            double portion = (double)counter/forOutput.size();
            PDF.insert(make_pair(rg1, portion));
            continue;
        }
        double rg2 = forOutput[i];
        if (rg1 != rg2)             //if not equal, we count another
        {
            double portion = (double)counter/forOutput.size();
            PDF.insert(make_pair(rg1, portion));
            counter = 0;
        }
    }
    
    multimap<double, double>::iterator mapdopos;
    //fprintf(fp, "#%lf\n", RHO_0); //rho_0
    for (mapdopos = PDF.begin(); mapdopos != PDF.end(); ++mapdopos) {
        
        fprintf(fp, "%lf %lf\n", mapdopos->first, mapdopos->second);        //rg, portion
    }
	fclose(fp);
    cout << "File Prg.txt is outputed." << endl;
}

void Field::DumpCCDF(int x, char *filename1, char *filename2)
{
	if (!filename1 || !filename2)
		return;
    
	FILE *fp;
    //output connectivity area distribution
    fp = fopen(filename1, "w");     //areaCCDF
    if(!fp) {
		fprintf(stderr, "open file %s to write failed!\n", filename1);
		return;
	}
    //analyze the connectivity area distribution
    set<GridArea*> chosenAreas;
    for (int i = 0; i < nGrid; ++i) {
        if (gridPool[i].humanNum > x) {     //humanNum > x
            chosenAreas.insert(&gridAreaPool[gridPool[i].areaID]);
        }
    }
    multimap<double, GridArea*, greater<double> > forOutput;
    set<GridArea*>::iterator gapos;
    for (gapos = chosenAreas.begin(); gapos != chosenAreas.end(); ++gapos)
    {
        GridArea *ga = *gapos;
        double area = ga->gridArea.size()*(gridPool[0].size*gridPool[0].size);    //here we assume all grids are squares.
        forOutput.insert(make_pair(area, ga));
    }
    int counter = 0;
    multimap<double, GridArea*, greater<double> >::iterator mapgapos = forOutput.begin();
    multimap<double, double> CCDF;              //area, portion
    while (mapgapos != forOutput.end())
    {
        double area1 = mapgapos->first;
        ++counter;
        ++mapgapos;
        if (mapgapos == forOutput.end()) {
            double portion = (double)counter/chosenAreas.size();
            CCDF.insert(make_pair(area1, portion));
            continue;
        }
        double area2 = mapgapos->first;
        if (area1 != area2)             //if not equal, we count another
        {
            double portion = (double)counter/chosenAreas.size();
            CCDF.insert(make_pair(area1, portion));
            //counter = 0;
        }
    }
    
    multimap<double, double>::iterator mapdopos;
    fprintf(fp, "#%lf\n", RHO_0); //rho_0
    for (mapdopos = CCDF.begin(); mapdopos != CCDF.end(); ++mapdopos) {
        
        fprintf(fp, "%lf %lf\n", mapdopos->first, mapdopos->second);        //area, portion
    }
	fclose(fp);
    
    ////////////////////////////////////////////////////////////////////////////////
    //output human density distribution
    fp = fopen(filename2, "w");
    if(!fp) {
		fprintf(stderr, "open file %s to write failed!\n", filename2);
		return;
	}
    
    //analyze the density distribution
    forOutput.clear();      //here, we store density, grid area pointer
    
    for (int i = 0; i < nGridAreas; ++i)
    {
        GridArea *ga = &gridAreaPool[i];
        forOutput.insert(make_pair(ga->humanDensity, ga));
    }
    
    counter = 0;
    mapgapos = forOutput.begin();
    CCDF.clear();                       //density, portion
    while (mapgapos != forOutput.end())
    {
        double area1 = mapgapos->first;
        ++counter;
        ++mapgapos;
        if (mapgapos == forOutput.end()) {
            double portion = (double)counter/forOutput.size();
            CCDF.insert(make_pair(area1, portion));
            continue;
        }
        double area2 = mapgapos->first;
        if (area1 != area2)             //if not equal, we count another
        {
            double portion = (double)counter/forOutput.size();
            CCDF.insert(make_pair(area1, portion));
            //counter = 0;
        }
    }
    fprintf(fp, "#%lf\n", RHO_0); //rho_0
    for (mapdopos = CCDF.begin(); mapdopos != CCDF.end(); ++mapdopos) {
        
        fprintf(fp, "%lf %lf\n", mapdopos->first, mapdopos->second); //area, portion
    }
	fclose(fp);
}

void Field::DumpDistanceMap(char *filename1)    //output distance map for each grid area
{
    if (!filename1)
		return;
    
	FILE *fp;
    fp = fopen(filename1, "w");     //areaCCDF
    if(!fp) {
		fprintf(stderr, "open file %s to write failed!\n", filename1);
		return;
	}
    
    for (int i = 0; i < nGridAreas; ++i) {
        GridArea *ga = &gridAreaPool[i];
        fprintf(fp, "%lf %lf\n", ga->averageDistanceToCenter, ga->humanDensity); //average distance, density
    }
    fclose(fp);
}

void Field::DumpGrid(char *filename)
{
	if (!filename)
		return;
    
	FILE *fp;
	fp = fopen(filename, "w");
	
    if(!fp) {
		fprintf(stderr, "open file %s to write failed!\n", filename);
		return;
	}
    
    for (int i  = 0; i < nGrid; ++i) {
        Grid *g = &gridPool[i];
        if (g->hasHuman) {
            fprintf(fp, "%d %lf %lf\n", i, g->location.x, g->location.y); //grid ID, coordinate x, coordinate y
        }
        
    }
	fclose(fp);
}
#if 0
void Field::DumpMobilityTopologyAnalysis(char *filename1, char *filename2, char *filename3)
{
    //analyze GridArea degree distribution, transport probability distribution, GridArea size distribution
    for (int i = 0; i < nHuman; ++i) {
        //Explore GridArea neighbors
        Human *h = &humanPool[i];
		if (h->jumpGridRecord.empty())
			continue;
        //cout << "Human " << i << endl;
        Grid *inFirstArea = NULL;
        Grid *inSecondArea = NULL;
        bool firstJump = true;                              //we do not count the initial stand at the center
        for (int j = 0; j < h->jumpGridRecord.size(); ++j) {
            
            Grid *preG = h->jumpGridRecord[j].first;
            Grid *currentG = h->jumpGridRecord[j].second;
            
            if (preG->areaID >= 0) {
                inFirstArea = preG;
            }
            else if (inFirstArea == NULL) {
                continue;
            }
             
            //inFirstArea = preG;
            //assert(inFirstArea->areaID >= 0);
            
            
            if (currentG->areaID >= 0) {       //it means this grid was touched once, but now it is empty, no human
                inSecondArea = currentG;
            }
            else {
                continue;
            }
            
            if (firstJump) {                                               //we do not count the initial point
                firstJump = false;
                ++gridAreaPool[inSecondArea->areaID].gridAreaVisitedTimes;  //gridAreaVisitedTimes determins the gridarea size
            }
            else {
                ++gridAreaPool[inFirstArea->areaID].gridAreaVisitedTimes;
                ++gridAreaPool[inSecondArea->areaID].gridAreaVisitedTimes;
            }
            
            //compare inFirstArea->areaID and inSecondArea->areaID, if == then they are in the same grid area
            if (inFirstArea->areaID == inSecondArea->areaID) {
                continue;
            }
            else {
                GridArea *cGA = &gridAreaPool[inSecondArea->areaID];
                map<GridArea*, int>::iterator gapos = gridAreaPool[inFirstArea->areaID].mobilityTopologyNeighbors.find(cGA);
                
                if (gapos == gridAreaPool[inFirstArea->areaID].mobilityTopologyNeighbors.end()) {      //not find, add new neighbor
                    gridAreaPool[inFirstArea->areaID].mobilityTopologyNeighbors.insert(make_pair(cGA, 1));
                    cGA->mobilityTopologyNeighbors.insert(make_pair(&gridAreaPool[inFirstArea->areaID], 1));
                }
                else {
                    ++gapos->second;                                                            //increase transport times
                }
            }
            inFirstArea = NULL;
            inSecondArea = NULL;
        }
    }
    //mobilityTopologyNeighbors.size() is GridArea degree; mobilityTopologyNeighbors -> second is for transport probability;
    //gridAreaVisitedTimes is for GridArea size
    if (!filename1 || !filename2 || !filename3)
		return;
    
	FILE *fp;
    
    //output CCDF files
    //1. output GridArea size distribution
    multiset<int, greater<int> > forOutput;
    for (int i = 0; i < nGridAreas; ++i) {
        int vT = gridAreaPool[i].gridAreaVisitedTimes;
        forOutput.insert(vT);
    }
    
    fp = fopen(filename1, "w");
	
    if(!fp) {
		fprintf(stderr, "open file %s to write failed!\n", filename1);
		return;
	}
    
    int counter = 0;
    multiset<int, greater<int> >::iterator setpos = forOutput.begin();
    multimap<int, double> CCDF;         //gridarea node size, portion
    
    multimap<int, double>::iterator ccdfpos;
    CCDF.clear();
    while (setpos != forOutput.end())
    {
        double size1 = *setpos;
        ++counter;
        ++setpos;
        if (setpos == forOutput.end()) {
            double portion = (double)counter/forOutput.size();
            CCDF.insert(make_pair(size1, portion));
            continue;
        }
        double size2 = *setpos;
        if (size1 != size2)             //if not equal, we count another
        {
            double portion = (double)counter/forOutput.size();
            CCDF.insert(make_pair(size1, portion));
            //counter = 0;
        }
    }
    for (ccdfpos = CCDF.begin(); ccdfpos != CCDF.end(); ++ccdfpos) {
        
        fprintf(fp, "%d %lf\n", ccdfpos->first, ccdfpos->second); //gridarea node size, portion
    }
    cout << "File gridAreaSizeCCDF is outputed" << endl;
	fclose(fp);
    
    //2. output GridArea degree distribution
    fp = fopen(filename2, "w");
    if(!fp) {
		fprintf(stderr, "open file %s to write failed!\n", filename2);
		return;
	}
    
    forOutput.clear();
    for (int i = 0; i < nGridAreas; ++i) {
        assert(!gridAreaPool[i].mobilityTopologyNeighbors.empty());
        int de = gridAreaPool[i].mobilityTopologyNeighbors.size();
        forOutput.insert(de);
    }
    
    counter = 0;
    setpos = forOutput.begin();
    CCDF.clear();
    while (setpos != forOutput.end())
    {
        double degree1 = *setpos;
        ++counter;
        ++setpos;
        if (setpos == forOutput.end()) {
            double portion = (double)counter/forOutput.size();
            CCDF.insert(make_pair(degree1, portion));
            continue;
        }
        double degree2 = *setpos;
        if (degree1 != degree2)             //if not equal, we count another
        {
            double portion = (double)counter/forOutput.size();
            CCDF.insert(make_pair(degree1, portion));
            //counter = 0;
        }
    }
    for (ccdfpos = CCDF.begin(); ccdfpos != CCDF.end(); ++ccdfpos) {
        
        fprintf(fp, "%d %lf\n", ccdfpos->first, ccdfpos->second); //gridarea degree, portion
    }
    cout << "File grid area DegreeCCDF is outputed" << endl;
	fclose(fp);
    
    //3. output transport probability distribution
    fp = fopen(filename3, "w");
    if(!fp) {
		fprintf(stderr, "open file %s to write failed!\n", filename3);
		return;
	}
    
    //forOutput.clear();
    int totalTransport = 0;
    map<GridArea*, int>::iterator mobilitypos;
    
    //transform into probability
    for (int i = 0; i < nGridAreas; ++i) {
        for (mobilitypos = gridAreaPool[i].mobilityTopologyNeighbors.begin(); mobilitypos != gridAreaPool[i].mobilityTopologyNeighbors.end(); ++mobilitypos) {
            
            totalTransport += mobilitypos->second;
            
        }
        for (mobilitypos = gridAreaPool[i].mobilityTopologyNeighbors.begin(); mobilitypos != gridAreaPool[i].mobilityTopologyNeighbors.end(); ++mobilitypos) {
            
            double prob = (double)mobilitypos->second/totalTransport;
            gridAreaPool[i].mobilityTransportProb.insert(make_pair(mobilitypos->first, prob));
            
        }
        
        totalTransport = 0;
    }
    
    multiset<double, greater<double> > forOutputProb;
    map<GridArea*, double>::iterator probpos;
    for (int i = 0; i < nGridAreas; ++i) {
        for (probpos = gridAreaPool[i].mobilityTransportProb.begin(); probpos != gridAreaPool[i].mobilityTransportProb.end(); ++probpos) {
            
            forOutputProb.insert(probpos->second);
            
        }
    }
    multimap<double, double> probCCDF;
    multiset<double, greater<double> >::iterator probsetpos;
    multimap<double, double>::iterator probccdfpos;
    counter = 0;
    probsetpos = forOutputProb.begin();
    //CCDF.clear();
    while (probsetpos != forOutputProb.end())
    {
        double prob1 = *probsetpos;
        ++counter;
        ++probsetpos;
        if (probsetpos == forOutputProb.end()) {
            double portion = (double)counter/forOutputProb.size();
            probCCDF.insert(make_pair(prob1, portion));
            continue;
        }
        double prob2 = *probsetpos;
        if (prob1 != prob2)             //if not equal, we count another
        {
            double portion = (double)counter/forOutputProb.size();
            probCCDF.insert(make_pair(prob1, portion));
            //counter = 0;
        }
    }
    for (probccdfpos = probCCDF.begin(); probccdfpos != probCCDF.end(); ++probccdfpos) {
        
        fprintf(fp, "%lf %lf\n", probccdfpos->first, probccdfpos->second); //transport probability, portion
    }
    cout << "File transport prob CCDF is outputed" << endl;
	fclose(fp);
}
#endif

void Field::DumpMobilityTopologyAnalysis(char *filename1, char *filename2, char *filename3, char *filename4, char *filename5)
{
    //analyze GridArea degree distribution, transport probability distribution, GridArea size distribution
    for (int i = 0; i < nHuman; ++i) {
        //Explore GridArea neighbors
        Human *h = &humanPool[i];
		if (h->jumpGridRecord.empty())
			continue;
        
        Grid *inFirstArea = NULL;
        Grid *inSecondArea = NULL;
        bool firstJump = true;                              //we do not count the initial stand at the center
        for (int j = 0; j < h->jumpGridRecord.size(); ++j) {
            
            Grid *preG = h->jumpGridRecord[j].first;
            Grid *currentG = h->jumpGridRecord[j].second;
            
            if (preG->areaID >= 0) {
                inFirstArea = preG;
            }
            else if (inFirstArea == NULL) {
                continue;
            }
            
            //inFirstArea = preG;
            //assert(inFirstArea->areaID >= 0);
            
            if (currentG->areaID >= 0) {       //it means this grid was touched once, but now it is empty, no human
                inSecondArea = currentG;
            }
            else {
                continue;
            }
            
            if (firstJump) {                                               //we do not count the initial point
                firstJump = false;
                ++gridAreaPool[inSecondArea->areaID].gridAreaVisitedTimes;  //gridAreaVisitedTimes determins the gridarea size
            }
            else {
                ++gridAreaPool[inFirstArea->areaID].gridAreaVisitedTimes;
                ++gridAreaPool[inSecondArea->areaID].gridAreaVisitedTimes;
            }
            
            //compare inFirstArea->areaID and inSecondArea->areaID, if == then they are in the same grid area
            if (inFirstArea->areaID == inSecondArea->areaID) {
                continue;
            }
            else {
                GridArea *cGA = &gridAreaPool[inSecondArea->areaID];
                map<GridArea*, int>::iterator gapos = gridAreaPool[inFirstArea->areaID].mobilityTopologyNeighbors.find(cGA);
                
                if (gapos == gridAreaPool[inFirstArea->areaID].mobilityTopologyNeighbors.end()) {      //not find, add new neighbor
                    gridAreaPool[inFirstArea->areaID].mobilityTopologyNeighbors.insert(make_pair(cGA, 1));
                    cGA->mobilityTopologyNeighbors.insert(make_pair(&gridAreaPool[inFirstArea->areaID], 1));
                }
                else {
                    ++gapos->second;                                                            //increase transport times
                }
            }
            inFirstArea = NULL;
            inSecondArea = NULL;
        }
    }
    //mobilityTopologyNeighbors.size() is GridArea degree; mobilityTopologyNeighbors -> second is for transport probability;
    //gridAreaVisitedTimes is for GridArea size
    if (!filename1 || !filename2 || !filename3 || !filename4 || !filename5)
		return;
    
	FILE *fp;
    
    //output CCDF files
    //1. output GridArea size distribution
    multiset<int, greater<int> > forOutput;
    int totalVT = 0;
    for (int i = 0; i < nGridAreas; ++i) {
        int vT = gridAreaPool[i].gridAreaVisitedTimes;
        totalVT += vT;
        forOutput.insert(vT);
    }
    
    fp = fopen(filename4, "w");     //output size.txt
	
    if(!fp) {
		fprintf(stderr, "open file %s to write failed!\n", filename4);
		return;
	}
    
    for (int i = 0; i < nGridAreas; ++i) {
        int vT = gridAreaPool[i].gridAreaVisitedTimes;
        double portion = (double)vT/totalVT;
        fprintf(fp, "%d %lf\n", i, portion); //grid ID, portion
    }
    cout << "File size.txt is outputed" << endl;
	fclose(fp);
    
    fp = fopen(filename1, "w");
	
    if(!fp) {
		fprintf(stderr, "open file %s to write failed!\n", filename1);
		return;
	}
    
    int counter = 0;
    multiset<int, greater<int> >::iterator setpos = forOutput.begin();
    multimap<int, double> CCDF;         //gridarea node size, portion
    
    multimap<int, double>::iterator ccdfpos;
    CCDF.clear();
    while (setpos != forOutput.end())
    {
        double size1 = *setpos;
        ++counter;
        ++setpos;
        if (setpos == forOutput.end()) {
            double portion = (double)counter/forOutput.size();
            CCDF.insert(make_pair(size1, portion));
            continue;
        }
        double size2 = *setpos;
        if (size1 != size2)             //if not equal, we count another
        {
            double portion = (double)counter/forOutput.size();
            CCDF.insert(make_pair(size1, portion));
        }
    }
    for (ccdfpos = CCDF.begin(); ccdfpos != CCDF.end(); ++ccdfpos) {
        
        fprintf(fp, "%d %lf\n", ccdfpos->first, ccdfpos->second); //gridarea node size, portion
    }
    cout << "File gridAreaSizeCCDF is outputed" << endl;
	fclose(fp);
    
    //2. output GridArea degree distribution
    fp = fopen(filename2, "w");
    if(!fp) {
		fprintf(stderr, "open file %s to write failed!\n", filename2);
		return;
	}
    
    forOutput.clear();
    for (int i = 0; i < nGridAreas; ++i) {
        assert(!gridAreaPool[i].mobilityTopologyNeighbors.empty());
        int de = gridAreaPool[i].mobilityTopologyNeighbors.size();
        forOutput.insert(de);
    }
    
    counter = 0;
    setpos = forOutput.begin();
    CCDF.clear();
    while (setpos != forOutput.end())
    {
        double degree1 = *setpos;
        ++counter;
        ++setpos;
        if (setpos == forOutput.end()) {
            double portion = (double)counter/forOutput.size();
            CCDF.insert(make_pair(degree1, portion));
            continue;
        }
        double degree2 = *setpos;
        if (degree1 != degree2)             //if not equal, we count another
        {
            double portion = (double)counter/forOutput.size();
            CCDF.insert(make_pair(degree1, portion));
            //counter = 0;
        }
    }
    for (ccdfpos = CCDF.begin(); ccdfpos != CCDF.end(); ++ccdfpos) {
        
        fprintf(fp, "%d %lf\n", ccdfpos->first, ccdfpos->second); //gridarea degree, portion
    }
    cout << "File grid area DegreeCCDF is outputed" << endl;
	fclose(fp);
    ////////////////////////////////////////
    //3. output transport probability distribution
    
    
    
    int totalTransport = 0;
    map<GridArea*, int>::iterator mobilitypos;
    
    //transform into probability
    for (int i = 0; i < nGridAreas; ++i) {
        for (mobilitypos = gridAreaPool[i].mobilityTopologyNeighbors.begin(); mobilitypos != gridAreaPool[i].mobilityTopologyNeighbors.end(); ++mobilitypos) {
            
            totalTransport += mobilitypos->second;
            
        }
    }
    
    fp = fopen(filename5, "w");     //output link.txt
    if(!fp) {
		fprintf(stderr, "open file %s to write failed!\n", filename5);
		return;
	}
    
    for (int i = 0; i < nGridAreas; ++i) {
        for (mobilitypos = gridAreaPool[i].mobilityTopologyNeighbors.begin(); mobilitypos != gridAreaPool[i].mobilityTopologyNeighbors.end(); ++mobilitypos) {
            
            double prob = (double)mobilitypos->second/totalTransport;
            gridAreaPool[i].mobilityTransportProb.insert(make_pair(mobilitypos->first, prob));  //neighbor area, prob
            fprintf(fp, "%d %d %lf\n", i, mobilitypos->first->index, prob); //grid ID, transport neighbor ID, portion
        }
    }
    
    cout << "File link.txt is outputed" << endl;
	fclose(fp);

    fp = fopen(filename3, "w");
    if(!fp) {
		fprintf(stderr, "open file %s to write failed!\n", filename3);
		return;
	}
    
    multiset<double, greater<double> > forOutputProb;
    map<GridArea*, double>::iterator probpos;
    for (int i = 0; i < nGridAreas; ++i) {
        for (probpos = gridAreaPool[i].mobilityTransportProb.begin(); probpos != gridAreaPool[i].mobilityTransportProb.end(); ++probpos) {
            
            forOutputProb.insert(probpos->second);
            
        }
    }
    multimap<double, double> probCCDF;
    multiset<double, greater<double> >::iterator probsetpos;
    multimap<double, double>::iterator probccdfpos;
    counter = 0;
    probsetpos = forOutputProb.begin();
    
    while (probsetpos != forOutputProb.end())
    {
        double prob1 = *probsetpos;
        ++counter;
        ++probsetpos;
        if (probsetpos == forOutputProb.end()) {
            double portion = (double)counter/forOutputProb.size();
            probCCDF.insert(make_pair(prob1, portion));
            continue;
        }
        double prob2 = *probsetpos;
        if (prob1 != prob2)             //if not equal, we count another
        {
            double portion = (double)counter/forOutputProb.size();
            probCCDF.insert(make_pair(prob1, portion));
            
        }
    }
    for (probccdfpos = probCCDF.begin(); probccdfpos != probCCDF.end(); ++probccdfpos) {
        
        fprintf(fp, "%lf %lf\n", probccdfpos->first, probccdfpos->second); //transport probability, portion
    }
    cout << "File transport prob CCDF is outputed" << endl;
	fclose(fp);
}

void Field::DumpDeltaR(char *filename)
{
    if (!filename)
		return;
    
	FILE *fp;
    fp = fopen(filename, "w");
    if(!fp) {
		fprintf(stderr, "open file %s to write failed!\n", filename);
		return;
	}
    
    multimap<double, double> probCCDF;
    multiset<double, greater<double> >::iterator deltaRsetpos;
    multimap<double, double>::iterator probccdfpos;
    int counter = 0;
    deltaRsetpos = deltaR.begin();
    
    while (deltaRsetpos != deltaR.end())
    {
        double prob1 = *deltaRsetpos;
        ++counter;
        ++deltaRsetpos;
        if (deltaRsetpos == deltaR.end()) {
            double portion = (double)counter/deltaR.size();
            probCCDF.insert(make_pair(prob1, portion));
            continue;
        }
        double prob2 = *deltaRsetpos;
        if (prob1 != prob2)             //if not equal, we count another
        {
            double portion = (double)counter/deltaR.size();
            probCCDF.insert(make_pair(prob1, portion));
            
        }
    }
    for (probccdfpos = probCCDF.begin(); probccdfpos != probCCDF.end(); ++probccdfpos) {
        
        fprintf(fp, "%lf %lf\n", probccdfpos->first, probccdfpos->second); //transport probability, portion
    }
    cout << "File deltaR is outputed" << endl;
	fclose(fp);
}
