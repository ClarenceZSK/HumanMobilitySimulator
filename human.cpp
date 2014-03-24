#include "human.h"

Human::Human()
{
	index = -1;
    field = NULL;
    grid = NULL;
    
    powerLawT = 0;
    timeCounter = 0;
    timePassed = 0;
    
    startIndex = 1;
    visitFrequency.clear();
    
    waitNow = false;
    readyToMove = false;
    waitingTimeSequence.push_back(make_pair(0, 1));
	jumpGridRecord.clear();
	normalizedFrequencyLocations.clear();
	userVisitedRecord.clear();
	visitedLocations.clear();
    rg_t.clear();
    maxRGT = 0;
    
	gridNum = -1;
    
    uRandR.setmin(0);
    uRandR.setmax(1);
    
    uRandT.setmin(0);
    uRandT.setmax(1);
    
    uRandM.setmin(0);
    uRandM.setmax(1);
    
    uRandV.setmin(0);
    uRandV.setmax(1);
    
    uRandA.setmin(0);
    uRandA.setmax(360);
}

Human::~Human()
{
}

void Human::init(int index_, Point location_)
{
	index = index_;
	currentMapPosition.position = location_;
    currentMapPosition.visitFrequency = 1;
    //waitingTimeSequence.push_back(make_pair(0, 1));
}
//move step unit version
#if 0
void Human::moveModel1(double randr, double randt, double angle, double moveP)
{
    //init moving
    if (visitedLocations.empty()) {
        visitedLocations.push_back(currentMapPosition);
    }
    
    //double r = uniRand.value();
    double x0 = 1;
    double x1 = 1000;   //range from 0 to 1000m
    double alpha = -1.55;
    double rho = 0.1;
    double gamma = -0.21;
    double powerLawR = pow(((pow(x1, alpha+1) - pow(x0, alpha+1))*randr + pow(x0, alpha+1)), 1/(alpha+1));  //jump size
    double Pnew = rho * pow((double)visitedLocations.size(), gamma); //Exploration probability
    
    //double t = uniRand.value();
    double t0 = 1;
    double t1 = 17*3600;     //range from 0 to 17 hours
    double beta = -1.8;
    double powerLawT = pow(((pow(t1, beta+1) - pow(t0, beta+1))*randt + pow(t0, beta+1)), 1/(beta+1));  //waiting time
    
    
    
    //double moveP = uniRand.value();
    double newX = currentMapPosition.position.x + powerLawR*cos(angle*PI/180);
    double newY = currentMapPosition.position.y + powerLawR*sin(angle*PI/180);
    if (moveP <= Pnew) {
        //Explore a new location
        Point tp(newX, newY);
        currentMapPosition.position = tp;
        currentMapPosition.visitFrequency = 1;
        visitedLocations.push_back(currentMapPosition);
    }
    else {
        //Go back to some places that the user visited
        /*
        Point tempLocation(newX, newY);
        //find the nearest visited location
        double minD = 0xffff;
        int target = 0;
        for (int i = 0; i < visitedLocations.size(); ++i) {
            double Dis = tempLocation.distance(visitedLocations[i].position);
            if (Dis < minD) {
                minD = Dis;
                target = i;
            }
        }
        visitedLocations[target].visitFrequency++;
        currentMapPosition = visitedLocations[target];
         */
        multimap<int, Maplocation>::iterator lvpos; //location visited pos
        int totalVisit = 0;
        for (lvpos = userVisitedRecord.begin(); lvpos != userVisitedRecord.end(); ++lvpos) {
            totalVisit += lvpos->first;
        }
        
        multimap<double, Maplocation> distribution;      //probability, location
        double preProb = 0;
        for (lvpos = userVisitedRecord.begin(); lvpos != userVisitedRecord.end(); ++lvpos) {
            double prob = (double)lvpos->first/totalVisit;
            distribution.insert(make_pair(preProb+prob, lvpos->second));
            preProb += prob;
        }
        
        //test
        multimap<double, Maplocation>::iterator dispos = distribution.end();
        --dispos;
        assert(dispos->first - 1 <= 0.00001);
        //end test
        //search the revisit location
        for (dispos = distribution.begin(); dispos != distribution.end(); ++dispos) {
            if (revisitP <= dispos->first) {   //this is the target
                currentMapPosition = dispos->second;
                //locate the corresponding location and update visited times
                for (lvpos = userVisitedRecord.begin(); lvpos != userVisitedRecord.end(); ++lvpos) {
                    if (lvpos->second == dispos->second) {
                        int value = lvpos->first;
                        Maplocation m = lvpos->second;
                        userVisitedRecord.erase(lvpos);
                        ++value;
                        userVisitedRecord.insert(make_pair(value, m));
                        break;
                    }
                }
                break;
            }
        }
        
    }
    waitingTimeSequence.push_back(make_pair(powerLawT, visitedLocations.size()));
}
#endif
//time unit version
void Human::moveModel1()      //(double randr, double randt, double angle, double moveP, double revisitP)
{
    //init moving
    if (userVisitedRecord.empty()) {
        userVisitedRecord.insert(make_pair(1, currentMapPosition));
		
    }
    
    //double powerLawT = 0;
    //double timeCounter = 0;
    
    //if it does not wait, compute wating time after moving
    if (!waitNow) {
        double t0 = 1;
        double t1 = 17*3600;     //range from 0 to 17 hours
        double beta = -1.8;
        double randt = uRandT.value();
        powerLawT = pow(((pow(t1, beta+1) - pow(t0, beta+1))*randt + pow(t0, beta+1)), 1/(beta+1));  //waiting time
        timeCounter = powerLawT;
        waitNow = true;
        
    }
    
    //if ready to move, compute jump distance. Otherwise, wait
    if (waitNow && !readyToMove) {
        timeCounter -= UPDATE_TIME_UNIT;
        if (timeCounter <= 0) {
            waitNow = false;
            readyToMove = true;
        }
        else
            waitNow = true;
    }
    else {
        readyToMove = false;
        double x0 = 1;
        double x1 = 1000;   //range from 0 to x1 meters
        double alpha = -1.55;
        double rho = 0.1;
        double gamma = -0.21;
        double randr = uRandR.value();
        double powerLawR = pow(((pow(x1, alpha+1) - pow(x0, alpha+1))*randr + pow(x0, alpha+1)), 1/(alpha+1));  //jump size
        double Pnew = rho * pow((double)userVisitedRecord.size(), gamma); //Exploration probability
        
        double moveP = uRandM.value();
        double angle = uRandA.value();
        if (moveP <= Pnew) {
            //Explore a new location
            readyToMove = false;
            double newX = currentMapPosition.position.x + powerLawR*cos(angle*PI/180);
            double newY = currentMapPosition.position.y + powerLawR*sin(angle*PI/180);
            Point tp(newX, newY);
            currentMapPosition.position = tp;
            currentMapPosition.visitFrequency = 1;
            userVisitedRecord.insert(make_pair(1, currentMapPosition));
            updateGridInfo();       //record every jump step and update all corresponding information
        }
        else {
            //Go back to some places that the user visited
            /*
            Point tempLocation(newX, newY);
            //find the nearest visited location
            double minD = 0xffff;
            int target = 0;
            for (int i = 0; i < visitedLocations.size(); ++i) {
                double Dis = tempLocation.distance(visitedLocations[i].position);
                if (Dis < minD) {
                    minD = Dis;
                    target = i;
                }
            }
            visitedLocations[target].visitFrequency++;
            currentMapPosition = visitedLocations[target];
             */
            multimap<int, Maplocation>::iterator lvpos; //location visited pos
            int totalVisit = 0;
            for (lvpos = userVisitedRecord.begin(); lvpos != userVisitedRecord.end(); ++lvpos) {
                totalVisit += lvpos->first;
            }
            
            multimap<double, Maplocation> distribution;      //probability, location
            double preProb = 0;
            for (lvpos = userVisitedRecord.begin(); lvpos != userVisitedRecord.end(); ++lvpos) {
                double prob = (double)lvpos->first/totalVisit;
                distribution.insert(make_pair(preProb+prob, lvpos->second));
                preProb += prob;
            }
            
            //test
            multimap<double, Maplocation>::iterator dispos = distribution.end();
            --dispos;
            assert(dispos->first - 1 <= 0.00001);
            //end test
            //search the revisit location
            double revisitP = uRandV.value();
            for (dispos = distribution.begin(); dispos != distribution.end(); ++dispos) {
                if (revisitP <= dispos->first) {   //this is the target
                    if (currentMapPosition != dispos->second) {
                        currentMapPosition = dispos->second;
                        readyToMove = false;
                    }
                    //locate the corresponding location and update visited times
                    for (lvpos = userVisitedRecord.begin(); lvpos != userVisitedRecord.end(); ++lvpos) {
                        if (lvpos->second == dispos->second) {
                            int value = lvpos->first;
                            Maplocation m = lvpos->second;
                            userVisitedRecord.erase(lvpos);
                            ++value;
                            userVisitedRecord.insert(make_pair(value, m));
                            break;
                        }
                    }
                    break;
                }
            }
            updateGridInfo();
        }
    }
    //waitingTimeSequence.push_back(make_pair(powerLawT, visitedLocations.size()));
}

void Human::moveModel2()          //(double randr, double randt, double angle, double moveP, double revisitP)
{
    //init moving
    if (userVisitedRecord.empty()) {
        userVisitedRecord.insert(make_pair(1, currentMapPosition));
    }
    
    
    //if it does not wait, compute wating time after moving
    if (!waitNow && !readyToMove) {
        double t0 = 1;
        //double t1 = 17*3600;     //range from 0 to 17 hours
		double t1 = 4000;
        double beta = -1.8;
        double randt = uRandT.value();
        powerLawT = pow(((pow(t1, beta+1) - pow(t0, beta+1))*randt + pow(t0, beta+1)), 1/(beta+1));  //waiting time
        timeCounter = powerLawT;
        waitNow = true;
        
    }
    
    //if ready to move, compute jump distance. Otherwise, wait
    if (waitNow && !readyToMove) {
        timeCounter -= UPDATE_TIME_UNIT;
        if (timeCounter <= 0) {
            waitNow = false;
            readyToMove = true;
        }
        else
            waitNow = true;
		
    }
    else {
		
        readyToMove = false;
        double x0 = 1;
        double x1 = 500;   //range from 0 to x1 meters
        double alpha = -1.55;       //n
        //alpha = -1.75;
        double rho = 0.1;
        double gamma = -0.21;
        //gamma = -0.2;
        double randr = uRandR.value();
        double powerLawR = pow(((pow(x1, alpha+1) - pow(x0, alpha+1))*randr + pow(x0, alpha+1)), 1/(alpha+1));  //jump size 1
        //field->deltaR.insert(powerLawR);
        
        
        double Pnew = rho * pow((double)userVisitedRecord.size(), gamma); //Exploration probability
		if (userVisitedRecord.empty())
		{
			cout << "userVisitedRecord is empty" << endl;
			exit(2);
		}
        double moveP = uRandM.value();
        double angle = uRandA.value();
        if (moveP <= Pnew) {
            //Explore a new location
            readyToMove = false;
            double newX = currentMapPosition.position.x + powerLawR*cos(angle*PI/180.0);
            double newY = currentMapPosition.position.y + powerLawR*sin(angle*PI/180.0);
            int loopTimes = 0;
            while ((newX < 0 || newY < 0 || newX > field->maxX || newY > field->maxY) && loopTimes < 4) {
                angle += 90;
                newX = currentMapPosition.position.x + powerLawR*cos(angle*PI/180.0);
                newY = currentMapPosition.position.y + powerLawR*sin(angle*PI/180.0);
                ++loopTimes;
            }
            
            if (newX < 0)
                newX = 0;
            if (newY < 0)
                newY = 0;
            if (newX > field->maxX)
                newX = field->maxX;
            if (newY > field->maxY)
                newY = field->maxY;
            
            Point tp(newX, newY);
            Point preLocation = currentMapPosition.position;
            currentMapPosition.position = tp;
            
            //gravity model
            int gNum = this->inGrid(GRID_SIZE);

			if (gNum < 0 || gNum >= 1000000)
			{
				cout << "Grid number is wrong" << endl;
				exit(2);
			}
            int rho_r = field->gridPool[gNum].humanNum;
            double Cr = rho_r+RHO_0;
			if (Cr <= 0)
				Cr = 1;
            
            powerLawR = powerLawR*pow(Cr, alpha);
            newX = preLocation.x + powerLawR*cos(angle*PI/180.0);
            newY = preLocation.y + powerLawR*sin(angle*PI/180.0);
            loopTimes = 0;
            while ((newX < 0 || newY < 0 || newX > field->maxX || newY > field->maxY) && loopTimes < 4) {
                angle += 90;
                newX = preLocation.x + powerLawR*cos(angle*PI/180.0);
                newY = preLocation.y + powerLawR*sin(angle*PI/180.0);
                ++loopTimes;
            }
            
            if (newX < 0)
                newX = 0;
            if (newY < 0)
                newY = 0;
            if (newX > field->maxX)
                newX = field->maxX;
            if (newY > field->maxY)
                newY = field->maxY;
             
            currentMapPosition.position.x = newX;
            currentMapPosition.position.y = newY;
            //////////////////////////////////////////////////////
            currentMapPosition.visitFrequency = 1;
            userVisitedRecord.insert(make_pair(1, currentMapPosition));
			
            updateGridInfo();
            
        }
        else {
            //Go back to some places that the user visited
            //Point tempLocation(newX, newY);
            /*
            //find the nearest visited location
            double minD = 0xffff;
            int target = 0;
            for (int i = 0; i < visitedLocations.size(); ++i) {
                double Dis = tempLocation.distance(visitedLocations[i].position);
                if (Dis < minD) {
                    minD = Dis;
                    target = i;
                }
            }
            visitedLocations[target].visitFrequency++;
            currentMapPosition = visitedLocations[target];
             */
			
            multimap<int, Maplocation>::iterator lvpos; //location visited pos
            int totalVisit = 0;
            for (lvpos = userVisitedRecord.begin(); lvpos != userVisitedRecord.end(); ++lvpos) {
                totalVisit += lvpos->first;
            }
            if (totalVisit == 0 || userVisitedRecord.empty())
            {
				cout << "totalVisit = 0" << endl;
				exit(2);
            }
            multimap<double, Maplocation> distribution;      //probability, location
            double preProb = 0;
            for (lvpos = userVisitedRecord.begin(); lvpos != userVisitedRecord.end(); ++lvpos) {
                double prob = (double)lvpos->first/totalVisit;
				Maplocation ml = lvpos->second;
                distribution.insert(make_pair(preProb+prob, ml));
                preProb += prob;
            }
			
            
            //test
            multimap<double, Maplocation>::iterator dispos = distribution.end();
            --dispos;
            if(distribution.empty() || dispos->first - 1 > 0.00001) {
				cout << "Final prob is " << dispos->first << endl;
				exit(2);
			}
            //end test
            //search the revisit location
            double revisitP = uRandV.value();
            for (dispos = distribution.begin(); dispos != distribution.end(); ++dispos) {
                if (revisitP <= dispos->first) {   //this is the target
                    if (currentMapPosition != dispos->second) {
                        currentMapPosition = dispos->second;
                        readyToMove = false;
                    
						//locate the corresponding location and update visited times
						for (lvpos = userVisitedRecord.begin(); lvpos != userVisitedRecord.end(); ++lvpos) {
							if (lvpos->second == dispos->second) {
								int value = lvpos->first;
								Maplocation m = lvpos->second;
								userVisitedRecord.erase(lvpos);
								++value;
								userVisitedRecord.insert(make_pair(value, m));
								break;
							}
						}
						break;
					}
                }
            }
			//if (TIME_UNIT_COUNTER == 11500  && this->index == 713168)
			//	cout << "Return" << endl;
            updateGridInfo();
        }
        
    }
    
}

void Human::moveModel2SingleTrack()       //(double randr, double randt, double angle, double moveP, double revisitP)
{
    //init moving
    if (visitedLocations.empty()) {
        visitedLocations.push_back(currentMapPosition);
    }
    
    //double r = uniRand.value();
    double x0 = 1;
    double x1 = 1000;   //range from 0 to 1000m
    double alpha = -1.55;
    //alpha = -1.75;
    double rho = 0.1;
    double gamma = -0.21;
    //gamma = -0.2;
    double randr = uRandR.value();
    double powerLawR = pow(((pow(x1, alpha+1) - pow(x0, alpha+1))*randr + pow(x0, alpha+1)), 1/(alpha+1));  //jump size
    double Pnew = rho * pow((double)visitedLocations.size(), gamma); //Exploration probability
    
    //double t = uniRand.value();
    double t0 = 1;
    double t1 = 17*3600;     //range from 0 to 17 hours
    double beta = -1.8;
    double randt = uRandT.value();
    powerLawT = pow(((pow(t1, beta+1) - pow(t0, beta+1))*randt + pow(t0, beta+1)), 1/(beta+1));  //waiting time
    
    
    
    //double moveP = uniRand.value();
    //double newX = currentMapPosition.position.x + powerLawR*cos(angle*PI/180);
    //double newY = currentMapPosition.position.y + powerLawR*sin(angle*PI/180);
    double moveP = uRandM.value();
    double angle = uRandA.value();
    if (moveP <= Pnew) {
        //Explore a new location
        double newX = currentMapPosition.position.x + powerLawR*cos(angle*PI/180.0);
        double newY = currentMapPosition.position.y + powerLawR*sin(angle*PI/180.0);
        int loopTimes = 0;
        while ((newX < 0 || newY < 0 || newX > field->maxX || newY > field->maxY) && loopTimes < 4) {
            angle += 90;
            newX = currentMapPosition.position.x + powerLawR*cos(angle*PI/180.0);
            newY = currentMapPosition.position.y + powerLawR*sin(angle*PI/180.0);
            ++loopTimes;
        }
        
        if (newX < 0)
            newX = 0;
        if (newY < 0)
            newY = 0;
        if (newX > field->maxX)
            newX = field->maxX;
        if (newY > field->maxY)
            newY = field->maxY;
        
        Point tp(newX, newY);
        Point preLocation = currentMapPosition.position;
        currentMapPosition.position = tp;
        
        //gravity model
        int gNum = this->inGrid(GRID_SIZE);
        
        int rho_r = field->gridPool[gNum].humanNum;
        double Cr = rho_r+RHO_0;
        if (Cr <= 0)
            Cr = 1;
        
        powerLawR = powerLawR*pow(Cr, alpha);
        newX = preLocation.x + powerLawR*cos(angle*PI/180.0);
        newY = preLocation.y + powerLawR*sin(angle*PI/180.0);
        loopTimes = 0;
        while ((newX < 0 || newY < 0 || newX > field->maxX || newY > field->maxY) && loopTimes < 4) {
            angle += 90;
            newX = preLocation.x + powerLawR*cos(angle*PI/180.0);
            newY = preLocation.y + powerLawR*sin(angle*PI/180.0);
            ++loopTimes;
        }
        
        if (newX < 0)
            newX = 0;
        if (newY < 0)
            newY = 0;
        if (newX > field->maxX)
            newX = field->maxX;
        if (newY > field->maxY)
            newY = field->maxY;
        
        currentMapPosition.position.x = newX;
        currentMapPosition.position.y = newY;
        //////////////////////////////////////////////////////
        currentMapPosition.visitFrequency = 1;
        userVisitedRecord.insert(make_pair(1, currentMapPosition));
        visitedLocations.push_back(currentMapPosition);
        
        updateGridInfo();
        
    }
    else {
        //Go back to some places that the user visited
        multimap<int, Maplocation>::iterator lvpos; //location visited pos
        int totalVisit = 0;
        for (lvpos = userVisitedRecord.begin(); lvpos != userVisitedRecord.end(); ++lvpos) {
            totalVisit += lvpos->first;
        }
        
        multimap<double, Maplocation> distribution;      //probability, location
        double preProb = 0;
        for (lvpos = userVisitedRecord.begin(); lvpos != userVisitedRecord.end(); ++lvpos) {
            double prob = (double)lvpos->first/totalVisit;
            distribution.insert(make_pair(preProb+prob, lvpos->second));
            preProb += prob;
        }
        
        //test
        multimap<double, Maplocation>::iterator dispos = distribution.end();
        --dispos;
        assert(dispos->first - 1 <= 0.00001);
        //end test
        //search the revisit location
        double revisitP = uRandV.value();
        for (dispos = distribution.begin(); dispos != distribution.end(); ++dispos) {
            if (revisitP <= dispos->first) {   //this is the target
                currentMapPosition = dispos->second;
                //for updating visitedLocations to output St and fk
                for (int i = 0; i < visitedLocations.size(); ++i)
                {
                    if (visitedLocations[i] == currentMapPosition)
                    {
                        visitedLocations[i].visitFrequency++;
                    }
                }
                //locate the corresponding location and update visited times
                for (lvpos = userVisitedRecord.begin(); lvpos != userVisitedRecord.end(); ++lvpos) {
                    if (lvpos->second == dispos->second) {
                        int value = lvpos->first;
                        Maplocation m = lvpos->second;
                        userVisitedRecord.erase(lvpos);
                        ++value;
                        userVisitedRecord.insert(make_pair(value, m));
                        break;
                    }
                }
                break;
            }
        }
        
        
    }
    timePassed += powerLawT;
    waitingTimeSequence.push_back(make_pair(timePassed, visitedLocations.size()));
}

void Human::normalizeLocations()
{
    double totalFreqeucy = 0;
    for (int i = 0; i < visitedLocations.size(); ++i) {
        totalFreqeucy += visitedLocations[i].visitFrequency;
    }
    
    for (int i = 0; i < visitedLocations.size(); ++i) {
        normalizedFrequencyLocations.insert(make_pair((double)visitedLocations[i].visitFrequency/totalFreqeucy, visitedLocations[i].position));
    }
    
}

int Human::inGrid(double gridSize)
{
	int i = (int)currentMapPosition.position.x/gridSize;
	int j = (int)currentMapPosition.position.y/gridSize;
	if (i < 0 || j < 0) {
		cout << "inGrid index wrong" << endl;
		exit(2);
	}
	if (i == field->maxX)
		i = field->maxX-1;
	if (j == field->maxY)
		j = field->maxY-1;
	int result = j*(field->maxX/gridSize)+i;
    return result;
}

void Human::updateGridInfo()
{
	Grid *preGrid = NULL;
	Grid *currentGrid = NULL;
    preGrid = &field->gridPool[gridNum];
	if (this->gridNum < 0 || this->gridNum >= 1000000)
	{
		cout << "Grid number is wrong" << endl;
		exit(2);
	}
    //field->gridPool[gridNum].humanNum--;            //decrease the human number of the previous grid
    //assert(field->gridPool[gridNum].humanNum >= 0);
    //if (field->gridPool[gridNum].humanNum == 0) {
    //    field->gridPool[gridNum].hasHuman = false;
    //}
    
    this->gridNum = this->inGrid(GRID_SIZE);        //locate to the current grid number
	if (this->gridNum < 0 || this->gridNum >= 1000000)
	{
		cout << "Grid number is wrong" << endl;
		exit(2);
	}
    
    currentGrid = &field->gridPool[gridNum];
    //field->gridPool[gridNum].humanNum++;            //add the human number of the current grid
    //if (!field->gridPool[gridNum].hasHuman) {       //if there was no human, change the state now
    //    field->gridPool[gridNum].hasHuman = true;
    //}
    grid = &field->gridPool[gridNum];               //update the grid pointer
    
    if (preGrid->index != currentGrid->index) {
		preGrid->humanNum--;
		if (preGrid->humanNum < 0) {
			cout << "humanNum is wrong" << endl;
			exit(2);
		}
		if (preGrid->humanNum == 0) {
			preGrid->hasHuman = false;
		}
		currentGrid->humanNum++;
		if (!currentGrid->hasHuman) {       //if there was no human, change the state now
			currentGrid->hasHuman = true;
		}
        jumpGridRecord.push_back(make_pair(preGrid, currentGrid));
	}
}

void Human::computeRGT()
{
    //double maxRGT = 0;
    for (int i = 0; i < waitingTimeSequence.size(); ++i) {
        double time_ = waitingTimeSequence[i].first;
        int visitIndex = waitingTimeSequence[i].second;     //the number of visited locations
        //apply the equation to calculate rg_t
        Point r_cm(0, 0);
        for (int j = 0; j < visitIndex; ++j) {
            Point p = visitedLocations[j].position;
            r_cm.x += p.x;
            r_cm.y += p.y;
        }
        r_cm.x = r_cm.x/visitIndex;
        r_cm.y = r_cm.y/visitIndex;
        
        double dist = 0;
        for (int j = 0; j < visitIndex; ++j) {
            Point p = visitedLocations[j].position;
            dist += p.distance2(r_cm);      //distance square
        }
        dist = dist/visitIndex;
        double rgt = sqrt(dist);
        rg_t.push_back(make_pair(time_, rgt));
        if (maxRGT < rgt) {
            maxRGT = rgt;
        }
    }
    //Normalize
    //for (int i = 0; i < rg_t.size(); ++i) {
    //    rg_t[i].second = rg_t[i].second/maxRGT;
    //}
}
