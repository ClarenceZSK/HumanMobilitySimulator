#include "geometry.h"
#include "common.h"
#include <algorithm>
#include "assert.h"


Point PointTransform(Point pt, double degree)
{
	double dist = pt.distance(Point(0,0));
	
	double newAng = atan(pt.y/pt.x) - (degree - PI/2.0);

	double newX = dist * cos(newAng);
	double newY = dist * sin(newAng);

	return Point(newX, newY);
}

Point PointTransformBack(Point pt, double degree)
{
	
	double dist = pt.distance(Point(0,0));
	double oldAng = atan2(pt.y, pt.x);
	if (oldAng < 0)
		oldAng += 2 * PI;
	
	double newAng = oldAng + (degree - PI/2.0);

	double newX = dist * cos(newAng);
	double newY = dist * sin(newAng);

	return Point(newX, newY);
	
}

int pnpoly_C(int npol, double *xp, double *yp, double x, double y)
{
	int i, j, c = 0;
	for (i = 0, j = npol-1; i < npol; j = i++) {
        if ((((yp[i]<=y) && (y<yp[j])) ||
			((yp[j]<=y) && (y<yp[i]))) &&
            (x < (xp[j] - xp[i]) * (y - yp[i]) / (yp[j] - yp[i]) + xp[i]))
			
			c = !c;
	}
	return c;
}


/////////////////////////////// Line segments /////////////////////////
LineSegment::LineSegment(Point start_, Point end_)
{
	start = start_;
	end = end_;
	a = end.y - start.y;
	b = -(end.x - start.x);
	c = -(start.x * end.y - end.x * start.y);
}

LineSegment::LineSegment(Point start_, double angle_)
{
	start = start_;
	
	end.x = start.x + 10 * cos(angle_);
	end.y = start.y + 10 * sin(angle_);
	
	a = end.y - start.y;
	b = -(end.x - start.x);
	c = -(start.x * end.y - end.x * start.y);
}

LineSegment::LineSegment()
{
	start = INVALID_POINT;
	end = INVALID_POINT;
}

LineSegment::~LineSegment()
{
}

void LineSegment::setEnds(Point start_, Point end_)
{
	start = start_;
	end = end_;
	a = end.y - start.y;
	b = -(end.x - start.x);
	c = -(start.x * end.y - end.x * start.y);
}

bool LineSegment::isHorizontal()
{
	return (start.y == end.y);
}

bool LineSegment::isVertical()
{
	return (start.x == end.x);
}

bool LineSegment::containPointOpen(Point pt)
{
	double zero = a * pt.x + b * pt.y + c;

	double tmp = DOUBLE_EQUAL(zero, 0);
	return (DOUBLE_EQUAL(zero, 0));
}

double LineSegment::distanceToPoint(Point pt)
{

	double x1 = start.x;
	double y1 = start.y;
	double x2 = end.x;
	double y2 = end.y;
	double x0 = pt.x;
	double y0 = pt.y;

	double res = (x2-x1)*(y1-y0) - (x1-x0)*(y2-y1);
	res = fabs(res);

	res /= sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1));

	return res;

	return 0;
}

Point LineSegment::projectionOfPoint(Point pt)
{
	if(containPointOpen(pt))
		return pt;

	double u = (pt.x - start.x) * (end.x - start.x) + (pt.y - start.y) * (end.y - start.y);
	u /= ((start.x - end.x) * (start.x - end.x) + (start.y - end.y) * (start.y - end.y));

	double x = start.x + u * (end.x - start.x);
	double y = start.y + u * (end.y - start.y);

	return Point(x,y);
}

bool LineSegment::overlap(LineSegment & line)
{
	
	if (start == line.start) {
		if (containPointClose(line.end) || line.containPointClose(end)) {
			return true;
		}
		else {
			return false;
		}
	}

	if (start == line.end) {
		if (containPointClose(line.start) || line.containPointClose(end)) {
			return true;
		}
		else {
			return false;
		}
	}

	
	if (end == line.start) {
		return (line.containPointClose(start) || containPointClose(line.end));
	}

	if (end == line.end) {
		return (line.containPointClose(start) || containPointClose(line.start));
	}
	

	return ((containPointClose(line.start) || containPointClose(line.end)) &&
		containPointOpen(Point::getMiddlePoint(line.start, line.end)));

}

bool LineSegment::containEnd(Point pt)
{
	return (start.approxEqual(pt) || end.approxEqual(pt));
}

bool LineSegment::containPointClose(Point pt)
{
	double minx = min(start.x, end.x);
	double maxx = max(start.x, end.x);
	double miny = min(start.y, end.y);
	double maxy = max(start.y, end.y);

	if (containPointOpen(pt) && 
		(pt.x > minx || DOUBLE_EQUAL(pt.x, minx)) && 
		(pt.x < maxx || DOUBLE_EQUAL(pt.x, maxx)) &&
		(pt.y > miny || DOUBLE_EQUAL(pt.y, miny)) && 
		(pt.y < maxy || DOUBLE_EQUAL(pt.y, maxy))) 
		return true;
	else
		return false;
}


bool LineSegment::containPointExclusiveEnds(Point pt)
{
	double minx = min(start.x, end.x);
	double maxx = max(start.x, end.x);
	double miny = min(start.y, end.y);
	double maxy = max(start.y, end.y);

	if (containPointOpen(pt) && 
		(pt.x > minx || DOUBLE_EQUAL(pt.x, minx)) && 
		(pt.x < maxx || DOUBLE_EQUAL(pt.x, maxx)) && 
		(pt.y > miny || DOUBLE_EQUAL(pt.y, miny)) && 
		(pt.y < maxy || DOUBLE_EQUAL(pt.y, maxy)) && 
		pt.approxInequal(start) && 
		pt.approxInequal(end)) 
		return true;
	else
		return false;
}

Point LineSegment::intersectionOpen(LineSegment &line)
{
	double xx = a * line.b - b * line.a;
	if (xx == 0) { //parallel lines! 
		return INVALID_POINT;
	}
	else 
		return Point(-(c*line.b-line.c*b)/xx, (c*line.a-line.c*a)/xx);
}

//including ends
Point LineSegment::intersectionClose(LineSegment &line)
{
	Point pt = intersectionOpen(line);
	if(line.containPointClose(pt) && this->containPointClose(pt))
		return pt;
	else
		return INVALID_POINT;
}


Point LineSegment::intersectionExclusiveEnds(LineSegment &line)
{

	double xx = a * line.b - b * line.a;
	if (xx == 0) { //parallel lines! 
		return INVALID_POINT;
	}
	//share an end
	else if (start == line.start || start == line.end || end == line.start
	|| end == line.end){
		return INVALID_POINT;
	}
	else {
		double x = -(c*line.b-line.c*b)/xx;
		double y = (c*line.a-line.c*a)/xx;

		Point interp(x,y);

		if (containPointExclusiveEnds(interp) &&
			line.containPointExclusiveEnds(interp)) 
			return interp;

		if (containPointExclusiveEnds(interp) && 
			(interp.approxEqual(line.start) || interp.approxEqual(line.end)))
			return interp;

		if (line.containPointExclusiveEnds(interp) && 
			(interp.approxEqual(start) || interp.approxEqual(end)))
			return interp;
		
	}
	return INVALID_POINT;
}

double LineSegment::angle()
{

  double x1 = start.x;
  double y1 = start.y;
  double x2 = end.x;
  double y2 = end.y;

  double line_len = sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1));

  double sin_theta, cos_theta;
  double theta;

  if(line_len == 0.0){
    fprintf(stderr, "the two nodes are the same\n");
    return -1.0;
  }

  sin_theta = (y2-y1)/line_len;
  cos_theta = (x2-x1)/line_len;

  theta = acos(cos_theta);
  
  if(sin_theta<0){
    theta = 2*PI - theta;
  }

  return theta;
}

double LineSegment::getX(double y)
{
	assert(a != 0);
	return - (b * y + c) / a;
}

double LineSegment::getY(double x)
{
	if (!b) {
		int stop = 1;
	}
	assert(b);

	return - (a * x + c) / b;
}

LineSegment LineSegment::getPerpendicularLineSegment(Point pt)
{
	assert(containPointOpen(pt));

	double ang = angle();
	return LineSegment(pt, ang + PI / 2);
}

vector<Point> LineSegment::getOnlinePointsGivenDistanceToOnlinePoint(Point pt, double dist)
{
	vector<Point> res;
	double ang = angle();
	Point p1, p2;
	p1.x = pt.x + dist * cos(ang);
	if (DOUBLE_EQUAL(p1.x, pt.x) || b == 0) {
		p1.y = pt.y + dist;
	}
	else {
		p1.y = getY(p1.x);
	}
	p2.x = pt.x - dist * cos(ang);
	if (DOUBLE_EQUAL(p2.x, pt.x) || b == 0) {
		p2.y = pt.y - dist;
	}
	else {
		p2.y = getY(p2.x);
	}
	res.push_back(p1);
	res.push_back(p2);
	return res;
}

Point LineSegment::getOnlinePoint(Point start, Point end, double dist)
{
	assert(containPointOpen(end));

	vector<Point> pts = getOnlinePointsGivenDistanceToOnlinePoint(start, dist);

	if (pts[0].distance(end) < pts[1].distance(end)) 
		return pts[0];
	else
		return pts[1];
}

bool LineSegment::betweenTwoAngle(double ang1, double ang2) //counter-clockwise!
{
	double curAng = angle();
	if (ang1 < ang2) {
		if (ang1 <= curAng && curAng <= ang2) {
			return true;
		}
	}
	else { //ang1 > ang2
		if ((ang1 <= curAng && curAng <= 2*PI) ||
			(0 <= curAng && curAng <= ang2)) {
			return true;
		}
	}

	return false;
}

int LineSegment::relativePosition(Point pt)
{
	if (containPointOpen(pt)) 
		return PS_ONLINE;

	double a = angle();
	double a1 = (LineSegment(start, pt)).angle();
				
	a1 -= a;
				
	if (a1 < 0) a1 += 2*PI;
	if (a1 < PI)
		return PS_TOLEFT;
	else
		return PS_TORIGHT;
}

void MyPolygon::print()
{
	for(unsigned int i = 0; i < vertices.size(); i++) {
		Point pt = vertices[i];
		printf("(%lf, %lf) \n", pt.x, pt.y);
	}
	printf("\n");
}

double MyPolygon::calcArea()
{
	int i, j;
	double area = 0;
	int N = vertices.size();
	
	for (i = 0; i < N; i++) {
		j = (i + 1) % N;
		area += vertices[i].x * vertices[j].y;
		area -= vertices[i].y * vertices[j].x;
	}
	
	area /= 2;
	return(area < 0 ? -area : area);

}

bool MyPolygon::isConvex()
{
	return (maxInnerAngle() < PI);
}

//does not include boundary
bool MyPolygon::containPointInside(Point p)
{

	//this is not needed theoretically, here it
	//is for accounting for the rounding error in calculations
	//of double type
	if (p == INVALID_POINT || containVertice(p)) {
		return false;
	}

	unsigned int i, j;
	int c = 0;
	for (i = 0, j = vertices.size()-1; i < vertices.size(); j = i++) {
		Point pi = vertices[i];
		Point pj = vertices[j];
        if ((((pi.y <= p.y) && (p.y < pj.y)) ||
			((pj.y <= p.y) && (p.y < pi.y))) &&
            (p.x < (pj.x - pi.x) * (p.y - pi.y) / (pj.y - pi.y) + pi.x))
			
			c = !c;
	}
	
	return (c == 1);	
}

bool MyPolygon::containPointIncludingBoundary(Point pt)
{
	return (containPointInside(pt) || containPointOnBoundary(pt));
}

double MyPolygon::innerAngle(Point pt)
{
	bool found = false;

	for(unsigned int i = 0; i < vertices.size(); i++) {
		if(pt == vertices[i]) {
			found = true;

			int start = (i - 1 + vertices.size()) % vertices.size();
			int end = (i + 1) % vertices.size();

			LineSegment e1(pt, vertices[start]);
			LineSegment e2(pt, vertices[end]);

			double angleDiff = e2.angle() - e1.angle();
			if (angleDiff <= 0) {
				angleDiff += 2 * PI;
			}

			return angleDiff;
		}
	}

	if (!found) {
		fprintf(stderr, "have not found the vertex for calulating inner angle\n");
		exit(0);
	}

	return -1;
}

int MyPolygon::maxInnerAngleVertex()
{
	double maxAngle = 0;
	double maxAngleIdx = 0;

	for(unsigned int i = 0; i < vertices.size(); i++) {
		double ang = innerAngle(vertices[i]);
		if (ang > maxAngle) {
			maxAngle = ang;
			maxAngleIdx = i;
		}
	}

	return (int)maxAngleIdx;
}

double MyPolygon::maxInnerAngle()
{
	double maxAngle = 0;

	for(unsigned int i = 0; i < vertices.size(); i++) {
		double ang = innerAngle(vertices[i]);
		if (ang > maxAngle) {
			maxAngle = ang;
		}
	}

	return maxAngle;
}

int MyPolygon::minInnterAngleVertex()
{
	assert(0);
	return 0;
}

void MyPolygon::getIncidentTwoEdges(int verIdx, LineSegment & e1, LineSegment & e2)
{
	assert(verIdx >= 0 && (unsigned) verIdx < vertices.size());
	
	int start = (verIdx - 1 + vertices.size()) % vertices.size();
	int end = (verIdx + 1) % vertices.size();

	e1.setEnds(vertices[start], vertices[verIdx]);
	e2.setEnds(vertices[verIdx], vertices[end]);
	
}

bool MyPolygon::containVertice(Point pt)
{
	//return	(find(vertices.begin(), vertices.end(), pt) != vertices.end());
	for(unsigned int i = 0; i < vertices.size(); i++) {
		//if (pt.approxEqual(vertices[i])) {
		if(pt == vertices[i]) {
			return true;
		}
	}

	return false;
}

bool MyPolygon::containEdge(Point p1, Point p2)
{
	for(unsigned int j = 0; j < vertices.size(); j++) {
		Point v1, v2;
		
		v1 = vertices[j];
		v2 = vertices[(j+1) % vertices.size()];
		
		if ((p1 == v1 && p2 == v2) || (p1 == v2 && p2 == v1))
			return true;
	}

	return false;
}

//including vertices
bool MyPolygon::containPointOnBoundary(Point pt)
{

	for(unsigned int i = 0; i < vertices.size(); i++) {
		int j = (i+1) % vertices.size();
		LineSegment ls(vertices[i], vertices[j]);
		if (ls.containPointClose(pt)) 
			return true;
	}

	return false;
}

vector<LineSegment> MyPolygon::getEdges() //clockwise
{
	vector<LineSegment> ls;

	for(unsigned int i = 0; i < vertices.size(); i++) {
		int j = (i+1) % vertices.size();
		ls.push_back(LineSegment(vertices[i], vertices[j]));
	}

	return ls;
}


Point MyPolygon::getNeighbor(Point pt, bool clockwise)
{
	vector<Point>::const_iterator iter = find(vertices.begin(), vertices.end(), pt);
	if (iter == vertices.end()) 
		return INVALID_POINT;

	if (clockwise) {
		iter++;
		if (iter == vertices.end()) 
			return *(vertices.begin());
		else
			return *iter;
	}
	else {
		if (iter == vertices.begin()) 
			return vertices[vertices.size() - 1];
		else {
			--iter;
			return *iter;
		}
	}
}

vector<Point> MyPolygon::getNeighborVertices(Point pt)
{
	vector<Point> neighbors;
	if (!containVertice(pt)) 
		return neighbors;

	for(unsigned int i = 0; i < vertices.size(); i++) {
		if (pt == vertices[i]) {
			if (i == 0) {
				neighbors.push_back(vertices[vertices.size()-1]);
				neighbors.push_back(vertices[1]);
			}
			else if (i == vertices.size() - 1) {				
				neighbors.push_back(vertices[i-1]);
				neighbors.push_back(vertices[0]);
			}
			else {				
				neighbors.push_back(vertices[i-1]);
				neighbors.push_back(vertices[i+1]);
			}
		}
	}

	return neighbors;
}

vector<pair<Point, Point>*> MyPolygon::getNeigbhorsPairs(Point pt)
{
	vector<pair<Point, Point>*> neighborPairs;

	pair<Point, Point> ngh;
	//neighbors.first = INVALID_POINT; //ccw ngh
	//neighbors.second = INVALID_POINT; //cw ngh


	/*
	if (pt == vertices[0]) {
		neighbors.first = vertices[1];
		neighbors.second = vertices[vertices.size()-2];
		return neighbors;
	}*/

	for(unsigned int i = 0; i < vertices.size(); i++) {
		if (vertices[i] == pt) {
			pair<Point, Point> *ngh = new pair<Point, Point>;
			ngh->first = (i == 0 ? vertices[vertices.size()-2] : vertices[i-1]);
			ngh->second = (i == vertices.size() - 1 ? vertices[1] : vertices[i+1]);
			neighborPairs.push_back(ngh);
		}
	}

	return neighborPairs;
}

//including intersection at vertices
//vector<Point> 
bool MyPolygon::intersectWithLineSeg(LineSegment line)
{
	if (containEdge(line.start, line.end))
		return false;

	for(unsigned int i = 0; i < vertices.size(); i++) {
		int j = (i+1) % vertices.size();
		LineSegment ls(vertices[i], vertices[j]);
		if (ls.intersectionClose(line) != INVALID_POINT ||
			line.containPointExclusiveEnds(vertices[i]) ||
			line.containPointExclusiveEnds(vertices[j])) {
			return true;
		}
	}

	return false;
}


bool MyPolygon::crossWithLine(LineSegment line)
{
	if (containEdge(line.start, line.end))
		return false;

	unsigned int i, j;
	vector<Point> inters;

	bool startPointInside = containPointInside(line.start);
	bool endPointInside = containPointInside(line.end);

	if (startPointInside && endPointInside) {
		fprintf(stderr, "both points are in the polygon!");
		exit(0);
	}
	
	if (startPointInside || endPointInside) 
		return true;
	
	for(i = 0; i < vertices.size(); i++) {
		j = (i+1) % vertices.size();
		LineSegment ls(vertices[i], vertices[j]);
		Point inter = ls.intersectionClose(line);
		
		if (containVertice(inter)) {
			if (find(inters.begin(), inters.end(), inter) == inters.end()) {
				inters.push_back(inter);
			}
		}
		else if (inter != INVALID_POINT/* && inter.approxInequal(line.start)
			&& inter.approxInequal(line.end)*/) {
			return true;
		}
	}

	if (inters.size() < 2)
		return false;

	for(i = 0; i < inters.size() - 1; i++) {
		Point midp((inters[i].x + inters[i+1].x)/2, (inters[i].y + inters[i+1].y)/2);
		if (containPointInside(midp)) {
			return true;
		}
	}

	return false;
}


bool MyPolygon::crossWithLine1(LineSegment line)
{
	int num = 200;
	for(int i = 0; i <= num; i++) {
		Point pt;
		pt.x = line.start.x + (line.end.x - line.start.x) * double(i) / double(num);
		pt.y = line.start.y + (line.end.y - line.start.y) * double(i) / double(num);

		if (containPointInside(pt)) {
			return true;
		}
	}

	return false;
}

bool MyPolygon::verticeInSight(Point watcher, Point ver)
{
	if (!containVertice(ver))
		return false;

	if (containPointOnBoundary(watcher)) {

		//is watcher on the edge incident on ver?
		for(unsigned int i = 0; i < vertices.size(); i++) {
			int j = (i+1) % vertices.size();
			if (vertices[i] == ver || vertices[j] == ver) {
				LineSegment line(vertices[i], vertices[j]);
				if (line.containPointClose(watcher)) 
					return true;
			}
		}		
		
		Point midPoint((watcher.x+ver.x)/2.0, (watcher.y+ver.y)/2.0);
		if (containPointInside(midPoint)) 
			return false;	
	}
	
	LineSegment ls(watcher, ver);
	
	for(unsigned int i = 0; i < vertices.size(); i++) {
		int j = (i+1) % vertices.size();
		if (vertices[i] == ver) {
			if (ls.containPointExclusiveEnds(vertices[j])) 
				return false;
			else 
				continue;
		}
		
		if (vertices[j] == ver) {
			if (ls.containPointExclusiveEnds(vertices[i])) 
				return false;
			else 
				continue;
		}
		
		LineSegment edge(vertices[i], vertices[j]);
		Point inter = edge.intersectionClose(ls);
		if (inter != INVALID_POINT && inter.approxInequal(watcher)) 
			return false;			
	}	

	return true;
}

bool MyPolygon::forwardWithRightHandRule(Point p1, Point p2)
{
	if (!containEdge(p1, p2)) {
		return true; //default
	}
	else {
		
		for(unsigned int j = 0; j < vertices.size(); j++) {
			Point v1, v2;
			
			v1 = vertices[j];
			v2 = vertices[(j+1) % vertices.size()];
			
			if ((p1 == v1 && p2 == v2) || (p1 == v2 && p2 == v1))
			{				
				if (p1 == v1 && p2 == v2)
					return true;
				else
					return false;
			}
		}
	}

	return true;
}

Point MyPolygon::findAnInternalPoint()
{
	if (vertices.size() == 3) {
		double midx = (vertices[0].x + vertices[1].x + 2 * vertices[2].x) / 4.0;
		double midy = (vertices[0].y + vertices[1].y + 2 * vertices[2].y) / 4.0;
		return Point(midx, midy);
	}
	else { // > 3
		for(unsigned int i = 0; i < vertices.size(); i++) {
			Point p1 = vertices[i];

			int j = (i + 2) % vertices.size();
			
			Point p2 = vertices[j];
			
			Point midp((p1.x+p2.x)/2.0, (p1.y+p2.y)/2.0);
			if (containPointInside(midp)) 
				return midp;
		}
	}

	return INVALID_POINT;

}

void MyPolygon::dividingByLine(LineSegment &ls, MyPolygon & leftPolygon, MyPolygon &rightPolygon)
{

	for(unsigned int i = 0; i < vertices.size(); i++) {
		Point ver = vertices[i];
		Point nextVer = (i == vertices.size() - 1) ?  vertices[0] : vertices[i+1];

		int pos = ls.relativePosition(ver);
		int nextPos = ls.relativePosition(nextVer);

		if (pos == PS_TOLEFT) {
			leftPolygon.vertices.push_back(ver);
		}
		else if (pos == PS_TORIGHT) {
			rightPolygon.vertices.push_back(ver);
		}
		else { //PS_ONLINE
			leftPolygon.vertices.push_back(ver);
			rightPolygon.vertices.push_back(ver);
			continue;
		}

		if (nextPos != PS_ONLINE && pos != nextPos) {
			LineSegment edge(ver, nextVer);
			Point inter = edge.intersectionOpen(ls);

			if (inter.approxInequal(nextVer)) {
				
				leftPolygon.vertices.push_back(inter);
				rightPolygon.vertices.push_back(inter);
			}
		}
	}
}

void MyPolygon::getContainingPartialPolygon(Point pt, LineSegment &ls, MyPolygon &polygon)
{
	polygon.vertices.clear();

	MyPolygon leftPolygon, rightPolygon;

	for(unsigned int i = 0; i < vertices.size(); i++) {
		Point ver = vertices[i];
		Point nextVer = (i == vertices.size() - 1) ?  vertices[0] : vertices[i+1];

		int pos = ls.relativePosition(ver);
		int nextPos = ls.relativePosition(nextVer);

		if (pos == PS_TOLEFT) {
			leftPolygon.vertices.push_back(ver);
		}
		else if (pos == PS_TORIGHT) {
			rightPolygon.vertices.push_back(ver);
		}
		else { //PS_ONLINE
			leftPolygon.vertices.push_back(ver);
			rightPolygon.vertices.push_back(ver);
			continue;
		}

		if (nextPos != PS_ONLINE && pos != nextPos) {
			LineSegment edge(ver, nextVer);
			Point inter = edge.intersectionOpen(ls);

			if (inter.approxInequal(nextVer)) {
				
				leftPolygon.vertices.push_back(inter);
				rightPolygon.vertices.push_back(inter);
			}
		}
	}

	if (leftPolygon.containPointIncludingBoundary(pt)) {
		polygon = leftPolygon;
	}
	else {
		polygon = rightPolygon;
	}
}

vector<Point> movePoints(vector<Point> pts, double dx, double dy)
{
	vector<Point> points;
	for(unsigned int i = 0; i < pts.size(); i++) {
		points.push_back(Point(pts[i].x + dx, pts[i].y + dy));
	}

	return points;
}

vector<Point> makeSquareHole(Point topleft, double width)
{
	vector<Point> points;

	points.push_back(topleft);
	points.push_back(Point(topleft.x + width, topleft.y));
	points.push_back(Point(topleft.x + width, topleft.y - width));
	points.push_back(Point(topleft.x, topleft.y - width));

	return points;
}


vector<Point> makeRectangleHole(Point topleft, double width, double height)
{
	vector<Point> points;

	points.push_back(topleft);
	points.push_back(Point(topleft.x + width, topleft.y));
	points.push_back(Point(topleft.x + width, topleft.y - height));
	points.push_back(Point(topleft.x, topleft.y - height));

	return points;
}



vector<Point> makeTriangle(Point p1, Point p2, Point p3)
{
	
	vector<Point> points;

	points.push_back(p1);
	points.push_back(p2);
	points.push_back(p3);

	return points;
}