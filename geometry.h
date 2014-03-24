#ifndef GEOMETRY_H_
#define GEOMETRY_H_


#include <vector>
#include "common.h"
#include "math.h"
#include "stdio.h"

using namespace std;

//does not include boundary!

extern int pnpoly_C(int npol, double *xp, double *yp, double x, double y);

class Point
{
public:
	Point() {x = 0; y = 0; z = 0;}
	Point(double x_, double y_, double z_) : x(x_), y(y_), z(z_){};
	Point(double x_, double y_) : x(x_), y(y_), z(0){};
	~Point(){};

	double X() {return x;}
	double Y() {return y;}
	double Z() {return z;}

	double x;
	double y;
	double z;

	inline double distance(Point P) { 
		double xDist = pow(P.x - x, 2);
		double yDist = pow(P.y - y, 2);
		double zDist = pow(P.z - z, 2);
		return sqrt(xDist + yDist + zDist);
	};
    
    inline double distance2(Point P) {
		double xDist = pow(P.x - x, 2);
		double yDist = pow(P.y - y, 2);
		double zDist = pow(P.z - z, 2);
		return xDist + yDist + zDist;
	};

	inline double distance(double X, double Y, double Z) { 		
		return distance(Point(X, Y, Z));
	};


	inline Point & operator = (const Point &right)
	{
		x = right.x;
		y = right.y;
		z = right.z;

		return *this;
	}

	inline bool operator == (const Point &right) const
	{
		return (x == right.x && y == right.y && z == right.z);
	}

	inline bool operator != (const Point &right) const
	{
		return (x != right.x || y != right.y || z != right.z);
	}

	inline bool approxInequal(const Point &right)
	{
		return (!DOUBLE_EQUAL(x, right.x) || !DOUBLE_EQUAL(y, right.y) || !DOUBLE_EQUAL(z, right.z));
	}

	inline bool approxEqual(const Point &right)
	{
		return (DOUBLE_EQUAL(x, right.x) && DOUBLE_EQUAL(y, right.y) && DOUBLE_EQUAL(z, right.z));
	}

	inline bool isValid()
	{
		return (x != IMPOSSIBLE_VALUE && y != IMPOSSIBLE_VALUE && z != IMPOSSIBLE_VALUE);
	}

	static Point getMiddlePoint(const Point & p1, const Point & p2) 
	{
		return Point((p1.x+p2.x)/2.0, (p1.y+p2.y)/2.0, (p1.z + p2.z)/2.0);
	}
};


//(directed) line segement
//(y2-y1)*x - (x2-x1)*y - (x1*y2-x2*y1) = 0
class LineSegment
{
public:
	LineSegment(Point start_, Point end_);
	LineSegment(Point start_, double angle_);
	LineSegment();
	~LineSegment();
	
	Point start;
	Point end;
	double a;  //end.y - start.y
	double b;  //-(end.x - start.x)
	double c;  //-(start.x * end.y - end.x * start.y)

	void setEnds(Point start_, Point end_);

	double getX(double y);
	double getY(double x);
	vector<Point> getOnlinePointsGivenDistanceToOnlinePoint(Point pt, double dist);
	Point getOnlinePoint(Point start, Point end, double dist);
	LineSegment getPerpendicularLineSegment(Point pt);
	bool isHorizontal();
	bool isVertical();
	double distanceToPoint(Point pt);
	Point projectionOfPoint(Point pt);
	bool overlap(LineSegment & line); //any overlap (other than endpoints) between the two lines?
	bool containEnd(Point pt); //pt is one of the ends
	bool containPointClose(Point pt); //including ends
	bool containPointExclusiveEnds(Point pt);
	bool containPointOpen(Point pt);
	Point intersectionClose(LineSegment &line); //including ends
	Point intersectionExclusiveEnds(LineSegment &line); 
	Point intersectionOpen(LineSegment &line);

	int relativePosition(Point pt);
	double angle();

	bool betweenTwoAngle(double ang1, double ang2); //right hand rule!

	bool operator == (const LineSegment &ls) 
	{
		return (a == ls.a && b == ls.b && c == ls.c);
	}

	bool isCollinearHorV(LineSegment line) {
		if (isHorizontal() && line.isHorizontal() && start.y == line.start.y) {
			return true;
		}
		else if (isVertical() && line.isVertical() && start.x == line.start.x) {
			return true;
		}
		else
			return false;
	}
};


//NOTE: vertices must be stored in clockwise direction!
class MyPolygon
{
public:

	MyPolygon(vector<Point> vertices_) {vertices = vertices_;}
	MyPolygon() {}
	~MyPolygon() {}

	void print();
	double calcArea();

	//return intersection points
	//vector<Point> 
	bool isConvex();
	bool intersectWithLineSeg(LineSegment line);
	bool crossWithLine1(LineSegment line);
	bool crossWithLine(LineSegment line);
	vector<LineSegment> getEdges(); //clockwise

	bool verticeInSight(Point watcher, Point ver);

	bool containVertice(Point pt);
	bool containEdge(Point pt1, Point pt2);
	bool containPointOnBoundary(Point pt);
	bool containPointInside(Point pt);
	bool containPointIncludingBoundary(Point pt);

	double innerAngle(Point verPt);
	double maxInnerAngle();
	int maxInnerAngleVertex();
	int minInnterAngleVertex();

	void getIncidentTwoEdges(int verIdx, LineSegment & e1, LineSegment & e2);

	Point getNeighbor(Point pt, bool clockwise = true);
	vector<Point> getNeighborVertices(Point pt);
	vector<pair<Point, Point>*> getNeigbhorsPairs(Point pt);

	bool forwardWithRightHandRule(Point pt1, Point pt2);
	
	Point findAnInternalPoint();

	//cut the polygon using a line, return the half polygon containing the point pt
	void dividingByLine(LineSegment &ls, MyPolygon & left, MyPolygon &right);
	void getContainingPartialPolygon(Point pt, LineSegment &ls, MyPolygon &polygon);

	bool operator == (const MyPolygon & poly) const
	{
		int i, j;
		if (vertices.size() != poly.vertices.size()) 
			return false;

		int n = vertices.size();

		int startidx = -1;
		for(i = 0; i < n; i++) {
			j = (i + 1) % n;

			if (vertices[0] == poly.vertices[i] &&
				vertices[1] == poly.vertices[j]) {
				startidx = i;
				break;
			}
		}

		if (startidx < 0) {
			return false;
		}
		else {
			for(i = 0; i < n; i++) {
				j = (i + startidx) % n;
				if (vertices[i] != poly.vertices[i]) {
					return false;
				}
			}		
		}

		/*
		for(int i = 0; i < vertices.size(); i++) {
			if (vertices[i] != poly.vertices[i]) {
				return false;
			}
		}
		*/
		return true;
	}
	
	bool operator != (const MyPolygon & poly) 
	{
		return !(*this == poly);
	}

	vector<Point> vertices;
};

vector<Point> movePoints(vector<Point> pts, double dx, double dy);
vector<Point> makeSquareHole(Point topleft, double width);
vector<Point> makeRectangleHole(Point topleft, double width, double height);
vector<Point> makeTriangle(Point p1, Point p2, Point p3);
	
Point PointTransform(Point pt, double degree);
Point PointTransformBack(Point pt, double degree);




#endif

