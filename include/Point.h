/*
 * Point.h
 *
 *  Created on: Jul 8, 2016
 *      Author: dbpc
 */

#ifndef POINT_H_
#define POINT_H_

class Point
{
  public:
	std::vector<double> coord; // coordinates for point in 3-dimensional space

	/** constructors */
	Point();
    Point(double x, double y, double z);
    Point(std::vector<double> coord);
    //Point(const Point &p);

    std::vector<double> getCoord() { return coord; };
    void setCoord(std::vector<double> coord) { this->coord = coord; };
    void   setX(double x_coord) { this->coord[0]=x_coord; };
    double getX() { return coord[0]; };
    void   setY(double y_coord) { this->coord[1]=y_coord; };
    double getY() { return coord[1]; };
    void   setZ(double z_coord) { this->coord[2]=z_coord; };
    double getZ() { return coord[2]; };

    std::string coordToString();
    bool equals(Point &p);
    int getBoxIndex(unsigned int level);
};




#endif /* POINT_H_ */
