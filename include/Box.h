/*
 * Box.h
 *
 *  Created on: Jul 13, 2016
 *      Author: dbpc
 */

#ifndef BOX_H_
#define BOX_H_

#include <complex>
#include <vector>
#include <iostream>

#include "Point.h"
//#include "FmmTree.h"

class Box
{
  public:
    //int DEFAULT_P = 16;
    //int DEFAULT_ORDER_OF_APPROXIMATION = 3;
    int DEFAULT_ORDER_OF_APPROXIMATION = 4;
    int DEFAULT_P = (DEFAULT_ORDER_OF_APPROXIMATION+1)*(DEFAULT_ORDER_OF_APPROXIMATION+1);

    // The default level of the box
    int DEFAULT_LEVEL=3;
    // The default index of the box for the boxes given level
    int DEFAULT_INDEX=0;

    int level;                             // refinement level of box
    int index;                             // cell index of box (or cell)
    int p;                                 // p is the index at which the series are truncated
    int order_of_approximation;            // abs_alpha - highest order used in series approximation

//    unsigned int nParticlesPerCell = 4;    // points in each cell at lowest refinement level (l = L)
//    int totalParticles;                     // total particles in domain (unit square)


    bool empty;                            // is box empty


    std::vector<std::vector<std::complex<double> > > c;
    std::vector<std::vector<std::complex<double> > > dtilde;
    std::vector<std::vector<std::complex<double> > > d;
    std::vector<std::vector<std::complex<double> > > dtilde_grad;
//    std::vector<std::vector<std::complex<double> > > d_grad;


    // The number of source points (and target points) in a box will depend
    // on the refinement level.  The lowest level l = L will have 4 particles
    // per box.  As we work up through the levels to level l = 0, the number
    // of particles will increase by a multiple of 4 when going from one level
    // up to the next (four children for each parent)

    std::vector<Point> x;  // source points
    std::vector<Point> y;  // target points

    // Creates a new instance of Node
    Box();
    Box(int level, int index, int p, int order_of_approximation);

    int       getLevel() { return level; };
    void      printLevel() { std::cout << "Box level is " << level << '\n'; };
    void      setLevel(int i) { this->level = i; };
    int       getIndex() { return index; };
    void      printIndex() { std::cout << "Box index is " << index << '\n'; };
    void      setIndex(int i) { this->index = i; };
    Point     getCenter();
    double    getSize() { return std::pow(2.0, -level); };

    void      setP(int p) { this->p = p;}; //this->c.resize(p); this->d.resize(p); this->dtilde.resize(p);};
	void      setOrderApprox(int order_of_approximation) { this->order_of_approximation = order_of_approximation; };
    bool      isEmpty() { return empty; };

    std::vector<std::vector<std::complex<double> > >                getC() {return c; };
//    void                               addToC(std::vector<double> &increment);
    void                               addToC(std::vector<std::vector<std::complex<double> > > &increment);
    std::string                        cToString();
    void                               printC();

    std::vector<std::vector<std::complex<double> > >                getD() {return d; };
    void                               addToD(std::vector<std::vector<std::complex<double> > > &increment);
    void                               addToDGrad(std::vector<std::vector<std::vector<std::complex<double> > > > &increment);
    std::string                        dToString();
    void                               printD();

    std::vector<std::vector<std::complex<double> > >                getDtilde() {return dtilde; };
    void                               addToDtilde(const std::vector<std::vector<std::complex<double> > > &increment);
    void                               addToDtildeGrad(std::vector<std::vector<std::vector<std::complex<double> > > > &increment);
    std::string                        dtildeToString();
    void                               printDtilde();


    void               addX(Point &p);
    int                getSizeX() { return this->x.size(); };
    std::vector<Point> getX() { return this->x; };
    void               printSizeX() { std::cout << "Box sizeX is " << x.size() << "\n"; };

    void               addY(Point &p);
    int                getSizeY() { return this->y.size(); };
    std::vector<Point> getY() { return this->y; };
    void               printSizeY() { std::cout << "Box sizeY is " << y.size() << "\n"; };

    std::string        toString();
    int                getParentIndex();
    void               getNeighborsIndex(std::vector<int> &neighbor_indexes);
    void               getParentsNeighborsIndex (std::vector<int> &parents_neighbor_indexes);

    void               getNeighborsE4Index (std::vector<int> &neighborE4_indexes);

    void               getChildrenIndex(std::vector<int> &children_indexes);


  private:

    void               getChildrenIndexOfBox(int levelOfBox, int indexOfBox, std::vector<int> &children_indexes_of_box);

};




#endif /* BOX_H_ */
