/*
 * Util.h
 *
 *  Created on: Jul 8, 2016
 *      Author: dbpc
 */

#ifndef UTIL_H_
#define UTIL_H_

class Util
{
  public:
	Util () {};
    int interleave(int x, int y, int z, int level);
    std::vector<double> uninterleave(int n, int L);
    int setbit(int n, int pos, int setto);
    int getbit(int n, int pos);

};




#endif /* UTIL_H_ */
