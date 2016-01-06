/*
 * testisi.h
 *
 *  Created on: Jul 19, 2015
 *      Author: filip
 */

#ifndef TESTISI_H_
#define TESTISI_H_

struct individual
{   unsigned chrom;		// 6/3/2011 carefull! max 32 bit per chromosome in this version
	double   fitness;
	int      xsite;
	int      parent[2];
};

struct bestever
{	unsigned chrom;
	double   fitness;
	int      generation;
};

void sga_parameters(int,int,int,int,int,int,double,double);
void set_sga_objfun(void (*objfun)(struct individual *critter));


#endif /* TESTISI_H_ */
