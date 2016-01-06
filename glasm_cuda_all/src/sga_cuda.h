/*
 * sga_cuda.h
 *
 *  Created on: Jun 27, 2015
 *      Author: filip
 */

#ifndef SGA_CUDA_H_
#define SGA_CUDA_H_


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
#ifdef DRAW_PNG
	void set_sga_drawfun(void (*drawfun)(int run, int gen, int popsize, struct individual *critter, struct bestever *bestfit));
#endif
bestever sga(void);

#endif /* SGA_CUDA_H_ */
