/*
 *  Copyright (C) 2014  Kristijan Lenac
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License along
 *  with this program; if not, write to the Free Software Foundation, Inc.,
 *  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */


//#include <iostream>
//using namespace std;


#include <stdlib.h>

#include "sga_cuda.h"
#include "utils.h"


//---------------------------------------------------------------------------
//------------------------ sga ----------------------------------------------
//---------------------------------------------------------------------------

namespace CUDA{

	// definizioni globali
	int nbitx;
	int nbity;
	int nbitrot;
	int popsize;
	int *d_popsize;
	int maxruns;
	int maxgen;
	double pcross;
	double pmutation;

	const int BITS_PER_BYTE = 8;// number of bits per byte on this machine
	//const int UINTSIZE = BITS_PER_BYTE*sizeof(unsigned int);// # of bits in unsigned

	int lchrom;
	int *d_lchrom;

	__device__ struct individual *oldpop;/* last generation of individuals */
	__device__ struct individual *newpop;/* next generation of individuals */
	__device__ struct bestever bestfit; /* fittest individual so far */
	__device__ double sumfitness;/* summed fitness for entire population */
	__device__ double mymax;/* maximum fitness of population */
	__device__ double avg;/* average fitness of population */
	__device__ double mymin;/* minumum fitness of population */

	int gen;/* current generation number */
	int *d_gen;
	int run;/* current run number */
	__device__ int nmutation;/* number of mutations */
	__device__ int ncross;/* number of crossovers */

	int *rand_numbers;
	int *d_rand_numbers;

	//ovo ce trebat promjenit. Najlakse je prepisat objfun na device
	void (*sga_objfun)(struct individual *critter);

	void set_sga_objfun(void (*objfun)(struct individual *critter))
	{
		sga_objfun=objfun;
	}
	//--------------------------------------------------------------

	void sga_parameters(int unbitx, int unbity, int unbitrot, int upopsize, int umaxruns, int umaxgen, double upcross, double upmutation)
	{
		nbitx=unbitx;
		nbity=unbity;
		nbitrot=unbitrot;
		lchrom=nbitx+nbity+nbitrot;

		popsize=upopsize;
		if((popsize%2) != 0) popsize++;
		maxruns=umaxruns;
		maxgen=umaxgen;
		pcross=upcross;
		pmutation=upmutation;

		cudaMalloc((void **)&d_popsize,sizeof(int));
		cudaMemcpy(d_popsize, &popsize, sizeof(int), cudaMemcpyHostToDevice);

		cudaMalloc((void **)&d_lchrom,sizeof(int));
		cudaMemcpy(d_lchrom, &lchrom, sizeof(int), cudaMemcpyHostToDevice);
	}

	__device__ static int *choices, nremain;
	__device__ static float *fraction;

	__device__ void d_selectMemory(int *popsize)
	{
		unsigned nbytes;

		nbytes = *popsize*sizeof(int);
		if((choices = (int *) malloc(nbytes)) == NULL) {} //Errore: Non posso allocare memoria dinamica per choices
		nbytes = *popsize*sizeof(float);
		if((fraction = (float *) malloc(nbytes)) == NULL) {} //Errore: Non posso allocare memoria dinamica per fraction
	}

	__global__ void d_initmalloc(int *popsize)
	{
		unsigned nbytes;

		nbytes = *popsize*sizeof(struct individual); // memory for old and new populations of individuals
		if((oldpop = (struct individual *) malloc(nbytes)) == NULL) {} //Errore: Non posso allocare memoria dinamica per oldpop
		if((newpop = (struct individual *) malloc(nbytes)) == NULL) {} //Errore: Non posso allocare memoria dinamica per newpop

		d_selectMemory(popsize);

		nmutation = 0;// initialize global counters/values
		ncross = 0;
		bestfit.fitness = 0.0;
		bestfit.generation = 0;
	}


	__device__ void d_sga_objfun(individual *critter)
	{
		critter->fitness = 4;
	}

	__global__ void d_initialize(int *lchrom, int* rand_num)
	{
		int j, k;// initialize population
		unsigned mask = 1;

		j = blockIdx.x;

		oldpop[j].chrom = 0;
		for(k = 0; k < *lchrom; k++)
		{
			oldpop[j].chrom = (oldpop[j].chrom<<1);
			if(rand_num[j]<499) oldpop[j].chrom = oldpop[j].chrom|mask;
		}
		oldpop[j].parent[0] = 0;// Initialize parent info
		oldpop[j].parent[1] = 0;
		oldpop[j].xsite = 0;

		d_sga_objfun(&(oldpop[j]));// Evaluate initial fitness
	}

	__global__ void d_statistics_1()
	{
		sumfitness = 0.0;
		mymin = oldpop[0].fitness;
		mymax = oldpop[0].fitness;
	}

	__global__ void d_statistics_2(int *gen)
	{
		struct individual *pop = oldpop;

		int j;
		j = blockIdx.x;

		__syncthreads();

		sumfitness = sumfitness + pop[j].fitness;// Accumulate
		if(pop[j].fitness > mymax) mymax = pop[j].fitness;// New maximum
		if(pop[j].fitness < mymin) mymin = pop[j].fitness;// New minimum
		if(pop[j].fitness > bestfit.fitness)// new global best-fit individual
		{
			bestfit.chrom = pop[j].chrom;
			bestfit.fitness= pop[j].fitness;
			bestfit.generation = *gen;
		}
	}

	__global__ void d_statistics_3(int *popsize)
	{
		avg = sumfitness/(*popsize);// Calculate average
	}

	void statistics_old()
	{
		d_statistics_1<<<1,1>>>();
		d_statistics_2<<<popsize,1>>>(d_gen);
		d_statistics_3<<<1,1>>>(d_popsize);
	}

	void randomize()
	{
		for(int i=0; i<popsize; i++)
		{
			rand_numbers[i] = _random(1000);
		}

		cudaMemcpy(d_rand_numbers, rand_numbers, sizeof(int), cudaMemcpyHostToDevice);
	}

	void initialize()
	{
		randomize();
		d_initialize<<<popsize,1>>>(d_lchrom, d_rand_numbers);
		statistics_old();
	}

	void initmalloc()
	{
		rand_numbers = (int *) malloc(popsize);

		cudaMalloc((void **)&d_gen,sizeof(int));
		cudaMalloc((void **)&d_rand_numbers,sizeof(int));

		d_initmalloc<<<1,1>>>(d_popsize);
	}

	void updateGen()
	{
		cudaMemcpy(d_gen, &gen, sizeof(int), cudaMemcpyHostToDevice);
	}

	void generation()
	{

	}

	bestever sga(void) // la funzione principale da chiamare dopo aver settato tutti i parametri del
	{
		initmalloc();

		for(run=1; run<=maxruns; run++)
		{
			initialize();

			for(gen=0; gen<maxgen; gen++)
			{
				generation();
			}
		}

		return bestfit;
	}

}
