#include <stdlib.h>
#include "utils.h"
#include "testisi.h"

namespace CUDA_TESTIS
{
	//global variables on host
	int nbitx;
	int nbity;
	int nbitrot;
	int popsize;
	int maxruns;
	int maxgen;
	double pcross;
	double pmutation;

	const int BITS_PER_BYTE = 8;// number of bits per byte on this machine

	int lchrom;

	struct individual *oldpop;/* last generation of individuals */
	struct individual *newpop;/* next generation of individuals */
	struct bestever bestfit; /* fittest individual so far */
	double sumfitness;/* summed fitness for entire population */
	double mymax;/* maximum fitness of population */
	double avg;/* average fitness of population */
	double mymin;/* minumum fitness of population */

	int gen;/* current generation number */
	int run;/* current run number */
	int nmutation;/* number of mutations */
	int ncross;/* number of crossovers */



	void (*sga_objfun)(struct individual *critter);		// pointer to application dependent objective function


	void set_sga_objfun(void (*objfun)(struct individual *critter))
	{
		sga_objfun=objfun;
	}
}
