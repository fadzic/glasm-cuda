/*
 * glasm_cuda.h
 *
 *  Created on: Jun 27, 2015
 *      Author: filip
 */

#ifndef GLASM_CUDA_H_
#define GLASM_CUDA_H_


#include <math.h>
#include "utils.h"
#include "cuda.h"
#include "cuda_runtime.h"
#include	 <device_launch_parameters.h>
#include "stdio.h"

#define CUDA_SAFE_CALL(call) do {cudaError_t err = call; if (cudaSuccess != err) { fprintf (stderr, "Cuda error in file '%s' in line %i : %s.", __FILE__, __LINE__, cudaGetErrorString(err) );exit(EXIT_FAILURE);}} while (0)


namespace GLASM_CUDA
{

void setSearchArea(double,double,double,double,double,double);
void setCenterOfSearchArea(position);
void setRefScan(const scan&);
void setNewScan(const scan&);
void setLookupTableParameters(int,int,double,double,double,double,double);
void setGeneticParameters(int,int,int,int,int,int,double,double);
void setDebugParameters(bool udraw_lookup,
		bool udraw_individual,
		bool udraw_generations,
		unsigned uverbose_level=1,
		double umap_size_x=15.0,
		double umap_size_y=15.0,
		std::string umap_filename="./cfg/bitmaps/empty_1000x1000.png",
		std::string ulookup_bitmap="map.png",
		std::string uvalid_bitmap="map.png");

void initialize_binary_lookup_table();
void initialize_gradient_lookup_table();
void initialize_lookup_table_from_bitmap(std::string,double,double,double);
void initialize_valid_table_from_bitmap(std::string);
void initialize_radial_gradient(double,double,double,double,double,double);
void initialize_invGray_table();
void delete_lookup_table();
void clean_up_cuda();

position glasm(double* outfitness);

}

#endif /* GLASM_CUDA_H_ */
