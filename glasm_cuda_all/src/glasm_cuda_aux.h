/*
 * glasm_cuda_aux.h
 *
 *  Created on: Jun 27, 2015
 *      Author: filip
 */

#ifndef GLASM_CUDA_AUX_H_
#define GLASM_CUDA_AUX_H_


#include "utils.h"
#include "glasm_cuda.h"


namespace GLASM_CUDA
{

	void init(const char* filename);                    // one time initialization. Reads and sets parameters.
	void pre_match(const scan &sref, const scan &snew); // before matching
	matching_result match();                            // matching
	void post_match();                                  // after matching
	void deinit();                                      // one time deinitialization

}

#endif /* GLASM_CUDA_AUX_H_ */
