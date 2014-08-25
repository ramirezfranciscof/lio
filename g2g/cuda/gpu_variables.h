#ifndef _GPU_VARIABLES_H
#define _GPU_VARIABLES_H

#include "../common.h"

/* grid */
__device__ uint gpu_atoms;
__device__ uint gpu_Iexch;
 #if FULL_DOUBLE
__device__ double3 gpu_atom_positions[MAX_ATOMS];
__device__ double gpu_normalization_factor;
 #else
__device__ float3 gpu_atom_positions[MAX_ATOMS];
__device__ float gpu_normalization_factor;
 #endif
 
 #endif
