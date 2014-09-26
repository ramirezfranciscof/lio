#ifndef _GPU_VARIABLES_H
#define _GPU_VARIABLES_H

#include "../common.h"

/* grid */
__device__ __constant__ uint gpu_atoms_[MAX_GPUS];
__device__ __constant__ uint gpu_Iexch_[MAX_GPUS];
#if FULL_DOUBLE
__device__ __constant__ double3 gpu_atom_positions_[MAX_GPUS][MAX_ATOMS];
__device__ __constant__ double gpu_normalization_factor_[MAX_GPUS];
#else
__device__ __constant__ float3 gpu_atom_positions_[MAX_GPUS][MAX_ATOMS];
__device__ __constant__ float gpu_normalization_factor_[MAX_GPUS];
#endif

#endif
