
/*
 * Funcion llamada para cada (i,j) en RMM, para calcular RMM(i,j) -> un thread por cada punto
 */

// TODO: esto desperdicia la mitad de los threads -> quizas se puede armar una grilla sin los bloques que no hagan nada

template<class scalar_type,bool check_pos>
__global__ void gpu_update_rmm(const scalar_type* __restrict__ factors, uint points, scalar_type* rmm, const scalar_type* __restrict__ function_values, uint m)
{
//    if (blockIdx.x * blockDim.x > blockIdx.y * blockDim.y) return;
    uint i, j, first_fi, first_fj;
    // There's more than one block; do the math to get position
    if (check_pos) {
        // Figure out where our block is in the lower triangle
        // We get 1D index k = blockIdx.x; column 1 value of k in row j is always k_1 = j*(j+1);
        // solving for j gives determinant 1+8*k_1
        // Thus, 1+8*k_1 must be square of (odd) integer; take sqrt(1+8*k) and get first odd integer below it -
        // get row and column from there
        uint n = sqrtf(1.0f+8.0f*blockIdx.x);
        n -= (1 - n % 2);
        uint block_j = (n - 1) / 2;
        uint block_i = blockIdx.x - (block_j + 1) * block_j / 2;

        first_fi = block_i*blockDim.x;
        first_fj = block_j*blockDim.y;
        i = first_fi + threadIdx.x;
        j = first_fj + threadIdx.y;
    // There's only one block; we don't need to take the square root
    } else {
        uint3 pos = index(blockDim, blockIdx, threadIdx);

        i = pos.x; // columna
        j = pos.y; // fila
        first_fi = blockIdx.x * blockDim.x;
        first_fj = blockIdx.y * blockDim.y;
    }

    bool valid_thread = (i < m && j < m && i <= j); // quiero triangulo inferior solamente TODO: sacar esto

    // calculate this rmm
    scalar_type rmm_local = 0.0f;

    __shared__ scalar_type functions_i_local[RMM_BLOCK_SIZE_XY][RMM_BLOCK_SIZE_XY+1];	// Fi[point][i]
    __shared__ scalar_type functions_j_local[RMM_BLOCK_SIZE_XY][RMM_BLOCK_SIZE_XY+1]; // Fj[point][j]
    __shared__ scalar_type factor_local[RMM_BLOCK_SIZE_XY * RMM_BLOCK_SIZE_XY];			// factor[point]    // TODO: esto seguramente tiene bank conflicts

    uint inc = RMM_BLOCK_SIZE_XY*RMM_BLOCK_SIZE_XY;
    for (uint point_base = 0; point_base < points; point_base += inc) {
        uint abs_threadIdx = threadIdx.y * blockDim.x + threadIdx.x;  // absolute threadId inside block

        __syncthreads();

        /* all threads load a point */
        if (point_base + abs_threadIdx < points)
            factor_local[abs_threadIdx] = factors[point_base + abs_threadIdx];

        uint last_point = point_base + inc;
#pragma unroll 16
        for (uint point = point_base; point < last_point; point += RMM_BLOCK_SIZE_XY) {
            if (point < points) {
                /* every RMM_BLOCK_SIZE_X iterations, Fi and Fj get filled with RMM_BLOCK_SIZE_Y functions,
                 * for RMM_BLOCK_SIZE_X different points */
                __syncthreads();
                bool ifexterno = point + threadIdx.x < points;
                bool ifinterno1 = (first_fi + threadIdx.y) < m;
                bool ifinterno2 = (first_fj + threadIdx.y) < m;
                bool valid1 = ifexterno * ifinterno1;
                bool valid2 = ifexterno * ifinterno2;
                scalar_type fi_times_factor =
                  function_values[valid1 * (COALESCED_DIMENSION(points) * (first_fi + threadIdx.y) + (point + threadIdx.x))] *
                  factor_local[valid1 * ((point-point_base) + threadIdx.x)];
                functions_i_local[threadIdx.x][threadIdx.y] =  valid1 * fi_times_factor;
              /* on blocks with some invalid threads, the computation contributes 0 to rmm_local (
               * this avoids instruction serialization) */
                functions_j_local[threadIdx.x][threadIdx.y] =
                  function_values[valid2 * (COALESCED_DIMENSION(points) * (first_fj + threadIdx.y) + (point + threadIdx.x))] *
                  valid2;

                __syncthreads();
                for (uint point_sub = 0; point_sub < RMM_BLOCK_SIZE_XY; point_sub++) {
                    rmm_local += functions_i_local[point_sub][threadIdx.x] * functions_j_local[point_sub][threadIdx.y];
                }
            }
        }
    }

    if (valid_thread) rmm[COALESCED_DIMENSION(m) * j + i] = rmm_local;
}
