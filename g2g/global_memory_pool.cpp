#include "global_memory_pool.h"
#include <cassert>

//TryAlloc return true if error, false if success
int GlobalMemoryPool::tryAlloc(size_t size) {
  if (!_init) init();
  int current_device; cudaGetDevice(&current_device);
  if (_freeGlobalMemory[current_device] < size)
      return 1;
  _freeGlobalMemory[current_device] -= size;
  return 0;
}

void GlobalMemoryPool::dealloc(size_t size) {
  if(!_init) init();
  int current_device; cudaGetDevice(&current_device);
  assert (_freeGlobalMemory[current_device] > size);
  _freeGlobalMemory[current_device] += size;
}

void GlobalMemoryPool::init(double free_global_memory) {
  //Aca tenemos que leer el GPU_global y restar un factor de tolerancia (1/4?)
#if !CPU_KERNELS
  int previous_device; cudaGetDevice(&previous_device);
  int gpu_count; cudaGetDeviceCount(&gpu_count);
  for (int i = 0; i < gpu_count; i++) {
    size_t free_memory, total_memory;
    cudaSetDevice(i);
    cudaGetMemoryInfo(free_memory, total_memory);
    double free_factor = free_global_memory;

    if (free_factor > 1.0f) free_factor = 1.0f;
    if (free_factor < 0.0f) free_factor = 0.0f;
    _freeFactor = free_factor;

    _freeGlobalMemory.push_back(static_cast<size_t>(static_cast<double>(free_memory)*_freeFactor));
    _totalGlobalMemory.push_back(total_memory);
  }
  cudaSetDevice(previous_device);
#else
  _totalGlobalMemory.push_back(0);
  _freeGlobalMemory.push_back(0);
  _freeFactor = 0.0f;
#endif
  _init = true;
}


std::vector<size_t> GlobalMemoryPool::_freeGlobalMemory;
std::vector<size_t> GlobalMemoryPool::_totalGlobalMemory;
float GlobalMemoryPool::_freeFactor = 0.8;
bool  GlobalMemoryPool::_init = false;
