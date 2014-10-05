#include <iostream>
#include <limits>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include "common.h"
#include "init.h"
#include "matrix.h"
#include "partition.h"
#include "timer.h"
using namespace std;

namespace G2G {
Partition partition;

ostream& operator<<(ostream& io, const Timers& t) {
#ifdef TIMINGS
  cout << "iteration: " << t.total << endl;
  cout << "rmm: " << t.rmm << " density: " << t.density << " pot: " << t.pot << " forces: " << t.forces << " resto: " << t.resto << " functions: " << t.functions << endl;
#endif
  return io;
}

/********************
 * PointGroup
 ********************/

template<class scalar_type>
void PointGroup<scalar_type>::get_rmm_input(HostMatrix<scalar_type>& rmm_input,
    FortranMatrix<double>& source) const {
  rmm_input.zero();
  for (uint i = 0, ii = 0; i < total_functions_simple(); i++) {
    uint inc_i = small_function_type(i);

    for (uint k = 0; k < inc_i; k++, ii++) {
      uint big_i = local2global_func[i] + k;
      for (uint j = 0, jj = 0; j < total_functions_simple(); j++) {
        uint inc_j = small_function_type(j);

        for (uint l = 0; l < inc_j; l++, jj++) {
          uint big_j = local2global_func[j] + l;
          if (big_i > big_j) continue;
          uint big_index = (big_i * fortran_vars.m - (big_i * (big_i - 1)) / 2) + (big_j - big_i);

          rmm_input(ii, jj) = (scalar_type)source.data[big_index];

          rmm_input(jj, ii) = rmm_input(ii, jj);
        }
      }
    }
  }
}

template<class scalar_type>
void PointGroup<scalar_type>::get_rmm_input(HostMatrix<scalar_type>& rmm_input) const {
  get_rmm_input(rmm_input, fortran_vars.rmm_input_ndens1);
}

template<class scalar_type>
void PointGroup<scalar_type>::get_rmm_input(HostMatrix<scalar_type>& rmm_input_a, HostMatrix<scalar_type>& rmm_input_b) const {
  get_rmm_input(rmm_input_a, fortran_vars.rmm_dens_a);
  get_rmm_input(rmm_input_b, fortran_vars.rmm_dens_b);
}

template<class scalar_type>
void PointGroup<scalar_type>::add_rmm_output(const HostMatrix<scalar_type>& rmm_output,
    FortranMatrix<double>& target ) const {
  for (uint i = 0, ii = 0; i < total_functions_simple(); i++) {
    uint inc_i = small_function_type(i);

    for (uint k = 0; k < inc_i; k++, ii++) {
      uint big_i = local2global_func[i] + k;
      for (uint j = 0, jj = 0; j < total_functions_simple(); j++) {
        uint inc_j = small_function_type(j);

        for (uint l = 0; l < inc_j; l++, jj++) {
          uint big_j = local2global_func[j] + l;
          if (big_i > big_j) continue;
          uint big_index = (big_i * fortran_vars.m - (big_i * (big_i - 1)) / 2) + (big_j - big_i);
          target(big_index) += (double)rmm_output(ii, jj);
        }
      }
    }
  }
}

template<class scalar_type>
void PointGroup<scalar_type>::add_rmm_output(const HostMatrix<scalar_type>& rmm_output) const {
  add_rmm_output(rmm_output, fortran_vars.rmm_output);
}

template<class scalar_type>
void PointGroup<scalar_type>::add_rmm_output(
    const HostMatrix<scalar_type>& rmm_output, HostMatrix<double>& target) const {
  for (uint i = 0, ii = 0; i < total_functions_simple(); i++) {
    uint inc_i = small_function_type(i);

    for (uint k = 0; k < inc_i; k++, ii++) {
      uint big_i = local2global_func[i] + k;
      for (uint j = 0, jj = 0; j < total_functions_simple(); j++) {
        uint inc_j = small_function_type(j);

        for (uint l = 0; l < inc_j; l++, jj++) {
          uint big_j = local2global_func[j] + l;
          if (big_i > big_j) continue;
          uint big_index = (big_i * fortran_vars.m - (big_i * (big_i - 1)) / 2) + (big_j - big_i);
          target(big_index) += (double)rmm_output(ii, jj);
        }
      }
    }
  }

}

template<class scalar_type>
void PointGroup<scalar_type>::add_rmm_output_a(const HostMatrix<scalar_type>& rmm_output) const {
  add_rmm_output(rmm_output, fortran_vars.rmm_output_a);
}

template<class scalar_type>
void PointGroup<scalar_type>::add_rmm_output_b(const HostMatrix<scalar_type>& rmm_output) const {
  add_rmm_output(rmm_output, fortran_vars.rmm_output_b);
}

template<class scalar_type>
void PointGroup<scalar_type>::add_rmm_open_output(const HostMatrix<scalar_type>& rmm_output_a,
    const HostMatrix<scalar_type>& rmm_output_b) const {
  add_rmm_output(rmm_output_a, fortran_vars.rmm_output_a);
  add_rmm_output(rmm_output_b, fortran_vars.rmm_output_b);
}

template<class scalar_type>
void PointGroup<scalar_type>::compute_nucleii_maps(void)
{
  if (total_functions_simple() != 0) {
    func2global_nuc.resize(total_functions_simple());
    for (uint i = 0; i < total_functions_simple(); i++) {
      func2global_nuc(i) = fortran_vars.nucleii(local2global_func[i]) - 1;
    }

    func2local_nuc.resize(total_functions());
    uint ii = 0;
    for (uint i = 0; i < total_functions_simple(); i++) {
      uint global_atom = func2global_nuc(i);
      uint local_atom = std::distance(local2global_nuc.begin(),
          std::find(local2global_nuc.begin(), local2global_nuc.end(), global_atom));
      uint inc = small_function_type(i);
      for (uint k = 0; k < inc; k++, ii++) func2local_nuc(ii) = local_atom;
    }
  }
}

template<class scalar_type>
void PointGroup<scalar_type>::add_point(const Point& p) {
  points.push_back(p);
  number_of_points++;
}

#define EXP_PREFACTOR 1.01057089636005 // (2 * pow(4, 1/3.0)) / M_PI

template<class scalar_type>
bool PointGroup<scalar_type>::is_significative(FunctionType type, double exponent, double coeff, double d2) {
  switch(type) {
    case FUNCTION_S:
      return (exponent * d2 < max_function_exponent-log(pow((2.*exponent/M_PI),3))/4);
    break;
    default:
    {
      double x = 1;
      double delta;
      double e = 0.1;
      double factor = pow((2.0*exponent/M_PI),3);
      factor = sqrt(factor*4.0*exponent) ;
      double norm = (type == FUNCTION_P ? sqrt(factor) : abs(factor)) ;
      do {
        double div = (type == FUNCTION_P ? log(x) : 2 * log(x));
        double x1 = sqrt((max_function_exponent - log(norm) + div) / exponent);
        delta = abs(x-x1);
        x = x1;
      } while (delta > e);
      return (sqrt(d2) < x);
    }
    break;
  }
}
template<class scalar_type>
bool PointGroup<scalar_type>::operator<(const PointGroup<scalar_type>& T) const{
    int my_cost = number_of_points * total_functions();
    int T_cost = T.number_of_points * T.total_functions();
    return my_cost < T_cost;
}
template<class scalar_type>
size_t PointGroup<scalar_type>::size_in_gpu() const
{
    uint total_cost=0;
    uint single_matrix_cost = COALESCED_DIMENSION(number_of_points) * total_functions();

    total_cost += single_matrix_cost;       //1 scalar_type functions
    if (fortran_vars.do_forces || fortran_vars.gga)
      total_cost += (single_matrix_cost*4); //4 vec_type gradient
    if (fortran_vars.gga)
      total_cost+= (single_matrix_cost*8);  //2*4 vec_type hessian
    return total_cost*sizeof(scalar_type);  // size in bytes according to precision
}

template<class scalar_type>
void PointGroup<scalar_type>::deallocate() {
#if !CPU_KERNELS
  if(inGlobal) {
    GlobalMemoryPool::dealloc(size_in_gpu());
    function_values.deallocate();
    gradient_values.deallocate();
    hessian_values.deallocate();
    hessian_values_transposed.deallocate();
    inGlobal = false;
  }
#endif
}
template<class scalar_type>
PointGroup<scalar_type>::~PointGroup<scalar_type>() {
  deallocate();
}

// Work stealing
void Partition::balance_load(vector<long long>& thread_duration, vector<vector<long long> >& work_duration) {
  int thread_round_counts = 0;
  const int max_thread_rounds = 5;
  // Hacemos multiples rondas de work stealing, por si tenemos varios threads que balancear
  while(thread_round_counts < max_thread_rounds) {
    thread_round_counts++;
    // Primero consigo el thread mas lento y el mas rapido
    long long min_time = numeric_limits<long long>::max();
    long long max_time = 0;
    int min_time_index = 0; int max_time_index = 0;
    for(int i = 0; i < thread_duration.size(); i++){
      if (thread_duration[i] > max_time) {
        max_time_index = i; max_time = thread_duration[i];
      }
      if (thread_duration[i] < min_time) {
        min_time_index = i; min_time = thread_duration[i];
      }
    }
    // Conseguimos el factor de correccion que vamos a usar para mover de un thread de
    // mas duracion a uno de menos duracion.
    const double time_correction = correction[max_time_index][min_time_index];
    // Agregamos un factor de correcion adicional producto de migrar el thread y recalcular las funciones.
    const double migration_factor = 7.0/6.0;
    int round_counts = 0;
    const int max_rounds = 15;
    // Si hay mas de un 2% de diferencia, migramos tareas
    while (max_time > min_time+(min_time/50) && round_counts < max_rounds) {
      round_counts++;
      long long delta = (max_time - min_time)/2;
      // Busco todos los trabajos del thread con max_time, el que mas cerca tenga
      // duracion como delta tiempo de los trabajos (/2 porque es un promedio) y se lo mando al
      // thread mas corto.
      // Multiplicamos por 7/6 para agregarle 1/6 de factor de penalidad por mover de placa.
      int best_index = 0;
      double best_index_delta =
        abs(work_duration[max_time_index][best_index] * time_correction * migration_factor - delta);
      for (int i = 1; i < work_duration[max_time_index].size(); i++) {
        double current_index_delta =
          abs((work_duration[max_time_index][i] * time_correction * migration_factor) - delta);
        if(current_index_delta < best_index_delta) {
          best_index = i; best_index_delta = current_index_delta;
        }
      }
      // best_index tiene el indice de trabajo del thread que mas tardo, que se podria
      // procesar en el otro thread
      int element_index = work[max_time_index][best_index];
      std::cout << "Voy a mover de " << max_time_index << " a " << min_time_index <<
        " un trabajo con duracion = " << work_duration[max_time_index][best_index] << std::endl;
      // Liberamos la memoria global del proceso que vayamos a migrar, por si se cachearon
      // esas funciones.
      if(element_index < cubes.size())
        cubes[element_index].deallocate();
      else
        spheres[element_index-cubes.size()].deallocate();
      // Actualizamos la migracion de cargas y las estimaciones, considerando las diferencias
      // de velocidad entre ambos.
      thread_duration[min_time_index] += work_duration[max_time_index][best_index] * time_correction * migration_factor;
      thread_duration[max_time_index] -= work_duration[max_time_index][best_index];

      min_time = thread_duration[min_time_index];
      max_time = thread_duration[max_time_index];

      work[min_time_index].push_back(element_index);
      work[max_time_index].erase(work[max_time_index].begin()+best_index);

      work_duration[min_time_index].push_back(work_duration[max_time_index][best_index] * time_correction * migration_factor);
      work_duration[max_time_index].erase(work_duration[max_time_index].begin()+best_index);
    }
  }
}

void Partition::solve(Timers& timers, bool compute_rmm,bool lda,bool compute_forces,
                      bool compute_energy, double* fort_energy_ptr, double* fort_forces_ptr, bool OPEN) {
  double cubes_energy = 0, spheres_energy = 0;
  double cubes_energy_i = 0, spheres_energy_i = 0;
  double cubes_energy_c = 0, spheres_energy_c = 0;
  double cubes_energy_c1 = 0, spheres_energy_c1 = 0;
  double cubes_energy_c2 = 0, spheres_energy_c2 = 0;

  int gpu_count;
  cudaGetDeviceCount(&gpu_count);
  if (gpu_count == 0) {
    std::cout << "Error: No GPU found" << std::endl;
    exit(1);
  }
  int total_threads = gpu_count;
  #ifndef _OPENMP
  total_threads = 1;
  gpu_count = 1;
  #else
  omp_set_num_threads(total_threads);
  #endif
  double energy_cubes[total_threads];
  double energy_spheres[total_threads];

  HostMatrix<double> fort_forces_ms[total_threads];

  if (compute_forces) {
      for(int i = 0; i < total_threads; i++) {
          fort_forces_ms[i].resize(fortran_vars.max_atoms, 3);
          fort_forces_ms[i].zero();
      }
  }

  HostMatrix<double> rmm_outputs[total_threads];
  if (compute_rmm) {
      for(int i = 0; i < total_threads; i++) {
          rmm_outputs[i].resize(fortran_vars.rmm_output.width, fortran_vars.rmm_output.height);
          rmm_outputs[i].zero();
      }
  }
  vector<long long> thread_duration(total_threads);
  vector<vector<long long> > work_duration(total_threads);
#pragma omp parallel shared(energy_spheres, energy_cubes)
  {
    int my_thread = 0;
    #ifdef _OPENMP
    my_thread = omp_get_thread_num();
    #endif
    energy_spheres[my_thread] = 0.0f;
    energy_cubes[my_thread] = 0.0f;
    if(cudaSetDevice(my_thread % gpu_count) != cudaSuccess)
      std::cout << "Error: can't set the device for thread " << my_thread << " using # " << gpu_count
        << " assigning device " << my_thread % gpu_count << std::endl;

    // Util este coeficiente para hacer pruebas de escalabilidad.
    double slow_coef = (my_thread==1)? (1) : 1;
    Timer t0;
    t0.start_and_sync();
    for(int i = 0; i < work[my_thread].size(); i++){
      int k = work[my_thread][i];
      Timer t1;
      t1.start_and_sync();
      if(k < cubes.size()) {
        cubes[k].solve(
            timers, compute_rmm,lda,compute_forces, compute_energy, energy_cubes[my_thread], cubes_energy_i,
            cubes_energy_c, cubes_energy_c1, cubes_energy_c2, fort_forces_ms[my_thread], rmm_outputs[my_thread], OPEN);
      }
      else {
        spheres[k - cubes.size()].solve(
            timers, compute_rmm,lda,compute_forces, compute_energy, energy_spheres[my_thread],
            spheres_energy_i, spheres_energy_c, spheres_energy_c1, spheres_energy_c2, fort_forces_ms[my_thread],
            rmm_outputs[my_thread], OPEN);
      }
      t1.stop_and_sync();
      work_duration[my_thread].push_back(slow_coef * (t1.getMicrosec() + 1000*1000*t1.getSec()));
    }
    t0.stop_and_sync();
    thread_duration[my_thread] = slow_coef*(t0.getMicrosec() + 1000*1000*t0.getSec());
    std::cout << "Workload " << my_thread << " " << thread_duration[my_thread] << std::endl;
  }

  if(total_threads > 1) {
    // Work stealing
    balance_load(thread_duration, work_duration);
  }

  if (compute_forces) {
      FortranMatrix<double> fort_forces_out(fort_forces_ptr, fortran_vars.atoms, 3, fortran_vars.max_atoms);
      for(int k = 0; k < total_threads; k++) {
          for(int i = 0; i < fortran_vars.atoms; i++) {
              for(int j = 0; j < 3; j++) {
                  fort_forces_out(i,j) += fort_forces_ms[k](i,j);
              }
          }
      }
  }

  if(compute_rmm) {
    for(int k = 0; k < total_threads; k++) {
      for(int i = 0; i < rmm_outputs[k].width; i++) {
        for(int j = 0; j < rmm_outputs[k].height; j++) {
          fortran_vars.rmm_output(i,j) += rmm_outputs[k](i,j);
        }
      }
    }
  }

  for(int i = 0; i< total_threads; i++) {
    cubes_energy += energy_cubes[i];
    spheres_energy += energy_spheres[i];
  }

  *fort_energy_ptr = cubes_energy + spheres_energy;
  if(*fort_energy_ptr != *fort_energy_ptr) {
      std::cout << "I see dead peaple " << std::endl;
#ifndef CPU_KERNELS
      cudaDeviceReset();
#endif
      exit(1);
   }

}



/**********************
 * Sphere
 **********************/
Sphere::Sphere(void) : atom(0), radius(0) { }
Sphere::Sphere(uint _atom, double _radius) : atom(_atom), radius(_radius) { }

/**********************
 * Cube
 **********************/

template class PointGroup<double>;
template class PointGroup<float>;
}
