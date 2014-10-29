#ifndef __CUBES_H__
#define __CUBES_H__

#include <vector>
#include <set>
#include <vector>
#include <algorithm>
#include <iostream>
#include "scalar_vector_types.h"
#include "timer.h"
#include "global_memory_pool.h"
#include "omp.h"

namespace G2G {
  struct Timers {
    Timer total, ciclos, rmm, density, forces, resto, pot, functions, density_derivs;
  };

  std::ostream& operator<<(std::ostream& io, const Timers& t);

/********************
 * Point information
 ********************/
struct Point {
	Point(uint _atom, uint _shell, uint _point, double3 _position, double _weight) :
		atom(_atom), shell(_shell), point(_point), position(_position), weight(_weight) {}

	uint atom, shell, point;
	double3 position;
	double weight;
};

enum FunctionType { FUNCTION_S = 1, FUNCTION_P = 3, FUNCTION_D = 6 };

template<class scalar_type>
class PointGroup {
  public:
    PointGroup(void) : number_of_points(0), s_functions(0), p_functions(0), d_functions(0), inGlobal(false) {  }
    virtual ~PointGroup(void);
    std::vector<Point> points;
    uint number_of_points;
    uint s_functions, p_functions, d_functions;

    G2G::HostMatrixUInt func2global_nuc; // size == total_functions_simple()
    G2G::HostMatrixUInt func2local_nuc; // size == total_functions()

    std::vector<uint> local2global_func; // size == total_functions_simple()
    std::vector<uint> local2global_nuc;  // size == total_nucleii()

    typedef vec_type<scalar_type,2> vec_type2;
    typedef vec_type<scalar_type,3> vec_type3;
    typedef vec_type<scalar_type,4> vec_type4;

    #if CPU_KERNELS
    G2G::HostMatrix<scalar_type> function_values;
    G2G::HostMatrix<vec_type3> gradient_values;
    G2G::HostMatrix<vec_type3> hessian_values;
    #else
    G2G::CudaMatrix<scalar_type> function_values;
    G2G::CudaMatrix<vec_type4> gradient_values;
    G2G::CudaMatrix<vec_type4> hessian_values_transposed;
    #endif

    long long cost() const;
    inline FunctionType small_function_type(uint f) const {
      if (f < s_functions) return FUNCTION_S;
      else if (f < s_functions + p_functions) return FUNCTION_P;
      else return FUNCTION_D;
    }
    //Las funciones totales, son totales del grupo, no las totales para todos los grupos.
    inline uint total_functions(void) const { return s_functions + p_functions * 3 + d_functions * 6; }
    inline uint total_functions_simple(void) const { return local2global_func.size(); } // == s_functions + p_functions + d_functions
    inline uint total_nucleii(void) const { return local2global_nuc.size(); }
    inline bool has_nucleii(uint atom) const { return (std::find(local2global_nuc.begin(), local2global_nuc.end(), atom) != local2global_nuc.end()); }

    void get_rmm_input(G2G::HostMatrix<scalar_type>& rmm_input, FortranMatrix<double>& source) const;
    void get_rmm_input(G2G::HostMatrix<scalar_type>& rmm_input) const;
    void get_rmm_input(G2G::HostMatrix<scalar_type>& rmm_input_a, G2G::HostMatrix<scalar_type>& rmm_input_b) const;

    void add_rmm_output(const G2G::HostMatrix<scalar_type>& rmm_output, FortranMatrix<double>& target) const;
    void add_rmm_output(const G2G::HostMatrix<scalar_type>& rmm_output) const;
    void add_rmm_output(const G2G::HostMatrix<scalar_type>& rmm_output, G2G::HostMatrix<double>& rmm_destination) const;
    void add_rmm_output_a(const G2G::HostMatrix<scalar_type>& rmm_output) const;
    void add_rmm_output_b(const G2G::HostMatrix<scalar_type>& rmm_output) const;
    void add_rmm_open_output(const G2G::HostMatrix<scalar_type>& rmm_output_a, const G2G::HostMatrix<scalar_type>& rmm_output_b) const;

    void compute_nucleii_maps(void);

    void add_point(const Point& p);
    void compute_weights(void);

    void compute_functions(bool forces, bool gga);
    void solve(Timers& timers, bool compute_rmm, bool lda, bool compute_forces, bool compute_energy,double&,double&,double&,double&,double&,HostMatrix<double>& fort_forces_ptr, HostMatrix<double>&, bool open);
    void solve_closed(Timers& timers, bool compute_rmm, bool lda, bool compute_forces, bool compute_energy,double&, HostMatrix<double>&, HostMatrix<double>&);
    void solve_opened(Timers& timers, bool compute_rmm, bool lda, bool compute_forces, bool compute_energy,double&,double&,double&,double&,double&,HostMatrix<double>&);

    bool is_significative(FunctionType, double exponent, double coeff, double d2);
    bool operator<(const PointGroup<scalar_type>& T) const;
    size_t size_in_gpu() const;
    void deallocate();

    virtual bool is_sphere(void) = 0;
    virtual bool is_cube(void) = 0;

    bool inGlobal;

};

// ===== Sphere Class =======//
#if FULL_DOUBLE
  typedef double base_scalar_type;
#else
  typedef float base_scalar_type;
#endif
class Sphere : public PointGroup<base_scalar_type> {
  public:
    Sphere(void);
    Sphere(uint _atom, double _radius);

    void assign_significative_functions(const std::vector<double>& min_exps, const std::vector<double>& min_coeff);
    bool is_sphere(void) { return true; }
    bool is_cube(void) { return false; }

    uint atom;
    double radius;
};

// ====== Cube Class ===========//
class Cube : public PointGroup<base_scalar_type> {
  public:
    void assign_significative_functions(const double3& cube_coord, const std::vector<double>& min_exps, const std::vector<double>& min_coeff);
    bool is_sphere(void) { return false; }
    bool is_cube(void) { return true; }

};

// =======Partition Class ========//

class Partition {
  public:
    void clear(void) {
      cubes.clear(); spheres.clear();
    }

    void generate_gpu_profile();
    void compute_work_partition();

    void regenerate(void);

    void solve(Timers& timers, bool compute_rmm,bool lda,bool compute_forces,
               bool compute_energy, double* fort_energy_ptr, double* fort_forces_ptr, bool OPEN);

    void balance_load(std::vector<long long>& thread_duration,
        std::vector<std::vector<long long> >& work_duration);

    void compute_functions(bool forces, bool gga)
    {
      Timer t1;
      t1.start_and_sync();
      for (std::vector<Cube>::iterator it = cubes.begin(); it != cubes.end(); ++it)
        it->compute_functions(forces, gga);
      for (std::vector<Sphere>::iterator it = spheres.begin(); it != spheres.end(); ++it)
        it->compute_functions(forces, gga);
      t1.stop_and_sync();
//      std::cout << "TIMER: funcs: " << t1 << std::endl;
    }

    std::vector<Cube> cubes;
    std::vector<Sphere> spheres;
    std::vector<std::vector<int> > work;
    std::vector<std::vector<double> > correction;
};

extern Partition partition;
}

#endif
