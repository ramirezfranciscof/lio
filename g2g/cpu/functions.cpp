#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include "../common.h"
//#include "../cuda_includes.h"
#include "../init.h"
#include "../matrix.h"
#include "../partition.h"
using namespace std;

namespace G2G {

inline int rang(int l, int i, int r)
{
    if(i < l) return 0;
    if(i >= r) return r - l;
    return i - l;
}

template<class scalar_type>
void PointGroup<scalar_type>::compute_functions(Timers & t, bool forces, bool gga)
{
  t.load_functions.start();
  if (this->inGlobal) return;
  if (!globalMemoryPool::tryAlloc(size_in_cpu())) {
    this->inGlobal = true;
  } else {
    this->inGlobal = false;
  }
  if (this->inGlobal) {
      function_values.set_permanent();
      function_values_transposed.set_permanent();

      gX.set_permanent();
      gY.set_permanent();
      gZ.set_permanent();

      hPX.set_permanent();
      hPY.set_permanent();
      hPZ.set_permanent();
      hIX.set_permanent();
      hIY.set_permanent();
      hIZ.set_permanent();
  }
  t.load_functions.pause();
  /* Load group functions */
  uint group_m = total_functions();

  t.create_matrices.start();
  function_values.resize(group_m, number_of_points);
  function_values_transposed.resize(number_of_points, group_m);

  gX.resize(group_m, number_of_points);
  gY.resize(group_m, number_of_points);
  gZ.resize(group_m, number_of_points);

  hPX.resize(group_m, number_of_points);
  hPY.resize(group_m, number_of_points);
  hPZ.resize(group_m, number_of_points);
  hIX.resize(group_m, number_of_points);
  hIY.resize(group_m, number_of_points);
  hIZ.resize(group_m, number_of_points);
  t.create_matrices.pause();

  const uint funcs = total_functions_simple();
  const int npoints = points.size();
  int point; uint i;
  
  t.calculate_matrices.start();
  #pragma omp parallel 
  {

  #pragma omp for collapse(2) schedule(static) private(point,i) nowait
  for(point = 0; point < npoints; point++) {
    for (i = 0; i < s_functions; i++) {
      const uint ii = 1 * rang(0, i, s_functions) +
                      3 * rang(s_functions, i, s_functions + p_functions) +
                      6 * rang(s_functions + p_functions, i, funcs);

      // compute exponential
      uint nuc = func2global_nuc(i);
      const double3 atompos = fortran_vars.atom_positions(nuc);
      const scalar_type vx = points[point].position.x - atompos.x;
      const scalar_type vy = points[point].position.y - atompos.y;
      const scalar_type vz = points[point].position.z - atompos.z;
      const scalar_type dist = vx * vx + vy * vy + vz * vz;

      scalar_type t = 0, tg = 0, th = 0;
      uint global_func = local2global_func[i];
      uint contractions = fortran_vars.contractions(global_func);

      #pragma vector always aligned
      for (uint contraction = 0; contraction < contractions; contraction++) {
        const scalar_type a = fortran_vars.a_values(global_func, contraction);
        const scalar_type c = fortran_vars.c_values(global_func, contraction);
        const scalar_type t0 = exp(-(a * dist)) * c;
        t += t0;
        tg += t0 * a;
        th += t0 * (a * a);
      }

      // compute s, p, d
      function_values(ii, point) = function_values_transposed(point, ii) = t;
      // gradient_values(ii, point) = vec_type3(v * (-2 * tg));
      gX(ii,point) = vx * -2 * tg;
      gY(ii,point) = vy * -2 * tg;
      gZ(ii,point) = vz * -2 * tg;
      //hessian_values(2 * ii + 0, point) = vec_type3((v * v) * 4 * th - 2 * tg); // Fxx, Fxy, Fxz
      hPX(ii, point) = vx * vx * 4 * th - 2 * tg;
      hPY(ii, point) = vy * vy * 4 * th - 2 * tg;
      hPZ(ii, point) = vz * vz * 4 * th - 2 * tg;
      //hessian_values(2 * ii + 1, point) = vec_type3(vxxy * vyzz * 4 * th); // Fxy, Fxz, Fyz
      hIX(ii, point) = vx * vy * 4 * th;
      hIY(ii, point) = vx * vz * 4 * th;
      hIZ(ii, point) = vy * vz * 4 * th;
    }
  }

  #pragma omp for collapse(2) schedule(static) private(point,i) nowait
  for(point = 0; point < npoints; point++) {
    for (i = s_functions; i < s_functions + p_functions; i++) {
      const uint ii = 1 * rang(0, i, s_functions) +
                      3 * rang(s_functions, i, s_functions + p_functions) +
                      6 * rang(s_functions + p_functions, i, funcs);

      // compute exponential
      uint nuc = func2global_nuc(i);
      const double3 atompos = fortran_vars.atom_positions(nuc);
      const scalar_type vx = points[point].position.x - atompos.x;
      const scalar_type vy = points[point].position.y - atompos.y;
      const scalar_type vz = points[point].position.z - atompos.z;
      const scalar_type dist = vx * vx + vy * vy + vz * vz;

      scalar_type t = 0, tg = 0, th = 0;
      uint global_func = local2global_func[i];
      uint contractions = fortran_vars.contractions(global_func);
      
      #pragma vector aligned always
      for (uint contraction = 0; contraction < contractions; contraction++) {
        const scalar_type a = fortran_vars.a_values(global_func, contraction);
        const scalar_type c = fortran_vars.c_values(global_func, contraction);
        const scalar_type t0 = exp(-(a * dist)) * c;
        t += t0;
        tg += t0 * a;
        th += t0 * (a * a);
      }

      function_values(ii + 0, point) = function_values_transposed(point, ii + 0) = vx * t;
      function_values(ii + 1, point) = function_values_transposed(point, ii + 1) = vy * t;
      function_values(ii + 2, point) = function_values_transposed(point, ii + 2) = vz * t;

      //gradient_values(ii + 0, point) = vec_type3(vec_type3(t, 0, 0) - v * 2 * tg * v.x());
      gX(ii, point) = t - vx * 2 * tg * vx;
      gY(ii, point) = 0 - vy * 2 * tg * vx;
      gZ(ii, point) = 0 - vz * 2 * tg * vx;
      //gradient_values(ii + 1, point) = vec_type3(vec_type3(0, t, 0) - v * 2 * tg * vy);
      gX(ii + 1, point) = 0 - vx * 2 * tg * vy;
      gY(ii + 1, point) = t - vy * 2 * tg * vy;
      gZ(ii + 1, point) = 0 - vz * 2 * tg * vy;
      //gradient_values(ii + 2, point) = vec_type3(vec_type3(0, 0, t) - v * 2 * tg * vz);
      gX(ii + 2, point) = 0 - vx * 2 * tg * vz;
      gY(ii + 2, point) = 0 - vy * 2 * tg * vz;
      gZ(ii + 2, point) = t - vz * 2 * tg * vz;
      //hessian_values(2 * (ii + 0) + 0, point) = vec_type3((v * v) *         4 * th * v.x() - vec_type3(6, 2, 2)     * tg * v.x());
      hPX(ii, point) = vx * vx * 4 * th * vx - 6 * tg * vx;
      hPY(ii, point) = vy * vy * 4 * th * vx - 2 * tg * vx;
      hPZ(ii, point) = vz * vz * 4 * th * vx - 2 * tg * vx;
      //hessian_values(2 * (ii + 0) + 1, point) = vec_type3((vxxy * vyzz) * 4 * th * v.x() - vec_type3(vy, vz, 0) * 2 * tg);
      hIX(ii, point) = vx * vy * 4 * th * vx - vy * 2 * tg;
      hIY(ii, point) = vx * vz * 4 * th * vx - vz * 2 * tg;
      hIZ(ii, point) = vy * vz * 4 * th * vx - 0 * 2 * tg;
      //hessian_values(2 * (ii + 1) + 0, point) = vec_type3((v * v) *        4 * th * vy - vec_type3(2, 6, 2)     * tg * vy);
      hPX(ii+1, point) = vx * vx * 4 * th * vy - 2 * tg * vy;
      hPY(ii+1, point) = vy * vy * 4 * th * vy - 6 * tg * vy;
      hPZ(ii+1, point) = vz * vz * 4 * th * vy - 2 * tg * vy;
      //hessian_values(2 * (ii + 1) + 1, point) = vec_type3((vxxy * vyzz) * 4 * th * vy - vec_type3(v.x(), 0, vz) * 2 * tg);
      hIX(ii+1, point) = vx * vy * 4 * th * vy - 2 * tg * vx;
      hIY(ii+1, point) = vx * vz * 4 * th * vy - 2 * tg * 0;
      hIZ(ii+1, point) = vy * vz * 4 * th * vy - 2 * tg * vz;
      //hessian_values(2 * (ii + 2) + 0, point) = vec_type3((v * v)       *   4 * th * vz - vec_type3(2, 2, 6)     * tg * vz);
      hPX(ii+2, point) = vx * vx * 4 * th * vz - 2 * tg * vz;
      hPY(ii+2, point) = vy * vy * 4 * th * vz - 2 * tg * vz;
      hPZ(ii+2, point) = vz * vz * 4 * th * vz - 6 * tg * vz;
//    hessian_values(2 * (ii + 2) + 1, point) = vec_type3((vxxy * vyzz) * 4 * th * vz - vec_type3(0, v.x(), vy) * 2 * tg);
      hIX(ii+2, point) = vx * vy * 4 * th * vz - 2 * tg * 0;
      hIY(ii+2, point) = vx * vz * 4 * th * vz - 2 * tg * vx;
      hIZ(ii+2, point) = vy * vz * 4 * th * vz - 2 * tg * vy;
    }
  }

  #pragma omp for collapse(2) schedule(static) private(point,i) nowait
  for(point = 0; point < npoints; point++) {
    for (i = s_functions + p_functions; i < funcs; i++) {
      const uint ii = 1 * rang(0, i, s_functions) +
                      3 * rang(s_functions, i, s_functions + p_functions) +
                      6 * rang(s_functions + p_functions, i, funcs);

      // compute exponential
      uint nuc = func2global_nuc(i);
      const double3 atompos = fortran_vars.atom_positions(nuc);
      const scalar_type vx = points[point].position.x - atompos.x;
      const scalar_type vy = points[point].position.y - atompos.y;
      const scalar_type vz = points[point].position.z - atompos.z;
      const scalar_type dist = vx * vx + vy * vy + vz * vz;

      scalar_type t = 0, tg = 0, th = 0;
      uint global_func = local2global_func[i];
      uint contractions = fortran_vars.contractions(global_func);

      #pragma vector aligned always
      for (uint contraction = 0; contraction < contractions; contraction++) {
        const scalar_type a = fortran_vars.a_values(global_func, contraction);
        const scalar_type c = fortran_vars.c_values(global_func, contraction);
        const scalar_type t0 = exp(-(a * dist)) * c;
        t += t0;
        tg += t0 * a;
        th += t0 * (a * a);
      }

      function_values(ii + 0, point) = function_values_transposed(point, ii + 0) = t * vx * vx * fortran_vars.normalization_factor;
      function_values(ii + 1, point) = function_values_transposed(point, ii + 1) = t * vy * vx;
      function_values(ii + 2, point) = function_values_transposed(point, ii + 2) = t * vy * vy * fortran_vars.normalization_factor;
      function_values(ii + 3, point) = function_values_transposed(point, ii + 3) = t * vz * vx;
      function_values(ii + 4, point) = function_values_transposed(point, ii + 4) = t * vz * vy;
      function_values(ii + 5, point) = function_values_transposed(point, ii + 5) = t * vz * vz * fortran_vars.normalization_factor;

      //gradient_values(ii + 0, point) = vec_type3((vec_type3(2 * v.x(), 0      , 0      ) * t - v * 2 * tg * v.x() * v.x()) * fortran_vars.normalization_factor);
      gX(ii, point) = (2 * vx * t - vx * 2 * tg * vx * vx) * fortran_vars.normalization_factor;
      gY(ii, point) = (0 * t - vy * 2 * tg * vx * vx) * fortran_vars.normalization_factor;
      gZ(ii, point) = (0 * t - vz * 2 * tg * vx * vx) * fortran_vars.normalization_factor;
      //gradient_values(ii + 1, point) = vec_type3(vec_type3(vy     , v.x()    , 0      ) * t - v * 2 * tg * vy * v.x());
      gX(ii+1, point) = vy * t - vx * 2 * tg * vy * vx;
      gY(ii+1, point) = vx * t - vy * 2 * tg * vy * vx;
      gZ(ii+1, point) = 0 * t - vz * 2 * tg * vy * vx;
      //gradient_values(ii + 2, point) = vec_type3((vec_type3(0      , 2 * vy, 0      ) * t - v * 2 * tg * vy * vy) * fortran_vars.normalization_factor);
      gX(ii+2, point) = (0 * t - vx * 2 * tg * vy * vy) * fortran_vars.normalization_factor;
      gY(ii+2, point) = (2 * vy * t - vy * 2 * tg * vy * vy) * fortran_vars.normalization_factor;
      gZ(ii+2, point) = (0 * t - vz * 2 * tg * vy * vy) * fortran_vars.normalization_factor;
      //gradient_values(ii + 3, point) = vec_type3( vec_type3(vz    , 0      , v.x()    ) * t - v * 2 * tg * vz * v.x());
      gX(ii+3, point) = vz * t - vx * 2 * tg * vz * vx;
      gY(ii+3, point) = 0 * t - vy * 2 * tg * vz * vx;
      gZ(ii+3, point) = vx * t - vz * 2 * tg * vz * vx;
      //gradient_values(ii + 4, point) = vec_type3(vec_type3(0       , vz    , vy    ) * t - v * 2 * tg * vz * vy);
      gX(ii+4, point) = 0 * t - vx * 2 * tg * vz * vy;
      gY(ii+4, point) = vz * t - vy * 2 * tg * vz * vy;
      gZ(ii+4, point) = vy * t - vz * 2 * tg * vz * vy;
      //gradient_values(ii + 5, point) = vec_type3((vec_type3(0      , 0      , 2 * vz) * t - v * 2 * tg * vz * vz) * fortran_vars.normalization_factor);
      gX(ii+5, point) = (0 * t - vx * 2 * tg * vz * vz) * fortran_vars.normalization_factor;
      gY(ii+5, point) = (0 * t - vy * 2 * tg * vz * vz) * fortran_vars.normalization_factor;
      gZ(ii+5, point) = (2 * vz * t - vz * 2 * tg * vz * vz) * fortran_vars.normalization_factor;

//    hessian_values(2 * (ii + 0) + 0, point) = vec_type3(((v * v)       * 4 * th * (v.x() * v.x()) - vec_type3(10, 2, 2) * tg * (v.x() * v.x())    + vec_type3(2 * t, 0    , 0)) * fortran_vars.normalization_factor);
      hPX(ii, point) = (vx * vx * 4 * th * vx * vx - 10 * tg * vx * vx + 2* t) * fortran_vars.normalization_factor;
      hPY(ii, point) = (vy * vy * 4 * th * vx * vx - 2 * tg * vx * vx + 0) * fortran_vars.normalization_factor;
      hPZ(ii, point) = (vz * vz * 4 * th * vx * vx - 2 * tg * vx * vx + 0) * fortran_vars.normalization_factor;
//    hessian_values(2 * (ii + 0) + 1, point) = vec_type3(((vxxy * vyzz) * 4 * th * (v.x() * v.x()) - vec_type3(4,  4, 0) * tg * (vxxy * vyzz)                                 ) * fortran_vars.normalization_factor);
      hIX(ii, point) = (vx * vy * 4 * th * vx * vx - 4 * tg * vx * vy) * fortran_vars.normalization_factor;
      hIY(ii, point) = (vx * vz * 4 * th * vx * vx - 4 * tg * vx * vz) * fortran_vars.normalization_factor;
      hIZ(ii, point) = (vy * vz * 4 * th * vx * vx - 0 * tg * vy * vz) * fortran_vars.normalization_factor;
//    hessian_values(2 * (ii + 1) + 0, point) = vec_type3(((v * v)       * 4 * th * (v.x() * vy) - vec_type3(6,  6, 2) * tg * (v.x() * vy)                                   ));
      hPX(ii+1, point) = vx * vx * 4 * th * vx * vy - 6 * tg * vx * vy;
      hPY(ii+1, point) = vy * vy * 4 * th * vx * vy - 6 * tg * vx * vy;
      hPZ(ii+1, point) = vz * vz * 4 * th * vx * vy - 2 * tg * vx * vy;
//    hessian_values(2 * (ii + 1) + 1, point) = vec_type3(((vxxy * vyzz) * 4 * th * (v.x() * vy) - vec_type3(2 * (v.x() * v.x() + vy * vy), 2 * vy * vz, 2 * v.x() * vz) * tg + vec_type3(t     , 0    , 0)));
      hIX(ii+1, point) = vx * vy * 4 * th * vx * vy - 2 * (vx * vx + vy * vy) * tg + t;
      hIY(ii+1, point) = vx * vz * 4 * th * vx * vy - 2 * vy * vz * tg + 0;
      hIZ(ii+1, point) = vy * vz * 4 * th * vx * vy - 2 * vx * vz * tg + 0;
//    hessian_values(2 * (ii + 2) + 0, point) = vec_type3(((v * v)       * 4 * th * (vy * vy) - vec_type3(2, 10, 2) * tg * (vy * vy) + vec_type3(0    , 2 * t, 0)) * fortran_vars.normalization_factor);
      hPX(ii+2, point) = (vx * vx * 4 * th * vy * vy - 2 * tg * vy * vy + 0) * fortran_vars.normalization_factor;
      hPY(ii+2, point) = (vy * vy * 4 * th * vy * vy - 10 * tg * vy * vy + 2 * t) * fortran_vars.normalization_factor;
      hPZ(ii+2, point) = (vz * vz * 4 * th * vy * vy - 2 * tg * vy * vy + 0) * fortran_vars.normalization_factor;
//    hessian_values(2 * (ii + 2) + 1, point) = vec_type3(((vxxy * vyzz) * 4 * th * (vy * vy) - vec_type3(4,  0, 4) * tg * (vxxy * vyzz)                                ) * fortran_vars.normalization_factor);
      hIX(ii+2, point) = (vx * vy * 4 * th * (vy * vy) - 4 * tg * vx * vy) * fortran_vars.normalization_factor;
      hIY(ii+2, point) = (vx * vz * 4 * th * (vy * vy) - 0 * tg * vx * vz) * fortran_vars.normalization_factor;
      hIZ(ii+2, point) = (vy * vz * 4 * th * (vy * vy) - 4 * tg * vy * vz) * fortran_vars.normalization_factor;
//    hessian_values(2 * (ii + 3) + 0, point) = vec_type3(((v * v)       * 4 * th * (v.x() * vz) - vec_type3(6,  2, 6) * tg * (v.x() * vz)                                  ));
      hPX(ii+3, point) = vx * vx * 4 * th * vx * vz - 6 * tg  * vx * vz;
      hPY(ii+3, point) = vy * vy * 4 * th * vx * vz - 2 * tg  * vx * vz;
      hPZ(ii+3, point) = vz * vz * 4 * th * vx * vz - 6 * tg  * vx * vz;
//    hessian_values(2 * (ii + 3) + 1, point) = vec_type3(((vxxy * vyzz) * 4 * th * (v.x() * vz) - vec_type3(2 * vy * vz, 2 * (v.x() * v.x() + vz * vz), 2 * v.x() * vy) * tg + vec_type3(0,      t,     0)));
      hIX(ii+3, point) = vx * vy * 4 * th * vx * vz - 2 * vy * vz * tg + 0;
      hIY(ii+3, point) = vx * vz * 4 * th * vx * vz - 2 * (vx * vx + vz * vz) * tg + t;
      hIZ(ii+3, point) = vy * vz * 4 * th * vx * vz - 2 * vx * vy * tg + 0;
//    hessian_values(2 * (ii + 4) + 0, point) = vec_type3(((v * v)       * 4 * th * (vy * vz) - vec_type3(2,  6, 6) * tg * (vy * vz)                                ));
      hPX(ii+4, point) = vx * vx * 4 * th * vy * vz - 2 * tg * vy * vz;
      hPY(ii+4, point) = vy * vy * 4 * th * vy * vz - 6 * tg * vy * vz;
      hPZ(ii+4, point) = vz * vz * 4 * th * vy * vz - 6 * tg * vy * vz;
//    hessian_values(2 * (ii + 4) + 1, point) = vec_type3(((vxxy * vyzz) * 4 * th * (vy * vz) - vec_type3(2 * v.x() * vz, 2 * v.x() * vy, 2 * (vy * vy + vz * vz)) * tg + vec_type3(0,      0,     t)));
      hIX(ii+4, point) = vx * vy * 4 * th * vy * vz - 2 * vx * vz * tg + 0;
      hIY(ii+4, point) = vx * vz * 4 * th * vy * vz - 2 * vx * vy * tg + 0;
      hIZ(ii+4, point) = vy * vz * 4 * th * vy * vz - 2 * (vy * vy + vz * vz) * tg + t;
//    hessian_values(2 * (ii + 5) + 0, point) = vec_type3(((v * v)       * 4 * th * (vz * vz) - vec_type3(2,  2, 10) * tg * (vz * vz) + vec_type3(0,      0, 2 * t)) * fortran_vars.normalization_factor);
      hPX(ii+5, point) =  (vx * vx * 4 * th * vz * vz - 2 * tg * vz * vz + 0) * fortran_vars.normalization_factor;
      hPY(ii+5, point) =  (vy * vy * 4 * th * vz * vz - 2 * tg * vz * vz + 0) * fortran_vars.normalization_factor;
      hPZ(ii+5, point) =  (vz * vz * 4 * th * vz * vz - 10 * tg * vz * vz + 2 * t) * fortran_vars.normalization_factor;
//    hessian_values(2 * (ii + 5) + 1, point) = vec_type3(((vxxy * vyzz) * 4 * th * (vz * vz) - vec_type3(0,  4, 4) * tg * (vxxy * vyzz)                                 ) * fortran_vars.normalization_factor);
      hIX(ii+5, point) =  (vx * vy * 4 * th * vz * vz - 0 * tg * vx * vy) * fortran_vars.normalization_factor;
      hIY(ii+5, point) =  (vx * vz * 4 * th * vz * vz - 4 * tg * vx * vz) * fortran_vars.normalization_factor;
      hIZ(ii+5, point) =  (vy * vz * 4 * th * vz * vz - 4 * tg * vy * vz) * fortran_vars.normalization_factor;
    }
  }
  
  }
  t.calculate_matrices.pause();
}

template class PointGroup<double>;
template class PointGroup<float>;
}
