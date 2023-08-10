#ifndef types_h
#define types_h
#include <assert.h>                                             // Some compilers don't have cassert
#include <complex>
#include "kahan.h"
#include "macros.h"
#include "namespace.h"
#include <stdint.h>
#include <vector>
#include "vec.h"

#include "ityr/ityr.hpp"

namespace EXAFMM_NAMESPACE {
  template <typename T>
  using global_span = ityr::global_span<T>;
  template <typename T>
  using raw_span = ityr::common::span<T>;
  template <typename T>
  using global_vec = ityr::global_vector<T>;

  inline constexpr std::size_t cutoff_body = 1024;
  inline constexpr std::size_t cutoff_cell = 128;

  // global vector options for collective allocation
  inline constexpr ityr::global_vector_options global_vec_coll_opts {
    .collective         = true,
    .parallel_construct = true,
    .parallel_destruct  = true,
    .cutoff_count       = cutoff_body,
  };

  // serial/parallel execution policies (similar to C++17 parallel STL)
  inline constexpr ityr::execution::sequenced_policy body_seq_policy {.checkout_count = cutoff_body};
  inline constexpr ityr::execution::sequenced_policy cell_seq_policy {.checkout_count = cutoff_cell};
  inline constexpr ityr::execution::parallel_policy body_par_policy {.cutoff_count = cutoff_body, .checkout_count = cutoff_body};
  inline constexpr ityr::execution::parallel_policy cell_par_policy {.cutoff_count = cutoff_cell, .checkout_count = cutoff_cell};

  // Basic type definitions
#if EXAFMM_SINGLE
  typedef float real_t;                                         //!< Floating point type is single precision
  const real_t EPS = 1e-8f;                                     //!< Single precision epsilon
#else
  typedef double real_t;                                        //!< Floating point type is double precision
  const real_t EPS = 1e-16;                                     //!< Double precision epsilon
#endif
  typedef std::complex<real_t> complex_t;                       //!< Complex type
  const complex_t I(0.,1.);                                     //!< Imaginary unit

  typedef vec<3,int> ivec3;                                     //!< Vector of 3 int types
  typedef vec<3,real_t> vec3;                                   //!< Vector of 3 real_t types
  typedef vec<4,real_t> vec4;                                   //!< Vector of 4 real_t types
  typedef vec<3,float> fvec3;                                   //!< Vector of 3 float types
  typedef vec<3,complex_t> cvec3;                               //!< Vector of 3 complex_t types

#if EXAFMM_USE_SIMD
  // SIMD vector types for AVX512, AVX, and SSE
  const int NSIMD = SIMD_BYTES / int(sizeof(real_t));           //!< SIMD vector length (SIMD_BYTES defined in macros.h)
  typedef vec<NSIMD,real_t> simdvec;                            //!< SIMD vector type
#endif

  // Kahan summation types (Achieves quasi-double precision using single precision types)
#if EXAFMM_USE_KAHAN
  typedef kahan<real_t> kreal_t;                                //!< Real type with Kahan summation
  typedef kahan<complex_t> kcomplex_t;                          //!< Complex type with Kahan summation
#if EXAFMM_USE_SIMD
  typedef kahan<simdvec> ksimdvec;                              //!< SIMD vector type with Kahan summation
#endif
#else
  typedef real_t kreal_t;                                       //!< Real type (dummy Kahan)
  typedef complex_t kcomplex_t;                                 //!< Complex type (dummy Kahan)
#if EXAFMM_USE_SIMD
  typedef simdvec ksimdvec;                                     //!< SIMD vector type (dummy Kahan)
#endif
#endif
  typedef vec<4,kreal_t> kvec4;                                 //!< Vector of 4 real types with Kahan summaiton
  typedef vec<4,kcomplex_t> kcvec4;                             //!< Vector of 4 complex types with Kahan summaiton

  //! Center and radius of bounding box
  struct Box {
    vec3   X;                                                   //!< Box center
    real_t R;                                                   //!< Box radius
  };

  //! Min & max bounds of bounding box
  struct Bounds {
    vec3 Xmin;                                                  //!< Minimum value of coordinates
    vec3 Xmax;                                                  //!< Maximum value of coordinates
  };

  //! Structure of aligned source for SIMD
  struct Source {                                               //!< Base components of source structure
    vec3      X;                                                //!< Position
#if EXAFMM_LAPLACE
    real_t    SRC;                                              //!< Scalar real values
#elif EXAFMM_HELMHOLTZ
    complex_t SRC;                                              //!< Scalar complex values
#elif EXAFMM_BIOTSAVART
    vec4      SRC;                                              //!< Vector real values
#endif
  };

  //! Structure of bodies
  struct Body : public Source {                                 //!< Base components of body structure
    int     IBODY;                                              //!< Initial body numbering for sorting back
    int     IRANK;                                              //!< Initial rank numbering for partitioning back
    int64_t ICELL;                                              //!< Cell index
    real_t  WEIGHT;                                             //!< Weight for partitioning
#if EXAFMM_LAPLACE
    kvec4   TRG;                                                //!< Scalar+vector3 real values
#elif EXAFMM_HELMHOLTZ
    kcvec4  TRG;                                                //!< Scalar+vector3 complex values
#elif EXAFMM_BIOTSAVART
    kvec4   TRG;                                                //!< Scalar+vector3 real values
#endif
  };

  typedef raw_span<Body> Bodies;                             //!< Vector of bodies
  typedef typename Bodies::iterator B_iter;                     //!< Iterator of body vector

  typedef global_span<Body> GBodies;
  typedef typename GBodies::iterator GB_iter;

  /*
#ifdef EXAFMM_PMAX
  const int Pmax = EXAFMM_PMAX;                                 //!< Max order of expansions
#else
  const int Pmax = 10;                                          //!< Max order of expansions
#endif
  const int Pmin = 4;                                           //!< Min order of expansions
  */

  //! Base components of cells
  struct CellBase {
    int IPARENT;                                                //!< Index of parent cell
    int ICHILD;                                                 //!< Index of first child cell
    int NCHILD;                                                 //!< Number of child cells
    int IBODY;                                                  //!< Index of first body
    int NBODY;                                                  //!< Number of descendant bodies
#if EXAFMM_COUNT_LIST
    int numP2P;                                                 //!< Size of P2P interaction list per cell
    int numM2L;                                                 //!< Size of M2L interaction list per cell
#endif
    uint64_t ICELL;                                             //!< Cell index
    real_t   WEIGHT;                                            //!< Weight for partitioning
    vec3     X;                                                 //!< Cell center
    real_t   R;                                                 //!< Cell radius
    GB_iter   BODY;                                              //!< Iterator of first body
  };
  //! Structure of cells
  struct Cell : public CellBase {
    global_vec<complex_t> M;                                   //!< Multipole expansion coefs
    global_vec<complex_t> L;                                   //!< Local expansion coefs
    using CellBase::operator=;
  };

  typedef raw_span<Cell> Cells;                              //!< Vector of cells
  typedef raw_span<CellBase> CellBases;                      //!< Vector of cell bases
  typedef typename Cells::iterator C_iter;                      //!< Iterator of cell vector
  typedef typename CellBases::iterator CB_iter;                 //!< Iterator of cell vector
                                                                //
  typedef global_span<Cell> GCells;
  typedef global_span<CellBase> GCellBases;
  typedef typename GCells::iterator GC_iter;
  typedef typename GCellBases::iterator GCB_iter;
}
#endif
