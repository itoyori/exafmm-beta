#ifndef up_down_pass_h
#define up_down_pass_h
#include "logger.h"
#include "namespace.h"
#include "types.h"

namespace EXAFMM_NAMESPACE {
  class UpDownPass {
  private:
    Kernel* kernel;                                            //!< Kernel class

  private:
    //! Post-order traversal for upward pass
    void postOrderTraversal(GC_iter C, GC_iter C0) const {
      auto cs = ityr::make_checkout(C, 1, ityr::checkout_mode::read);
      int ichild = cs[0].ICHILD;
      int nchild = cs[0].NCHILD;
      cs.checkin();

      ityr::for_each(
          ityr::execution::par,
          ityr::count_iterator<int>(0),
          ityr::count_iterator<int>(nchild),
          [=, *this](int i) {
            postOrderTraversal(C0 + ichild + i, C0);
          });

      auto cs2 = ityr::make_checkout(C, 1, ityr::checkout_mode::read);

      if(nchild==0) {
        kernel->P2M(&cs2[0]);                           // P2M kernel
      } else {                                                    // If not leaf cell
        auto csj = ityr::make_checkout(C0 + ichild, nchild, ityr::checkout_mode::read);
        kernel->M2M(&cs2[0], csj.data());                                      //  M2M kernel
      }                                                         // End if for non leaf cell
    };

    //! Pre-order traversal for downward pass
    void preOrderTraversal(GC_iter C, GC_iter C0) const {
      auto cs = ityr::make_checkout(C, 1, ityr::checkout_mode::read);
      auto& C_ = cs[0];

      auto csj = ityr::make_checkout(C0 + C_.IPARENT, 1, ityr::checkout_mode::read);
      auto& Cj0_ = csj[0];

      kernel->L2L(&C_, &Cj0_);                                        //  L2L kernel

      if (C_.NCHILD==0) {                                       //  If leaf cell
        kernel->L2P(&C_);                                          //  L2P kernel
      }                                                         // End if for leaf cell
#if 0
#if EXAFMM_USE_WEIGHT
      C_iter CP = C0 + C->IPARENT;                              // Parent cell
      C->WEIGHT += CP->WEIGHT;                                  // Add parent's weight
      if (C->NCHILD==0) {                                       // If leaf cell
        for (B_iter B=C->BODY; B!=C->BODY+C->NBODY; B++) {      //  Loop over bodies in cell
          B->WEIGHT += C->WEIGHT;                               //   Add cell weights to bodies
        }                                                       //  End loop over bodies in cell
      }                                                         // End if for leaf cell
#endif
#endif

      int ichild = C_.ICHILD;
      int nchild = C_.NCHILD;
      cs.checkin();
      csj.checkin();

      ityr::for_each(
          ityr::execution::par,
          ityr::count_iterator<int>(0),
          ityr::count_iterator<int>(nchild),
          [=, *this](int i) {
            preOrderTraversal(C0 + ichild + i, C0);
          });
    };

  public:
    //! Constructor
    UpDownPass(Kernel* _kernel) : kernel(_kernel) {}           // Initialize variables

    //! Upward pass (P2M, M2M)
    void upwardPass(GCells cells) {
      if (ityr::is_master()) {
        logger::startTimer("Upward pass");                        // Start timer
      }
      ityr::root_exec([=, *this] {
        if (!cells.empty()) {                                     // If cell vector is not empty
          GC_iter C0 = cells.begin();                              //  Set iterator of target root cell
          ityr::for_each(
              cell_par_policy,
              ityr::make_global_iterator(cells.begin(), ityr::checkout_mode::read_write),
              ityr::make_global_iterator(cells.end()  , ityr::checkout_mode::read_write),
              [=, *this](Cell& c) {
                c.M.resize(kernel->NTERM, 0.0);                       //   Allocate & initialize M coefs
                c.L.resize(kernel->NTERM, 0.0);                       //   Allocate & initialize L coefs
              });
          postOrderTraversal(C0, C0);                             //  Start post-order traversal from root
        }                                                         // End if for empty cell vector
      });
      if (ityr::is_master()) {
        logger::stopTimer("Upward pass");                         // Stop timer
      }
    }

    //! Downward pass (L2L, L2P)
    void downwardPass(GCells cells) {
      if (ityr::is_master()) {
        logger::startTimer("Downward pass");                      // Start timer
      }
      ityr::root_exec([=, *this] {
        if (!cells.empty()) {                                     // If cell vector is not empty
          GC_iter C0 = cells.begin();                              //  Root cell
          auto cs = ityr::make_checkout(C0, 1, ityr::checkout_mode::read);
          int ichild = cs[0].ICHILD;
          int nchild = cs[0].NCHILD;
          if (nchild == 0) {                                 //  If root is the only cell
            kernel->L2P(&cs[0]);                                       //   L2P kernel
          }                                                       //  End if root is the only cell
          cs.checkin();
          ityr::for_each(
              ityr::execution::par,
              ityr::count_iterator<int>(0),
              ityr::count_iterator<int>(nchild),
              [=, *this](int i) {
                preOrderTraversal(C0 + ichild + i, C0);                            //   Start pre-order traversal from root
              });
        }                                                         // End if for empty cell vector
      });
      if (ityr::is_master()) {
        logger::stopTimer("Downward pass");                       // Stop timer
      }
    }

    //! Get dipole of entire system
    vec3 getDipole(Bodies & bodies, vec3 X0) {
      vec3 dipole = 0;                                          // Initialize dipole correction
#if EXAFMM_LAPLACE
      for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {     // Loop over bodies
	dipole += (B->X - X0) * std::real(complex_t(B->SRC));   //  Calcuate dipole of the whole system
      }                                                         // End loop over bodies
#endif
      return dipole;                                            // Return dipole
    }

    //! Dipole correction
    void dipoleCorrection(Bodies & bodies, vec3 dipole, int numBodies, vec3 cycle) {
#if EXAFMM_LAPLACE
      real_t coef = 4 * M_PI / (3 * cycle[0] * cycle[1] * cycle[2]);// Precalcualte constant
      for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {     // Loop over bodies
	B->TRG[0] -= coef * norm(dipole) / numBodies / B->SRC;  //  Dipole correction for potential
	for (int d=0; d!=3; d++) {                              //  Loop over dimensions
	  B->TRG[d+1] -= coef * dipole[d];                      //   Dipole correction for forces
	}                                                       //  End loop over dimensions
      }                                                         // End loop over bodies
#endif
    }
  };
}
#endif
