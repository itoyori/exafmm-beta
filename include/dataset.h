#ifndef dataset_h
#define dataset_h
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include "namespace.h"
#include <sstream>
#include "types.h"

namespace EXAFMM_NAMESPACE {
  class Dataset {                                               // Class for datasets
  private:
    long filePosition;                                          //!< Position of file stream

    //! Split range and return partial range
    void splitRange(int & begin, int & end, int iSplit, int numSplit) const {
      assert(end > begin);                                      // Check that size > 0
      int size = end - begin;                                   // Size of range
      int increment = size / numSplit;                          // Increment of splitting
      int remainder = size % numSplit;                          // Remainder of splitting
      begin += iSplit * increment + std::min(iSplit,remainder); // Increment the begin counter
      end = begin + increment;                                  // Increment the end counter
      if (remainder > iSplit) end++;                            // Adjust the end counter for remainder
    }

#if 0
    //! Uniform distribution on [-1,1]^3 lattice
    Bodies lattice(int numBodies, int mpirank, int mpisize) {
      int nx = int(std::pow(numBodies*mpisize, 1./3));          // Number of points in x direction
      int ny = nx;                                              // Number of points in y direction
      int nz = nx;                                              // Number of points in z direction
      int begin = 0;                                            // Begin index in z direction
      int end = nz;                                             // End index in z direction
      splitRange(begin, end, mpirank, mpisize);                 // Split range in z direction
      int numLattice = nx * ny * (end - begin);                 // Total number of lattice points
      Bodies bodies(numLattice);                                // Initialize bodies
      B_iter B = bodies.begin();                                // Initialize body iterator
      for (int ix=0; ix<nx; ix++) {                             // Loop over x direction
	for (int iy=0; iy<ny; iy++) {                           //  Loop over y direction
	  for (int iz=begin; iz<end; iz++, B++) {               //   Loop over z direction
	    B->X[0] = (ix / real_t(nx-1)) * 2 - 1;              //    x coordinate
	    B->X[1] = (iy / real_t(ny-1)) * 2 - 1;              //    y coordinate
	    B->X[2] = (iz / real_t(nz-1)) * 2 - 1;              //    z coordinate
	  }                                                     //   End loop over z direction
	}                                                       //  End loop over y direction
      }                                                         // End loop over x direction
      return bodies;                                            // Return bodies
    }
#endif

    // TODO: support other particle distribution

    //! Random distribution in [-1,1]^3 cube
    void cube(GBodies bodies, int seed, int numSplit) const {
      for (int i=0; i<numSplit; i++, seed++) {                  // Loop over partitions (if there are any)
	int begin = 0;                                          //  Begin index of bodies
	int end = bodies.size();                                //  End index of bodies
	splitRange(begin, end, i, numSplit);                    //  Split range of bodies
	srand48(seed);                                          //  Set seed for random number generator

        ityr::for_each(
            body_seq_policy,
            ityr::make_global_iterator(bodies.begin() + begin, ityr::checkout_mode::write),
            ityr::make_global_iterator(bodies.begin() + end  , ityr::checkout_mode::write),
            [&](Body& B) {
              for (int d=0; d<3; d++) {                             //   Loop over dimension
                B.X[d] = drand48() * 2 * M_PI - M_PI;              //    Initialize coordinates
              }                                                     //   End loop over dimension
            });
      }                                                         // End loop over partitions
    }

#if 0
    //! Random distribution on r = 1 sphere
    Bodies sphere(int numBodies, int seed, int numSplit) {
      Bodies bodies(numBodies);                                 // Initialize bodies
      for (int i=0; i<numSplit; i++, seed++) {                  // Loop over partitions (if there are any)
	int begin = 0;                                          //  Begin index of bodies
	int end = bodies.size();                                //  End index of bodies
	splitRange(begin, end, i, numSplit);                    //  Split range of bodies
	srand48(seed);                                          //  Set seed for random number generator
	for (B_iter B=bodies.begin()+begin; B!=bodies.begin()+end; B++) {// Loop over bodies
	  for (int d=0; d<3; d++) {                             //   Loop over dimension
	    B->X[d] = drand48() * 2 - 1;                        //    Initialize coordinates
	  }                                                     //   End loop over dimension
	  real_t r = std::sqrt(norm(B->X));                     //   Distance from center
	  for (int d=0; d<3; d++) {                             //   Loop over dimension
	    B->X[d] *= M_PI / r;                                //    Normalize coordinates
	  }                                                     //   End loop over dimension
	}                                                       //  End loop over bodies
      }                                                         // End loop over partitions
      return bodies;                                            // Return bodies
    }

    //! Random distribution on one octant of a r = 1 sphere
    Bodies octant(int numBodies, int seed, int numSplit) {
      Bodies bodies(numBodies);                                 // Initialize bodies
      for (int i=0; i<numSplit; i++, seed++) {                  // Loop over partitions (if there are any)
	int begin = 0;                                          //  Begin index of bodies
	int end = bodies.size();                                //  End index of bodies
	splitRange(begin, end, i, numSplit);                    //  Split range of bodies
	srand48(seed);                                          //  Set seed for random number generator
	for (B_iter B=bodies.begin()+begin; B!=bodies.begin()+end; B++) {// Loop over bodies
	  real_t theta = drand48() * M_PI * 0.5;                //   Polar angle [0,pi/2]
	  real_t phi = drand48() * M_PI * 0.5;                  //   Azimuthal angle [0,pi/2]
	  B->X[0] = 2 * M_PI * std::sin(theta) * std::cos(phi) - M_PI;// x coordinate
	  B->X[1] = 2 * M_PI * std::sin(theta) * std::sin(phi) - M_PI;// y coordinate
	  B->X[2] = 2 * M_PI * std::cos(theta) - M_PI;          //    z coordinate
	}                                                       //  End loop over bodies
      }                                                         // End loop over partitions
      return bodies;                                            // Return bodies
    }

    //! Plummer distribution in a r = M_PI/2 sphere
    Bodies plummer(int numBodies, int seed, int numSplit) {
      Bodies bodies(numBodies);                                 // Initialize bodies
      for (int i=0; i<numSplit; i++, seed++) {                  // Loop over partitions (if there are any)
	int begin = 0;                                          //  Begin index of bodies
	int end = bodies.size();                                //  End index of bodies
	splitRange(begin, end, i, numSplit);                    //  Split range of bodies
	srand48(seed);                                          //  Set seed for random number generator
	B_iter B=bodies.begin()+begin;                          //  Body begin iterator
	while (B != bodies.begin()+end) {                       //  While body iterator is within range
	  real_t X1 = drand48();                                //   First random number
	  real_t X2 = drand48();                                //   Second random number
	  real_t X3 = drand48();                                //   Third random number
	  real_t R = 1.0 / sqrt( (pow(X1, -2.0 / 3.0) - 1.0) ); //   Radius
	  if (R < 100.0) {                                      //   If radius is less than 100
	    real_t Z = (1.0 - 2.0 * X2) * R;                    //    z component
	    real_t X = sqrt(R * R - Z * Z) * std::cos(2.0 * M_PI * X3);// x component
	    real_t Y = sqrt(R * R - Z * Z) * std::sin(2.0 * M_PI * X3);// y component
	    real_t scale = 3.0 * M_PI / 16.0;                   //    Scaling factor
	    X *= scale; Y *= scale; Z *= scale;                 //    Scale coordinates
	    B->X[0] = X;                                        //    Assign x coordinate to body
	    B->X[1] = Y;                                        //    Assign y coordinate to body
	    B->X[2] = Z;                                        //    Assign z coordinate to body
	    B++;                                                //    Increment body iterator
	  }                                                     //   End if for bodies within range
	}                                                       //  End while loop over bodies
      }                                                         // End loop over partitions
      return bodies;                                            // Return bodies
    }
#endif

  public:
    Dataset() : filePosition(0) {}                              // Constructor

    // TODO: parallelize
    //! Initialize source values
    void initSource(GBodies bodies, int seed, int numSplit) const {
      for (int i=0; i<numSplit; i++, seed++) {                  // Loop over partitions (if there are any)
	int begin = 0;                                          //  Begin index of bodies
	int end = bodies.size();                                //  End index of bodies
	splitRange(begin, end, i, numSplit);                    //  Split range of bodies
	srand48(seed);                                          //  Set seed for random number generator

#if EXAFMM_LAPLACE

	real_t average = 0;                                     //  Initialize average charge
        ityr::for_each(
            body_seq_policy,
            ityr::make_global_iterator(bodies.begin() + begin, ityr::checkout_mode::read_write),
            ityr::make_global_iterator(bodies.begin() + end  , ityr::checkout_mode::read_write),
            [&](Body& B) {
              B.SRC = drand48() - .5;                              //   Initialize charge
              average += B.SRC;                                    //   Accumulate average
            });

        average /= (end - begin);

        ityr::for_each(
            body_par_policy,
            ityr::make_global_iterator(bodies.begin() + begin, ityr::checkout_mode::read_write),
            ityr::make_global_iterator(bodies.begin() + end  , ityr::checkout_mode::read_write),
            [=](Body& B) {
              B.SRC -= average;
            });

#elif EXAFMM_HELMHOLTZ

        ityr::for_each(
            body_par_policy,
            ityr::make_global_iterator(bodies.begin() + begin, ityr::checkout_mode::read_write),
            ityr::make_global_iterator(bodies.begin() + end  , ityr::checkout_mode::read_write),
            [=](Body& B) {
              B.SRC = B.X[0] + I * B.X[1];
            });

#elif EXAFMM_BIOTSAVART

        ityr::for_each(
            body_seq_policy,
            ityr::make_global_iterator(bodies.begin() + begin, ityr::checkout_mode::read_write),
            ityr::make_global_iterator(bodies.begin() + end  , ityr::checkout_mode::read_write),
            [&](Body& B) {
	      for (int d=0; d<3; d++) {                             //   Loop over dimensions
	        B.SRC[d] = drand48() / bodies.size();              //    Initialize source
	      }                                                     //   End loop over dimensions
	      B.SRC[3] = powf(bodies.size() * numSplit, -1./3) * 2 * M_PI * 0.01; // Initialize core radius
            });

#endif
      }                                                         // End loop over partitions
    }

    //! Read source values from file
    void readSources(Bodies & bodies, int mpirank) {
      std::stringstream name;                                   // File name
      name << "source" << std::setfill('0') << std::setw(4)     // Set format
	   << mpirank << ".dat";                                // Create file name
      std::ifstream file(name.str().c_str(),std::ios::in);      // Open file
      file.seekg(filePosition);                                 // Set position in file
      for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {     // Loop over bodies
	file >> B->X[0];                                        //  Read data for x coordinates
	file >> B->X[1];                                        //  Read data for y coordinates
	file >> B->X[2];                                        //  Read data for z coordinates
#if EXAFMM_BIOTSAVART
	file >> B->SRC[0];                                      //  Read data for x source
	file >> B->SRC[1];                                      //  Read data for y source
	file >> B->SRC[2];                                      //  Read data for z source
	file >> B->SRC[3];                                      //  Read data for core radius
#else
	file >> B->SRC;                                         //  Read data for charge
#endif
      }                                                         // End loop over bodies
      filePosition = file.tellg();                              // Get position in file
      file.close();                                             // Close file
    }

    //! Write source values to file
    void writeSources(Bodies & bodies, int mpirank) {
      std::stringstream name;                                   // File name
      name << "source" << std::setfill('0') << std::setw(4)     // Set format
	   << mpirank << ".dat";                                // Create file name
      std::ofstream file(name.str().c_str(),std::ios::out);     // Open file
      for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {     // Loop over bodies
	file << B->X[0] << std::endl;                           //  Write data for x coordinates
	file << B->X[1] << std::endl;                           //  Write data for y coordinates
	file << B->X[2] << std::endl;                           //  Write data for z coordinates
#if EXAFMM_BIOTSAVART
	file << B->SRC[0] << std::endl;                         //  Write data for x source
	file << B->SRC[1] << std::endl;                         //  Write data for y source
	file << B->SRC[2] << std::endl;                         //  Write data for z source
	file << B->SRC[3] << std::endl;                         //  Write data for core radius
#else
	file << B->SRC << std::endl;                            //  Write data for charge
#endif
      }                                                         // End loop over bodies
      file.close();                                             // Close file
    }

    //! Initialize target values
    void initTarget(GBodies bodies) const {
      ityr::for_each(
          body_par_policy,
          ityr::make_global_iterator(bodies.begin(), ityr::checkout_mode::read_write),
          ityr::make_global_iterator(bodies.end()  , ityr::checkout_mode::read_write),
          ityr::count_iterator<int>(0),
          [=](auto&& B, int i) {
            B.TRG = 0;                                             //  Clear target values
            B.IBODY = i;                            //  Initial body numbering
            B.ICELL = 0;                                           //  Initial cell index
            B.WEIGHT = 1;                                          //  Initial weight
          });
    }

    //! Initialize dsitribution, source & target value of bodies
    void initBodies(GBodies bodies, const char * distribution,
	 	    int mpirank=0, int mpisize=1, int numSplit=1) const {
      switch (distribution[0]) {                                // Switch between data distribution type
#if 0
      case 'l':                                                 // Case for lattice
	bodies = lattice(numBodies,mpirank,mpisize);            //  Uniform distribution on [-1,1]^3 lattice
	break;                                                  // End case for lattice
#endif
      case 'c':                                                 // Case for cube
	cube(bodies,mpirank,numSplit);              //  Random distribution in [-1,1]^3 cube
	break;                                                  // End case for cube
#if 0
      case 's':                                                 // Case for sphere
	bodies = sphere(numBodies,mpirank,numSplit);            //  Random distribution on surface of r = 1 sphere
	break;                                                  // End case for sphere
      case 'o':                                                 // Case for octant
	bodies = octant(numBodies,mpirank,numSplit);            //  Random distribution on octant of a r = 1 sphere
	break;                                                  // End case for octant
      case 'p':                                                 // Case plummer
	bodies = plummer(numBodies,mpirank,numSplit);           //  Plummer distribution in a r = M_PI/2 sphere
	break;                                                  // End case for plummer
#endif
      default:                                                  // If none of the above
	fprintf(stderr, "Unknown data distribution %s\n", distribution);// Print error message
      }                                                         // End switch between data distribution type
      initSource(bodies,mpirank,numSplit);                      // Initialize source values
      initTarget(bodies);                                       // Initialize target values
    }

    //! Read target values from file
    void readTargets(Bodies & bodies, int mpirank) {
      std::stringstream name;                                   // File name
      name << "target" << std::setfill('0') << std::setw(4)     // Set format
	   << mpirank << ".dat";                                // Create file name
      std::ifstream file(name.str().c_str(),std::ios::in);      // Open file
      file.seekg(filePosition);                                 // Set position in file
      for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {     // Loop over bodies
	file >> B->TRG[0];                                      //  Read data for potential
	file >> B->TRG[1];                                      //  Read data for x acceleration
	file >> B->TRG[2];                                      //  Read data for y acceleration
	file >> B->TRG[3];                                      //  Read data for z acceleration
      }                                                         // End loop over bodies
      filePosition = file.tellg();                              // Get position in file
      file.close();                                             // Close file
    }

    //! Write target values to file
    void writeTargets(Bodies & bodies, int mpirank) {
      std::stringstream name;                                   // File name
      name << "target" << std::setfill('0') << std::setw(4)     // Set format
	   << mpirank << ".dat";                                // Create file name
      std::ofstream file(name.str().c_str(),std::ios::out);     // Open file
      for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {     // Loop over bodies
	file << B->TRG[0] << std::endl;                         //  Write data for potential
	file << B->TRG[1] << std::endl;                         //  Write data for x acceleration
	file << B->TRG[2] << std::endl;                         //  Write data for y acceleration
	file << B->TRG[3] << std::endl;                         //  Write data for z acceleration
      }                                                         // End loop over bodies
      file.close();                                             // Close file
    }

    //! Downsize target bodies by even sampling
    global_vec<Body> sampleBodies(GBodies bodies, int numTargets) {
      if (numTargets < int(bodies.size())) {                    // If target size is smaller than current
        global_vec<Body> sampled(numTargets);
        int stride = bodies.size() / numTargets;                //  Stride of sampling
        for (int i=0; i<numTargets; i++) {                      //  Loop over target samples
          sampled[i] = bodies[i*stride];                         //   Sample targets
        }                                                       //  End loop over target samples
        return sampled;
      } else {
        global_vec<Body> sampled(bodies.begin(), bodies.end());
        return sampled;
      }
    }
  };
}
#endif
