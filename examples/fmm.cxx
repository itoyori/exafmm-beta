#include "args.h"
#include "bound_box.h"
#include "build_tree.h"
#include "dataset.h"
#include "kernel.h"
#include "logger.h"
#include "namespace.h"
#include "traversal.h"
#include "up_down_pass.h"
#include "verify.h"
using namespace EXAFMM_NAMESPACE;

// ugly hack of making them global variables to have the same virtual address,
// so that all Itoyori threads can access them via raw pointers
Kernel kernel;
Traversal traversal;

void run_fmm(const Args& args) {
  // SPMD region

  const vec3 cycle = 2 * M_PI;
  const real_t eps2 = 0.0;
  const complex_t wavek = complex_t(10.,1.) / real_t(2 * M_PI);

  global_vec<Body> bodies_vec(global_vec_coll_opts);
  global_vec<Body> jbodies_vec(global_vec_coll_opts);
  global_vec<Body> buffer_vec(global_vec_coll_opts);
  global_vec<Cell> cells_vec(global_vec_coll_opts);

  global_vec<BuildTree::BinaryTreeNode> binary_node_vec(global_vec_coll_opts);

  GBodies bodies, jbodies, buffer;
  BoundBox boundBox;
  Bounds bounds;
  BuildTree buildTree(args.ncrit, args.nspawn);
  GCells cells, jcells;
  Dataset data;

  // call constructor
  kernel = {args.P, eps2, wavek};
  traversal = {&kernel, args.theta, args.nspawn, args.images, args.path};

  UpDownPass upDownPass(&kernel);
  Verify verify(args.path);
  num_threads(args.threads);

  verify.verbose = args.verbose;
  logger::verbose = args.verbose;
  logger::path = args.path;

  bodies_vec.resize(args.numBodies);
  bodies = GBodies{bodies_vec.begin(), bodies_vec.end()};

  buffer_vec.resize(args.numBodies);
  buffer = GBodies{buffer_vec.begin(), buffer_vec.end()};

  if (ityr::is_master()) {
    logger::printTitle("FMM Parameters");
    args.print(logger::stringLength);
  }

  if (args.IneJ) {
#if 0
    for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {
      B->X[0] += M_PI;
      B->X[0] *= 0.5;
    }
    jbodies = data.initBodies(args.numBodies, args.distribution, 1);
    for (B_iter B=jbodies.begin(); B!=jbodies.end(); B++) {
      B->X[0] -= M_PI;
      B->X[0] *= 0.5;
    }
#endif
    std::cout << "IneJ unimplemented" << std::endl;
    abort();
  }

  bool pass = true;
  bool isTime = false;
  for (int t=0; t<args.repeat; t++) {
    ityr::root_exec([=] {
      data.initBodies(bodies, args.distribution, 0);
    });

    if (ityr::is_master()) {
      logger::printTitle("FMM Profiling");
      logger::startTimer("Total FMM");
      logger::startPAPI();
      logger::startDAG();
    }
    int numIteration = 1;
    if (isTime) numIteration = 10;
    for (int it=0; it<numIteration; it++) {
      if (ityr::is_master()) {
        std::stringstream title;
        title << "Time average loop " << it;
        logger::printTitle(title.str());
      }

      bounds = boundBox.getBounds(bodies);

      if (args.IneJ) {
#if 0
        bounds = boundBox.getBounds(jbodies, bounds);
#endif
      }

      buildTree.buildTree(bodies, buffer, bounds, cells_vec, binary_node_vec);
      cells = GCells{cells_vec.begin(), cells_vec.end()};

      upDownPass.upwardPass(cells);
      traversal.initListCount(cells);
      traversal.initWeight(cells);

      if (args.IneJ) {
#if 0
        jcells = buildTree.buildTree(jbodies, buffer, bounds);
        upDownPass.upwardPass(jcells);
        traversal.traverse(cells, jcells, cycle, args.dual);
#endif
      } else {
        traversal.traverse(cells, cells, cycle, args.dual);

        jbodies_vec = bodies_vec;
        jbodies = GBodies{jbodies_vec.begin(), jbodies_vec.end()};
      }

      upDownPass.downwardPass(cells);
    }

    auto [should_break, p] = ityr::root_exec([&]() {
      logger::printTitle("Total runtime");
      logger::stopDAG();
      logger::stopPAPI();
      double totalFMM = logger::stopTimer("Total FMM");
      totalFMM /= numIteration;
      logger::resetTimer("Total FMM");
      if (args.write) {
        logger::writeTime();
      }
      traversal.writeList(cells, 0);

      bool should_break = false;
      bool pass_ = true;

      if (!isTime) {
        const int numTargets = 100;

        global_vec<Body> bodies_sampled_vec = data.sampleBodies(bodies, numTargets);
        GBodies bodies_sampled = GBodies{bodies_sampled_vec.begin(), bodies_sampled_vec.end()};

        global_vec<Body> bodies2_vec = bodies_sampled_vec;
        GBodies bodies2 = GBodies{bodies2_vec.begin(), bodies2_vec.end()};

        data.initTarget(bodies_sampled);
        logger::startTimer("Total Direct");
        traversal.direct(bodies_sampled, jbodies, cycle);
        logger::stopTimer("Total Direct");

        auto bs = ityr::make_checkout(bodies_sampled.data(), bodies_sampled.size(), ityr::checkout_mode::read_write);
        auto b2 = ityr::make_checkout(bodies2.data()       , bodies2.size()       , ityr::checkout_mode::read_write);

        Bodies bodies_sampled_(bs.data(), bs.size());
        Bodies bodies2_(b2.data(), b2.size());

        double potDif = verify.getDifScalar(bodies_sampled_, bodies2_);
        double potNrm = verify.getNrmScalar(bodies_sampled_);
        double accDif = verify.getDifVector(bodies_sampled_, bodies2_);
        double accNrm = verify.getNrmVector(bodies_sampled_);
        double potRel = std::sqrt(potDif/potNrm);
        double accRel = std::sqrt(accDif/accNrm);
        logger::printTitle("FMM vs. direct");
        verify.print("Rel. L2 Error (pot)",potRel);
        verify.print("Rel. L2 Error (acc)",accRel);

        buildTree.printTreeData(cells);
        traversal.printTraversalData();
        logger::printPAPI();

        pass_ = verify.regression(args.getKey(), isTime, t, potRel, accRel);

        if (pass_) {
          if (verify.verbose) std::cout << "passed accuracy regression at t: " << t << std::endl;
          if (args.accuracy) should_break = true;
        }
      } else {
        pass_ = verify.regression(args.getKey(), isTime, t, totalFMM);
        if (pass_) {
          if (verify.verbose) std::cout << "passed time regression at t: " << t << std::endl;
          should_break = true;
        }
      }
      return std::make_tuple(should_break, pass_);
    });
    pass = p;

    if (!isTime && pass && !args.accuracy) {
      t = -1;
      isTime = true;
    }

    if (should_break) break;

    ityr::root_exec([=] {
      data.initTarget(bodies);
    });
  }

  if (!pass) {
    if (verify.verbose) {
      if(!isTime) std::cout << "failed accuracy regression" << std::endl;
      else std::cout << "failed time regression" << std::endl;
    }
    abort();
  }
}

int main(int argc, char** argv) {
  Args args(argc, argv);
  ityr::init();

  ityr::common::profiler::event_initializer<prof_event_user_M2L>        ITYR_ANON_VAR;
  ityr::common::profiler::event_initializer<prof_event_user_P2P>        ITYR_ANON_VAR;
  ityr::common::profiler::event_initializer<prof_event_user_M2L_kernel> ITYR_ANON_VAR;
  ityr::common::profiler::event_initializer<prof_event_user_P2P_kernel> ITYR_ANON_VAR;

  if (ityr::is_master()) {
    setlocale(LC_NUMERIC, "en_US.UTF-8");
    printf("=============================================================\n"
           "[ExaFMM]\n"
           "# of processes:               %d\n"
           "-------------------------------------------------------------\n",
           ityr::n_ranks());

    printf("[Compile Options]\n");
    ityr::print_compile_options();
    printf("-------------------------------------------------------------\n");
    printf("[Runtime Options]\n");
    ityr::print_runtime_options();
    printf("=============================================================\n");
    printf("PID of the main worker: %d\n", getpid());
    printf("\n");
    fflush(stdout);
  }

  run_fmm(args);

  ityr::fini();
  return 0;
}
