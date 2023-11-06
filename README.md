# exafmm-beta

A distributed-memory implementation of ExaFMM with Itoyori.

The original exafmm-beta repository: https://github.com/exafmm/exafmm-beta

The current vertion only implements the `laplace`, `helmholtz`, and `biotsavart` kernels and the `cube` particle distribution with Itoyori.
The Itoyori library is already included as a Git submodule at `itoyori/`.
The changes from the original exafmm-beta implementation can be checked in [this commit](https://github.com/itoyori/exafmm-beta/commit/035558d7cc73f94e3ce2f9efd1412232a7077058).

Build:
```sh
git clone --recursive https://github.com/itoyori/exafmm-beta.git
cd exafmm-beta
autoreconf --force --install
./configure CXXFLAGS="-O3 -g -DNDEBUG" --enable-mpi --disable-simd
make
```

- `--disable-simd` is set to avoid some SIMD-related compile errors (depending on environments), but this is not related to Itoyori.
- `-DNDEBUG` is set for faster performance of Itoyori, as it has many assertions on the critical path. Remove `-DNDEBUG` for debugging.

Run:
```sh
# Set the noncollective memory size allocated for each process at program startup
export ITYR_ORI_NONCOLL_ALLOCATOR_SIZE=268435456 # 256 MB/process

# Run exafmm with the laplace kernel
mpirun [MPI options...] -n <nprocs> setarch $(uname -m) --addr-no-randomize ./examples/laplace -v -D -n <n_bodies> [params...]
```

Check [Itoyori's documentation](https://github.com/itoyori/itoyori) for more information, such as how to properly choose MPI implementations for better performance.

In ExaFMM, data for each cell are dynamically allocated with the noncollective allocation policy.
The noncollective memory region is preallocated at program startup by default, and the option `ITYR_ORI_NONCOLL_ALLOCATOR_SIZE` is used to set the maximum allocation size for noncollective memory.
Therefore, this option value should be sufficiently large to accommodate dynamically allocated memory.

The following error suggests that the value of `ITYR_ORI_NONCOLL_ALLOCATOR_SIZE` is too small:
```
[ityr::common::allocator] Could not allocate memory for malloc_local()
```

## Profiling Build

Configure:
```sh
./configure CXXFLAGS="-O3 -g -DNDEBUG -DITYR_PROFILER_MODE=stats -DITYR_ITO_DAG_PROF=workspan -DITYR_ORI_CACHE_PROF=stats" --enable-mpi --disable-simd
```

- The following three profilers are enabled in the above example:
    - `ITYR_PROFILER_MODE=stats` enables profiling for Itoyori's internal events and user-defined events (M2L, P2P).
    - `ITYR_ITO_DAG_PROF=workspan` enables profiling for work and span of the DAG to calculate parallelism.
    - `ITYR_ORI_CACHE_PROF=stats` enables profiling for software caching for global memory access.
- Each profiler can be independently enabled/disabled
- As profiling overhead is not negligible, they should be configured at compile time.

Note: Please run `make clean` after reconfiguraion.

## Enable Almost Deterministic Work Stealing (ADWS)

Configure:
```sh
./configure CXXFLAGS="-O3 -g -DNDEBUG -DITYR_ITO_SCHEDULER=adws -DITYR_ORI_DEFAULT_MEM_MAPPER=block_adws" --enable-mpi --disable-simd
```

- The following compile options should be set to enable ADWS:
    - `ITYR_ITO_SCHEDULER=adws` chooses ADWS as the thread scheduler.
    - `ITYR_ORI_DEFAULT_MEM_MAPPER=block_adws` sets the default memory distribution policy to *ADWS-aware* block distribution.

## Build on Fujitsu A64FX systems

Configure on login node:
```sh
BOOST_ROOT=<path_to_boost>
./configure \
  CXXFLAGS="-Nclang -O3 -g -DNDEBUG -DITYR_ALLOCATOR_USE_BOOST=1 -I${BOOST_ROOT}/include" \
  LDFLAGS="-L${BOOST_ROOT}/lib -Wl,-R${BOOST_ROOT}/lib" \
  LIBS="-lboost_container -ltofucom" \
  --enable-mpi --disable-simd --host=aarch64-linux-gnu
```

- `--host=aarch64-linux-gnu` is needed for cross compilation.
- `mpiFCCpx` should be used for compilation.
- Use clang mode (`-Nclang`) for better performance.
- Boost library is needed because Fujitsu compiler does not support `<memory_resource>` (C++17 feature).
    - Set the Boost installation path to `BOOST_ROOT`.
    - Set `-DITYR_ALLOCATOR_USE_BOOST=1` for Itoyori to use Boost's memory allocator.
    - Specify the boost path with the `-I` and `-L` options (`-Wl,-R` is for setting RPATH).
    - Only Boost containers (`-lboost_container`) is required.
- `-ltofucom` is needed to link to the uTofu library because Itoyori partly uses the low-level uTofu communication layer
