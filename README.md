# exafmm-beta

A distributed-memory implementation of ExaFMM with Itoyori.

The original exafmm-beta repository: https://github.com/exafmm/exafmm-beta

The current vertion only implements the `laplace`, `helmholtz`, and `biotsavart` kernels and the `cube` particle distribution with Itoyori.

Build:
```sh
git clone --recursive https://github.com/itoyori/exafmm-beta.git
autoreconf --force --install
./configure CXXFLAGS="-O3 -g" --enable-mpi --disable-simd
make
```

Note that Itoyori is already included as a Git submodule at `itoyori/`.

Run:
```sh
# Set the noncollective memory size allocated for each process at program startup
export ITYR_ORI_NONCOLL_ALLOCATOR_SIZE=268435456 # 256 MB/process

mpirun <MPI options...> -n <nprocs> setarch $(uname -m) --addr-no-randomize ./examples/laplace -v -D -n <n_bodies>
```

Check [Itoyori's documentation](https://github.com/itoyori/itoyori) for more information, such as how to properly choose MPI implementations for better performance.

In ExaFMM, data for each cell are dynamically allocated with the noncollective allocation policy.
The noncollective memory region is preallocated at program startup by default, and the option `ITYR_ORI_NONCOLL_ALLOCATOR_SIZE` is used to set the maximum allocation size for noncollective memory.
Therefore, this option value should be sufficiently large to accommodate dynamically allocated memory.

The following error suggests that the value of `ITYR_ORI_NONCOLL_ALLOCATOR_SIZE` is too small:
```
[ityr::common::allocator] Could not allocate memory for malloc_local()
```
