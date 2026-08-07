# Cluster installation: ASPECT with FastScape-Fortran glacial erosion

These two branches are intended to be installed together:

- `Minerallo/fastscapelib-fortran:aspect_tester_glacial_erosion`
- `Minerallo/aspect:aspect_posthack26_all_prs`

The FastScape branch is based on Anne Glerum's `aspect_tester` branch and adds
prescribed-field glacial erosion, cumulative-erosion checkpoint restoration,
regression tests, and glacial examples. The ASPECT branch contains the matching
Fortran calls and parameters, the TVD selector, restart serialization, the
empty-surface MPI fix, and the other post-hack ASPECT features.

## 1. Load one compatible compiler/MPI stack

Use the same C, C++, Fortran, and MPI family that was used to build deal.II.
Example module names vary by cluster:

```sh
module purge
module load gcc cmake openmpi
module load deal.II
```

Do not mix an Intel-built deal.II with GNU-built ASPECT/FastScape, or different
OpenMPI installations.

Choose explicit locations in your project or scratch space:

```sh
export SOFTWARE_ROOT=/path/to/your/software
export FASTSCAPE_SOURCE=$SOFTWARE_ROOT/fastscapelib-fortran
export FASTSCAPE_BUILD=$SOFTWARE_ROOT/build-fastscape-glacial
export ASPECT_SOURCE=$SOFTWARE_ROOT/aspect
export ASPECT_BUILD=$SOFTWARE_ROOT/build-aspect-fastscape-glacial
export DEAL_II_DIR=/path/to/deal.II/lib/cmake/deal.II
```

## 2. Build the matching FastScape-Fortran branch

```sh
git clone --branch aspect_tester_glacial_erosion \
  https://github.com/Minerallo/fastscapelib-fortran.git \
  "$FASTSCAPE_SOURCE"

cmake -S "$FASTSCAPE_SOURCE" -B "$FASTSCAPE_BUILD" \
  -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_TESTING=ON \
  -DBUILD_EXAMPLES=ON \
  -DCMAKE_POSITION_INDEPENDENT_CODE=ON

cmake --build "$FASTSCAPE_BUILD" --parallel 8
ctest --test-dir "$FASTSCAPE_BUILD" --output-on-failure
```

Keep the build directory. ASPECT reads its `Makefile` to locate the FastScape
source and version, and links `libfastscapelib_fortran` from this directory.

## 3. Build the matching ASPECT branch

```sh
git clone --branch aspect_posthack26_all_prs \
  https://github.com/Minerallo/aspect.git \
  "$ASPECT_SOURCE"

cmake -S "$ASPECT_SOURCE" -B "$ASPECT_BUILD" \
  -DCMAKE_BUILD_TYPE=Release \
  -DDEAL_II_DIR="$DEAL_II_DIR" \
  -DASPECT_WITH_FASTSCAPE=ON \
  -DFASTSCAPE_DIR="$FASTSCAPE_BUILD"

cmake --build "$ASPECT_BUILD" --target aspect --parallel 8
```

During configuration, verify that CMake reports `ASPECT_WITH_FASTSCAPE = ON`,
finds FastScape version `2.9.1-devASPECT`, and finds the FastScape library.

## 4. Verify the coupling

```sh
cd "$ASPECT_SOURCE/cookbooks/fastscape_eroding_box"
srun -n 2 "$ASPECT_BUILD/aspect" fastscape_glacial_smoke.prm
```

Use `mpirun -np 2` instead of `srun -n 2` on clusters without Slurm. After
this small model succeeds, submit a short reduced-resolution version of the
production model before increasing refinement and MPI rank count.

## 5. Model parameters

The Fortran plugin is selected with `top : fastscape` (not `fastscapecc`). Its
new controls are:

```prm
subsection Mesh deformation
  set Mesh deformation boundary indicators = top : fastscape

  subsection Fastscape
    set Advection scheme = TVD

    subsection Glacial erosion
      set Enable = true
      set Sliding velocity exponent = 1
      set Ice thickness scale = 150
      set Minimum ice thickness = 5
      # Define the three spatial/time functions as in the smoke cookbook.
    end
  end
end
```

This initial glacial law consumes prescribed erodibility, ice thickness, and
basal sliding velocity fields. It does not solve ice dynamics or glacial
sediment transport. The cumulative erosion diagnostic is preserved across
ASPECT checkpoints by this paired branch.
