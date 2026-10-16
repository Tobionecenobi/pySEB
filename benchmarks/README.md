# Profiling pySEB

Run the reproducible workloads from the repository root:

```bash
OPENBLAS_NUM_THREADS=1 .venv-sympy/bin/python benchmarks/profile_pyseb.py --output /tmp/pyseb-profile
.venv-sympy/bin/python -m pstats /tmp/pyseb-profile/python.prof
```

After rebuilding with a change, compare against that baseline:

```bash
OPENBLAS_NUM_THREADS=1 .venv-sympy/bin/python benchmarks/profile_pyseb.py \
  --output /tmp/pyseb-profile-after --compare /tmp/pyseb-profile/results.json
```

The script measures warmed native vector evaluation for analytical, connected,
atomic-cloud, and hybrid models, plus cloud radius of gyration. It writes median
and minimum timings and numerical results to `results.json`. Model construction
and the first evaluation are excluded. Check results across builds as well as
timings; timings depend on the compiler, machine load, and model size. These
are diagnostic workloads, not a performance threshold enforced by CI.

Python cProfile identifies time in the bindings but cannot resolve C++ internals.
For a focused native profile of atomic-cloud radius of gyration, build the core
libraries first and then compile the small native driver with gprof instrumentation:

```bash
g++ -O3 -g -pg -std=c++14 -DSPECIALFUNCTIONS_IMPL -I. -Iseb -Iseb-symbolic \
  -c seb/Subunits/AtomicSAXSCloud.cpp -o /tmp/pyseb-cloud-profile.o
g++ -O3 -pg -std=c++14 -DSPECIALFUNCTIONS_IMPL -I. -Iseb -Iseb-symbolic \
  benchmarks/profile_atomic_rg.cpp /tmp/pyseb-cloud-profile.o \
  build-python/libseb.a build-python/libseb-symbolic.a \
  -lgsl -lgslcblas -lpthread -o /tmp/pyseb-cloud-profile
(cd /tmp && ./pyseb-cloud-profile)
gprof /tmp/pyseb-cloud-profile /tmp/gmon.out > /tmp/pyseb-cloud-gprof.txt
```

These commands instrument the cloud implementation, not the whole library;
they use the portable backend and Linux compiler/linker conventions. Adjust
paths and dependencies for another build configuration. Use ordinary release
build timings to assess speed: instrumentation can change optimization and
runtime, and short sampled profiles have limited timing precision.

## First optimization

Baseline: `7c3b5c41`, Release `-O3 -DNDEBUG`, Python 3.11, 101 q values,
seven timed repeats after warm-up. Native profiling showed 26,316,800 calls to
`amplitudeMoment` for 100 Rg evaluations of 512 points. Both atoms' moments were
recomputed for every pair, although they depend only on the atom and current
solvent options. Locally reusing each atom's beta reduces this to 153,600 calls
(including the existing total-beta and contrast checks). Pair distances and
summation order remain unchanged. No approximation or persistent result cache
is introduced; each call uses the current solvent and hydration settings.

Measured medians on this machine (seven repeats, warmed geometry):

| Operation | Before | After | Speedup |
|---|---:|---:|---:|
| Hydrated 128-point cloud Rg squared | 0.703 ms | 0.029 ms | 24.1x |
| Hydrated 512-point cloud Rg squared | 12.057 ms | 0.234 ms | 51.6x |

All seven workload outputs matched the baseline at `rtol=1e-12, atol=1e-12`.
The 42 atomic-cloud, hydration-reference, forward-matching, evaluation, and
validation tests passed, including a new low-q check after signed-contrast,
solvent-density, and hydration changes. The after timings were taken without
the test suite running concurrently. This change targets Rg calculations;
timing differences in unmodified intensity workloads are not evidence that
their algorithms became faster. No full fitting speedup is claimed.

For context, after-change vector intensity medians were 2.50 ms for a sphere,
106.90 ms for eight connected Gaussian polymers, 179.39 ms for a 512-point
hydrated cloud, and 5.04 ms for the small hybrid, all on 101 q values.

The remaining cloud intensity calculation still evaluates the direct Debye pair
sum for each q. The connected analytical example also remains relatively costly.
Further work should profile those native paths before changing algorithms;
histogram acceleration would introduce a separately controlled approximation.

## Large analytical structures

Measure chains and sphere-with-polymer-arm stars without symbolic conversion:

```bash
OPENBLAS_NUM_THREADS=1 .venv-sympy/bin/python benchmarks/profile_analytical.py \
  --sizes 8 16 32 --q-count 31 --repeats 3 --output /tmp/pyseb-analytical
```

Pass `--compare /path/to/earlier/results.json` after rebuilding to check every
curve against an earlier run at `rtol=1e-12, atol=1e-12`. The star's size counts
arms (there is also one sphere); the chain's size counts polymer subunits.
Construction is measured separately. These models deliberately use shared
parameter tags, so increasing size changes topology rather than parameter count.

For native traversal profiling, the `profile_analytical.cpp` driver evaluates a
32-segment chain. Instrument World and Subunit separately to avoid rebuilding
all dependencies with profiling enabled:

```bash
g++ -O3 -g -pg -std=c++14 -DSPECIALFUNCTIONS_IMPL -I. -Iseb -Iseb-symbolic \
  -Ithird_party -c seb/World.cpp -o /tmp/pyseb-profile-world.o
g++ -O3 -g -pg -std=c++14 -DSPECIALFUNCTIONS_IMPL -I. -Iseb -Iseb-symbolic \
  -c seb/Subunit.cpp -o /tmp/pyseb-profile-subunit.o
g++ -O3 -pg -std=c++14 -DSPECIALFUNCTIONS_IMPL -I. -Iseb -Iseb-symbolic \
  benchmarks/profile_analytical.cpp /tmp/pyseb-profile-world.o \
  /tmp/pyseb-profile-subunit.o build-python/libseb.a build-python/libseb-symbolic.a \
  -lgsl -lgslcblas -lpthread -o /tmp/pyseb-profile-analytical
(cd /tmp && ./pyseb-profile-analytical 32)
gprof -b /tmp/pyseb-profile-analytical /tmp/gmon.out > /tmp/pyseb-analytical-gprof.txt
```

The baseline profile made 61,504 `findpath` calls for four 31-q evaluations:
496 distinct subunit pairs, with an identical search repeated at every q.
`findpath`, `getNeighbors` and `isVisited` together accounted for roughly 45%
of sampled self time in this partially instrumented binary. Additional search
costs appear in string handling, maps, and allocations. These percentages are
specific to this workload, not a whole-library performance claim.

Vector intensity, amplitude, and phase evaluation now reuse topology paths
within each call. Numeric values are still evaluated at each q and parameter
set. The cache is discarded on return; it does not require invalidation across
World edits or introduce shared mutable evaluation state. Storage is bounded
at 4,096 paths and 262,144 reference entries; uncached paths use the original
search when that budget is exhausted. Pair enumeration and long-chain phase
products remain, so this does not remove all growth with structure size.

A separate exploratory test found that the Python symbolic conversion of the
32-segment chain exceeds the string parser's nesting limit. The native vector
benchmark avoids that conversion. Reusing or compiling a complete symbolic
expression is therefore not yet a reliable general workaround for large
structures; that conversion limitation remains separate work.

The after profile reduced `findpath` calls from 61,504 to 1,984, exactly one
search per pair per vector call. All six native benchmark curves matched the
saved baseline at `rtol=1e-12, atol=1e-12`; all 54 relevant Python tests passed.
A new regression compares vector intensity, amplitude, and phase against scalar
calls for repeated, negative, and zero q after parameter and topology changes.

For a same-build timing cross-check, add `--scalar-reference` to the command.
This times uncached scalar traversal as well as optimized vector traversal;
the scalar loop follows the old vector algorithm but includes extra Python
call overhead. Warmed medians for three repeats on 31 q values were:

| Model | Uncached scalar loop | Optimized vector | Ratio |
|---|---:|---:|---:|
| 8 polymer segments | 29.34 ms | 19.48 ms | 1.51x |
| 16 polymer segments | 254.43 ms | 112.74 ms | 2.26x |
| 32 polymer segments | 2616.18 ms | 774.65 ms | 3.38x |
| Sphere + 8 arms | 83.66 ms | 39.79 ms | 2.10x |
| Sphere + 16 arms | 585.75 ms | 161.27 ms | 3.63x |
| Sphere + 32 arms | 7631.87 ms | 754.50 ms | 10.12x |

These are same-machine observations, not universal speedup guarantees. The
original before-build vector measurements were 3128.65 ms for the 32-segment
chain and 10403.04 ms for the 32-arm star; separate-run timing variation is why
the same-build scalar comparison is also recorded. No full optimizer timing
was measured. This optimization applies automatically to native vector calls;
it does not require a generated kernel or a change of symbolic backend.
