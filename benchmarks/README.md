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
at 4,096 paths and 4,096 neighbor lists, sharing a budget of 262,144 reference
entries; searches fall back to uncached work when that budget is exhausted.
Pair enumeration and long-chain phase products remain, so this does not remove
all growth with structure size.

The baseline Python symbolic conversion of a 32-segment chain exceeded the
string parser's nesting limit. The numerical benchmark avoids that conversion.
The subsequent symbolic-generation work below removes that failure for the
portable backend by transferring its expression tree directly to SymPy.

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

## Analytical equation generation

Regression tests were added before optimization. The original code failed the
32-segment Python conversion test by returning one unparsed symbol; the smaller
chain, nested, signed-contrast, distributed-reference, and world-edit checks
passed. The final tests also cover a 300-level native expression, expression
modes, parameter tracking, zero-q limits, and shared versus independent sampled
attachment points. Tests compare mathematical results with numerical traversal
and a Gaussian-chain closed form, rather than fixing a printed equation layout.

Generation now memoizes intensities, attachment amplitudes, and phase factors
within one call, including depth and expression mode in the keys. Long sums and
products are built as balanced trees. Path searches retain their original
neighbor order while using a hash set for visited-reference membership. The
portable Python bridge transfers native nodes directly to SymPy, batches
associative sums/products, and reuses converted symbols, constants, and shared
subexpressions. All caches are local to one generation/conversion call. Other
symbolic backends retain their existing conversion path.

Run the public Python benchmark:

```bash
OPENBLAS_NUM_THREADS=1 .venv-sympy/bin/python benchmarks/profile_symbolic.py \
  --output /tmp/pyseb-symbolic
OPENBLAS_NUM_THREADS=1 .venv-sympy/bin/python benchmarks/profile_symbolic.py \
  --sizes 64 --output /tmp/pyseb-symbolic-64
```

Each successful conversion is evaluated and compared with native traversal.
The report explicitly records failed conversions instead of treating their
shorter runtimes as valid measurements. Warmed medians of three repeats on
this machine, before (`7f1c2f2d`) and after:

| Equation | Before: generation + Python conversion | After | Status |
|---|---:|---:|---|
| 8-segment chain | 28.20 ms | 2.24 ms | Equivalent |
| 16-segment chain | 134.28 ms | 8.19 ms | Equivalent |
| 32-segment chain | Failed conversion | 63.98 ms | Now valid |
| Sphere + 8 arms | 24.15 ms | 2.47 ms | Equivalent |
| Sphere + 16 arms | 90.57 ms | 10.63 ms | Equivalent |
| Sphere + 32 arms | Failed conversion | 76.25 ms | Now valid |

The additional 64-segment chain and 64-arm star both converted successfully
and matched native traversal; generation plus Python conversion took 736 ms
and 1515 ms respectively. No before/after speedup is claimed for those cases.

To separate native construction and native string rendering from Python costs:

```bash
g++ -O3 -std=c++14 -DSPECIALFUNCTIONS_IMPL -I. -Iseb -Iseb-symbolic \
  benchmarks/profile_symbolic.cpp build-python/libseb.a build-python/libseb-symbolic.a \
  -lgsl -lgslcblas -lpthread -o /tmp/pyseb-profile-symbolic
/tmp/pyseb-profile-symbolic 32
```

For 32 segments, native generation changed from 79.88 ms to 52.14 ms (about
1.5x faster); native rendering remained around 9 ms. The generated expression
at q=0.1 agreed to approximately 2e-15 relative precision. The larger Python
speedups therefore include improvements in conversion and SymPy construction,
not just the C++ traversal. No scattering approximation or backend switch was
introduced. Timings are observations on this machine, not performance guarantees.

Validation: 123 Python tests passed with one opt-in network test skipped;
all 35 C++ structure/symbolic tests passed in the GiNaC-enabled build, including
new equation-versus-numerical comparisons across both available backends.

## Callgrind and heap-allocation profiling

A second investigation used Valgrind 3.19 Callgrind and DHAT on Linux, starting
from `76880adf`. Callgrind markers in `profile_symbolic.cpp` isolate the four
native generation calls from model construction, rendering, and evaluation.
DHAT measures the complete native driver, including those other phases. These
are different measurement scopes; allocation totals below are not attributed
solely to generation.

The baseline generation profile recorded 1,822,989,986 instructions, with
`findpath` and its callees accounting for about 88%. DHAT identified two path
copying allocation sites totaling about 122 MB across the workload, plus eager
GSL integration workspaces for subunits that never perform quadrature.

The changes preserve seed and neighbor order while using predecessor indices
to reconstruct only a successful breadth-first-search path. Adjacency lists
are reused within the existing call-local context; visited-reference filtering
is still specific to each search. Iterating links by reference avoids copying
string pairs. Cached paths and adjacency share the bounded reference budget.
Integration options are still validated immediately, while workspaces are
allocated on the first nonempty integral and reused until method or workspace
size changes. No scattering approximation or relaxed tolerance is involved.

Compile the marked driver with Valgrind headers available (ordinary builds
without `PYSEB_CALLGRIND` do not need Valgrind):

```bash
g++ -O3 -g -std=c++14 -DSPECIALFUNCTIONS_IMPL -DPYSEB_CALLGRIND \
  -I. -Iseb -Iseb-symbolic benchmarks/profile_symbolic.cpp \
  build-python/libseb.a build-python/libseb-symbolic.a \
  -lgsl -lgslcblas -lpthread -o /tmp/pyseb-symbolic-callgrind
valgrind --vgdb=no --tool=callgrind --instr-atstart=no \
  --callgrind-out-file=/tmp/pyseb-symbolic.callgrind /tmp/pyseb-symbolic-callgrind 32
callgrind_annotate --inclusive=yes /tmp/pyseb-symbolic.callgrind
```

For heap profiling, compile the ordinary native driver as above in the
previous section, then run:

```bash
valgrind --vgdb=no --tool=dhat --dhat-out-file=/tmp/pyseb-symbolic.dhat \
  /tmp/pyseb-profile-symbolic 32
```

Valgrind was downloaded and extracted under `/tmp` for this investigation;
it is not a package dependency or a system installation requirement. Instrumented
runtimes are not used to claim speedups. Ordinary release-build warmed medians
(three repeats, on this machine) changed as follows:

| Measurement | Before | After |
|---|---:|---:|
| Native generation, 32-segment chain | 51.05 ms | 13.20 ms |
| Native generation, 64-segment chain | 627.20 ms | 139.07 ms |
| Generation + Python conversion, 32-arm star | 78.56 ms | 21.33 ms |
| Generation + Python conversion, 64-arm star | 978.14 ms | 155.37 ms |
| Generation instructions, Callgrind, 32 segments | 1,822,989,986 | 472,333,479 |
| Total allocated bytes, DHAT complete driver | 325,642,163 | 139,271,419 |
| Allocation count, DHAT complete driver | 4,884,219 | 1,924,700 |
| Peak live allocated bytes, DHAT complete driver | 43,303,750 | 16,678,982 |

Native string-rendering times were essentially unchanged, and both native
benchmarks returned exactly the same recorded value at q=0.1 and the same
output length. All six 16/32/64-segment/arm Python benchmark curves matched
the saved baseline at `rtol=1e-12, atol=1e-12`. DHAT peak heap describes live
requested allocations, not process RSS. No full fitting speedup is claimed.

Validation: the full Python suite ran 123 tests successfully with its opt-in
network test skipped; all 67 C++ structure, symbolic, numerical, and integration
tests passed. New tests pin endpoint/reversal/branch path ordering and check
integration-option changes before first use and after integration, zero-width
bounds, reversed bounds, and rejection of invalid options.

## Numerical evaluation with Valgrind

`profile_numerical.cpp` isolates complete, warmed vector intensity evaluations
for a Gaussian chain, hydrated atomic cloud, and cloud–sphere hybrid. It uses
31 q points from 0.001 to 0.5. Pass the case, number of segments/atoms, and
repeat count. The cloud cases retain excluded solvent and hydration; the hybrid
retains interference. These synthetic workloads do not cover every subunit or
quadrature method.

Build an ordinary release driver:

```bash
g++ -O3 -std=c++14 -DSPECIALFUNCTIONS_IMPL -I. -Iseb -Iseb-symbolic \
  benchmarks/profile_numerical.cpp build-python/libseb.a build-python/libseb-symbolic.a \
  -lgsl -lgslcblas -lpthread -o /tmp/pyseb-profile-numerical
/tmp/pyseb-profile-numerical chain 16 3
/tmp/pyseb-profile-numerical cloud 128 3
/tmp/pyseb-profile-numerical hybrid 128 3
```

For Callgrind, compile with `-DPYSEB_CALLGRIND` and the Valgrind headers
(add `-I/path/to/valgrind/include` if necessary). Use a different output binary:

```bash
valgrind --vgdb=no --tool=callgrind --instr-atstart=no \
  --callgrind-out-file=/tmp/numerical.callgrind /tmp/pyseb-profile-numerical-callgrind chain 16 1
callgrind_annotate --inclusive=yes /tmp/numerical.callgrind
valgrind --vgdb=no --tool=dhat --dhat-out-file=/tmp/numerical.dhat \
  /tmp/pyseb-profile-numerical chain 16 1
valgrind --vgdb=no --tool=memcheck --leak-check=full \
  --show-leak-kinds=definite,indirect --errors-for-leak-kinds=definite,indirect \
  --error-exitcode=99 /tmp/pyseb-profile-numerical hybrid 128 1
```

Measured with Valgrind 3.19.0, portable Release backend, including the preceding
uncommitted path-search and lazy-workspace improvements:

| Workload | Evaluation-only instructions | Ordinary median time per curve |
| --- | ---: | ---: |
| 16-segment Gaussian chain | 1,177,809,989 | 116.12 ms |
| Hydrated cloud, 128 atoms | 34,544,546 | 2.85 ms |
| Hydrated cloud–sphere hybrid, 128 atoms | 56,247,752 | 4.68 ms |
| Hydrated cloud, 512 atoms | not profiled | 43.52 ms |

Ordinary times are medians of five runs, each measuring three curves after
warm-up, with no concurrent profiling. They are baselines, not speedup claims.
Callgrind measures one curve after warm-up, excluding construction.

Findings:

- The chain spends 47.20% of recorded instructions in
  `evaluateNumericExpression`, including substitution and evaluation of
  temporary expression trees. `CachedPath` accounts for 1.59% including its
  callees. These inclusive percentages overlap with parent frames.
- DHAT reports 202,521,436 cumulative allocated bytes across 3,677,796 blocks
  for the whole chain process (construction, warm-up, one measured curve,
  destruction). Peak live requested heap is only 3,100,004 bytes and exit heap
  is zero. This is allocation churn, not a leak; expression cloning and
  substitution appear prominently in allocation stacks.
- In the 128-atom cloud, sine implementation frames account for 52.66% of
  instructions and the `sinc` wrapper another 11.10% (exclusive counts).
  This reflects the exact Debye pair sum. The 512-atom timing also shows its
  approximately quadratic scaling.
- The hybrid evaluates atomic form factors for both cloud intensity and its
  reference amplitude: the Cromer–Mann evaluation instruction count doubles
  relative to the isolated cloud. Reuse within one evaluation is a potential
  optimization, but must respect solvent, hydration, and parameter changes.
- Memcheck reports zero errors and zero heap blocks at exit for the tested
  hydrated hybrid. This result is limited to that workload.

The next candidates are evaluating analytical expressions without creating
substituted trees, and reusing repeated per-q amplitudes/phases within one
curve evaluation. Both require parity and parameter-mutation tests. Cloud pair
sums need separate attention; histogram approximations would change the
calculation and are not introduced by this profiling work. No additional core
optimization was made in this numerical profiling pass. Keep raw profiles in
`/tmp`, outside version control.

## Direct numerical expression evaluation

The portable backend now evaluates immutable expression trees with parameter
lookups rather than constructing parameter-substituted copies. Subunit numeric
evaluation uses this path; other backends retain their substitution fallback.
The existing arithmetic and special-function dispatch are shared with ordinary
`eval()`. Parameter-map precedence for `q` is retained, and no mutable cache or
cross-call parameter state is introduced.

For the same 16-segment, 31-q numerical chain workload above:

| Measurement | Before | After |
| --- | ---: | ---: |
| Evaluation-only Callgrind instructions | 1,177,809,989 | 711,826,723 |
| Whole-process DHAT allocated bytes | 202,521,436 | 61,019,276 |
| Whole-process DHAT allocated blocks | 3,677,796 | 1,018,030 |
| Peak live requested heap | 3,100,004 bytes | 3,100,004 bytes |

This is 39.6% fewer recorded instructions and 69.9% fewer cumulative allocated
bytes, with an identical measured intensity checksum. DHAT includes startup,
construction, warm-up, measurement and destruction; exit heap remains zero.
Peak heap does not decrease because the transient allocation churn was not the
peak-memory driver. These figures describe the portable backend, not GiNaC.

Direct evaluation is checked against substitution for individual arithmetic and
special-function operations in `Tests/seb_symbolic_test.cpp`, including repeated
parameter changes, missing parameters, and immutability. The numerical Python
regression tests additionally check signed cloud and hybrid intensities against
an independent Debye sum. The existing tiny-nonzero-q sphere cancellation defect
remains an expected failure; this optimization does not change its formula.

Ordinary release timings, five interleaved before/after runs with three curves
per run after warm-up (no concurrent build, tests or Valgrind):

| Workload | Before | After | Speedup |
| --- | ---: | ---: | ---: |
| 16-segment Gaussian chain | 108.89 ms | 71.15 ms | 1.53x |
| 32-segment Gaussian chain | 898.24 ms | 548.71 ms | 1.64x |
| Hydrated cloud, 128 atoms | 2.57 ms | 2.65 ms | 0.97x |
| Hydrated cloud–sphere hybrid, 128 atoms | 4.57 ms | 4.14 ms | 1.10x |

All before/after workload checksums match exactly. The cloud's approximately
3% timing variation is not evidence of an improvement; its pair loop is
unchanged. Results are local observations, not guaranteed speedup factors.
Validation: all 130 Python tests completed successfully (one optional skip and
one existing expected failure), plus 69 C++ tests across symbolic, structure,
and numerical-subunit/integration suites. Per-q reuse and pair-loop changes
remain separate optimization candidates.

## Per-sample reference-value reuse

Numerical vector evaluation now reuses leaf reference amplitudes and phase
factors within each q sample. Each cache is capped at 4096 entries and falls
back to evaluation when full. All numeric entries are cleared before the next
sample, including repeated q values, and the context is discarded at the end
of the call. Thus external callback state, parameter updates and changes to
cloud solvent/hydration settings are picked up on subsequent calls. Reference
names include subunit identity; equal parameter tags do not merge distinct
subunits. Scalar evaluation remains an uncached reference path.

Link scans also use references to stored links instead of copying both strings
for every candidate. No summation order, quadrature tolerance, cloud pair loop,
or hydration model is changed. The callback regression in
`test_numerical_optimization.py` checks external-state and parameter changes,
repeated/unordered q values, and successful evaluation following an error.

As with other scattering evaluators, a callback must return a consistent
mathematical value for fixed q/parameters during one evaluation. Stateful
callbacks may change between evaluations; callback invocation counts are not
an API guarantee.

Release timings against commit `929a8a40`, using five interleaved before/after
runs with three measured 31-q curves per run after warm-up and no concurrent
builds or profiling:

| Workload | Before | After | Speedup |
| --- | ---: | ---: | ---: |
| 16-segment Gaussian chain | 68.34 ms | 19.03 ms | 3.59x |
| 32-segment Gaussian chain | 565.93 ms | 113.03 ms | 5.01x |
| Hydrated cloud, 128 atoms | 3.28 ms | 3.43 ms | 0.95x |
| Hydrated cloud–sphere hybrid, 128 atoms | 5.42 ms | 5.48 ms | 0.99x |

All before/after checksums match exactly. Cloud and small-hybrid timings show
no meaningful improvement; their repeated-reference work is limited. These
measurements describe these particular workloads, not a general speedup claim.

Callgrind records 194,821,031 instructions for the warmed 16-segment curve,
compared with 711,826,723 after direct expression evaluation (72.6% fewer).
DHAT whole-process allocation drops from 61,019,276 bytes in 1,018,030 blocks
to 28,607,536 bytes in 330,698 blocks. Peak live requested heap remains
3,100,004 bytes, and exit heap is zero. Instrumented runtimes are not used for
speedup claims.

Validation completed: 131 Python tests (one optional skip and the existing
expected tiny-q sphere failure) and 69 C++ tests passed. The focused eight-test
numerical optimization suite also passed with that same expected failure.
