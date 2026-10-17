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

## Call-local atomic form factors

Atomic-cloud intensity and reference-amplitude calculations now compute the
common hydrogen and hydration-water form factors once per operation at each q,
reusing them across atoms. Geometry/options validation occurs once before the
internal atom loop. The public single-atom function retains validation and uses
its own fresh context. Shared factors are evaluated lazily, so unused or
inactive hydration points do not trigger water normalization. No factors are
stored on the cloud or shared between calls, and per-atom vacuum coefficients,
occupancy, temperature factors, solvent displacement and hydrogen distances
remain independent. The exact pair loop and compensated summation are unchanged.

Release measurements against `07dc5bad`: seven interleaved before/after runs,
31 q points, with 100 measured curves per run for eight atoms, 30 for 32 atoms,
10 for 128 atoms and three for 512 atoms. Construction/warm-up are excluded;
no builds, tests or profilers ran concurrently with timing.

| Workload | Before | After | Speedup |
| --- | ---: | ---: | ---: |
| Hydrated cloud, 8 atoms | 0.0493 ms | 0.0314 ms | 1.57x |
| Hydrated cloud, 32 atoms | 0.3079 ms | 0.2178 ms | 1.41x |
| Hydrated cloud, 128 atoms | 3.1284 ms | 3.0032 ms | 1.04x |
| Hydrated cloud, 512 atoms | 45.2525 ms | 45.0645 ms | 1.00x |
| Hydrated cloud–sphere hybrid, 8 atoms | 1.1180 ms | 1.0676 ms | 1.05x |
| Hydrated cloud–sphere hybrid, 128 atoms | 4.7311 ms | 4.1478 ms | 1.14x |

All workload checksums match exactly. Gains are strongest where atom amplitude
work is appreciable; the large-cloud pair sum is unaffected. Small differences
near 1.00x should not be treated as evidence of a speedup.

Callgrind's warmed 32-atom curve records 3,422,732 instructions before and
2,373,847 after (30.6% fewer). Profiles and timing JSON are kept under `/tmp`.
The new regression compares complete-cloud intensity and reference amplitude
against individually evaluated atoms, exercising every hydration-water model,
solvent-density changes, occupancy, temperature factors, implicit hydrogens,
and repeated/near-zero q. Existing independent Debye and hybrid tests remain.

Validation: 132 Python tests completed successfully (one optional skip and the
existing expected tiny-q sphere failure); all 69 C++ tests passed. Memcheck on
the hydrated 128-atom cloud–sphere hybrid reports zero errors and zero heap
blocks at exit. The focused nine-test numerical regression suite also passed
with the same expected failure.

## Preparing exported analytical equations for repeated evaluation

The public `evaluate_expression` helper reparses, substitutes parameter values,
and generates a NumPy function on each call. The optional explicit preparation
path performs code generation once and accepts new scalar bindings each time:

```python
# Keep the equation available for inspection/export; prepare it once for a fit.
equation = world.FormFactor_Unnormalized('particle')
evaluate = pyseb.compile_expression(equation)
curve = evaluate(parameters, q)
# An optimizer may now change parameters without regenerating numerical code.
next_curve = evaluate(updated_parameters, q)
```

`compile_expression` generates reusable NumPy code with SymPy common
subexpression elimination. It is not a machine-code C++ kernel or a replacement
for native evaluation. Every non-q free symbol is rebound on every call;
`evaluate.parameter_names` lists required bindings. Extra dictionary entries
are ignored, and the q argument supplies q. Outputs preserve scalar/array/empty
shapes, including constant expressions. Missing/non-finite inputs and non-finite
outputs raise errors. SymPy's generated code must only be prepared from trusted
expressions. SciPy is optional and extends special-function support.

Unresolved integrals, derivatives, sums/products, and unknown symbolic functions
are rejected. Use native World evaluation for quadrature, numerical clouds and
hybrids containing numerical subunits. Generic equations carry no structural
metadata to resolve removable q=0 singularities; a non-finite result is reported
rather than silently replaced. The existing one-shot helper is unchanged.

Reproduce the measurement, including generation, one-time preparation, repeated
parameter changes, and comparison with native vector results:

```bash
OPENBLAS_NUM_THREADS=1 .venv-sympy/bin/python \
  benchmarks/profile_equation_evaluation.py --output /tmp/pyseb-equations.json
```

Observed medians for seven changed-parameter evaluations, 101 q values from
0.03 to 0.5, portable backend, after warming optional-library imports:

| Analytical model | Existing helper | Prepared evaluator | One-time preparation |
| --- | ---: | ---: | ---: |
| Gaussian chain, 8 segments | 69.59 ms | 0.149 ms | 28.14 ms |
| Gaussian chain, 16 segments | 66.49 ms | 0.161 ms | 32.74 ms |
| Sphere + 8 Gaussian arms | 14.44 ms | 0.135 ms | 10.85 ms |
| Sphere + 16 Gaussian arms | 11.54 ms | 0.091 ms | 7.13 ms |

These local timings reflect reused equations with shared polymer parameters,
not general performance guarantees. Generation is recorded separately in the
JSON and is not included in preparation or curve times. Every changed-parameter
curve agrees with native evaluation at rtol=2e-8, atol=1e-8. Tests additionally
cover scalar/empty/multidimensional shapes, constants, special functions,
missing/rebound parameters, invalid inputs, and preparation without SciPy.

Validation: all eight prepared-expression tests pass; the complete Python suite
runs 140 tests successfully with one optional skip and the existing expected
tiny-q sphere failure. No C++ engine changes are introduced by this helper.

## Large-cloud q-blocking experiment (not retained)

An experiment evaluated exact Debye sums in blocks of 4, 8 and 16 q values,
sharing geometry reads while retaining each q's original pair order, sinc
function and compensated summation. Both q-major and atom-major amplitude
layouts were tested. A prototype's apparent improvement did not translate into
consistent gains with the production Release flags (`-O3 -DNDEBUG -fPIC`).

Nine interleaved runs of five warmed 31-q curves, alternating scalar/vector
execution order, gave these results for the eight-q, atom-major variant:

| Workload | Scalar reference | Experimental q blocks | Speedup |
| --- | ---: | ---: | ---: |
| Synthetic grid, 512 atoms | 37.71 ms | 42.04 ms | 0.90x |
| 1TIU model 1, 688 heavy atoms | 85.93 ms | 85.36 ms | 1.01x |
| Same protein linked to a sphere | 88.93 ms | 87.04 ms | 1.02x |

Every returned curve value matched the scalar reference exactly. The synthetic
regression and negligible protein gain did not justify keeping the engine/API
changes, so they were removed. This is a negative experiment, not a delivered
speedup. It does not rule out SIMD sine implementations or parallel evaluation,
which would need separate accuracy, portability and threading work.

`profile_cloud_batch.cpp` retains a baseline workload for future experiments.
It loads the checked-in protein file explicitly as model 1, without merging NMR
conformations, and includes implicit hydrogens, solvent subtraction and hydration.
Construction/SASA and warm-up are outside the measurement. Run from the repo root:

```bash
g++ -O3 -std=c++14 -DSPECIALFUNCTIONS_IMPL -I. -Iseb -Iseb-symbolic \
  benchmarks/profile_cloud_batch.cpp build-python/libseb.a build-python/libseb-symbolic.a \
  -lgsl -lgslcblas -lpthread -o /tmp/pyseb-cloud-baseline
/tmp/pyseb-cloud-baseline scalar protein 5
/tmp/pyseb-cloud-baseline vector protein 5
/tmp/pyseb-cloud-baseline vector protein 5 hybrid
/tmp/pyseb-cloud-baseline vector 512 5
```

These commands measure the retained native evaluator, not the discarded q-block
implementation. The experimental patch and timing JSON are under `/tmp` for
this session. The new large-cloud test checks scalar/vector equality around
potential block boundaries, normalization, solvent/hydration updates, signed
contrasts, invalid q, custom callbacks and hybrid interference.

After removing the experiment and rebuilding, all 141 Python tests complete
successfully (one optional skip and the existing expected tiny-q sphere
failure). The focused ten-test numerical regression suite passes with that
same expected failure. Baseline scalar/vector curves agree exactly for both
512 synthetic points and 688 protein atoms, including their sphere hybrids.

## Independently parameterized analytical structures

`profile_symbolic_stages.cpp` embeds Python solely to call the production
converter, separating native equation generation from SymPy conversion without
adding profiling hooks to the public API. Models have different Rg and signed
beta parameters per polymer. It checks intensity, normalized reference amplitude
and phase against native evaluation at q=0.03, 0.1 and 0.3. Nested cases contain
two distinct polymers per block; no parameters are shared between blocks.

Profiling identified repeated BFS searches and linear link scans as significant
native-generation costs. Evaluation contexts now retain ordered BFS predecessor
trees for repeatedly queried sources and a sorted index of actual links. The
first two queries from a source keep the original early-exit search; subsequent
queries can reuse a full search. Seed order, neighbor order, reversed searches,
coalesced empty paths and error fallback remain unchanged. Public scalar path
queries still use the original algorithm.

Search reuse is limited to 64 sources and 4096 nodes per search. Its stored
reference/object strings share the existing 262144-entry path/adjacency budget.
Link indexing is limited to worlds with at most 4096 links. All caches belong
to one generation/vector-evaluation call; larger requests fall back to the
existing search, and later World edits get a fresh context. No expression
algebra, scattering weights or expression representation is changed.

Native intensity-generation timings against `99a5b4f7`, portable Release
backend, median of seven calls after warm-up, with no concurrent builds/tests:

| Structure | Before | After | Speedup |
| --- | ---: | ---: | ---: |
| Chain, 32 independent polymers | 10.48 ms | 6.54 ms | 1.60x |
| Chain, 64 independent polymers | 109.45 ms | 40.62 ms | 2.69x |
| Sphere + 32 independent arms | 17.20 ms | 8.78 ms | 1.96x |
| Sphere + 64 independent arms | 160.69 ms | 38.35 ms | 4.19x |
| 16 nested blocks / 32 polymers | 4.97 ms | 3.61 ms | 1.38x |
| 32 nested blocks / 64 polymers | 34.00 ms | 19.67 ms | 1.73x |

Amplitude generation also improved on larger models (e.g. 5.68 to 2.35 ms for
64 arms); phase generation was roughly unchanged. These are native-stage
measurements, not end-to-end Python speedups. First-call SymPy canonicalization
can still dominate the public call. Conversion is measured separately and is
not optimized by this change. Extra converter memoization and shared-root
experiments were not retained because their timings were inconsistent.

Build and run from the repository root (match python3-config to the Python
version in the selected environment):

```bash
g++ -O3 -DNDEBUG -std=c++14 -DSPECIALFUNCTIONS_IMPL \
  -I. -Iseb -Iseb-symbolic -Ipyseb/bindings \
  $(.venv-sympy/bin/python -m pybind11 --includes) \
  benchmarks/profile_symbolic_stages.cpp pyseb/bindings/bindingsSymbolic.cpp \
  build-python/libseb.a build-python/libseb-symbolic.a \
  -lgsl -lgslcblas -lpthread $(python3-config --embed --ldflags) \
  -o /tmp/pyseb-symbolic-stages
OPENBLAS_NUM_THREADS=1 /tmp/pyseb-symbolic-stages chain 64 \
  /home/tobi/pySEB/.venv-sympy/lib/python3.11/site-packages 7
# Also use 'star 64' and 'nested 32' for comparable large structures.
```

JSON lines record both first-call and warmed timings. Compile with
`PYSEB_CALLGRIND` and the Valgrind include directory to isolate native generation
using Callgrind's `--instr-atstart=no`; instrumented elapsed times must not be
used as speedup claims. Reports are stored under `/tmp/pyseb-symbolic-final-*`.
The new tests exercise independently parameterized nested blocks, native
operation/operand order, topology edits after repeated searches, and a 70-arm
star beyond the 64-source cache limit using an independent pair-count reference.

Validation: the full Python suite passes all 144 tests with one optional skip
and the existing expected tiny-q sphere failure. All 69 C++ tests across the
GiNaC-enabled symbolic, structure and numerical/integration suites pass. No
converter-cache or expression-root-sharing changes remain in the implementation.

## SymPy conversion and downstream costs

The stage driver also measures CSE-enabled `lambdify` preparation, first plain
text and LaTeX rendering, output byte counts, and warmed 31-q evaluation. These
stages follow conversion and are timed separately. Downstream timings are
single observations; conversion reports a first call and a median of repeated
calls. Every operation is checked against native evaluation at three q values.
The JSON records the installed SymPy version and cache capacity. An optional
sixth argument `profile` instruments only the first conversion with cProfile
and writes `/tmp/pyseb-conversion-<topology>-<size>-<operation>.pstats`:

```bash
# Build /tmp/pyseb-symbolic-stages using the command above.
OPENBLAS_NUM_THREADS=1 SYMPY_CACHE_SIZE=1000 /tmp/pyseb-symbolic-stages \
  chain 64 /home/tobi/pySEB/.venv-sympy/lib/python3.11/site-packages 3 profile
.venv-sympy/bin/python -c \
  "import pstats; pstats.Stats('/tmp/pyseb-conversion-chain-64-intensity.pstats').strip_dirs().sort_stats('cumulative').print_stats(25)"
```

Instrumented first-call times are explicitly marked and must not be compared
with ordinary timings. Optional argument `simplify` instead measures general
SymPy simplification; use small cases (e.g. `chain 2`) because large-expression
simplification has no practical runtime bound. A value of -1 means it was not
requested. Profiling and simplification are opt-in and do not affect the library.

A cProfile run of the independently parameterized 64-polymer intensity found
8.60 million calls (including recursive calls). SymPy `Mul.flatten` accounted
for roughly 79% of the instrumented public call's cumulative time. Factor
sorting/comparison and exponent gathering were prominent underneath it;
these inclusive costs overlap. This is canonical construction of scalar
products, rather than an explicit call to general `simplify` in our converter.

Fresh-process observations with SymPy 1.14.0 on the portable Release backend, with no concurrent
builds, tests or profiling during the measurements:

| 64-polymer chain intensity | SymPy cache 1000 | SymPy cache 10000 |
| --- | ---: | ---: |
| Native generation, warmed median | 35.1 ms | 35.9 ms |
| First conversion | 1637 ms | 1466 ms |
| Repeated conversion, median of three | 1468 ms | 61.3 ms |
| CSE evaluator preparation | 9.56 s | 7.60 s |
| Plain text rendering | 6.45 s | 5.11 s |
| LaTeX rendering | 8.05 s | 7.46 s |
| Warmed 31-q evaluation, mean of 100 | 2.94 ms | 2.83 ms |

Both runs produced exactly the same recorded validation values and output
lengths: 1,083,165 plain-text bytes and 1,350,717 LaTeX bytes for the intensity.
Cache size is a process-start SymPy setting, not a pySEB option; to repeat the
second run, set `SYMPY_CACHE_SIZE=10000` before launching Python/the driver.
A larger cache helps repeated extraction substantially in this case, but it
does not eliminate first-time construction, retains more cached objects, and
is not automatically enabled by pySEB. The small first-call/downstream timing
differences are not established speedups from a replicated experiment.
For comparison, the default-cache 32-polymer chain converted in 324 ms on its
first call and 12.2 ms on repeated calls: warmed timings alone conceal the
cache-capacity threshold. These benchmarks use independent parameter tags;
shared-tag chains are much easier for SymPy to combine.

An experimental converter combined exponential factors using the exact scalar
identity `exp(a)*exp(b) = exp(a+b)` before constructing each product. It was
rejected: the 64-polymer case became much slower and was interrupted rather
than accepted as an extraction speedup. Both that converter and its rebuilt
Python extension were restored to the committed implementation. No unevaluated
operations or process-global cache changes were retained. The experimental
patch and raw profiling/timing reports are session artifacts under `/tmp`.

The retained change is measurement coverage. First-time conversion remains a
bottleneck; generating a more compact, factorized native equation is a candidate
for future work, with independent checks of algebra, parameter changes,
preparation, rendering and numerical evaluation. Deferred canonicalization by
itself would not establish an end-to-end improvement.

Validation: the ten symbolic-generation regression tests pass after restoring
the original converter. The extended driver compiles and validates intensity,
amplitude and phase on chain, star and nested examples, including its optional
cProfile and small-case simplification modes. No scattering-engine code changed.

## OpenUSD export stages and exact-output optimization

The native `profile_usd.cpp` driver separates model construction from warmed
exports. It covers a 64-cylinder analytical chain, a 16-cylinder chain with
readable layout, a 512-point Debye-cloud/sphere hybrid with and without an
envelope, and a 4096-atom/sphere hybrid. The default exporter quality settings
are retained: seed 42, 96 curve samples, 48 surface samples, 64 orientation
trials, four relaxation sweeps and envelope resolution 24. Cloud coordinates,
mixed Debye radii and signed weights are deterministic; no external data,
scattering evaluations or fitting are included in the timed exports.

`build_usd_profile.py` creates an instrumented copy of the selected exporter
under `/tmp`. No profiling callbacks, timers or public API additions enter the
production build. Exact source anchors fail if the exporter changes. The copied
exporter defines the USD symbols, so linking the regular portable static library
does not extract its USD object. Build the normal portable Release core first,
then run from the repository root:

```bash
git show 98ba5ca2:seb/USDExport.cpp > /tmp/pyseb-usd-baseline.cpp
python3 benchmarks/build_usd_profile.py \
  --source /tmp/pyseb-usd-baseline.cpp --build-dir /tmp/pyseb-usd-before
python3 benchmarks/build_usd_profile.py --build-dir /tmp/pyseb-usd-after
python3 benchmarks/profile_usd.py \
  --before /tmp/pyseb-usd-before/profile-usd \
  --after /tmp/pyseb-usd-after/profile-usd \
  --output /tmp/pyseb-usd-comparison --repeats 5
```

The comparison alternates before/after order over five fresh processes per
version and case, with one warm-up and one measured export per process. It
checks SHA-256 equality for every output, within and between versions, and
saves individual samples, stage medians and representative scenes. Run without
concurrent builds, tests or other benchmarks. The native driver also accepts
individual cases, for example:

```bash
/tmp/pyseb-usd-after/profile-usd analytical 32 readable points 3 /tmp/chain.usda
/tmp/pyseb-usd-after/profile-usd debye 512 random envelope 3 /tmp/cloud.usda
```

The non-overlapping stages are preparation (occurrence traversal, parameter
binding, analytical geometry sampling and proxy construction), layout,
Debye-envelope meshing, remaining geometry/USDA serialization, file writing
(including temporary-file handling and rename), and cleanup. Reference sampling
performed lazily during placement is included in layout. Envelope meshing runs
inside serialization; its timer is subtracted from that stage. These boundaries
measure implementation stages, not hypothetical independent geometry and
formatting pipelines. Filesystem writes are warmed/buffered, not forced to
stable storage. Independent stage medians need not sum to the median total.

The measured bottleneck was numeric formatting during serialization: the old
exporter constructed a classic-locale, 17-digit string stream and temporary
string for each coordinate and floating-point field. It now validates each
number and writes it directly into the export-local stream, configured once
with exactly the same locale and precision. The final file transfer reads that
stream's buffer instead of allocating a second full USDA string. No mutable
formatting state is shared between exports. Layout algorithms, random draws,
sampling quality, meshes, instancing, references and atomic metadata are unchanged.

Portable Release measurements (g++ `-O3 -DNDEBUG -fPIC`, five interleaved
fresh-process samples per version; medians on this machine):

| Scene | Bytes | Before | After | Speedup |
| --- | ---: | ---: | ---: | ---: |
| 64-cylinder chain, random layout | 35,382,510 | 1.227 s | 0.717 s | 1.71x |
| 16-cylinder chain, readable layout | 8,845,178 | 0.907 s | 0.714 s | 1.27x |
| 512-point Debye hybrid, envelope | 683,178 | 37.9 ms | 26.5 ms | 1.43x |
| 512-point Debye hybrid, points only | 309,372 | 15.3 ms | 8.60 ms | 1.78x |
| 4096-atom hybrid | 1,032,161 | 53.0 ms | 27.3 ms | 1.94x |

Serialization of the 64-cylinder chain fell from 1.053 s to 0.519 s. In the
readable chain, layout stayed essentially unchanged (480 ms versus 481 ms),
while serialization fell from 355 ms to 182 ms. Debye-envelope meshing also
stayed essentially unchanged (10.1 ms versus 9.92 ms). The observed improvement
comes primarily from serialization, not from faster placement or surface
construction. Local timings vary; these are not universal speedup guarantees.

On Linux the driver reports process peak RSS in KiB, including construction,
warm-up and measured exports. Removing the document copy did not materially
reduce peak RSS in these cases: approximately 101 versus 97 MiB for the large
chain, and 266 MiB for both atomic-hybrid runs. Buffer-growth peaks and cloud
construction can dominate this process-wide measure. This is not evidence of
a large peak-memory improvement or a streaming exporter; the complete USDA
buffer is still retained until writing finishes.

All five scene files are byte-identical before and after the optimization.
They also load with OpenUSD 0.26.8; mesh face counts/indices and point-instancer
prototype targets/array lengths were checked. All 17 Python USD-export tests
pass, including existing nested hierarchy, deterministic placement, readable
clearance and attachment checks. The new regression covers 17-digit/scientific
formatting, signed zero, rejection of non-finite fields without replacing an
existing scene, and recovery after that error. Generated scenes, profiles and
raw JSON remain under `/tmp`, outside version control.

The complete Python suite also passes: 145 tests, with the existing opt-in
network skip and known expected tiny-q sphere failure. The portable Python
extension was rebuilt with the retained exporter changes before validation.

## Unevaluated SymPy construction experiment (not retained)

Two converter variants were tested without changing native generation or
special-function mappings:

1. Construct every collected `Add` and `Mul` with `evaluate=False`.
2. Keep normal construction except for `Mul` with more than eight collected
   arguments, where `evaluate=False` is used.

The existing `profile_symbolic_stages.cpp` driver measures conversion separately
from CSE-enabled evaluator preparation, plain text/LaTeX printing, and warmed
31-q evaluation. Both variants were compiled against the same restored
`3f969b1c` core library and run on independent-parameter chains, stars and
nested structures with signed contrasts. Each operation is checked against
native evaluation at q=0.03, 0.1 and 0.3. Recorded baseline/experimental q=0.3
values also agree within relative/absolute 1e-10.

Screening observations with SymPy 1.14.0 and cache capacity 1000, one fresh
process per version/case, alternating version order between cases, with three
warmed conversion calls and 100 curve evaluations. No builds, tests or other
profilers ran concurrently. These are individual observations, not replicated
speedup guarantees. This first screen uses the stage driver's original
`dummify=False` evaluator setting:

| Intensity case | First conversion, normal | First conversion, all unevaluated | Compiled curve, normal | Compiled curve, all unevaluated |
| --- | ---: | ---: | ---: | ---: |
| Chain, 32 polymers | 295 ms | 83.4 ms | 0.793 ms | 2.719 ms |
| Sphere + 32 arms | 181 ms | 75.7 ms | 0.684 ms | 1.101 ms |
| 16 nested blocks / 32 polymers | 132 ms | 72.6 ms | 0.435 ms | 0.849 ms |
| Chain, 64 polymers | 1274 ms | 264 ms | 2.583 ms | 19.378 ms |

For the 64-polymer intensity, evaluator preparation increased from 6.47 to
8.14 seconds. Plain-text rendering was essentially unchanged (4.40 versus
4.32 seconds), as was LaTeX rendering (6.13 versus 6.13 seconds). The text grew
from 1,083,165 to 1,106,033 bytes. Thus the faster initial transfer left more
work for preparation and repeated evaluation, rather than compacting the
expression. CSE and `docstring_limit=0` remained enabled throughout.

The all-unevaluated variant also failed two existing structural/canonical-form
assertions: a sum of 300 identical symbols stayed a repeated sum rather than
`300*x`, and a constant expression stayed `1 + 1 + 2*1*1*1` rather than `4`.
These are changes in representation, not demonstrated numerical errors. The
other eight symbolic-generation tests passed.

The long-product-only variant passed all ten symbolic-generation tests, but
still made the 64-polymer evaluator slower: first conversion 462 ms,
preparation 8.46 seconds and repeated curves 19.59 ms. The 32-arm star was
roughly unchanged, while the smaller chain and nested cases also had slower
evaluation. No beneficial general default was established.

Neither variant was retained. The converter and rebuilt Python extension use
the original evaluated constructors. Reproduce the broad experiment by replacing
only the converter's collected `Add`/`Mul` constructor call with its
`evaluate=False` equivalent, then compiling two stage-driver binaries using the
command above, keeping their core libraries identical. The restricted variant
applies that keyword only to products with more than eight arguments. The
experimental source copies, patches, raw JSON and test logs are under
`/tmp/pyseb-unevaluated-*`; no experimental setting or new dependency is exposed
in pySEB.

A follow-up uses the same `dummify=True`, CSE and docstring settings as the
public `compile_expression` helper's numerical kernel. Parameter renaming did
not recover the performance: for the 64-polymer intensity, ordinary conversion
plus kernel preparation took 1.388 + 9.453 seconds, versus 0.270 + 21.644 seconds
with all sums/products unevaluated. Warmed 31-q kernels took 2.413 versus
20.663 ms. These kernel timings exclude the helper's parameter/q validation
wrapper. This confirms the regression is relevant to the helper's preparation
settings as well as the original stage-driver setting. The restricted variant
also regressed with parameter renaming: preparation took 20.46 seconds and
curves 21.99 ms.

The stage driver now accepts optional sixth argument `dummify` to reproduce
that preparation policy, and records its value in JSON. Existing default,
`profile`, and `simplify` modes are unchanged:

```bash
OPENBLAS_NUM_THREADS=1 SYMPY_CACHE_SIZE=1000 /tmp/pyseb-symbolic-stages \
  chain 64 /home/tobi/pySEB/.venv-sympy/lib/python3.11/site-packages 3 dummify
```

After restoring and rebuilding the production converter, all ten symbolic
regressions and all eight compiled-evaluator regressions pass. The retained
changes are benchmark policy coverage and these findings, not a production
conversion optimization.

## Prepared sasmodels fitting: combined intensity and amplitude

The runtime adapter can now request the complete unnormalized intensity and
signed reference amplitude together. C++ reuses its traversal paths and leaf
amplitudes within each q sample, and returns numeric buffers instead of Python
float lists. Intensity is still evaluated independently; this does not replace
it with the square of the mean amplitude. Complete parameter bindings still
run on every distribution sample, and no state is cached across fitting calls.
The adapter uses this combined path when an interaction-enabled model needs F
and I. The existing intensity-only adapter path is retained: its initial
buffer-only experiment did not establish a useful improvement.

`profile_sasmodels_fitting.py` measures already-prepared runtime models while
changing radius on every call. Both hybrid and interacting AOT cases include
Pinhole1D resolution on a 101-point measured grid (108 evaluation points).
Declared radius distributions use 21 samples for the hybrid and 11 for AOT.
Model construction, symbolic preparation, and kernel creation are excluded.
`--fit` additionally fits AOT radius and shell thickness to noiseless synthetic
data, starting at (25, 8) Angstrom and recovering (29, 6.5) Angstrom. This is a
performance/correctness fixture, not an experimental-data fit.

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
SAS_DLL_PATH=/tmp/pyseb-fitting-kernels \
python benchmarks/profile_sasmodels_fitting.py \
  "$PWD" /tmp/pyseb-fitting.json --fit
```

The first positional argument selects the source/build tree and requires its
`build-python/pyseb` extension. This permits the same driver to measure a
separately built older snapshot. The script verifies the selected native module
path, saves changing-parameter curves and timing samples, and checks parameter
recovery. `--profile` separately profiles repeated evaluations after timing.
It needs the optional sasmodels extra and SciPy, and adds no dependencies to
ordinary imports.

Local portable Release comparison against `3971dd52`, GCC 12.2 with the same
-O3/LTO settings, Python 3.11, NumPy 2.4.2, SciPy 1.17.1 and sasmodels 1.1.0:

| Prepared runtime workload | Before | Combined evaluation | Speedup |
| --- | ---: | ---: | ---: |
| Interacting AOT curve, radius distribution of 11 | 193.3 ms | 154.4 ms | 1.25x |
| Complete two-parameter SciPy fit, same distribution and resolution | 3.607 s | 2.788 s | 1.29x |

These are medians of three fresh-process medians, with three timed samples per
process and alternating old/new process order. No builds, tests or other
profilers ran concurrently. All 18 fits recovered identical parameters and
reported six function evaluations. Changing-parameter curves agreed at relative
tolerance 2e-12, with zero absolute tolerance. A separate five-round comparison
of old/new adapter callbacks sharing the exact same native binary measured
1.35x for the distributed interacting AOT curve. Local run times vary; these
measurements do not establish a universal speedup for every particle.

This targets repeated **runtime** fitting evaluations. Compiled analytical
kernels and their SymPy preparation are unchanged. No speedup is claimed for
intensity-only hybrids. Numeric-buffer regressions cover strided/converted
inputs, owned output storage, empty and invalid grids, signed contrast,
contrast matching, hydration changes and cloud/analytical interference.
Generated curves, profiles, logs and raw measurements remain under `/tmp`.
