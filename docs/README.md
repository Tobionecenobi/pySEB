# Building the website

The GitHub Actions [documentation workflow](../.github/workflows/documentation.yml)
contains the complete build environment, including Python dependencies and
compiled C++ tutorial programs. Local documentation environments also need
`jupyter-cache`:

```bash
python -m pip install jupyter-cache
quarto render docs
```

Activate the environment containing pySEB and the tutorial dependencies first.
If Quarto should use a specific environment, set `QUARTO_PYTHON` to its Python
executable. The maintained C++ examples must be built before rendering pages
that run them; the workflow shows the CMake configuration and build commands.

For a prose change to one page, render only that page:

```bash
quarto render docs/tutorials/first-model.qmd
```

## What is cached

Notebook execution results are reused across renders. Quarto still rebuilds
the HTML, prose, navigation, and styles. Python cell changes, scientific input
files, YAML models, native/Python calculation code, C++ examples, helper
modules, installed package versions, Quarto version, and detected GSL/compiler
versions invalidate the execution cache.

This invalidation is deliberately conservative: a calculation change refreshes
all notebooks, keeping cached results and generated downloads in one consistent
calculation context. Prose and navigation edits reuse execution results.

The pre-render hook checks those inputs for both local and CI builds.
The post-render hook records generated downloads. If any recorded download is
missing, notebook results are refreshed so scenes are regenerated rather than
leaving broken links. GitHub Actions persists both the notebook caches and
`tutorials/generated` between compatible builds. Cache directories and
generated scene downloads are not committed.

The first build, or a change to calculation code, inputs, or the environment,
still runs every example. Later prose-only builds reuse those results, including the expensive F127
fits and paired model/unfolded scenes.

## Force a fresh calculation

```bash
quarto render docs --cache-refresh
```

The Documentation workflow also provides a **Recompute every documentation
example** checkbox for manual runs. Use a refresh when intentionally checking
full reproducibility or investigating a suspected stale result. Inputs outside
the checkout and manually changed native build flags require an explicit refresh.
