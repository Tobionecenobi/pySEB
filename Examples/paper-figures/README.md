# Original paper figure programs

These C++ programs generate scattering curves for figures 3 and 10–13 in the
SEB paper. They are preserved as paper examples; they are not an automatic test
suite. Each source file is a standalone program using the SEB C++ library.

| Filename prefix | Structure |
| --- | --- |
| `fig3_` | Diblock reference-point choices |
| `fig10_` | Linear chains |
| `fig11_` | Dendrimers |
| `fig12_` | Surface attachments |
| `fig13_` | Chain of diblock stars |

The published images are in [`docs/images/paper-figures/`](../../docs/images/paper-figures/).
The current executable Python walkthroughs are the Figure 10–13 tutorials under
[`docs/tutorials/`](../../docs/tutorials/). See the repository's publication guide
for citation details.

Programs write curve files into their working directory. Run them from a
separate output directory to keep generated data out of the source tree.
