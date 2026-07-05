# Cubature from Rational Approximation

This package accompanies the manuscript `Cubature from rational approximation`.
It contains the MATLAB scripts, paper source, compiled paper PDF, and generated
figure PDFs needed to reproduce the numerical experiments.

## Contents

- `matlab/` - all MATLAB source files for the experiments.
- `paper/` - LaTeX source, SIAM class file, compiled PDF, and figure PDFs.
- `docs/` - project notes and prototype result summary.
- `logs/` - the most recent paper build log included with this package.

No external data files are required.

## Dependencies

- MATLAB with `integral2`, `exportgraphics`, and standard plotting support.
- Chebfun, for `jacpts` and `legpts`.
- An `aaa.m` implementation. Chebfun's `aaa` is sufficient; if an `aaa`
  supporting the optional `'sign'` flag is available, `cub_rule.m` uses it.
- A LaTeX installation with `pdflatex`, to rebuild the manuscript.

If dependencies are not already on the MATLAB path, either set
`CHEBFUN_PATH` / `AAA_PATH` before starting MATLAB, or copy
`matlab/local_paths_template.m` to `matlab/local_paths.m` and edit the paths.
`cub_config.m` also checks a few common relative and home-directory locations
for Chebfun.

## Quick Start

In MATLAB:

```matlab
cd matlab
run_reproduce
```

This runs every experiment used in the paper through `run_freeze.m`, regenerates
the figure PDFs in `matlab/`, and copies the paper figures into
`paper/figures/`.

For a faster smoke test:

```matlab
cd matlab
cub1_disk
cub2_ellipse
```
