# Installation:

Analysator is now a package! Clone the repository locally

```bash
git clone https://github.com/fmihpc/analysator.git
```
and install via `pip` from the cloned path. Dependency handling is via pip. Using the editable flag is recommended for getting updates via `git pull origin`:

```bash
pip install --editable ./analysator
```

For installing the correct VTK depencies, especially BVTKnodes support, use either of
```bash
pip install --editable ./analysator[vtk]
```
```bash
pip install --editable ./analysator[bvtk]
```
respectively.

## Backward compatibility

The packaged analysator is (should be!) backwards-compatible with previous analysator versions: `import pytools` will import analysator onto the pytools name, and emits a deprecation warning. `export PYTHONPATH=[path-to-analysator-repo]` installation as before is also supported. If `--editable` flag is not available, regular `pip install` should work, but in that case using `export PYTHONPATH=[path-to-analysator]` is recommended if you wish to contribute to Analysator (and who wouldn't!).

# Using Analysator:
```python
import analysator as pt  # Import Analysator

# Navigating functions:
pt.calculations.pitch_angles? #press [Enter]
pt.vlsvfile.VlsvReader? #press [Enter]
pt.plot.plot_colormap? #press [Enter]
pt.plot.plot_vdf? #press [Enter]
```

```bash
# For non-interactive mode (also when no X is available):
# set the environment variable PTNONINTERACTIVE to any value before launching
# python/ipython. If in interactive mode, experimental non-blocking
# windows via matplotlib.pyplot.ion() are in use.
#################################################
export PTNONINTERACTIVE=1

# For selecting a backend manually (if Agg or TkAgg is not available)
#################################################
export PTBACKEND=Qt5Agg

# For disabling full LaTeX formatted output (if texlive is not installed)
#################################################
export PTNOLATEX=1
```

# Outputs

For setting the default output directory (default: $HOME/Plots)
```
export PTOUTPUTDIR=/proj/USERNAME/Plots/
```
If using a jupyter notebook on the UH hub system, you can activate interactive plot windows with `%matplotlib ipympl` or `%matplotlib notebook` before importing pytools

Examples and instructions for batch scripting (on CSC's system) are found in
`examples/generate_panel.py` and `examples/generate_movie.sh`

Analysator is using logging. It is controlled by setting the enviroment variable ANALYSATOR_LOG_LEVEL
Supported: DEBUG, INFO (default), WARNING, ERROR, CRITICAL

For example, disable INFO prints via:
```bash
export ANALYSATOR_LOG_LEVEL='WARNING'
```

# About

For more information visit the wiki: https://github.com/fmihpc/analysator/wiki

Analysator reference: https://fmihpc.github.io/analysator/

For citations, use the DOI https://doi.org/10.5281/zenodo.4462514 or the ready-made button the the right in GitHub!
