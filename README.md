Analysator is a set of Python tools and scripts used to read and analyze [Vlasiator](https://github.com/fmihpc/vlasiator) output files (VLSV). Vlasiator and its related tools are developed at [University of Helsinki](https://www.helsinki.fi/en/researchgroups/vlasiator), and are open source.

# Installation:

Analysator is now a package!

The packaged analysator is (should be!) backwards-compatible with previous analysator versions: `export PYTHONPATH=[path-to-analysator-repo]` installation as before is still supported. `import pytools` will import analysator onto the pytools name, and emits a deprecation warning.


## 1. Via `git` and `pip`

1. Clone the analysator repository:
```bash
git clone https://github.com/fmihpc/analysator.git
```

2. Install via `pip` from the cloned path, so `pip` will handle dependencies automatically. Using the editable flag is recommended for getting updates via `git pull origin`:

```bash
pip install --editable ./analysator
```
Nb. you are giving pip the root folder of your cloned repository, not some package on PyPi. `--editable` will have your Python environment point at your Analysator repository. You can run this command from within a virtual environment, which is recommended.

Here's a full example of setting up a virtual environment with analysator:
```bash
git clone https://github.com/fmihpc/analysator.git
python -m venv my_analysating_environment
source my_analysating_environment/bin/activate
pip install --editable ./analysator
```

### Updating

With `--editable`, you can update via `git pull` as usual in your cloned analysator repo:
```bash
cd ./analysator
git pull origin
``` 
If `--editable` flag is not available, regular `pip install` should work, but in that case using `export PYTHONPATH=[path-to-analysator]` is recommended if you wish to contribute to Analysator (and who wouldn't!).

### VTK

For installing the correct VTK depencies, especially BVTKnodes support, use either of
```bash
pip install --editable ./analysator[vtk]
```
```bash
pip install --editable ./analysator[bvtk]
```
respectively.

## 2. Via `git` and `PYTHONPATH`

This requires you to handle dependency installation manually (as before, or see the shortcut below).

1. Clone the analysator repository:
```bash
git clone https://github.com/fmihpc/analysator.git
```

2. Add to your `.bashrc` or similar
```bash
export PYTHONPATH=[path-to-analysator]
```

Contained example for Linux installation to user $HOME folder and .bashrc:
```bash
cd
git clone https://github.com/fmihpc/analysator.git
echo "export PYTHONPATH=$PYTHONPATH:$HOME/analysator" >> $HOME/.bashrc
```

### Updating
You can update via `git pull` as usual in your cloned analysator repo:
```bash
cd ./analysator
git pull origin
``` 

### Dependencies shortcut
See `requirements.txt` - `pip install -r ./analysator/requirements.txt` will install all of them.



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

For citations, use the DOI https://doi.org/10.5281/zenodo.4462514 or the ready-made button on the right in GitHub!
