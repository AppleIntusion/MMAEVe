# MMAEVe - Membranes, Micelles, And Even Vesicles

[![DOI](https://zenodo.org/badge/485386956.svg)](https://zenodo.org/doi/10.5281/zenodo.11105899)

MMAEVe is a simple python package for creating complex biomolecular systems. It can be used to distribute biomolecules about different surfaces, remove overlap between constructed system elements, write pdb files suitable for use as initial starting structures for Amber and Gromacs, and write Gromacs topology files. It is a simple but powerful tool that can be used to systematically generate complex structures.

## Installation

### Requirements

MMAEVe will work with python version >=3.5. Additionally, it requires the following packages:

```bash
numpy
scipy
pandas
```

### Installation

I intened to make MMAEVe available pip or conda at some point. For now just clone the repo then add it to your PYTHONPATH.

```bash
git clone https://github.com/AppleIntusion/MMAEVe.git
export PYTHONPATH="${PYTHONPATH}":some_path/MMAEVe
```

Then it can be imported.

```python
import MMAEVe as mav
```

## Usage

A comprehensive [tutorial](tutorial/README.md) is provided. It serves as an introduction to how MMAEVe works and showcase of systems that it can be used to create.

## Example Systems

MMAEVe can be used to easily create complex systems. See the [tutorial](tutorial/README.md) for usage. A couple examples of systems that can be created are shown below.

### Membrane-Bridging
![](tutorial/images/vesicle_bi_a2t.png)

### Covid Virion
![](tutorial/images/covid_virion.png)

## Citation

Please cite using the citation under "About" if you use MMAEVe or incorporate it into a project.
