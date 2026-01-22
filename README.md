# MultiNEAs

## Numerical Tools for Near-Earth Asteroid Dynamics and Population

[![version](https://img.shields.io/pypi/v/multineas?color=blue)](https://pypi.org/project/multineas/)
[![downloads](https://img.shields.io/pypi/dw/multineas)](https://pypi.org/project/multineas/)
[![license](https://img.shields.io/pypi/l/multineas)](https://pypi.org/project/multineas/)
[![pythonver](https://img.shields.io/pypi/pyversions/multineas)](https://pypi.org/project/multineas/)

`MultiNEAs` is a `Python` package designed to provide numerical tools for studying the dynamics and population of Near-Earth Asteroids (NEAs). The package offers a comprehensive suite of utilities for orbital calculations, statistical analysis, and visualization of NEA populations.

## Features

- **Orbital Dynamics**: Tools for computing and analyzing near-earth asteroid orbits
- **Population Studies**: Statistical methods for studying NEA populations
- **Visualization**: Plotting and visualization utilities for asteroid data
- **Data Management**: Efficient handling of asteroid catalogs and databases

## Installation

### From PyPI

`MultiNEAs` will be available on PyPI at https://pypi.org/project/multineas/. Once published, you can install it with:

```bash
pip install -U multineas
```

### From Sources

You can also install from the [GitHub repository](https://github.com/seap-udea/MultiNEAs):

```bash
git clone https://github.com/seap-udea/MultiNEAs
cd MultiNEAs
pip install .
```

For development, use an editable installation:

```bash
cd MultiNEAs
pip install -e .
```

### In Google Colab

If you use Google Colab, you can install `MultiNEAs` by executing:

```python
!pip install -U multineas
```

## Quick Start

Getting started with `MultiNEAs` is straightforward. Import the package:

```python
import multineas as mn
```

> **NOTE**: If you are working in Google Colab, load the matplotlib backend before producing plots:
>
> ```python
> %matplotlib inline
> ```

## Documentation

Full API documentation will be available soon.

## Examples

Working examples and tutorials will be added as the package develops.

## Contributing

We welcome contributions! If you're interested in contributing to MultiNEAs, please:

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Submit a pull request

## Citation

If you use `MultiNEAs` in your research, please cite:

```bibtex
@software{multineas2026,
  author = {Zuluaga, Jorge I. and Agudelo, Juanita A.},
  title = {MultiNEAs: Numerical tools for near-earth asteroid dynamics and population},
  year = {2026},
  url = {https://github.com/seap-udea/MultiNEAs}
}
```

## License

This project is licensed under the GNU Affero General Public License v3.0 (AGPL-3.0) - see the [LICENSE](LICENSE) file for details.

## Authors

- **Jorge I. Zuluaga** - jorge.zuluaga@udea.edu.co
- **Juanita A. Agudelo** - juanita.agudelo@udea.edu.co

## Acknowledgments

This package is being developed at the Solar, Earth and Planetary Physics Group (SEAP) at Universidad de Antioquia, Medellín, Colombia.

## What's New

For a detailed list of changes and new features, see [WHATSNEW.md](WHATSNEW.md).

---

MultiNEAs (C) 2026 - Jorge I. Zuluaga and Juanita A. Agudelo
