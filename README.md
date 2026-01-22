# MultiNEAs

[![version](https://img.shields.io/pypi/v/multineas?color=blue)](https://pypi.org/project/multineas/)
[![license](https://img.shields.io/badge/License-AGPL%20v3-blue.svg)](https://github.com/seap-udea/MultiNEAs/blob/master/LICENSE)
[![pythonver](https://img.shields.io/pypi/pyversions/multineas)](https://pypi.org/project/multineas/)
[![Powered by SpiceyPy](https://img.shields.io/badge/Powered%20by-SpiceyPy-blue)](https://github.com/AndrewAnnex/SpiceyPy)
<!-- [![downloads](https://img.shields.io/pypi/dw/multineas)](https://pypi.org/project/multineas/) -->

<p></p>
<div align="center">
  <a href="https://multineas.readthedocs.io/">
  <img src="https://raw.githubusercontent.com/seap-udea/MultiNEAs/master/docs/MultiNEAs-logo.webp" alt="MultiNEAs Logo" width="600"/>
  </a>
</div>
<p></p>

`MultiNEAs` is a `Python` package designed to provide numerical tools for studying the dynamics and population of Near-Earth Asteroids (NEAs). The package offers a comprehensive suite of utilities for orbital calculations, statistical analysis, and visualization of NEA populations.

## Features

- **Orbital Dynamics**: Tools for computing and analyzing near-earth asteroid orbits.
- **Population Studies**: Statistical methods for studying NEA populations.
- **Visualization**: Plotting and visualization utilities for asteroid data.
- **Data Management**: Efficient handling of asteroid catalogs and databases.

## Documentation

Full API documentation is available at [https://multineas.readthedocs.io](https://multineas.readthedocs.io).

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

## Examples

Working examples and tutorials will be added as the package develops.

## Citation

The numerical tools and codes provided in this package have been developed and tested over several years of scientific research. Several publications have been produced using the code included in this package.

This is a list of the publications where we have tested the code included in this package:

- **Conditions for visual and high-resolution bistatic radar observations of Apophis in 2029.**
  *Agustín Vallejo, Jorge I. Zuluaga, Germán Chaparro*
  Monthly Notices of the Royal Astronomical Society, 518(3), 4438-4448.
  [DOI: 10.1093/mnras/stac3046](https://doi.org/10.1093/mnras/stac3046) | [arXiv: 2201.12205](https://arxiv.org/abs/2201.12205)

- **Location, orbit and energy of a meteoroid impacting the moon during the Lunar Eclipse of January 21, 2019.**
  *Jorge I. Zuluaga, Matipon Tangmatitham, Pablo A. Cuartas-Restrepo, Jonathan Ospina, Fritz Pichardo, Sergio A. Lopez, Karls Pena, J. Mauricio Gaviria-Posada*
  Monthly Notices of the Royal Astronomical Society, Volume 492, Issue 3, March 2020, Pages 3666–3673
  [DOI: 10.1093/mnras/stz3531](https://doi.org/10.1093/mnras/stz3531) | [arXiv: 1901.09573](https://arxiv.org/abs/1901.09573)

- **Can we predict the impact conditions of metre-sized meteoroids?**
  *Jorge I. Zuluaga, Pablo A. Cuartas-Restrepo, Jhonatan Ospina, Mario Sucerquia*
  Monthly Notices of the Royal Astronomical Society: Letters, 486(1), L69-L73.
  [DOI: 10.1093/mnrasl/slz060](https://doi.org/10.1093/mnrasl/slz060) | [arXiv: 1904.12807](https://arxiv.org/abs/1904.12807)

- **Speed Thresholds for Hyperbolic Meteors: The Case of the 2014 January 8 CNEOS Meteor.**
  *Jorge I. Zuluaga*
  Research Notes of the AAS, 3(5), 68.
  [DOI: 10.3847/2515-5172/ab1de3](https://doi.org/10.3847/2515-5172/ab1de3)

- **Towards a theoretical determination of the geographical probability distribution of meteoroid impacts on Earth.**
  *Jorge I. Zuluaga, Mario Sucerquia*
  Monthly Notices of the Royal Astronomical Society, 477(2), 1970-1983.
  [DOI: 10.1093/mnras/sty702](https://doi.org/10.1093/mnras/sty702) | [arXiv: 1801.05720](https://arxiv.org/abs/1801.05720)

- **A General Method for Assessing the Origin of Interstellar Small Bodies: The Case of 1I/2017 U1 ('Oumuamua).**
  *Jorge I. Zuluaga, Oscar Sanchez-Hernandez, Mario Sucerquia, Ignacio Ferrin*
  The Astronomical Journal, 155(6), 236.
  [DOI: 10.3847/1538-3881/aabd7c](https://doi.org/10.3847/1538-3881/aabd7c) | [arXiv: 1711.09397](https://arxiv.org/abs/1711.09397)

- **The orbit of the Chelyabinsk event impactor as reconstructed from amateur and public footage.**
  *Jorge I. Zuluaga, Ignacio Ferrin, Stefan Geens*
  arXiv Repository
  [DOI: 10.48550/arXiv.1303.1796](https://doi.org/10.48550/arXiv.1303.1796) | [arXiv: 1303.1796](https://arxiv.org/abs/1303.1796)

- **A preliminary reconstruction of the orbit of the Chelyabinsk Meteoroid.**
  *Jorge I. Zuluaga, Ignacio Ferrin*
  arXiv Repository
  [DOI: 10.48550/arXiv.1302.5377](https://doi.org/10.48550/arXiv.1302.5377) | [arXiv: 1302.5377](https://arxiv.org/abs/1302.5377) 

If you use `MultiNEAs` in your research, please cite:

```bibtex
@software{multineas2026,
  author = {Zuluaga, Jorge I. and Agudelo, Juanita A.},
  title = {MultiNEAs: Numerical tools for near-earth asteroid dynamics and population},
  year = {2026},
  url = {https://github.com/seap-udea/MultiNEAs}
}
```

## What's New

For a detailed list of changes and new features, see [WHATSNEW.md](WHATSNEW.md).

## Authors and Licensing

This project is developed by the Solar, Earth and Planetary Physics Group (SEAP) at Universidad de Antioquia, Medellín, Colombia. The main developers are:

- **Jorge I. Zuluaga** - jorge.zuluaga@udea.edu.co
- **Juanita A. Agudelo** - juanita.agudelo@udea.edu.co

This project is licensed under the GNU Affero General Public License v3.0 (AGPL-3.0) - see the [LICENSE](LICENSE) file for details.

## Contributing

We welcome contributions! If you're interested in contributing to MultiNEAs, please:

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Submit a pull request

Please read the [CONTRIBUTING.md](CONTRIBUTING.md) file for more information.

## File attribution

Most of this file was vibe coded by the authors using [Gemini 3 Pro](https://gemini.google.com/pro) in [Antigravity](https://antigravity.google/).

---

MultiNEAs (C) 2026 - Jorge I. Zuluaga and Juanita A. Agudelo
