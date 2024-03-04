<p align="center">
  <picture>
      <source media="(prefers-color-scheme: dark)" srcset="assets/alecs_logo_dark.svg?sanitize=true" height=130>
      <source media="(prefers-color-scheme: light)" srcset="assets/alecs_logo_light.svg?sanitize=true" height=130>
      <img alt="Database logo" src="">
  </picture>
</p>

This repository presents the **A**strochemistry **L**ow-energy **e**lectron **C**ross-**S**ection ([ALeCS](astrobrandt.github.io/ALeCS/)) database, the result of an international collaboration between astronomers and chemists. The database includes electron-impact cross section data for a wide range of molecules of astrochemical interest, ranging from simple diatomics such as carbon monoxide (CO) to complex organics (CH3OH, methanol) to prebiotics (NH2CHO, formamide).

[![doi](https://img.shields.io/badge/doi-10.48550/arXiv.2310.10739-blue?logo=DOI&logoColor=white)](https://ui.adsabs.harvard.edu/abs/2023arXiv231010739G/abstract)

![version](https://img.shields.io/badge/version-2024.02-orange)

[astrobrandt.github.io/ALeCS/]: astrobrandt.github.io/ALeCS/

## Initial Release
The initial release of the database, described in [Gaches et al. (2024)](https://ui.adsabs.harvard.edu/abs/2023arXiv231010739G/abstract), includes total electron-impact ionization cross sections for over 200 neutral molecules and a few ions spanning electron energies from 10 eV to 5 keV. The cross sections were calculated following the [Binary-encounter Bethe](https://ui.adsabs.harvard.edu/abs/1994PhRvA..50.3954K/abstract) (BEB) formalism using optimize molecule geometries and electronic structures computed using the [Gaussian16](https://www.gaussian.com/) and [MolPro](https://www.molpro.net/) codes. The initial release also includes computed ionization potentials for most of the molecules, computed at two different levels of theory: CCSD(T)/aug-cc-pVTZ+CAM-B3LYP/aug-cc-pVTZ and CCSD(T)/CBS.

> [!NOTE]
> Evaluations of ionization potentials can vary drastically between each other, and for many molecules experimental values vary as well. We suggest those interested in the ionization potentials to also consult the [NIST Chemistry WebBook](https://webbook.nist.gov/chemistry/ie-ser/).

In the database we also include the optimized geometries and computed electron orbitals. The orbitals can be used to recompute the BEB cross sections in applications as desired, with the format in `NIST_orbitals/` matching that in the [NIST database](https://physics.nist.gov/PhysRefData/Ionization/molTable.html).

> [!NOTE]
> Using the orbitals directly is preferable for those wishing to use the cross sections outside the energy range in the raw data, to avoid extrapolation issues. We have added a notebook in the `scripts/` folder to demonstrate this.

Finally, we include chemical networks for the ionization chemistry in both the [UMIST](http://udfa.ajmarkwick.net/index.php) and [KIDA](https://kida.astrochem-tools.org/) formats.

## Database Format
The database folders are described as follows:
- `NIST_orbitals/`: This folder stores the molecule orbitals in the same format as the NIST database. It has all the information needed to compute the BEB cross sections.
  * This file is formatted with columns as follows: `Orbital B(eV) U(eV) N Q`, where `Orbital` is the orbital number, `B(eV)` is the electron binding energy in eV, `U(eV)` is the average kinetic energy of the electron, `N` is the number of electrons in the orbital and `Q` is the dipole constant, but in this data is kept to be 1. This is done to be consistent with the format of the NIST database.
- `full_orbitals/`: Here we store the full information for the molecular orbital calculation. This includes the different orbital energy level calculations from the optimization and subsequent population analysis and from electron propagator theory (EPT). This data is presented for transparency, and we recommend using the NIST_orbital format for calculations.
- `geoms/`: This folder stores [pdb format](https://en.wikipedia.org/wiki/Protein_Data_Bank_(file_format)) files with the optimized geometries at the HF, MP2 and CCSD(T) levels of theory.
  * The MP2 and CCSD(T) geometries are saved as PDB files. The HF geometries are saved as XYZ files. These were checked and readable by GaussView and Avogadro (versions 1 and 2).
- `ion_xs/`: This folder contains the total electron-impact ionization cross sections at HF, MP2 and CCSD(T) levels of theory. The folders contain different subsets of the molecules, with some overlap between them.
  * These files have two or three columns. The first column is always the energy in eV. The second column is the BEB cross section in units of $a_0^2$ where $a_0$ is the Bohr radius. The third column, if it is there, is the damped BEB cross section in units of $a_0^2$.
  * In the `recommended/` folder there is the full sample of cross sections where for each species we give the recommended cross sections. These files are all two column containing **only** the BEB cross sections.
- `ips/`: There are several text files here containing the calculated ionization potentials at CAM-B3LYP/aug-cc-pVQZ and CCSD(T)/CBS levels of theory and the recommended values from the NIST Chemistry WebBook.
  * There is a file for each subset, with the first column being the chemical name and the second the ionization potential in eV.


## Collaboration
Our current collaboration is
- Brandt Gaches (PI), Cosmic Origins Fellow, Chalmers University of Technology, Sweden
- Stefano Bovino, Associate Professor, Sapienza University of Rome, Italy
- Giulia Bovolenta, PhD Student, Universidad de Concepción, Chile
- Prasanta Gorai, Cosmic Origins Fellow, Chalmers University of Technology & Gothenburg University, Sweden
- Tommaso Grassi, Scientist, Max Planck Institute for Extraterrestrial Physics, Germany
- David Heathcote, Postdoctoral Researcher, University of Oxford, United Kingdom
- Marco Padovani, Scientist, INAF-Observatorio Astrofisico di Arcetri, Italy
- Claire Vallance, Professor, University of Oxford, United Kingdom
- Stefan Vogt-Geisse, Professor, Universidad de Concepción, Chile

The ALeCS database is an ongoing live project and open to community involvement. If you wish to get involved in this work, please feel free to send the PI an email. We have a broad interst in expanding the database to include a wide range of physical and energetic chemical processes relevant for astrochemical modeling.

## Citing the database
We strongly believe that astrophysical and astrochemical modeling should cite the intrinsic data that goes into the chemical networks and model. If you use the data within this database, we request that you cite the database as
```
@ARTICLE{Gaches2024,
       author = {{Gaches}, Brandt A.~L. and {Grassi}, Tommaso and {Vogt-Geisse}, Stefan and {Bovolenta}, Giulia M. and {Vallance}, Claire and {Heathcote}, David and {Padovani}, Marco and {Bovino}, Stefano and {Gorai}, Prasanta},
        title = "{The Astrochemistry Low-energy Electron Cross-Section (ALeCS) database I. Semi-empirical electron-impact ionization cross-section calculations and ionization rates}",
      journal = {arXiv e-prints},
     keywords = {Astrophysics - Astrophysics of Galaxies, Astrophysics - Instrumentation and Methods for Astrophysics, Physics - Chemical Physics},
         year = 2023,
        month = oct,
          eid = {arXiv:2310.10739},
        pages = {arXiv:2310.10739},
          doi = {10.48550/arXiv.2310.10739},
archivePrefix = {arXiv},
       eprint = {2310.10739},
 primaryClass = {astro-ph.GA},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2023arXiv231010739G},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
```
If you use the Hartree-Fock computed data from [Heathcote & Valance (2018)](https://dx.doi.org/10.1088/1361-6455/aadd42) and [Zhou et al. 2019](https://ui.adsabs.harvard.edu/abs/2019MolPh.117.3066Z), we request you cite the papers:
```
@article{Heathcote2018,
    doi = {10.1088/1361-6455/aadd42},
    url = {https://dx.doi.org/10.1088/1361-6455/aadd42},
    year = {2018},
    month = {sep},
    publisher = {IOP Publishing},
    volume = {51},
    number = {19},
    pages = {195203},
    author = {David Heathcote and Claire Vallance},
    title = {Total electron ionization cross-sections for neutral molecules relevant to astrochemistry},
    journal = {Journal of Physics B: Atomic, Molecular and Optical Physics}
}
```
```
@ARTICLE{Zhou2019,
       author = {{Zhou}, Weiwei and {Wilkinson}, Lorna and {Lee}, Jason W.~L. and {Heathcote}, David and {Vallance}, Claire},
        title = "{Total electron ionization cross-sections for molecules of astrochemical interest}",
      journal = {Molecular Physics},
         year = 2019,
        month = nov,
       volume = {117},
       number = {21},
        pages = {3066-3075},
          doi = {10.1080/00268976.2019.1583389},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2019MolPh.117.3066Z},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
```
Finally, if you use any data from subsequent releases, we request that you cite the paper relevant for that release.

## Support
We acknowledge support from Chalmers University of Technology and the Chalmers Initiative on Cosmic Origins. The MP2 and CAM-B3LYP calculations were performed on the Vera computing facility managed by the Chalmers Centre for Computational Science and Engineering.
