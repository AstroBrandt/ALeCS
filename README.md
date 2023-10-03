<p align="center">
  <picture>
      <source media="(prefers-color-scheme: dark)" srcset="assets/alecs_logo_dark.svg?sanitize=true" height=130>
      <source media="(prefers-color-scheme: light)" srcset="assets/alecs_logo_light.svg?sanitize=true" height=130>
      <img alt="Database logo" src="">
  </picture>
</p>

This repository presents the **A**strochemistry **L**ow-energy **e**lectron **C**ross **S**ection ([ALeCS](alecs.brandt-gaches.space)) database, the result of an international collaboration between astronomers and chemists. The database includes electron-impact cross section data for a wide range of molecules of astrochemical interest, ranging from simple diatomics such as carbon monoxide (CO) to complex organics (CH3OH, methanol) to prebiotics (NH2CHO, formamide).

![doi](https://img.shields.io/badge/doi-10.0000-blue?logo=DOI&logoColor=white)

![version](https://img.shields.io/badge/version-0.1-orange)

[alecs.brandt-gaches.space]: https://alecs.brandt-gaches.space

## Initial Release
The initial release of the database, described in [Gaches et al. (2023c)](), includes total electron-impact ionization cross sections for over 200 neutral molecules spanning electron energies from 10 eV to 5 keV. The cross sections were calculated following the [Binary-encounter Bethe](https://ui.adsabs.harvard.edu/abs/1994PhRvA..50.3954K/abstract) (BEB) formalism using optimize molecule geometries and electronic structures computed using the [Gaussian16](https://www.gaussian.com/) and [MolPro](https://www.molpro.net/) codes. The initial release also includes computed ionization potentials for most of the molecules, computed at two different levels of theory: CCSD(T)/aug-cc-pVTZ+CAM-B3LYP/aug-cc-pVTZ and CCSD(T)/CBS.

> [!NOTE]
> Evaluations of ionization potentials can vary drastically between each other, and for many molecules experimental values vary as well. We suggest those interested in the ionization potentials to also consult the [NIST Chemistry WebBook](https://webbook.nist.gov/chemistry/ie-ser/).

In the database we also include the optimized geometries and computed electron orbitals. The orbitals can be used to recompute the BEB cross sections in applications as desired, with the format in `NIST_orbitals/` matching that in the [NIST database](https://physics.nist.gov/PhysRefData/Ionization/molTable.html).

> [!NOTE]
> Using the orbitals directly is preferable for those wishing to use the cross sections outside the energy range in the raw data, to avoid extrapolation issues.

Finally, we include chemical networks for the ionization chemistry in both the [UMIST](http://udfa.ajmarkwick.net/index.php) and [KIDA](https://kida.astrochem-tools.org/) formats. For each format, there are two possible networks, one using a Voyager-like spectrum (**L**) and one using a "High" proton spectrum to reproduce diffuse gas measurements (**H**).

## Database Format
The database folders are described as follows:
- `NIST_orbitals/`: This folder stores the molecule orbitals in the same format as the NIST database. It has all the information needed to compute the BEB cross sections.
- `full_orbitals/`: Here we store the full information for the molecular orbital calculation. This includes the different orbital energy level calculations from the optimization and subsequent population analysis and from electron propagator theory (EPT). This data is presented for transparency, and we recommend using the NIST_orbital format for calculations.
- `geoms/`: This folder stores [pdb format](https://en.wikipedia.org/wiki/Protein_Data_Bank_(file_format)) files with the optimized geometries at the HF, MP2 and CCSD(T) levels of theory.
- `ion_xs/`: This folder contains the total electron-impact ionization cross sections at HF, MP2 and CCSD(T) levels of theory. The folders contain different subsets of the molecules, with some overlap between them.
- `ips/`: There are several text files here containing the calculated ionization potentials at CAM-B3LYP/aug-cc-pVQZ and CCSD(T)/CBS levels of theory and the recommended values from the NIST Chemistry WebBook.


## Collaboration
Our current collaboration is
- Brandt Gaches (PI), Cosmic Origins Fellow, Chalmers University of Technology, Sweden
- Stefano Bovino, Associate Professor, Italy
- Prasanta Gorai, Cosmic Origins Fellow, Chalmers University of Technology, Sweden
- Tommaso Grassi, Scientist, Max Planck Institute for Extraterrestrial Physics, Germany
- Marco Mardovani, Scientist, INAF-Observatorio Astrofisico di Arcetri, Italy
- Stefan Geisse, Professor, Universidad de Concepción, Chile
- Giulia Bovolenta, PhD Student, Universidad de Concepción, Chile
- Claire Vallance, Professor, University of Oxford, United Kingdom
- David Heathcote, Postdoctoral Researcher, University of Oxford, United Kingdom

The ALeCS database is an ongoing live project and open to community involvement. If you wish to get involved in this work, please feel free to send the PI an email. We have a broad interst in expanding the database to include a wide range of physical and energetic chemical processes relevant for astrochemical modeling.

## Citing the database
We strongly believe that astrophysical and astrochemical modeling should cite the intrinsic data that goes into the chemical networks and model. If you use the data within this database, we request that you cite the database as
```
@article{Gaches2023c,
author = {Gaches, Brandt and Bovino, Stefano and Gorai, Prasanta and Grassi, Tommaso and Padovani, Marco and Geisse, Stefan and Bovolenta, Giulia and Vallance, Claire and Heathcote, David},
doi = {10.0000/00000},
journal = {A&A},
month = oct,
number = {1},
pages = {1},
title = {{The Astrochemistry Low-energy Electron Cross-Section (ALeCS) database I. Semi-empirical cross-section calculations and ionization rates}},
volume = {1},
year = {2023}
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
    journal = {Journal of Physics B: Atomic, Molecular and Optical Physics}}
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
