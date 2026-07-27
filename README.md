# KerrSpinningFluxes

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19609875.svg)](https://doi.org/10.5281/zenodo.19609875)

This is a Wolfram Mathematica package for gravitational-wave fluxes from spinning particles orbiting a Kerr black hole. These codes are based on and used in the following papers:
- Viktor Skoupý, _A new approach to the calculation of extreme-mass-ratio inspirals with a spinning secondary_, [arXiv:2603.13482](https://arxiv.org/abs/2603.13482)
- Lisa V. Drummond, Scott A. Hughes, Viktor Skoupý, Philip Lynch, Gabriel Andres Piovano, _Shifted-geodesic approximation for spinning-body gravitational wave fluxes_  [arXiv:2603.12189](https://arxiv.org/abs/2603.12189)
- Viktor Skoupý, Gabriel A. Piovano, Vojtěch Witzany, _Spherical inspirals of spinning bodies into Kerr black holes_, [arXiv:2506.20726](https://arxiv.org/abs/2506.20726)
- Viktor Skoupý, Georgios Lukes-Gerakopoulos, Lisa V. Drummond, Scott A. Hughes, _Asymptotic gravitational-wave fluxes from a spinning test body on generic orbits around a Kerr black hole_, [arXiv:2303.16798](https://arxiv.org/abs/2303.16798)

Examples with the calculation of the analytical and numerical trajectory and the fluxes can be found in `Example.nb`.

## Requirements

### For nonlinearized fluxes

- [KerrGeodesics](https://bhptoolkit.org/KerrGeodesics/)
- [Teukolsky](https://bhptoolkit.org/Teukolsky/)
- [SpinWeightedSpheroidalHarmonics](https://bhptoolkit.org/SpinWeightedSpheroidalHarmonics/)

### For linearized fluxes

- modified [KerrGeodesics](https://github.com/vskoupy/KerrGeodesics)
- modified SpinWeightedSpheroidalHarmonics by Gabriel A. Piovano
- LinearizedTeukolskyEquations by Gabriel A. Piovano
- Functions_for_spherical_orbits_fixed_turning_point_on_average.nb by Gabriel A. Piovano

## Citing

When using this package, please add the following citations:

- Skoupý, V. (2026). vskoupy/KerrSpinningFluxes: KerrSpinningFluxes 1.0.1 (v1.0.1). Zenodo. https://doi.org/10.5281/zenodo.19609875
- Skoupý, V., “A new approach to the calculation of extreme-mass-ratio inspirals with a spinning secondary”, _ArXiv e-prints_, Art. no. arXiv:2603.13482, 2026. [doi:10.48550/arXiv.2603.13482](https://doi.org/10.48550/arXiv.2603.13482).

```
@software{skoupy_2026_19609875,
  author       = {Skoupý, Viktor},
  title        = {vskoupy/KerrSpinningFluxes: KerrSpinningFluxes
                   1.0.1
                  },
  month        = apr,
  year         = 2026,
  publisher    = {Zenodo},
  version      = {v1.0.1},
  doi          = {10.5281/zenodo.19609875},
  url          = {https://doi.org/10.5281/zenodo.19609875},
  swhid        = {swh:1:dir:fb9afa4a59b77b086e89f190f13f113304bf623a
                   ;origin=https://doi.org/10.5281/zenodo.18377281;vi
                   sit=swh:1:snp:76e60e5f3f16b47baf16a9a2c090d7d22580
                   bf51;anchor=swh:1:rel:b27348f8eeaed35e0201ea2e1615
                   747446fbffb6;path=vskoupy-
                   KerrSpinningFluxes-3764133
                  },
}

@article{Skoupy:2026ewu,
    author = "Skoup{\'y}, Viktor",
    title = "{A new approach to the calculation of extreme-mass-ratio inspirals with a spinning secondary}",
    eprint = "2603.13482",
    archivePrefix = "arXiv",
    primaryClass = "gr-qc",
    month = "3",
    year = "2026"
}
```

## Contributors

- Viktor Skoupý
