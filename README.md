# KerrSpinningFluxes

This is a Wolfram Mathematica package for gravitational-wave fluxes from spinning particles orbiting a Kerr black hole. These codes are based on and used in the following papers:
- Viktor Skoupý, _A new approach to the calculation of extreme-mass-ratio inspirals with a spinning secondary_, [arXiv:2603.13482](https://arxiv.org/abs/2603.13482)
- Lisa V. Drummond, Scott A. Hughes, Viktor Skoupý, Philip Lynch, Gabriel Andres Piovano, _Shifted-geodesic approximation for spinning-body gravitational wave fluxes_  [arXiv:2603.12189](https://arxiv.org/abs/2603.12189)
- Viktor Skoupý, Gabriel A. Piovano, Vojtěch Witzany, _Spherical inspirals of spinning bodies into Kerr black holes_, [arXiv:2506.20726](https://arxiv.org/abs/2506.20726)
- Viktor Skoupý, Georgios Lukes-Gerakopoulos, Lisa V. Drummond, Scott A. Hughes, _Asymptotic gravitational-wave fluxes from a spinning test body on generic orbits around a Kerr black hole_, [arXiv:2303.16798](https://arxiv.org/abs/2303.16798)

Examples with the calculation of the analytical and numerical trajectory and the fluxes can be found in 'Example.nb'.

## Requirements

### For nonlinearized fluxes

- [KerrGeodesics](https://bhptoolkit.org/KerrGeodesics/)
- [Teukolsky](https://bhptoolkit.org/Teukolsky/)
- [SpinWeightedSpheroidalHarmonics](https://bhptoolkit.org/SpinWeightedSpheroidalHarmonics/)

### For linearized fluxes

- modified KerrGeodesics
- modified SpinWeightedSpheroidalHarmonics by Gabriel A. Piovano
- LinearizedTeukolskyEquations by Gabriel A. Piovano
- Functions_for_spherical_orbits_fixed_turning_point_on_average.nb by Gabriel A. Piovano

## Contributors

- Viktor Skoupý
