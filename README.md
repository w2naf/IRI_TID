# IRI_TID
Perturb PyIRI with TIDs

Nathaniel Frissell
February 2024

This code can run PyIRI to create a 3D array of ionospheric electron densities. These can then be perturbed with sinusoidal model traveling ionospheric disturbances. 2D slice profiles can then be generated from the perturbed 3D array. These 2D arrays can then be fed into HF raytracing programs like PyLap/PHaRLAP for HF propagation modeling.

The script ```scripts/fhe_tid_2018.py``` was used to generate the ionosphere in Figure 2b of the Frissell et al. (2024) Geophysical Research Letters manuscript on Gravity Waves, Medium Scale Traveling Ionospheric Disturbances, and Large Scale Traveling Ionospheric Disturbances.
