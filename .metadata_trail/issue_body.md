### -> submitter ORCID (or name)

https://orcid.org/0000-0001-9424-2315

### -> slug

Lu-2026-ALEIB-FreeSurface

### -> license

CC-BY-4.0

### -> alternative license URL

_No response_

### -> model category

model published in study, community benchmark

### -> model status

completed

### -> associated publication DOI

https://doi.org/10.5194/gmd-19-5191-2026

### -> model creators

https://orcid.org/0000-0001-9424-2315

### -> title

A novel ALE scheme with the internal boundary for true free surface simulation in geodynamic models

### -> description

We propose a novel scheme for modelling the true free surface within the finite element method (FEM), which integrates the “sticky air” approach into the ALE scheme. This method employs an internal boundary to accurately represent the free surface, referred to as ALE-IB. We implement this scheme for free surface simulations in the geodynamic codes Underworld 2 and Underworld 3.

### -> abstract

The accurate simulation of Earth's surface is essential for understanding lithospheric and mantle dynamics, especially in processes such as subduction and surface deformation. Traditional top boundary conditions, such as free-slip or no-slip, do not fully capture the complex interactions occurring at the surface. The commonly used “Sticky Air” method, while practical, suffers from several limitations, including increased computational cost and marker fluctuation issues. Additionally, free surface numerical fluctuations, known as the “drunken sailor instability”, are characteristic of all free surface simulations, including true Lagrangian free surface treatments and Arbitrary Lagrangian–Eulerian (ALE) methods. In this study, we propose a novel scheme within the finite element framework that integrates the “Sticky Air” concept into an ALE formulation by employing an internal boundary to simulate a true free surface, referred to as the ALE-IB. This approach effectively addresses the limitations of existing methods, notably by reducing marker fluctuation issues and enhancing numerical stability. Moreover, it maintains a true surface in the computational domain that can be further reshaped by surface processes such as erosion and deposition, and provides a foundational scheme for further coupling framework of tectonic modelling and landscape evolution modelling. We detail the theoretical formulation, implementation strategies, and validation through a series of numerical experiments. The results demonstrate that our method achieves higher accuracy and broader applicability compared to conventional techniques. Ultimately, this framework provides a more realistic and robust tool for geodynamic modelling of the Earth's free surface.

### -> scientific keywords

free surface, ALE, internal boundary

### -> funder

https://ror.org/05mmh0f86, DP240102450

### -> model embargo?

_No response_

### -> include model code ?

- [x] yes

### -> model code/inputs DOI

https://github.com/NengLu/ALEIB_FreeSurface

### -> model code/inputs notes

Models are set up with Python scripts. Models require underworld2 or nderworld3 to be run.

### -> include model output data?

- [x] yes

### -> data creators

https://orcid.org/0000-0001-9424-2315

### -> model output data DOI

_No response_

### -> model output data notes

output is primarily h5 files, ~13.1 GB of data.

### -> model output data size

~13.1 GB

### -> software framework DOI/URI

https://zenodo.org/records/15128361

### -> software framework source repository

https://github.com/underworldcode/underworld2

### -> name of primary software framework (e.g. Underworld, ASPECT, Badlands, OpenFOAM)

Underworld2, Underworld3

### -> software framework authors

_No response_

### -> software & algorithm keywords

python, finite element, particle in cell

### -> computer URI/DOI

https://nci.org.au/

### -> add landing page image and caption

<img width="2385" height="1907" alt="Image" src="https://github.com/user-attachments/assets/04273c52-defb-4497-8a3c-75eb2ce17821" />
Classification of methods used for simulating a free surface (indicated by the magenta line). Colored points represent markers for different materials. Methods include: (a) ALE scheme, (b) ALE scheme with the internal boundary (ALE-IB) and the “sticky air” method, and (c) Eulerian scheme with the “sticky air” method. The real interface refers to the actual free surface, while the virtual interface represents the surface obtained from numerical modelling.

### -> add an animation (if relevant)

_No response_

### -> add a graphic abstract figure (if relevant)

_No response_

### -> add a model setup figure (if relevant)

<img width="2384" height="1847" alt="Image" src="https://github.com/user-attachments/assets/5b75bf45-c51a-49c7-bf5f-22a77013bf8e" />
Model setup for (a) viscous relaxation of sinusoidal topography, (b) Rayleigh–Taylor instability, (c) delamination, (d) rising sphere, and (e) subduction. Ω indicates different material domains, with Ω0 specifically representing the air domain – used only in Eulerian and ALE-IB schemes. The dashed line marks the free surface, and the stars denote tracer locations used in some experiments.

### -> add a description of your model setup

1. Topography relaxation
The loading of the Earth's surface can be described as the initial periodic surface displacement of an isoviscous fluid within the infinite half-space ([Turcotte and Schubert](https://gmd.copernicus.org/articles/19/5191/2026/#bib1.bibx44), [2002](https://gmd.copernicus.org/articles/19/5191/2026/#bib1.bibx44)).

2. Rayleigh–Taylor instability
The Rayleigh–Taylor instability model is adapted from [Kaus et al.](https://gmd.copernicus.org/articles/19/5191/2026/#bib1.bibx25) ([2010](https://gmd.copernicus.org/articles/19/5191/2026/#bib1.bibx25)) and [Duretz et al.](https://gmd.copernicus.org/articles/19/5191/2026/#bib1.bibx16) ([2011](https://gmd.copernicus.org/articles/19/5191/2026/#bib1.bibx16)). A dense and more viscous layer is sinking through a less dense fluid.

3. Delamination
This experiment builds upon the models developed in [Beall et al.](https://gmd.copernicus.org/articles/19/5191/2026/#bib1.bibx4) ([2017](https://gmd.copernicus.org/articles/19/5191/2026/#bib1.bibx4)) to examine conditions leading to triggered dripping and lithospheric delamination.

4. Rising Sphere
The rising sphere model is adapted from Case 2 in [Crameri et al.](https://gmd.copernicus.org/articles/19/5191/2026/#bib1.bibx11) ([2012](https://gmd.copernicus.org/articles/19/5191/2026/#bib1.bibx11)) for validating the sticky air approach.

5. Subduction
The subduction model is modified from [Crameri et al.](https://gmd.copernicus.org/articles/19/5191/2026/#bib1.bibx12) ([2017](https://gmd.copernicus.org/articles/19/5191/2026/#bib1.bibx12)). It is a thermo-mechanical model designed to simulate the subduction of a visco-plastic slab into the mantle and generate realistic topography signals.

### Please provide any feedback on the model submission process?

_No response_