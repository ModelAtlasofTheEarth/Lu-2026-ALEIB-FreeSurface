# New [M@TE](https://mate.science/)! model: 
 _we have provided a summary of your model as a starting point for the README, feel free to edit_
## Section 1: Summary of your model   

**Model Submitter:**  

Neng Lu ([0000-0001-9424-2315](https://orcid.org/0000-0001-9424-2315))

**Model Creator(s):**  

- Neng Lu ([0000-0001-9424-2315](https://orcid.org/0000-0001-9424-2315))  
  
**Model slug:**  

`Lu-2026-ALEIB-FreeSurface` 

(this will be the name of the model repository when created) 

**Model name:**  

_A novel ALE scheme with the internal boundary for true free surface simulation in geodynamic models_  

**License:**  

[Creative Commons Attribution 4.0 International]( https://creativecommons.org/licenses/by/4.0/legalcode.txt)

**Model Category:**  

- model published in study   
- community benchmark   
  
**Model Status:**  

- completed   
  
**Associated Publication title:**  

_[A novel ALE scheme with the internal boundary for true free surface simulation in geodynamic models](https://doi.org/10.5194/gmd-19-5191-2026)_ 

**Short description:**  

We propose a novel scheme for modelling the true free surface within the finite element method (FEM), which integrates the “sticky air” approach into the ALE scheme. This method employs an internal boundary to accurately represent the free surface, referred to as ALE-IB. We implement this scheme for free surface simulations in the geodynamic codes Underworld 2 and Underworld 3.

**Abstract:**  

The accurate simulation of Earth's surface is essential for understanding lithospheric and mantle dynamics, especially in processes such as subduction and surface deformation. Traditional top boundary conditions, such as free-slip or no-slip, do not fully capture the complex interactions occurring at the surface. The commonly used “Sticky Air” method, while practical, suffers from several limitations, including increased computational cost and marker fluctuation issues. Additionally, free surface numerical fluctuations, known as the “drunken sailor instability”, are characteristic of all free surface simulations, including true Lagrangian free surface treatments and Arbitrary Lagrangian–Eulerian (ALE) methods. In this study, we propose a novel scheme within the finite element framework that integrates the “Sticky Air” concept into an ALE formulation by employing an internal boundary to simulate a true free surface, referred to as the ALE-IB. This approach effectively addresses the limitations of existing methods, notably by reducing marker fluctuation issues and enhancing numerical stability. Moreover, it maintains a true surface in the computational domain that can be further reshaped by surface processes such as erosion and deposition, and provides a foundational scheme for further coupling framework of tectonic modelling and landscape evolution modelling. We detail the theoretical formulation, implementation strategies, and validation through a series of numerical experiments. The results demonstrate that our method achieves higher accuracy and broader applicability compared to conventional techniques. Ultimately, this framework provides a more realistic and robust tool for geodynamic modelling of the Earth's free surface.

**Scientific Keywords:**  

- free surface   
- ALE   
- internal boundary   
  
**Funder(s):**  
- ARC (https://ror.org/05mmh0f86)  
  
## Section 2: your model code, output data  

**No embargo on model contents requested** 

**Include model code:**   

True 

**Model code existing URL/DOI:**   

https://github.com/NengLu/ALEIB_FreeSurface 

**Model code notes:**   

Models are set up with Python scripts. Models require underworld2 or nderworld3 to be run. 

**Include model output data:**   

True 

**Model output data notes:**   

output is primarily h5 files, ~13.1 GB of data. 

## Section 3: software framework and compute details   
**Software Framework DOI/URL:**  

Found software: _[Underworld2, Underworld3](https://doi.org/10.5281/zenodo.15128361)_ 

**Software Repository:**   

https://github.com/underworldcode/underworld2 

**Name of primary software framework:**  

Underworld2, Underworld3 

**Software framework authors:**  
-  Romain Beucher ([0000-0003-3891-5444](0000-0003-3891-5444))  
-  Julian Giordani ([0000-0003-4515-9296](0000-0003-4515-9296))  
-  Louis Moresi ([0000-0003-3685-174X](0000-0003-3685-174X))  
-  John Mansour ([0000-0001-5865-1664](0000-0001-5865-1664))  
-  Owen Kaluza ([0000-0001-6303-5671](0000-0001-6303-5671))  
-  Mirko Velic   
-  Rebecca Farrington ([0000-0002-2594-6965](0000-0002-2594-6965))  
-  Steve Quenette ([0000-0002-0368-7341](0000-0002-0368-7341))  
-  Adam Beall ([0000-0002-7182-1864](0000-0002-7182-1864))  
-  Dan Sandiford ([0000-0002-2207-6837](0000-0002-2207-6837))  
-  Luke Mondy ([0000-0001-7779-509X](0000-0001-7779-509X))  
-  Claire Mallard ([0000-0003-2595-2414](0000-0003-2595-2414))  
-  Patrice Rey ([0000-0002-1767-8593](0000-0002-1767-8593))  
-  Guillaume Duclaux ([0000-0002-9512-7252](0000-0002-9512-7252))  
-  Arijit Laik ([0000-0002-3484-7985](0000-0002-3484-7985))  
-  Sara Morón ([0000-0002-1270-4377](0000-0002-1270-4377))  
-  Adam Beall ([0000-0002-7182-1864](0000-0002-7182-1864))  
-  Ben Knight ([0000-0001-7919-2575](0000-0001-7919-2575))  
-  Neng Lu ([0000-0001-9424-2315](0000-0001-9424-2315))  
  
**Software & algorithm keywords:**  

- python   
- finite element   
- particle in cell   
  
## Section 4: web material (for mate.science)   
**Landing page image:**  

Filename: [graphics/Image.png](https://github.com/user-attachments/assets/04273c52-defb-4497-8a3c-75eb2ce17821)  
Caption: Classification of methods used for simulating a free surface (indicated by the magenta line). Colored points represent markers for different materials. Methods include: (a) ALE scheme, (b) ALE scheme with the internal boundary (ALE-IB) and the “sticky air” method, and (c) Eulerian scheme with the “sticky air” method. The real interface refers to the actual free surface, while the virtual interface represents the surface obtained from numerical modelling.  
  
**Animation:**  

Filename: [None]()  
  
**Graphic abstract:**  

Filename: [None]()  
  
**Model setup figure:**  

Filename: [graphics/Image.png](https://github.com/user-attachments/assets/5b75bf45-c51a-49c7-bf5f-22a77013bf8e)  
Caption: Model setup for (a) viscous relaxation of sinusoidal topography, (b) Rayleigh–Taylor instability, (c) delamination, (d) rising sphere, and (e) subduction. Ω indicates different material domains, with Ω0 specifically representing the air domain – used only in Eulerian and ALE-IB schemes. The dashed line marks the free surface, and the stars denote tracer locations used in some experiments.  
Description:  1. Topography relaxation
The loading of the Earth's surface can be described as the initial periodic surface displacement of an isoviscous fluid within the infinite half-space ([Turcotte and Schubert](https://gmd.copernicus.org/articles/19/5191/2026/#bib1.bibx44), [2002](https://gmd.copernicus.org/articles/19/5191/2026/#bib1.bibx44)).

2. Rayleigh–Taylor instability
The Rayleigh–Taylor instability model is adapted from [Kaus et al.](https://gmd.copernicus.org/articles/19/5191/2026/#bib1.bibx25) ([2010](https://gmd.copernicus.org/articles/19/5191/2026/#bib1.bibx25)) and [Duretz et al.](https://gmd.copernicus.org/articles/19/5191/2026/#bib1.bibx16) ([2011](https://gmd.copernicus.org/articles/19/5191/2026/#bib1.bibx16)). A dense and more viscous layer is sinking through a less dense fluid.

3. Delamination
This experiment builds upon the models developed in [Beall et al.](https://gmd.copernicus.org/articles/19/5191/2026/#bib1.bibx4) ([2017](https://gmd.copernicus.org/articles/19/5191/2026/#bib1.bibx4)) to examine conditions leading to triggered dripping and lithospheric delamination.

4. Rising Sphere
The rising sphere model is adapted from Case 2 in [Crameri et al.](https://gmd.copernicus.org/articles/19/5191/2026/#bib1.bibx11) ([2012](https://gmd.copernicus.org/articles/19/5191/2026/#bib1.bibx11)) for validating the sticky air approach.

5. Subduction
The subduction model is modified from [Crameri et al.](https://gmd.copernicus.org/articles/19/5191/2026/#bib1.bibx12) ([2017](https://gmd.copernicus.org/articles/19/5191/2026/#bib1.bibx12)). It is a thermo-mechanical model designed to simulate the subduction of a visco-plastic slab into the mantle and generate realistic topography signals.

  
