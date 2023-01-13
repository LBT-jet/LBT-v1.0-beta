# LBT
A linear Boltzmann transport model for jet quenching \
Authors: Tan Luo, Shanshan Cao, Yayun He and Xin-Nian Wang

## Introduction
The linear Boltzmann transport (LBT) model is aimed to simulate the interaction between high energtic partons and QCD matter in high energy heavy-ion collision. It employs the Boltzmann transport equation with QCD matrix amptitude to describe the multiple scatterings. Particularly, it places great emphasis on jet-induced medium response by tracing the propagation of not only jet shower partons but also the recoil partons. It also includes the back reaction, called "negative particles", to conserve the total energy and momentum during each scattering. For the current stage, it only considers the interactions between an energetic jet parton (E > 2 GeV) and a thermal parton in local equilibrium in the medium, and any interactions between jet partons are neglected. 

The LBT model has been widely used to calculate hadron and jet suppresion, anisotropy and correlations with $\gamma/Z$ boson. It has been validated to probe the light/heavy flavor parton-medium interactions by comfronting the experimental data. 

The LBT project was initiated and conducted by Prof. Xin-Nian Wang, a senior scientist at LBNL, and implemented with $t-$ channel elastic scattering and medium-induced gluon radiation by his students or collaborators Han-Li Lin and Yan Zhu written in Fortran. To include the complete elastic scattering channels and make the code more standard and modern, it was re-written by Tan Luo and Yayun He. Afterwards, it was updated by Tan Luo and Shanshan Cao in order to describe heavy flavor transport.

## Installation
It's easy to install the LBT codes on one's computer or the server. The codes have integrated with PYTHIA 8 and CERN ROOT, so one has to install them at first. I have also placed a clean version free from them, see examples. Since the size of radiation table excesses the limit for uploading, anyone who wants to run the LBT codes can contact the authors.

1. Download or git clone the repository to one's directory.
2. cd codes
3. make

Then the executive $\textit{lbt}$ can be run on the computer.


## Parameters
The parameter settings that a new user can adjust is in the Codes/LBTinput.txt. Currently, it's only for a static and uniform medium. For a more realistic ideal/viscous hydro background, one's has to modify the interface in Codes/LBT.cpp. 

## Run an example
Here is a simple example for a single jet shower parton propagation in a uniform and static medium.

cd ./examples/codes
make
./lbt 0

