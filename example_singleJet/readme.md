This interface in main.cpp is well written so that parton showers generated from PYTHIA 8 can directly propagate in the LBT model. One can analyze the output for jet or parton observables in heavy-ion collisions. The interface and analysis examples for dijet and gamma/Z-triggered jet will also be uploaded in a few weeks.

We add smooth hydro profiles for Pb+Pb 5.02 TeV 0-10% collisions in this project. Hydro profiles for other center-of-mass energy and centrality bins will be uploaded in a few weeks.

It's simple to run the LBT codes. Here are the steps:

1 make clean && make
2 modify parameters in m01.cmnd for pp collisions and LBT.in for AA collisions, and make sure the number of events is the same for both pp and AA
modify the path for scattering rate tables and hydro profiles
3 ./lbt
The output files are pp.txt for pp collisions, positive.txt and negative.txt for AA collisions.

To analyze the final particles, calculate the pT spectrum for jets as an example in dNdpT.cc

1 modify the path of the final particles, the number of events in dNdpT.cc
2 make clean && make
3 ./dNdpT
One should note that the back reaction in negative.txt should be subtracted during jet reconstruction.

For the use of the LBT model, please kindly cite our papers:
1 LBT elastic and recoil: https://arxiv.org/abs/1503.03313
2 LBT radiation and applications: https://arxiv.org/abs/2306.13742

For more detailed calculations of single jet suppression and jet anisotropy flow in the LBT, one can refer
1 single jet: https://arxiv.org/abs/1809.02525
2 jet anisotropy: https://arxiv.org/abs/2201.08408

The hydro profiles are generated within CLVisc hydrodynamic model, one can refer
1 https://arxiv.org/abs/1411.7767
2 https://arxiv.org/abs/1205.5019
3 https://arxiv.org/abs/1802.04449
