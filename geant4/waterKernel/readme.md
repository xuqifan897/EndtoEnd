# Introduction
The purpose of this subproject is to calculate the point dose kernel in water. That is to calculate the dose distribution assuming that the photon interacts with the water molecule at a fixed point. This kernel is used to convolve with the Terma distribution to dose distribution for mono-energetic photons. For poly-energetic photons, approximations should be made both in the the mass attenuation coefficient $\mu$ and the point dose kernel.

# First step
The first step is to create a simple case, in which mono-energetic photons are generated and interact with water molecules. No extra functionalities are required.