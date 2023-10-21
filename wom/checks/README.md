# ppc geometry
As a sanity check it is sufficient to simulate only a small section of the IceCube detector. The vanilla simulation geometry comprises 4 WOMs on two strings spaced 1 m apart in the horizontal and 1 m apart in the vertical 
(z axis). Even though the WOM does not have DOM-like flashers, ppc does not know that and will assume that the flasher is the same as for the DOM. To simplify the geometry we modify the position of the flasher to the center of 
the module. This is done in pro.cxx by setting FLZ=0.0f, FLR=0.0f, FLB=0.0f above the statement if(p.ka>999.f). We set the zenith angle of the flasher with FZCR, the azimuthal angle with FDLR assuming an angle defined by the
positve x-direction counter-clockwise and FWID for the width of the flasher, here set to FDLR = 0.001. Furthermore, we set the absorption and scattering coefficients close to zero (1E-10) for all ice layers in icemodel.dat. 

Finally, we run ppc in command line and read out the hits with a dedicated function later with a data reader. For the different checks we use modified simulation setups depending on which quantity we are interested in checking.

For all but the effective area sanity check we run the simulation with the flasher set to 300 nm. This is done by modifying the om.wv_140.0 file which is saved in the respected directory. The om.wv file only contains three points [299,300,301] nm where the effective area is the absolute effective area from the original file at 300 nm but weighted with a pdf centerred around 300 nm with the entries [0.25,0.5,0.25]. 

We further adapt PPC to print out the hit position and direction in cartesian coordinates x,y,z which are accessible in r[0], r[1], r[2] and nx, ny, nz in the f2k.cxx file. Additionally, to include the PMT selection probability when assigning a hit to a PMT we implement the function WOM_z_dependence at the very top of f2k.cxx which is called later in the code to select a PMT based on rejection sampling. The parameters of the PMT selection probability function, the WOM height and the slope are currently defined as constants but eventually would be read in from input files (om.c

We set the effective scattering and absorption coefficient in icemodel.dat to 1E-10 (nearly zero) and turn of birefringence in the command line input (BFRM=2).
icemodel.dat is the unmodified file including scatting and absorption
icemodel_no_sca_abs.dat is the icemodel file where sca + abs are set to 1E-10 (disabled). We create the file using manipulate_dat_file.py
