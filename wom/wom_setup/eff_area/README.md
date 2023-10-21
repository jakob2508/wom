# ppc_wom/wom_setup/eff_area

# eff_area_wom_air.csv
Contains the wavelength [nm] and the effective area [cm^2] of the WOM tube embedded in air as reported in the WOM paper (arXiv:2112.12258).
The data has been obtained by digitizing Fig. 13.


# eff_area.py
However, in ice the effective area is reduced almost by a factor of 10 as transmission from the vessel 
back into the ice is more likely to happen. We therefore scale the effective area such that the module reaches a maximum effective area of 18cm^2.

Additionally, because the WOM tube is only sensitive on the inner tube (height = 76 cm) on not on the full height of the module's outer vessel (height = 130 cm)
we shrink the sensor geometry to 76 cm but upscale the effective area by a geometry scaling factor of 130/76 to ensure that we collect the same amount of photons.

This script plots is the rescaled effective area and the Cherenkov weighted, rescaled effective area.

# om.wv_140.0
The rescaled effective area is then parsed into the file om.wv_140.0 which contains the wavelength and effective area of the WOM in ice in steps of 1 nm.
