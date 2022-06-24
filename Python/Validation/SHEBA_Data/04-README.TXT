README.TXT
06 Nov 2003
Snow and ice thickness and ocean heat flux at the SHEBA site Pittsburgh.
Compiled by Hendrik Huwald, huwald@env.ethz.ch

i) Time convention:
-------------------

Julian day (jday) is defined as follows:
Reference date (JD 0.0) = January 1, 1997, 12:00 am,
i.e. January 1, 1997, 01:00 am = JD 0.0417
and  January 1, 1997, 11:00 pm = JD 0.9583
JD 303. = October 31, 1997

ii) Snow / ice thickness: (hs-hi_thickness.dat)
-------------------------

The initial snow/ice interface is at z=0 on the vertical axis z.
hs denotes the level of the snow surface,
hf denotes the level of the snow/ice interface 
   (or the ice surface if no snow is present)
hi denotes the level of the ice base

All values are given in meters.
Thickness may be calculated by adding levels.

a) the ice base is the mean of thickness gauges no.53, 69 and 71 
   at the Pittsburgh site, shifted upward by 0.1166m.
   Initial ice thickness on JD 303. is 1.94m
b) the ice surface is the level measured at gauge no.69 at the Pittsburgh site
   including a few corrections, i.e. ice growth during melt period removed
c) the level of the snow surface is determined from an analysis of the
   internal temperature profile (change of gradient)

All observations are linearly interpolated on a 1hr time step

iii) Ocean heat flux (Fo): (ocean_heat_flux.dat)
--------------------------

The Fo provided was determined from observations at the Pittsburgh site.
Fo is calculated as a residual of the basal energy balance using

1) the energy of melt at the ice base
2) the conductive heat flux at the ice base

The resulting ocean heat flux was smoothed with a 30-days running mean.
The residual method can result in negative values for Fo.
In these (few) cases, Fo was set equal to zero.

1) energy of melt

- a best estimate ice base (hs-hi_thickness.dat) is used
- growth/melt rate at time step i is calculated from the difference
  hi(i)-hi(i-1), where hi is the level of the ice base
- basal ice salinity: 4ppt
- the latent heat is computed according to Bitz and Lipscomb, 1999,
  resulting in a value of 0.87*L0, 
  where L0 is the latent heat of fusion of pure ice

2) conductive heat flux

- "reference level": 0.3m above the ice base
- linear interpolation of the observed internal ice temperature at the levels:
  a) ice base, b) ice base + 0.3m 
- approximation of the temperature gradient at the ice base using a linear 
  (first order) difference scheme
- the thermal conductivity is calculated according Bitz and Lipscomb, 1999,
  resulting in a value of 1.74 W/m/K (basal ice: -1.6C, 4ppt).
