# ll2cmap
 Program to convert a global latitude-longitude ARL formatted
 meteorology file for HYSPLIT to an ARL formatted file on a conformal
 projection, either polar stereographic, Mercator, or Lambert. The 
 converted file will be a geographic extract of
 the orginal, containing the minimum number of variables required
 to run HYSPLIT. Additional variables may be added to the output
 as a command line argument. The vertical coordinate of the output
 data will be unchanged unless the pressure to sigma interpolation
 option is selected, available only for pressure level input data.
 The remapping uses bi-linear interpolation of the meteorology
 variables unless the nearest neighbor option is selected.
