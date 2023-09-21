#### Code for calculating specific surface area from WP4C

# S = specific surface area (m^2/kg)
# p = density of water (1000 kg/m^2)
# k = hamaker constant (-6 * 10^-20 J)
# y = water potential (J/kg) 
  # note - water potential measured as MPa on WP4C
  # 1 J/kg = 1 kPa at p = 1000 kg/m^2
# = water content (g/g)

## Enter values here

#Mass of Tray

tray_g = 24.71

#Tray + Wet Sample

tray_wet_sed_g = 26.02

#Tray + Dry Sample

tray_dry_sed_g = 25.92


## water content

w = (tray_wet_sed_g - tray_dry_sed_g)/(tray_dry_sed_g - tray_g)

## water potential - enter 

y_MPa = -201.06

y_kPa = y_MPa * 1000

## Do not change these constants
k = as.numeric(-6*(10^-20))
p = 1000


##update to a function at some point

S = (w/((k/(6 * pi* (p * y_kPa)))^(1/3))*p)

