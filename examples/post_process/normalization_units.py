from scipy import constants as const

value,unit,precision=const.physical_constants["speed of light in vacuum"]
c=value

value,unit,precision=const.physical_constants["electron mass"]
e_m=value

value,unit,precision=const.physical_constants["elementary charge"]
e_e=value

value,unit,precision=const.physical_constants["electric constant"]
ephi0=value

#> the reference frequency: s^-1
omiga0 = c/1.0

#> normalized temerature unit: eV
norm_temp = (e_m*c*c)/e_e

#> normalized density unit: m^-3
norm_n = ephi0*e_m*omiga0*omiga0/(e_e*e_e)

#> normalized length unit: m
norm_l = c/omiga0

#> normalized time unit: s
norm_t = 1.0/omiga0

#> normalized electric field unit:
norm_E = e_m*c*omiga0/e_e

#> normalized potential unit:
norm_U = norm_E*norm_l



if __name__ == '__main__':
	print "norm_U = ", norm_U
	print "norm_n = ", norm_n
	print "norm_E = ", norm_E
	print "norm_temp = ", norm_temp


##print const.physical_constants.viewkeys()
