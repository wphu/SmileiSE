from scipy import constants as const

value,unit,precision=const.physical_constants["speed of light in vacuum"]
c=value

value,unit,precision=const.physical_constants["electron mass"]
e_m=value

value,unit,precision=const.physical_constants["elementary charge"]
e_e=value

value,unit,precision=const.physical_constants["electric constant"]
ephi0=value

omiga0=c/1e-5

x=e_e/(e_m*c*c)
y=ephi0*e_m*omiga0*omiga0/(e_e*e_e)


print x
print c
print e_m
print e_e
print omiga0

print 1e19/y
##print const.physical_constants.viewkeys()
