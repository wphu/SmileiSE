import math

mass_e = 9.109382616e-31
mass_D1 = 2.0 * 1.67262158e-27
T_e = 24.8
T_D1 = 65.5

const_pi = 3.1415926

potential_debye = -0.5 * T_e * math.log( ( 2.0 * const_pi * mass_e / mass_D1 ) * ( 1.0 + T_D1 / T_e ) )
potential_total = potential_debye + 0.7 * T_e

print "potential_debye = ", potential_debye
print "potential_total", potential_total