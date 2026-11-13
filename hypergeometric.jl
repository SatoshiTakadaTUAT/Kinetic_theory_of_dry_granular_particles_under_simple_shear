import Pkg; Pkg.add("HypergeometricFunctions")

using SpecialFunctions
using HypergeometricFunctions

const phi_f = 0.49
const phi_J = 0.639

##################################################
# Radial distribution function
##################################################

function g0(phi)
	
	if phi <= phi_f
		return (1 - phi/2) / (1 - phi)^3
	else
		return g0(phi_f) * (phi_J - phi_f) / (phi_J - phi)
	end
end

##################################################
# Hypergeometric functions
##################################################
function F(s, t, u, theta)

	a = ((s+1)/2, (t+1)/2)
	b = ((s+t+u+3)/4, (s+t+u+5)/4)
	z = -1/(16*theta)

	return pFq(a, b, z)
end

function G(s, t, u, theta)
	
	a = (1/2, (s+2)/2, (t+2)/2)
	b = (3/2, (s+t+u+5)/4, (s+t+u+7)/4)
	z = -1/(16*theta)

	return pFq(a, b, z)
end

##################################################
# Auxiliary functions
##################################################
function H1(s, t, u, theta)

	coeff = 1 / (3+s+t+u)
	val = (1 + s)*F(s+2, t, u, theta) - (1 + u)*F(s, t, u+2, theta)

	return coeff*val
end

function H2(s, t, u, theta)
	
	coeff = 1 / (5+s+t+u)
	val   = (2+s)*G(s+2, t, u, theta) - (1+u)*G(s, t, u+2, theta)
	
	return coeff*val
end

function H3(s, t, u, theta)
	
	coeff = 1 / ((3+s+t+u)*(5+s+t+u))
	val   = (1+s)*(3+s)*F(s+4, t, u, theta) - 2*(1+s)*(1+u)*F(s+2, t, u+2, theta) + (1+u)*(3+u)*F(s, t, u+4, theta)
	
	return coeff*val
end

function H4(s, t, u, theta)
	
	coeff = 1 / ((3+s+t+u)*(5+s+t+u))
	val   = (1+s)*(1+t)*F(s+2, t+2, u, theta) - (1+s)*(1+u)*F(s+2, t, u+2, theta) - (1+t)*(1+u)*F(s, t+2, u+2, theta) + (1+u)*(3+u)*F(s, t, u+4, theta)
	
	return coeff*val
end

##################################################
# Lambda_aa
##################################################
function Lambda_aa0(e, theta)
	
	coeff1 = -2*sqrt(2*pi) / (15*theta)
	val1   = F(2, 2, 0, theta) + G(1, 1, 0, theta) + (1/(14*theta))*G(3, 3, 0, theta)
	
	coeff2 = (1-e)*4*sqrt(2*pi)
	val2   = F(0, 0, 0, theta) + (1/(60*theta))*F(2, 2, 0, theta) + (1/(20*theta))*G(1, 1, 0, theta) + (1/(840*theta^2))*G(3, 3, 0, theta)
	
	return coeff1*val1 + coeff2*val2
end

function Lambda_aa1(e, theta)
	
	coeff1 = -2*sqrt(2*pi) / (15*theta)
	val1   = H2(1, 1, 0, theta)
	
	coeff2 = (1-e)*4*sqrt(2*pi)
	val2   = (3/2)*H1(0, 0, 0, theta) + (1/(20*theta))*H2(1, 1, 0, theta)
	
	return coeff1*val1 + coeff2*val2
end

function Lambda_aa2(e, theta)
	
	coeff = -2*sqrt(2)*pi*(1-3*e) / (15*sqrt(theta))
	
	return coeff
end

function Lambda_aa3(e, theta)
	
	coeff1 = sqrt(2*pi) / (60*theta)
	val1   = H3(2, 2, 0, theta)
	
	coeff2 = (1-e)*3*sqrt(2*pi) / 4
	val2   = H3(0, 0, 0, theta)
	
	return coeff1*val1 + coeff2*val2
end

function Lambda_aa4(e, theta)
	
	coeff1 = sqrt(2*pi) / (30*theta)
	val1   = H4(2, 2, 0, theta)
	
	coeff2 = (1-e)*3*sqrt(2*pi) / 2
	val2   = H4(0, 0, 0, theta)
	
	return coeff1*val1 + coeff2*val2
end

function Lambda_aa5(e, theta)
	
	coeff1 = sqrt(2*pi) / (105*theta)
	val1   = F(4, 4, 0, theta)
	
	coeff2 = (1-e)*sqrt(2*pi) / 5
	val2   = F(2, 2, 0, theta)
	
	return coeff1*val1 + coeff2*val2
end

function Lambda_aa(phi, e, theta, Pi_xx, Pi_yy, Pi_xy)

	coeff = (3*sqrt(2) / pi)*(1+e)*phi*g0(phi)*theta*sqrt(theta)
	val   = Lambda_aa0(e, theta) + Lambda_aa1(e, theta)*(Pi_xx + Pi_yy) + Lambda_aa2(e, theta)*Pi_xy + Lambda_aa3(e, theta)*(Pi_xx^2 + Pi_yy^2) + Lambda_aa4(e, theta)*Pi_xx*Pi_yy + Lambda_aa5(e, theta)*Pi_xy^2
	
	return coeff*val
end

##################################################
# Lambda_xy
##################################################
function Lambda_xy0(e, theta)
	
	coeff = -sqrt(2)*pi / (15*sqrt(theta))
	val   = (1 - 3*e) - (1/(14*theta))*(1 + e)
	
	return coeff*val
end

function Lambda_xy1(e, theta)
	
	coeff = -2*sqrt(2)*pi / (105*sqrt(theta)) * (4 - 3*e)
	
	return coeff
end

function Lambda_xy2(e, theta)
	
	coeff1 = sqrt(2*pi)
	val1   = (8/3)*F(2, 0, 0, theta) - (16/15)*F(2, 2, 0, theta) + (4/(35*theta))*G(3, 1, 0, theta) - (12/(105*theta))*G(3, 3, 0, theta)
	
	coeff2 = (1-e)*sqrt(2*pi)
	val2   = (4/5)*F(2, 2, 0, theta) + (2/(35*theta))*G(3, 3, 0, theta)
	
	return coeff1*val1 + coeff2*val2
end

function Lambda_xy3(e, theta)
	
	coeff1 = sqrt(2*pi)
	val1   = (2/15)*F(2, 2, 0, theta) + (1/3)*H1(2, 0, 0, theta) + (1/3)*H1(0, 2, 0, theta) - (8/15)*H1(2, 2, 0, theta) + (1/(105*theta))*H1(4, 4, 0, theta)
	
	coeff2 = (1 - e)*sqrt(2*pi) / 5
	val2   = H1(2, 2, 0, theta)
	
	return coeff1*val1 + coeff2*val2
end

function Lambda_xy(phi, e, theta, Pi_xx, Pi_yy, Pi_xy)

	coeff = (3*sqrt(2) / pi)*(1+e)*phi*g0(phi)*theta*sqrt(theta)
	val   = Lambda_xy0(e, theta) + Lambda_xy1(e, theta)*(Pi_xx + Pi_yy) + Lambda_xy2(e, theta)*Pi_xy + Lambda_xy3(e, theta)*(Pi_xx + Pi_yy)*Pi_xy
	
	return coeff*val
end

##################################################
# Lambda_xx
##################################################
function Lambda_xx0(e, theta)
	
	coeff1 = -2*sqrt(2*pi) / (35*theta)
	val1   = F(4, 2, 0, theta) + G(3, 1, 0, theta) + (5/(66*theta))*G(5, 3, 0, theta)
	
	coeff2 = (1 - e)*4*sqrt(2*pi)/3
	val2   = F(2, 0, 0, theta) + (3/(140*theta))*F(4, 2, 0, theta) + (9/(140*theta))*G(3, 1, 0, theta) + (1/(616*theta^2))*G(5, 3, 0, theta)
	
	return coeff1*val1 + coeff2*val2
end

function Lambda_xx1(e, theta)
	
	coeff1 = sqrt(2*pi)
	val1   = (8/3)*F(2, 0, 0, theta) + (4/(35*theta))*G(3, 1, 0, theta) - (8/3)*H1(2, 0, 0, theta) - (6/(35*theta))*H2(3, 1, 0, theta)
	
	coeff2 = (1 - e)*2*sqrt(2*pi)
	val2   = H1(2, 0, 0, theta) + (3/(70*theta))*H2(3, 1, 0, theta)
	
	return coeff1*val1 + coeff2*val2
end

function Lambda_xx2(e, theta)
	
	coeff1 = -sqrt(2*pi)
	val1   = (8/3)*H1(0, 2, 0, theta) + (6/(35*theta))*H2(1, 3, 0, theta)
	
	coeff2 = (1 - e)*2*sqrt(2*pi)
	val2   = H1(0, 2, 0, theta) + (3/(70*theta))*H2(1, 3, 0, theta)
	
	return coeff1*val1 + coeff2*val2
end

function Lambda_xx3(e, theta)
	
	coeff = -2*sqrt(2)*pi / (105*sqrt(theta)) * (5 - 9*e)
	
	return coeff
end

function Lambda_xx4(e, theta)
	coeff1 = 2*sqrt(2*pi) / 3
	val1   = H1(2, 0, 0, theta) - H3(2, 0, 0, theta) + (3/(280*theta))*H3(4, 2, 0, theta)
	
	coeff2 = (1 - e)*sqrt(2*pi) / 4
	val2   = H3(2, 0, 0, theta)
	
	return coeff1*val1 + coeff2*val2
end

function Lambda_xx5(e, theta)
	
	coeff1 = 2*sqrt(2*pi) / 3
	val1   = -H3(0, 2, 0, theta) + (3/(280*theta))*H3(2, 4, 0, theta)
	
	coeff2 = (1 - e)*sqrt(2*pi) / 4
	val2   = H3(0, 2, 0, theta)
	
	return coeff1*val1 + coeff2*val2
end

function Lambda_xx6(e, theta)
	
	coeff1 = 2*sqrt(2*pi) / 3
	val1   = H1(0, 2, 0, theta) - 2*H4(2, 0, 0, theta) + (3/(140*theta))*H4(4, 2, 0, theta)
	
	coeff2 = (1 - e)*sqrt(2*pi) / 2
	val2   = H4(2, 0, 0, theta)
	
	return coeff1*val1 + coeff2*val2
end

function Lambda_xx7(e, theta)
	
	coeff1 = sqrt(2*pi)
	val1   = (4/15)*F(2, 2, 0, theta) - (8/35)*F(4, 2, 0, theta) + (1/(231*theta))*F(6, 4, 0, theta)
	
	coeff2 = (1 - e)*3*sqrt(2*pi) / 35
	val2   = F(4, 2, 0, theta)
	
	return coeff1*val1 + coeff2*val2
end

function Lambda_xx(phi, e, theta, Pi_xx, Pi_yy, Pi_xy)

	coeff = (3*sqrt(2) / pi)*(1+e)*phi*g0(phi)*theta*sqrt(theta)
	val   = Lambda_xx0(e, theta) + Lambda_xx1(e, theta)*Pi_xx + Lambda_xx2(e, theta)*Pi_yy + Lambda_xx3(e, theta)*Pi_xy + Lambda_xx4(e, theta)*Pi_xx^2 + Lambda_xx5(e, theta)*Pi_yy^2 + Lambda_xx6(e, theta)*Pi_xx*Pi_yy + Lambda_xx7(e, theta)*Pi_xy^2
	
	return coeff*val
end

##################################################
# Lambda_yy
##################################################
function Lambda_yy0(e, theta)
	return Lambda_xx0(e, theta)
end

function Lambda_yy1(e, theta)
	return Lambda_xx2(e, theta)
end

function Lambda_yy2(e, theta)
	return Lambda_xx1(e, theta)
end

function Lambda_yy3(e, theta)
	return Lambda_xx3(e, theta)
end

function Lambda_yy4(e, theta)
	return Lambda_xx5(e, theta)
end

function Lambda_yy5(e, theta)
	return Lambda_xx4(e, theta)
end

function Lambda_yy6(e, theta)
	return Lambda_xx6(e, theta)
end

function Lambda_yy7(e, theta)
	return Lambda_xx7(e, theta)
end

function Lambda_yy(phi, e, theta, Pi_xx, Pi_yy, Pi_xy)

	coeff = (3*sqrt(2) / pi)*(1+e)*phi*g0(phi)*theta*sqrt(theta)
	val   = Lambda_yy0(e, theta) + Lambda_yy1(e, theta)*Pi_xx + Lambda_yy2(e, theta)*Pi_yy + Lambda_yy3(e, theta)*Pi_xy + Lambda_yy4(e, theta)*Pi_xx^2 + Lambda_yy5(e, theta)*Pi_yy^2 + Lambda_yy6(e, theta)*Pi_xx*Pi_yy + Lambda_yy7(e, theta)*Pi_xy^2
	
	return coeff*val
end

##################################################
# Pi^c
##################################################
function Pi_c0(theta)

	coeff = 2*pi / 3
	val   = 1 + 1/(30*theta)
	
	return coeff*val
end

function Pi_c1(theta)

	coeff = -(4*sqrt(pi)) / (45*sqrt(theta))
	val   = G(1, 1, 0, theta)
	
	return coeff*val
end

function Pi_c2(theta)

	coeff = sqrt(pi) / (45*sqrt(theta))
	val   = H1(2, 2, 0, theta)
	
	return coeff*val
end

function Pi_c(phi, e, theta, Pi_xx, Pi_yy, Pi_xy)

	coeff = (3/pi)*(1+e)*phi*g0(phi)
	val   = Pi_c0(theta) + Pi_c1(theta)*Pi_xy + Pi_c2(theta)*(Pi_xx + Pi_yy)*Pi_xy
	
	return coeff*val
end

##################################################
# Pi^c_xy
##################################################
function Pi_c_xy0(theta)

	coeff = -2*sqrt(pi) / (15*sqrt(theta))
	val   = F(2, 2, 0, theta) + G(1, 1, 0, theta) + (1/(14*theta))*G(3, 3, 0, theta)
	
	return coeff*val
end

function Pi_c_xy1(theta)

	coeff = -2*sqrt(pi) / (15*sqrt(theta))
	val   = H2(1, 1, 0, theta)
	
	return coeff*val
end

function Pi_c_xy2(theta)

	return 4*pi / 15
end

function Pi_c_xy3(theta)

	coeff = sqrt(pi) / (60*sqrt(theta))
	val   = H3(2, 2, 0, theta)
	
	return coeff*val
end

function Pi_c_xy4(theta)

	coeff = sqrt(pi) / (30*sqrt(theta))
	val   = H4(2, 2, 0, theta)
	
	return coeff*val
end

function Pi_c_xy5(theta)

	coeff = sqrt(pi) / (105*sqrt(theta))
	val   = F(4, 4, 0, theta)
	
	return coeff*val
end

function Pi_c_xy(phi, e, theta, Pi_xx, Pi_yy, Pi_xy)

	coeff = (3/pi)*(1+e)*phi*g0(phi)
	val   = Pi_c_xy0(theta) + Pi_c_xy1(theta)*(Pi_xx + Pi_yy) + Pi_c_xy2(theta)*Pi_xy + Pi_c_xy3(theta)*(Pi_xx^2 + Pi_yy^2) + Pi_c_xy4(theta)*Pi_xx*Pi_yy + Pi_c_xy5(theta)*Pi_xy^2
	
	return coeff*val
end

##################################################
# Pi^c_xx
##################################################
function Pi_c_xx0(theta)

	return (2*pi/3)*(1 + 3/(70*theta))
end

function Pi_c_xx1(theta)

	return 4*pi/15
end

function Pi_c_xx2(theta)

	coeff = -4*sqrt(pi) / (35*sqrt(theta))
	val   = G(3, 1, 0, theta)
	
	return coeff*val
end

function Pi_c_xx3(theta)

	coeff = sqrt(pi) / (35*sqrt(theta))
	val   = H1(4, 2, 0, theta)
	
	return coeff*val
end

function Pi_c_xx4(theta)

	coeff = sqrt(pi) / (35*sqrt(theta))
	val   = H1(2, 4, 0, theta)
	
	return coeff*val
end

function Pi_c_xx(phi, e, theta, Pi_xx, Pi_yy, Pi_xy)

	coeff = (3/pi)*(1+e)*phi*g0(phi)
	val   = Pi_c_xx0(theta) + Pi_c_xx1(theta)*Pi_xx + Pi_c_xx2(theta)*Pi_xy + Pi_c_xx3(theta)*Pi_xx*Pi_xy + Pi_c_xx4(theta)*Pi_yy*Pi_xy
	
	return coeff*val
end

##################################################
# Pi^c_yy
##################################################
function Pi_c_yy0(theta)

	return Pi_c_xx0(theta)
end

function Pi_c_yy1(theta)

	return Pi_c_xx1(theta)
end

function Pi_c_yy2(theta)

	return Pi_c_xx2(theta)
end

function Pi_c_yy3(theta)

	return Pi_c_xx4(theta)
end

function Pi_c_yy4(theta)

	return Pi_c_xx3(theta)
end

function Pi_c_yy(phi, e, theta, Pi_xx, Pi_yy, Pi_xy)

	coeff = (3/pi)*(1+e)*phi*g0(phi)
	val   = Pi_c_yy0(theta) + Pi_c_yy1(theta)*Pi_yy + Pi_c_yy2(theta)*Pi_xy + Pi_c_yy3(theta)*Pi_xx*Pi_xy + Pi_c_yy4(theta)*Pi_yy*Pi_xy
	
	return coeff*val
end
