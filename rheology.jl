##################################################
# eta
##################################################
function cal_eta(phi, e, theta, Pi_xx, Pi_yy, Pi_xy)

	return -theta*(Pi_xy + Pi_c_xy(phi, e, theta, Pi_xx, Pi_yy, Pi_xy))
end

##################################################
# N1
##################################################
function cal_N1(phi, e, theta, Pi_xx, Pi_yy, Pi_xy)

	denom = Pi_xx + Pi_c_xx(phi, e, theta, Pi_xx, Pi_yy, Pi_xy) - (Pi_yy + Pi_c_yy(phi, e, theta, Pi_xx, Pi_yy, Pi_xy))
	
	return denom / (1 + Pi_c(phi, e, theta, Pi_xx, Pi_yy, Pi_xy))
end

##################################################
# N2
##################################################
function cal_N2(phi, e, theta, Pi_xx, Pi_yy, Pi_xy)

	Pi_zz   = -(Pi_xx + Pi_yy)
	Pi_c_zz = 3*Pi_c(phi, e, theta, Pi_xx, Pi_yy, Pi_xy) - Pi_c_xx(phi, e, theta, Pi_xx, Pi_yy, Pi_xy) - Pi_c_yy(phi, e, theta, Pi_xx, Pi_yy, Pi_xy)
	
	denom = Pi_yy + Pi_c_yy(phi, e, theta, Pi_xx, Pi_yy, Pi_xy) - (Pi_zz + Pi_c_zz)
	
	return denom / (1 + Pi_c(phi, e, theta, Pi_xx, Pi_yy, Pi_xy))
end

##################################################
# mu
##################################################
function cal_mu(phi, e, theta, Pi_xx, Pi_yy, Pi_xy)

	return -(Pi_xy + Pi_c_xy(phi, e, theta, Pi_xx, Pi_yy, Pi_xy)) / (1 + Pi_c(phi, e, theta, Pi_xx, Pi_yy, Pi_xy))
end

##################################################
# angle
##################################################
function cal_angle(phi, e, theta, Pi_xx, Pi_yy, Pi_xy)

	val1 = Pi_xx + Pi_c_xx(phi, e, theta, Pi_xx, Pi_yy, Pi_xy) - (Pi_yy + Pi_c_yy(phi, e, theta, Pi_xx, Pi_yy, Pi_xy))
	val2 = 2*(Pi_xy + Pi_c_xy(phi, e, theta, Pi_xx, Pi_yy, Pi_xy))
	return -(1/2)*atan(val2, val1)
end