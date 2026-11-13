using Printf

threshold = 1e-10
step_interval = 1e4

function Bagnold(phi, e)

	theta_Bag     = 5*pi*(2+e) / (432*(1-e)*(1+e)^2*(3-e)^2*phi^2*g0(phi)^2)
	Del_theta_Bag = 25*pi*(2+e) / (432*(1+e)^2*(3-e)^3*phi^2*g0(phi)^2)
	del_theta_Bag = Del_theta_Bag
	Pi_xy_Bag     = -sqrt(5*(1-e)*(2+e) / (3*(3-e)^2))
	
	return theta_Bag, Del_theta_Bag, del_theta_Bag, Pi_xy_Bag
end

function cal_RHS(phi, e, theta, Del_theta, del_theta, Pi_xy)
	
	Pi_xx = (    Del_theta/3 + del_theta/3) / theta
	Pi_yy = (- 2*Del_theta/3 + del_theta/3) / theta
	
	Lam_aa = Lambda_aa(phi, e, theta, Pi_xx, Pi_yy, Pi_xy)
	Lam_xx = Lambda_xx(phi, e, theta, Pi_xx, Pi_yy, Pi_xy)
	Lam_yy = Lambda_yy(phi, e, theta, Pi_xx, Pi_yy, Pi_xy)
	Lam_xy = Lambda_xy(phi, e, theta, Pi_xx, Pi_yy, Pi_xy)
	
	RHS1 = -(2/3)*theta*Pi_xy - (1/3)*Lam_aa
	RHS2 = -2*theta*Pi_xy - (Lam_xx - Lam_yy)
	RHS3 = -2*theta*Pi_xy - (2*Lam_xx + Lam_yy - Lam_aa)
	RHS4 = (-RHS1*Pi_xy - (theta - (2/3)*Del_theta  + (1/3)*del_theta) - Lam_xy) / theta
	
	return RHS1, RHS2, RHS3, RHS4
end

function cal_evol(phi, e)

	println("$(phi), $(e)")
	
	theta, Del_theta, del_theta, Pi_xy = Bagnold(phi, e)
	theta_old = 0
	step = 0

	while true
	
		step += 1
		
		RHS1, RHS2, RHS3, RHS4 = cal_RHS(phi, e, theta, Del_theta, del_theta, Pi_xy)
		RHS_MAX = max(abs(RHS1), abs(RHS2), abs(RHS3), abs(RHS4))
		RHS_MAX = 10^ceil(log10(RHS_MAX))
		
		#Base.prompt("Stop")
		
		dt = min(1e-2, 1e-2 / RHS_MAX)
		
		theta     += RHS1*dt
		Del_theta += RHS2*dt
		del_theta += RHS3*dt
		Pi_xy     += RHS4*dt
		
		if step % step_interval == 0
		
			error = abs(theta - theta_old)/theta
			
			#@printf("%.8e\t%.3e\n", theta, error)
			
			if error < threshold
				break
			end
			
			theta_old = theta
		end
		
		error = abs(theta - theta_old)/theta
	end
	
	return theta, Del_theta, del_theta, Pi_xy
end