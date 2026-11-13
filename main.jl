include("hypergeometric.jl")
include("time_evol.jl")
include("rheology.jl")
using Printf

##################################################
# parameter
##################################################
phi_list1 = LinRange(0.01,0.5,50)
e_list1 = [0.7, 0.80, 0.90, 0.95]#[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99]

phi_list2 = [0.40]#[0.01, 0.10, 0.30, 0.40, 0.45]
e_list2 = LinRange(0, 0.99, 100)

##################################################
# run
##################################################
function run1(phi_list, e_list)

	# initial
	theta0, Del_theta0, del_theta0, Pi_xy0 = Bagnold(phi_list1[1], e_list1[1])
	
	for e in e_list
		println("$(e)")
		filename = "e_$(@sprintf("%.2f", e)).txt"
		
		open(filename, "w") do io
			
			for phi in phi_list
				println("$(phi)")
				
				theta, Del_theta, del_theta, Pi_xy = cal_evol(phi, e)
				
				Pi_xx  = (    Del_theta/3 + del_theta/3) / theta
				Pi_yy  = (- 2*Del_theta/3 + del_theta/3) / theta
				eta    = cal_eta(phi, e, theta, Pi_xx, Pi_yy, Pi_xy)
				N1     = cal_N1(phi, e, theta, Pi_xx, Pi_yy, Pi_xy)
				N2     = cal_N2(phi, e, theta, Pi_xx, Pi_yy, Pi_xy)
				mu     = cal_mu(phi, e, theta, Pi_xx, Pi_yy, Pi_xy)
				angle  = cal_angle(phi, e, theta, Pi_xx, Pi_yy, Pi_xy)
				
				@printf(io, "%.2f\t%s\t%s\t%s\t%s\t%s\t%s\n",
					phi, theta, eta, N1, N2, mu, angle)
			end
		end
	end
end

function run2(phi_list, e_list)

	# initial
	theta0, Del_theta0, del_theta0, Pi_xy0 = Bagnold(phi_list[1], e_list[1])
	
	for phi in phi_list
		println("$(phi)")
		filename = "phi_$(@sprintf("%.2f", phi)).txt"
		
		open(filename, "w") do io
			
			for e in e_list
				println("$(e)")
				
				theta, Del_theta, del_theta, Pi_xy = cal_evol(phi, e)
				
				Pi_xx  = (    Del_theta/3 + del_theta/3) / theta
				Pi_yy  = (- 2*Del_theta/3 + del_theta/3) / theta
				eta    = cal_eta(phi, e, theta, Pi_xx, Pi_yy, Pi_xy)
				N1     = cal_N1(phi, e, theta, Pi_xx, Pi_yy, Pi_xy)
				N2     = cal_N2(phi, e, theta, Pi_xx, Pi_yy, Pi_xy)
				mu     = cal_mu(phi, e, theta, Pi_xx, Pi_yy, Pi_xy)
				angle  = cal_angle(phi, e, theta, Pi_xx, Pi_yy, Pi_xy)
				
				@printf(io, "%.2f\t%s\t%s\t%s\t%s\t%s\t%s\n",
					e, theta, eta, N1, N2, mu, angle)
			end
		end
	end
end

##################################################
# main
##################################################
run1(phi_list1, e_list1)
run2(phi_list2, e_list2)