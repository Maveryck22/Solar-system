# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# This code was created to run with Julia V0.6.4. Although the logic can be applied with any ohter
# language, such as C++, scilab, matlab, etc.
# The forward-Euler method is used to calculate position (x,y,z), velocities (U,V,W) and acelerat-
# ion (Ax, Ay, Az) of every particle at the next step of time.
# ------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------- Created by Maveryck Garzon

# NOTE: you need a folder called Data where the information will be store as csv files.



# Initial conditions
#	Sun	Merc.	Venus	Earth	Mars	Jup.	Satur.	Uran.	Nept.	Plut.
Rp = [	0.0	0.387	0.723	1	1.524	5.203	9.54	19.2	30.1	39.4]	# Orbit's radio
Up = [	0.0	-78	-70	-59	-47	-26	-18.8	-13.2	-10.5	-9]	# Speed
Mp = [	331e3	0.06	0.82	1	0.1	318	95.1	14.6	17.2	0.002]	# Mass (relative to Earth's)
Ip = [	0.4	0.2	0.2	0.2	0.2	0.3	0.3	0.3	0.2	0.1]	# I use this for plotting

Obj = 200 # Number of objects in the system. This includes the Sun, 9 planets and many asteroids at the asteroids belt

X = zeros(SharedArray{Float64}(Obj));
Y = zeros(SharedArray{Float64}(Obj));
Z = zeros(SharedArray{Float64}(Obj));
U = zeros(SharedArray{Float64}(Obj));
V = zeros(SharedArray{Float64}(Obj));
W = zeros(SharedArray{Float64}(Obj));
M = zeros(SharedArray{Float64}(Obj));
I = zeros(SharedArray{Float64}(Obj));

G = 1.0e-02 		 # Notice all parameters are relative to earth's ones. So G must be fitted.
Kp= 4.0*pi^2.0/(G*Mp[1]) # Kepler's constant that dependens mainly on the mass of the Sun
			 # T^2=K*R^3 Kepler's law asumes a circular orbit which is accurate enough

for i in 1:10
	   R = Rp[i]
#	each planet is located at a random position along its orbit
	X[i] = rand(0.0:0.01:R)*rand(-1.0:2.0:1.0)
	Y[i] = sqrt(R^2.0 - X[i]^2.0)*rand(-1.0:2.0:1.0)
	Z[i] = 0.0

#	 vel = abs(Up[i])		# The speed can be predetermined
	 vel = (2.0*pi)/sqrt(Kp*R)	# Speed calculated with the Kepler's constant
	U[i] = -vel*Y[i]/R
	V[i] =  vel*X[i]/R
	W[i] = 0.0
	M[i] = Mp[i]
	I[i] = Ip[i]
end

for i in 11:Obj 			# All these objects are asteroids at the asteroids belt 
	   R = rand(2.6:0.001:2.73)	# between Mars's and Jupiter's orbits
	X[i] = rand(0.0:0.01:R)*rand(-1.0:2.0:1.0)
	Y[i] = sqrt(R^2.0 - X[i]^2.0)*rand(-1.0:2.0:1.0)
	Z[i] = 0.0

	 vel = (2.0*pi)/sqrt(Kp*R)	# Speed calculated with the Kepler's constant
	U[i] = -vel*Y[i]/R
	V[i] =  vel*X[i]/R
	W[i] = 0.0
	M[i] = rand(0.01:0.0001:0.04)	# A random mass is asigned, it must be very small so Jupiter doesn't affect them
	I[i] = 0.05			# I use this number to tell paraview what color asign to each object 
					# 	so it's not really relevant in the calculations
end

# The Sun can be moving in every direction:
U[1] = 0.0
V[1] = 0.0
W[1] = 0.01

Xc = copy(X)
Yc = copy(Y)
Zc = copy(Z)
Uc = copy(U)
Vc = copy(V)
Wc = copy(W)

# The initial conditions are stored in csv file, which can be easily viewed with Paraview
	open("Data/Data0.csv","w") do Data;
	   write(Data, "x_coor, y_coor, z_coor, Ind \n");
	end
	  for i in 1:Obj
	    open("Data/Data0.csv","a") do Data;
	    vel = sqrt(U[i]^2+V[i]^2)
	    write(Data, "$(Xc[i]), $(Yc[i]), $(Zc[i]), $(I[i])\n");
	    end;
	  end;

# ------------------------------------------------------------------------------------------------
# ------------------------------ Here comes the math part ----------------------------------------
# ------------------------------------------------------------------------------------------------

# Initial Parameters:
 N = 500;	# Number of time steps
 t = 0.0;	# Initial time
dt = 0.0001;	# Time step

for k in 1:N;
  t = t + dt

	for i in 1:Obj
		Ax = 0.0;	# Acceleration in X direction
		Ay = 0.0;	# Acceleration in Y direction
		Az = 0.0;	# Acceleration in Z direction

	     for j in 1:Obj
		if i==j
		  # Nothing happens
		else
		 r = sqrt(  (X[j]-X[i])^(2) + (Y[j]-Y[i])^(2) + (Z[j]-Z[i])^(2) )

		     if r < 0.3 # This is to avoid r=0 and a numerical error
		        r = 0.3
		     end

		   A = G*M[j]/(r^2)	     # Newton's Gravitational law, where:
		  Ax = Ax + A*(X[j]-X[i])/r  # 	r is the distance between 2 objects
		  Ay = Ay + A*(Y[j]-Y[i])/r  #	A is the total acceleration
		  Az = Az + A*(Z[j]-Z[i])/r  #	G is the gravitational constant
		end
	     end

		Uc[i] = U[i] + Ax*dt;	# Speed X in the next time step is calculated
		Vc[i] = V[i] + Ay*dt;	# Speed Y in the next time step is calculated
		Wc[i] = W[i] + Az*dt;	# Speed Z in the next time step is calculated
		Xc[i] = X[i] + Uc[i]*dt;	# Position X in the next time step
		Yc[i] = Y[i] + Vc[i]*dt;	# Position Y in the next time step
		Zc[i] = Z[i] + Wc[i]*dt;	# Position Z in the next time step
	end

# All positions and velocities are actualized so the next timestep can start
X = Xc
Y = Yc
Z = Zc
U = Uc
V = Vc
W = Wc

# Every Data file contains the position of that body at every 10 time-steps
	if k%10==0

	println("Progress: ", 100*k/N, "%");
	open("Data/Data$k.csv","w") do Data;
	   write(Data, "x_coor, y_coor, z_coor, Ind \n");
	end
	  for i in 1:Obj
	    open("Data/Data$k.csv","a") do Data;
	    write(Data, "$(Xc[i]), $(Yc[i]), $(Zc[i]), $(I[i])\n");
	    end;
	  end;
	end

end;

