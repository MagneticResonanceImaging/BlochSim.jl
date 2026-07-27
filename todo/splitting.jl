#=
splitting.jl
Bloch simulation with relaxation effects using symmetric operator splitting
06.08.2019
Christina Graf c.graf@tugraz.at
=#

function bloch_symmetric_strang_splitting(u, v, w, d, T1_k, T2_k, relax_k)
	xdis = d["xdis"]
	dt = d["dt"]
	M0c = d["M0c"]
	M0 = d["M0"]
	Nx = Int(d["Nx"])
	gamma = d["gamma"]

	Nu = size(u, 1)
	for j = 1:Nu
		if u[j] == 0
			u[j] = 10^-14
		end
		if v[j] == 0
			v[j] = 10^-14
		end
	end
	#u[u==0]=10^-14
	#v[v==0]=10^-14

	gadt = gamma * dt / 2
	B = (gadt * transpose(u - 1im * v) ) .* ones(Nx, 1)
	K = gadt * xdis * w'
	phi = -(abs.(B).^2 + K.^2).^(1 / 2)
	D = [exp(-1 / T2_k * relax_k * dt) 0 0;0 exp(-1 / T2_k * relax_k * dt) 0;0 0 exp(-1 / T1_k * relax_k * dt)]
	b = [0;0;M0c] - [0;0;M0c * exp(-1 / T1_k * relax_k * dt)]

	M = zeros(3, Nx, Nu + 1)

	for z = 1:Nx
		Mn = M0[:,z]
		M[:,:,1] = M0
		for n = 1:Nu # time loop
			phi_j = phi[z,n]
			B1x = real(B[z,n])
			B1y = imag(B[z,n])
			Gx = K[z,n]
			n_j = [B1x; B1y; Gx] ./ abs(phi_j)
			for aux = 1:3
				if isnan(n_j[aux])
					n_j[aux] = 0
				end
			end
			n1 = n_j[1]
			n2 = n_j[2]
			n3 = n_j[3]
			cs = cos(phi_j)
			si = sin(phi_j)

			Bd = [n1^2 * (1 - cs) + cs n1 * n2 * (1 - cs) - n3 * si n1 * n3 * (1 - cs) + n2 * si;
				n2 * n1 * (1 - cs) + n3 * si n2^2 * (1 - cs) + cs n2 * n3 * (1 - cs) - n1 * si;
				n3 * n1 * (1 - cs) - n2 * si n3 * n2 * (1 - cs) + n1 * si n3^2 * (1 - cs) + cs]


			Mrot = Bd * Mn

			Mt = D * Mrot + b

			Mn = Bd * Mt
			M[:,z,n + 1] = Mn
		end
	end
	return M
end


function bloch_symmetric_strang_splitting_vectorized(u, v, w, d, T1_k, T2_k, relax_k)
	xdis = d["xdis"]
	dt = d["dt"]
	M0c = d["M0c"]
	M0 = d["M0"]
	Nx = Int(d["Nx"])
	gamma = d["gamma"]

	Nu = size(u, 1)
	for j = 1:Nu
		if u[j] == 0
			u[j] = 10^-14
		end
		if v[j] == 0
			v[j] = 10^-14
		end
	end


	gadt = gamma * dt / 2
	B = (gadt * transpose(u - 1im * v)) .* ones(Nx, 1)
	K = gadt * xdis * w'
	phi = -(abs.(B).^2 + K.^2).^(1 / 2)
	D = [exp(-1 / T2_k * relax_k * dt) 0 0;0 exp(-1 / T2_k * relax_k * dt) 0;0 0 exp(-1 / T1_k * relax_k * dt)]
	b = [0;0;M0c] - [0;0;M0c * exp(-1 / T1_k * relax_k * dt)]

	cs = cos.(phi)
	si = sin.(phi)
	n1 = real(B) ./ abs.(phi)
	n2 = imag(B) ./ abs.(phi)
	n3 = K ./ abs.(phi)
	n1[isnan.(n1)] .= 1
	n2[isnan.(n2)] .= 0
	n3[isnan.(n3)] .= 0
	aux = ones(Nx, Nu)
	Bd1 = n1 .* n1 .* (aux - cs) .+ cs
	Bd2 = n1 .* n2 .* (aux - cs) .- n3 .* si
	Bd3 = n1 .* n3 .* (aux - cs) .+ n2 .* si
	Bd4 = n2 .* n1 .* (aux - cs) .+ n3 .* si
	Bd5 = n2 .* n2 .* (aux - cs) .+ cs
	Bd6 = n2 .* n3 .* (aux - cs) .- n1 .* si
	Bd7 = n3 .* n1 .* (aux - cs) .- n2 .* si
	Bd8 = n3 .* n2 .* (aux - cs) .+ n1 .* si
	Bd9 = n3 .* n3 .* (aux - cs) .+ cs

	M = zeros(3, Nx, Nu + 1)
	M[:,:,1] = M0
	Mt = M0
	Mrot = zeros(3, Nx)
	for n = 1:Nu
		Mrot[1,:] = Bd1[:,n]' .* Mt[1,:]' .+ Bd2[:,n]' .* Mt[2,:]' .+ Bd3[:,n]' .* Mt[3,:]'
		Mrot[2,:] = Bd4[:,n]' .* Mt[1,:]' .+ Bd5[:,n]' .* Mt[2,:]' .+ Bd6[:,n]' .* Mt[3,:]'
		Mrot[3,:] = Bd7[:,n]' .* Mt[1,:]' .+ Bd8[:,n]' .* Mt[2,:]' .+ Bd9[:,n]' .* Mt[3,:]'


		Mt = D * Mrot + b .* ones(3, Nx)

		Mrot[1,:] = Bd1[:,n]' .* Mt[1,:]' .+ Bd2[:,n]' .* Mt[2,:]' .+ Bd3[:,n]' .* Mt[3,:]'
		Mrot[2,:] = Bd4[:,n]' .* Mt[1,:]' .+ Bd5[:,n]' .* Mt[2,:]' .+ Bd6[:,n]' .* Mt[3,:]'
		Mrot[3,:] = Bd7[:,n]' .* Mt[1,:]' .+ Bd8[:,n]' .* Mt[2,:]' .+ Bd9[:,n]' .* Mt[3,:]'


		Mt = Mrot;
		M[:,:,n + 1] = Mrot
	end
	return M
end

using Pkg
using LinearAlgebra
using Plots

Pkg.add("MAT")
using MAT


vars = matread("Beispiel1.mat")
u = vars["u"]
v = vars["v"]
w = vars["w"]
d = vars["d"]


d["B1c"] = 57.8704 # for peak B1 of 12.5 muT
d["dt"] = (8.6409e-7) * 2 * 10^3 # for peak B1 of 12.5 muT and pi/2 in ms
B1c = d["B1c"]
u = u * B1c * 10^-3 # in mT
v = v * B1c * 10^-3 # in mT
d["G3"] = 18 / 2 # mT/m, fwhm 2mm
w = w * d["G3"] # in mT/m
d["xdis"] = range(-0.005, length = 101, stop = 0.005) # in m

d["T1"] = [10^-9 1331 400 832 1420]
d["T2"] = [10^-9 110 5 79.6 31.7]
d["relax"] = [0 1 1 1 1]
example = ["no relax", "grey matter", "tendons", "white matter", "muscle"]

# Init check (valuable in the future where people are using their data and not the test case above)
#Do RF and Gs exist?
if ~@isdefined u #uly way of implementing, however && does not work??
	if ~@isdefined v
		if ~@isdefined w
			println("No RF and Gs found, assumed to be zero with length 100")
			u=zeros(100,1)
			v=zeros(100,1)
			w=zeros(100,1)
		end
	end
end

# same length of RF and Gs?
if length(u)==length(v) && length(u)==length(w)
	# continue
else
	println("input dimensions not correct")
	exit()
end

# dt defined?
if ~haskey(d, :"dt")
	println("no time step length defined, assume 1")
	d["dt"]=1
end

# xdis defined?
if ~haskey(d, :"xdis")
	println("xdis not defined, assume (-0.005, 0.005)m")
	d["xdis"] = range(-0.005, length = 101, stop = 0.005) # in m
end
# T1 and T2 defined?
if ~haskey(d, :"T1") && haskey(d, :"T2")
	println("T1 and T2 not defined, assumed to be infinity")
	d["T1"] = 10^9
	d["T2"] = 10^9
	d["relax"]=0
end
if haskey(d, :"T1") && ~haskey(d, :"T2")
	println("T2 not defined, assumed to be infinity")
	d["T2"] = 10^9
	d["relax"]=0
end
if haskey(d, :"T2") && ~haskey(d, :"T1")
	println("T1 not defined, assumed to be infinity")
	d["T1"] = 10^9
	d["relax"]=0
end


k = 2 # k from 1 to 5 corresponding to the 5 examples
T1=d["T1"]
T2=d["T2"]
relax=d["relax"]
T1_k = T1[k]
T2_k = T2[k]
relax_k = relax[k]
M = bloch_symmetric_strang_splitting_vectorized(u, v, w, d, T1_k, T2_k, relax_k)
plot(d["xdis"] * 1000, [M[1,:,end] M[2,:,end] M[3,:,end]], title = example[k], label = ["x" "y" "z"], lw = 3)
