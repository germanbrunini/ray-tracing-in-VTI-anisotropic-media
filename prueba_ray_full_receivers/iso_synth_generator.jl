function iso_synth_generator(dt :: Float64,
    ns    :: Int64,
    nr    :: Int64,
    R     :: Array{Float64,1},
    t0s   :: Float64,
    S     :: Array{Float64,1},
    f0    :: Float64,
    vel   :: Float64, # equivalent isotropic phase velocity
    vᵣₑ   :: Float64, # anisotropic phase velocity on receiver
    vₛₒ   :: Float64, # anisotropic phase velocity on source
    ρᵣₑ   :: Float64, # density on receiver
    ρₛₒ   :: Float64, # density on source
    M     :: Array{Float64,2},
    phase :: String)::Array{Float64,2};


# Microseismic Monitoring. Vladimir Gretchka. 2017. Pag.109.
# Altoughtthe results expresed in equation 3.80 is obtained for homogeneos
# isotropic media, it turns out to be valid universally, that is, for arbitrary
# anisotropy and heterogeneity (e.g. Cerveny 2001, Chapter 2; also chapetr4)

    # Change units
    R   =   R/1e3;   # receiver coordinates
    S   =   S/1e3;   # Source coordinates
    vel = vel/1e3;   # P-velocity
	vᵣₑ = vᵣₑ/1e3;
	vₛₒ = vₛₒ/1e3;
    ρᵣₑ = ρᵣₑ/1e3;
    ρₛₒ = ρₛₒ/1e3;


    vₛₒ²  = vₛₒ^2               ; # eq. 4.93 Aki-Rich & eq.4.108 Gretchka Book
	rw    = Ricker(dt=dt, f0=f0);
    nrw   = length(rw)          ;
    rw    = reshape(rw,(1,nrw)) ;
    nf    = 4 * nextpow(2,ns)   ;
    dw    = 2.0 * pi/(nf*dt)    ;
    rwpad = hcat(rw,zeros(nf-nrw)');
    RW    = fft(rwpad)             ;

    U     = zeros(Float64, (ns, 3*nr));

    # U_spherical = zeros(Float64, (ns, 3*nr));
	# # changing coordinates from (x,y,z), to (r,θ,ϕ)
    # x_to_sph = [sind(θ)*cosd(ϕ), -sind(ϕ),  cosd(θ)*cosd(ϕ)];
    # y_to_sph = [sind(θ)*sind(ϕ),  cosd(ϕ),  cosd(θ)*sind(ϕ)];
    # z_to_sph = [cosd(θ)        ,   0.0   , -sind(θ)        ];
    if phase=="pp"
        @inbounds for j = 1:nr

            SR = -(S - R[:,j]);
            r  = norm(SR,2);
            v  = SR./r;
			# eq. 4.93 Aki-Richards and eq. 4.108 Gretchka Book
			cp = 1.0/(4.0 * pi * sqrt(ρₛₒ*ρᵣₑ * vₛₒ*vᵣₑ) * r * vₛₒ²);
            #cp = 1.0/(4.0*pi*rho*vel_3*r); #homogeneous case
			ap = v*cp*v'*M*v;
            tp = r/vel;

            Uj  = zeros(Complex{Float64}, nf, 3);
            Ujp = zeros(Complex{Float64}, nf, 3);
            nf2 = floor(Int, nf/2);
            for k = 2:nf2
                w    = (k-1)*dw;
                imw  = im*w;
                imwR = imw*RW[k];
                tsp  = t0s + tp;
                Ujp[k, :]     = ap*exp(-imw*tsp);
                Uj[k, :]      = imwR * (Ujp[k,:]);
                Uj[nf+2-k, :] = conj(Uj[k, :]);
            end
            uj = real(ifft(Uj, 1))[1:ns, :];
            jx = 3*(j-1) + 1;
            jy = jx + 1     ;
            jz = jx + 2     ;
            U[:, jx:jz] = uj;

            # Ux_spherical = rotation1(U[:,jx],x_to_sph,ns)
            # Uy_spherical = rotation1(U[:,jy],y_to_sph,ns)
            # Uz_spherical = rotation1(U[:,jz],z_to_sph,ns)
            # U_spherical  = Ux_spherical + Uy_spherical + Uz_spherical;
			# writedlm("1_pp_change.txt",U)
			# writedlm("1_pp_spheri.txt",U_spherical);
    	end
		return U
        # return U_spherical[:,1];
    end

    if phase=="sh"
        @inbounds for j = 1:nr
            SR = -(S - R[:,j]);
            r  = norm(SR,2);
            v  = SR./r;
			# eq. 4.93 Aki-Richards and eq. 4.108 Gretchka Book
			cs = 1.0/(4.0 * pi * sqrt(ρₛₒ*ρᵣₑ * vₛₒ*vᵣₑ) * r * vₛₒ²);
            #cs = 1.0/(4.0*pi*rho*vel_3*r);
            as = cs*(M*v - v*v'*M*v);
            ts = r/vel;
            Uj  = zeros(Complex{Float64}, nf, 3);
            Ujs = zeros(Complex{Float64}, nf, 3);
            nf2 = floor(Int, nf/2);
            for k = 2:nf2
                w    = (k-1)*dw;
                imw  = im*w;
                imwR = imw*RW[k];
                tss  = t0s + ts;
                Ujs[k, :] = as*exp(-imw*tss);
                Uj[k, :]  = imwR * (Ujs[k,:]);
                Uj[nf+2-k, :] = conj(Uj[k, :]);
            end
            uj = real(ifft(Uj, 1))[1:ns, :];
            jx = 3*(j-1) + 1;
			jy = jx + 1     ;
            jz = jx + 2;
            U[:, jx:jz] = uj;

			# Ux_spherical = rotation(U[:,jx],x_to_sph,ns)
			# Uy_spherical = rotation(U[:,jy],y_to_sph,ns)
			# Uz_spherical = rotation(U[:,jz],z_to_sph,ns)
			# U_spherical  = Ux_spherical + Uy_spherical + Uz_spherical;
			# writedlm("1_sh_change.txt",U);
			# writedlm("1_sh_spheri.txt",U_spherical);
		end
		# return U_spherical[:,2];
		return U;
    end

    if phase=="sv"
        @inbounds for j = 1:nr
            SR = -(S - R[:,j]);
            r  = norm(SR,2);
            v  = SR./r;
			# eq. 4.93 Aki-Richards and eq. 4.108 Gretchka Book
			cs = 1.0/(4.0 * pi * sqrt(ρₛₒ*ρᵣₑ * vₛₒ*vᵣₑ) * r * vₛₒ²);
            #cs = 1.0/(4.0*pi*rho*vel_3*r);
            as = cs*(M*v - v*v'*M*v);
            ts = r/vel;
            Uj  = zeros(Complex{Float64}, nf, 3);
            Ujs = zeros(Complex{Float64}, nf, 3);
            nf2 = floor(Int, nf/2);
            for k = 2:nf2
                w    = (k-1)*dw;
                imw  = im*w;
                imwR = imw*RW[k];
                tss  = t0s + ts;
                Ujs[k, :] = as*exp(-imw*tss);
                Uj[k, :]  = imwR * (Ujs[k,:]);
                Uj[nf+2-k, :] = conj(Uj[k, :]);
            end
            uj = real(ifft(Uj, 1))[1:ns, :];
            jx = 3*(j-1) + 1;
			jy = jx + 1     ;
            jz = jx + 2;
            U[:, jx:jz] = uj;

			# Ux_spherical = rotation(U[:,jx],x_to_sph,ns)
            # Uy_spherical = rotation(U[:,jy],y_to_sph,ns)
            # Uz_spherical = rotation(U[:,jz],z_to_sph,ns)
            # U_spherical  = Ux_spherical + Uy_spherical + Uz_spherical;
			# writedlm("1_sv_change.txt",U);
			# writedlm("1_sv_spheri.txt",U_spherical);
    	end
        # return U_spherical[:,3];
        return U;
    end

end


# function rotation1(U,card_to_sph,ns)
#     U_shp = zeros(Float64,(ns,3))
#     for i=1:ns
#         U_shp[i,:] = U[i]*card_to_sph;
#     end
# return U_shp;
# end


"""
Ricker(; <keyword arguments>)
Create a Ricker wavelet.
# Keyword arguments
* `dt::Real=0.002`: sampling interval in secs.
* `f0::Real=20.0`: central frequency in Hz.
# Examples
```julia
julia> w = Ricker(); plot(w);
julia> w = Ricker(dt=0.004, f0=20); plot(w);
```
# Reference
Sheriff, Robert, 2002, Encyclopedic Dictionary of Applied Geophysics, fourth
ed.: Society of Exploration Geophysicists. Geophysical Reference Series No. 13.
"""

function Ricker(;dt::Real=0.002, f0::Real=20.0)::Array{Float64,1}
    nw = 2.0/(f0*dt);
    nc = floor(Int, nw/2);
    t  = dt*collect(-nc:1:nc);
    p = [f0];

    # writedlm("dRick",dRick.(t,p));
    # writedlm("Rick",Rick.(t,p));
    return Rick.(t,p);
end

# Ricker
Rick(t,p)  = (1.0-2.0*(pi*p[1]*t)^2)*exp(-(pi*p[1]*t)^2);
