
function aniso_ray_tracing(
   z::Array{Float64,1},
   vp::Array{Float64,1},
   vs::Array{Float64,1},
   rho::Array{Float64,1},
   an::Array{Float64,2},
   s::Array{Float64,2},
   r::Array{Float64,2},
   geom)

   dt  = geom[1]
   ns  = geom[2]
   t0s = geom[3]
   f0  = geom[4]
   M   = geom[5]

   # in  z,vp,vs,an,s,r
   # out tp_out,tsh_out,tsv_out
   n_s = size(s,1);
   n_r = size(r,1);
   n_z = length(z);
   e       = zeros(Float64,n_z + 1)    ; # modelos de espesores
   dp_out  = zeros(Float64,n_s,n_r,n_z);
   dsh_out = zeros(Float64,n_s,n_r,n_z);
   dsv_out = zeros(Float64,n_s,n_r,n_z);
   n_int   = zeros(Int64,n_s,n_r)      ;
   z_ind   = zeros(Int64,n_z)          ; # indice de las discontinuidades

   tp_out  = zeros(Float64,n_s,n_r);
   tsh_out = zeros(Float64,n_s,n_r);
   tsv_out = zeros(Float64,n_s,n_r);

   cp_out  = zeros(Float64,n_s,n_r,n_z,3);
   csh_out = zeros(Float64,n_s,n_r,n_z,3);
   csv_out = zeros(Float64,n_s,n_r,n_z,3);

   C = Thomsen_to_Cij(an,vp,vs,rho);

   signal  = Matrix{Matrix{Float64}}(undef, n_s, n_r);
   sigpp   = Matrix{Matrix{Float64}}(undef, n_s, n_r);
   sigsh   = Matrix{Matrix{Float64}}(undef, n_s, n_r);
   sigsv   = Matrix{Matrix{Float64}}(undef, n_s, n_r);
   amps    = Array{Array{Float64}}(undef, n_s, n_r);
   Upp_xyz = zeros(Float64,ns,3);
   Ush_xyz = zeros(Float64,ns,3);
   Usv_xyz = zeros(Float64,ns,3);

   # para cada par fuente-receptor busco los tiempos de viaje
   for i = 1:n_s
      for j = 1:n_r

         Upp_xyz = zeros(Float64,ns,3);
         Ush_xyz = zeros(Float64,ns,3);
         Usv_xyz = zeros(Float64,ns,3);

         source = s[i,:];
         receiv = r[j,:];

         Z_s  = source[3];
         Z_r  = receiv[3];
         direct = direction(Z_s,Z_r); # string with ray direct. ("up" or "down")

         maxZ = max(Z_s,Z_r);
         minZ = min(Z_s,Z_r);
         # averiguo n_int[i,j]: nro de discontinuidades entre cada fuente y cada receptor
         n_int[i,j] = count(k->(minZ < k < maxZ),z);
         n_ij       = n_int[i,j]; # store value to be used any moment

         # si n_int[i,j]=0 la fuente y el receptor están en la misma capa
         # y saco el tiempo de viaje directamente..
         # si n_int[i,j] no es cero hay que estimar los rayos iterativamente

         if (n_ij == 0)
            # me fijo cual es la capa donde están la fuente y el receptor
            disc     = collect(1:n_z); # numero de discontinuidades
            # de todas las discontinuidades que estan mas profundas que maxZ,
            # me quedo con la primera (la menor)
            z_ind[1] = disc[findall(x->(x>maxZ),z)[1]]
            tp_out[i,j] ,
            tsv_out[i,j],
            tsh_out[i,j] = straight_ray(
            vp[z_ind[1]],   # la velocidad P de dicha capa
            vs[z_ind[1]],   # la velocidad S de dicha capa
            an[z_ind[1],:], # la anisotropia de dicha capa
            s[i,:],         # la fuente i
            r[j,:]);        # el receptor j

            pp_ray,sh_ray,sv_ray =  trayectories(source,receiv)          ;
            pp_ang,sh_ang,sv_ang =        angles(source,receiv)          ;

            pp_gsf,sh_gsf,sv_gsf     = gsfactor(pp_ray,sh_ray,sv_ray)   ;
            pp_diramp,pp_dirphase,ID = dir_amplitudes_an(C,pp_ang,sh_ang,sv_ang,receiv,source,z,rho,"pp",direct);
            sh_diramp,sh_dirphase,ID = dir_amplitudes_an(C,pp_ang,sh_ang,sv_ang,receiv,source,z,rho,"sh",direct);
            sv_diramp,sv_dirphase,ID = dir_amplitudes_an(C,pp_ang,sh_ang,sv_ang,receiv,source,z,rho,"sv",direct);

            pp_noramp = coef_normalization(rho,pp_ray,pp_ang,vp,vs,an,ID,"pp"); # eq. 4.109 "Gretcka Book".
            sh_noramp = coef_normalization(rho,sh_ray,sh_ang,vp,vs,an,ID,"sh"); # eq. 4.109 "Gretcka Book".
            sv_noramp = coef_normalization(rho,sv_ray,sv_ang,vp,vs,an,ID,"sv"); # eq. 4.109 "Gretcka Book".

            #pp_diramp,pp_dirphase,pp_traj,ID = dir_amplitudes(pp_ang,receiv,source,z,vp,vs,rho,1);
            #sh_diramp,sh_dirphase,sh_traj,ID = dir_amplitudes(sh_ang,receiv,source,z,vp,vs,rho,7);
            #sv_diramp,sv_dirphase,sv_traj,ID = dir_amplitudes(sv_ang,receiv,source,z,vp,vs,rho,2);
            T_ppdiramp = prod((pp_diramp.*pp_noramp)) * pp_gsf;
            T_shdiramp = prod((sh_diramp.*sh_noramp)) * sh_gsf;
            T_svdiramp = prod((sv_diramp.*sv_noramp)) * sv_gsf;

            # -----------------------------------------------------------------
            # Calculate wavelets
            # -----------------------------------------------------------------

            # Theory:
            # I   - I will consider a source in the last refraction point (LRP),
            # therefore, the direction of incidence on the receiver point (RP)
            # will be the correct.
            # II  - The wavelet needs to appear at calculated times, so I will
            # calculate the equivalent isot. vel. between LRP and RP.
            #     Vᵢₛₒ = Δ(LRP-RP)/Δt.
            # III - I will run isot. synthetic generator with this geometry.
            #    a) signal will appear at correct times.
            #    b) Source mechanism remains unaltered in all trayectory.
            #    c) X, Y and Z amplitud relations (normalized) are dominated by :
            #             Rᵢⱼₖ tensor  (dominated by dir of incidence, see "I -").
            #             Mᵢⱼ  tensor. (recall, signal is normalized, M₀= 1.0).

            # cpp_c: same as cpp but for currect rec_sou signal. Row dimention varies
            # depending on i-source and j-receiver i.e. number of intersections

            LRP_pp       = source  ; # last refraction point
            LRP_sh       = source  ; # last refraction point
            LRP_sv       = source  ; # last refraction point
            RP           = receiv  ; # receiver point
            LRP_RP_pp    = RP-LRP_pp         ; # vector from LRP_pp to RP.
            LRP_RP_sh    = RP-LRP_sh         ; # vector from LRP_sh to RP.
            LRP_RP_sv    = RP-LRP_sv         ; # vector from LRP_sv to RP.
            Δ_LRP_RP_pp  = norm(LRP_RP_pp)   ; # distance from LRP_pp to RP.
            Δ_LRP_RP_sh  = norm(LRP_RP_sh)   ; # distance from LRP_sh to RP.
            Δ_LRP_RP_sv  = norm(LRP_RP_sv)   ; # distance from LRP_sv to RP.
            t_pp =  tp_out[i,j]              ; # p-phase travel time
            t_sh = tsh_out[i,j]              ; # sh-phase travel time
            t_sv = tsv_out[i,j]              ; # sv-phase travel time
            vpp_i = Δ_LRP_RP_pp/t_pp    ; # equivalent vp isot. vel
            vsh_i = Δ_LRP_RP_sh/t_sh    ; # equivalent vsh isot. vel
            vsv_i = Δ_LRP_RP_sv/t_sv    ; # equivalent vsv isot. vel

            # figure 4.20 aki-richards. Figure 4.10 aki-richards
            α_pp = take_off_angle(pp_ang)   ; # t-off ang. Where the ray comes from
            α_sh = take_off_angle(sh_ang)   ; # t-off ang.
            α_sv = take_off_angle(sv_ang)   ; # t-off ang.

            β_pp = take_in_angle(pp_ang)   ; # t-off ang. Where the ray comes from
            β_sh = take_in_angle(sh_ang)   ; # t-off ang.
            β_sv = take_in_angle(sv_ang)   ; # t-off ang.

            ϕ_pp = take_off_azimuth(LRP_RP_pp) ; # t_off azimuth
            ϕ_sh = take_off_azimuth(LRP_RP_sh) ; # t_off azimuth
            ϕ_sv = take_off_azimuth(LRP_RP_sv) ; # t_off azimuth
            # eq.4.88 Aki-Richards: Projections for P- SH- and SV- rays
            pd_pp,ϕd_pp = wave_directionals(α_pp,ϕ_pp);
            pd_sh,ϕd_sh = wave_directionals(α_sh,ϕ_sv);
            pd_sv,ϕd_sv = wave_directionals(α_sh,ϕ_sv);

            idᵣ = ID[1]; # id of receiver layer
            idₛ = ID[2]; # id of source layer

            ρᵣ   = rho[idᵣ];
            ρₛ   = rho[idₛ];
            vppᵣ = an_vel(deg2rad(α_pp),vp[idᵣ],vs[idᵣ],an[idᵣ,:],"vp");
            vshᵣ = an_vel(deg2rad(α_sh),vp[idᵣ],vs[idᵣ],an[idᵣ,:],"sh");
            vsvᵣ = an_vel(deg2rad(α_sv),vp[idᵣ],vs[idᵣ],an[idᵣ,:],"sv");
            vppₛ = an_vel(deg2rad(β_pp),vp[idₛ],vs[idₛ],an[idₛ,:],"vp");
            vshₛ = an_vel(deg2rad(β_sh),vp[idₛ],vs[idₛ],an[idₛ,:],"sh");
            vsvₛ = an_vel(deg2rad(β_sv),vp[idₛ],vs[idₛ],an[idₛ,:],"sv");

            Upp_xyz = iso_synth_generator(dt,ns,1,RP,t0s,LRP_pp,f0,vpp_i,vppᵣ,vppₛ,ρᵣ,ρₛ,M,"pp");
            Ush_xyz = iso_synth_generator(dt,ns,1,RP,t0s,LRP_sh,f0,vsh_i,vshᵣ,vshₛ,ρᵣ,ρₛ,M,"sh");
            Usv_xyz = iso_synth_generator(dt,ns,1,RP,t0s,LRP_sv,f0,vsv_i,vsvᵣ,vsvₛ,ρᵣ,ρₛ,M,"sv");

            # normfact = normdata(Upp_xyz + Ush_xyz + Usv_xyz)[2];
            # Upp_xyz = Upp_xyz/normfact;
            # Ush_xyz = Ush_xyz/normfact;
            # Usv_xyz = Usv_xyz/normfact;

            # theory:
            # for each ray (say pp(x,y,z), sh(x,y,z) or sv(x,y,z)), I can
            # rotate to a spherical coordinate system with axis on radial direction (r)
            # on the direction of the ray where only P lives, and a perpendicular plane,
            # where only S lives, with directions (i:inc, ϕ:azimuth).
            # If I am rotating pp(x,y,z) to pp(r,ϕ,i), I can discard pp(i), pp(ϕ), as their
            # enery is zero.
            # If I am rotating sh(x,y,z) to sh(r,ϕ,i), I can discard sh(r), sh(i), as their
            # enery is zero.
            # If I am rotating sv(x,y,z) to sv(r,ϕ,i), I can discard sv(r), sv(ϕ), as their
            # enery is zero.

            # rotate from xyz to spherical acocrding to their ray direction
            # extract radial mechanism (where pp lives)
            Upp_r = rotation(Upp_xyz,α_pp,ϕ_pp,direct,"xyz2rϕi")[:,1];
            # extract azimuthal mechanism (where sh lives)
            Ush_ϕ = rotation(Ush_xyz,α_sh,ϕ_sh,direct,"xyz2rϕi")[:,2];
            # extract inclination mechanism (where sv lives)
            Usv_i = rotation(Usv_xyz,α_sv,ϕ_sv,direct,"xyz2rϕi")[:,3];

            # input their calculated amplitudes.
            Upp_amp,Ush_amp,Usv_amp = rescale(Upp_r,Ush_ϕ,Usv_i,T_ppdiramp,T_shdiramp,T_svdiramp);

            # rotate back
            Upp_xyz = rotation(Upp_amp,α_pp,ϕ_pp,direct,"rϕi2xyz");
            # extract azimuthal mechanism (where sh lives)
            Ush_xyz = rotation(Ush_amp,α_sh,ϕ_sh,direct,"rϕi2xyz");
            # extract inclination mechanism (where sv lives)
            Usv_xyz = rotation(Usv_amp,α_sv,ϕ_sv,direct,"rϕi2xyz");

            final_signal = Upp_xyz + Ush_xyz + Usv_xyz;
         else

            # !==============================================================
            # ! el rayo refractado anisotropo se calcula en dos partes:
            # !
            # ! 1) se calcula el caso de modelo isotropo
            # ! 2) usando el resultado de 1) como modelo incial se calcula
            # !    el rayo anisotropo
            # !
            # !==============================================================

            # indexo las discontinuidades entre la fuente y el receptor.
            disc          = collect(1:n_z);
            zcond         = findall(x->(minZ < x < maxZ),z);
            z_ind[1:n_ij] = disc[zcond];

            # defino el modelo de  espesores
            e[1]      = z[z_ind[1]]      - minZ              ;
            e[2:n_ij] = z[z_ind[2:n_ij]] - z[z_ind[1:n_ij-1]];
            e[n_ij+1] = maxZ             - z[z_ind[n_ij]]    ;

            # calculo los tiempos de arribo para el rayo refrectado,
            # tambien se sacan las distancias desde la fuente hacia
            # cada una de las intersecciones, en la linea que une
            # la fuente y el receptor
            dp_out[i,j,1:n_ij] ,
            dsv_out[i,j,1:n_ij],
            dsh_out[i,j,1:n_ij],
            tp_out[i,j]        , # PP trav.time from S to R
            tsv_out[i,j]       , # SV trav.time from S to R
            tsh_out[i,j]       , # SH trav.time from S to R
            tpp_mem            , # Δt[i] for P : sum(Δt[i]) = tp_out
            tsv_mem            , # Δt[i] for SV: sum(Δt[i]) = tsv_out
            tsh_mem = refracted_ray(
            e[1:n_ij+1]                 , # modelo completo de espesores
            vp[z_ind[1]:z_ind[n_ij]+1]  , # modelo de velocidades P
            vs[z_ind[1]:z_ind[n_ij]+1]  , # modelo de velocidades S
            an[z_ind[1]:z_ind[n_ij]+1,:], # Modelo de anisotropias
            s[i,:]                      , # coordenadas de la fuente i
            r[j,:]                      , # coordenadas del receptor j
            dp_out[i,j,1:n_ij] ,
            dsv_out[i,j,1:n_ij],
            dsh_out[i,j,1:n_ij])

            # -----------------------------------------------------------------
            # guardo la coord xyz de cada segmento del rayo
            # -----------------------------------------------------------------

            # si la fuente está por debajo de los receptores
            # invierto el indice de z_ind para tener la coordenada z
            # en el orden correcto
            if(Z_s > Z_r)
               z_ind[1:n_ij] = z_ind[n_ij:-1:1];
            end

            # guardo la coordenada xyz de cada interseccion de cada rayo
            # con cada discontinuidad
            normxy_rs = norm(source[1:2]-receiv[1:2])
            cos_a = (r[j,1]-s[i,1])/normxy_rs;
            sin_a = (r[j,2]-s[i,2])/normxy_rs;
            # rayo p
            cp_out[i,j,1:n_ij,1] = dp_out[i,j,1:n_ij]*cos_a;   # coord x
            cp_out[i,j,1:n_ij,2] = dp_out[i,j,1:n_ij]*sin_a;   # coord y
            cp_out[i,j,1:n_ij,3] = z[z_ind[1:n_ij]]        ;   # coord z
            # rayo sh
            csh_out[i,j,1:n_ij,1] = dsh_out[i,j,1:n_ij]*cos_a;  # coord x
            csh_out[i,j,1:n_ij,2] = dsh_out[i,j,1:n_ij]*sin_a;  # coord y
            csh_out[i,j,1:n_ij,3] = z[z_ind[1:n_ij]]         ;  # coord z
            # rayo sv
            csv_out[i,j,1:n_ij,1] = dsv_out[i,j,1:n_ij]*cos_a;  # coord x
            csv_out[i,j,1:n_ij,2] = dsv_out[i,j,1:n_ij]*sin_a;  # coord y
            csv_out[i,j,1:n_ij,3] = z[z_ind[1:n_ij]]         ;  # coord z

            cpp =  cp_out[i,j,:,:];
            csh = csh_out[i,j,:,:];
            csv = csv_out[i,j,:,:];

            # -----------------------------------------------------------------
            # Calculate trayectories & angles
            # -----------------------------------------------------------------
            pp_ray,sh_ray,sv_ray = trayectories(source,receiv,cpp,csh,csv,n_ij);
            pp_ang,sh_ang,sv_ang =       angles(source,receiv,cpp,csh,csv,n_ij);
            # -----------------------------------------------------------------
            # Calculate geometrical spreading  & amplitudes
            # -----------------------------------------------------------------
            pp_gsf,sh_gsf,sv_gsf     = gsfactor(pp_ray,sh_ray,sv_ray)                                           ;
            pp_diramp,pp_dirphase,ID = dir_amplitudes_an(C,pp_ang,sh_ang,sv_ang,receiv,source,z,rho,"pp",direct);
            sh_diramp,sh_dirphase,ID = dir_amplitudes_an(C,pp_ang,sh_ang,sv_ang,receiv,source,z,rho,"sh",direct);
            sv_diramp,sv_dirphase,ID = dir_amplitudes_an(C,pp_ang,sh_ang,sv_ang,receiv,source,z,rho,"sv",direct);

            pp_noramp = coef_normalization(rho,pp_ray,pp_ang,vp,vs,an,ID,"pp"); # eq. 4.109 "Gretcka Book".
            sh_noramp = coef_normalization(rho,sh_ray,sh_ang,vp,vs,an,ID,"sh"); # eq. 4.109 "Gretcka Book".
            sv_noramp = coef_normalization(rho,sv_ray,sv_ang,vp,vs,an,ID,"sv"); # eq. 4.109 "Gretcka Book".

            #pp_diramp,pp_dirphase,pp_traj,ID = dir_amplitudes(pp_ang,receiv,source,z,vp,vs,rho,1);
            #sh_diramp,sh_dirphase,sh_traj,ID = dir_amplitudes(sh_ang,receiv,source,z,vp,vs,rho,7);
            #sv_diramp,sv_dirphase,sv_traj,ID = dir_amplitudes(sv_ang,receiv,source,z,vp,vs,rho,2);
            T_ppdiramp = prod((pp_diramp.*pp_noramp)) * pp_gsf;
            T_shdiramp = prod((sh_diramp.*sh_noramp)) * sh_gsf;
            T_svdiramp = prod((sv_diramp.*sv_noramp)) * sv_gsf;
            # -----------------------------------------------------------------
            # Calculate wavelets
            # -----------------------------------------------------------------

            # Theory:
            # I   - I will consider a source in the last refraction point (LRP),
            # therefore, the direction of incidence on the receiver point (RP)
            # will be the correct.
            # II  - The wavelet needs to appear at calculated times, so I will
            # calculate the equivalent isot. vel. between LRP and RP.
            #     Vᵢₛₒ = Δ(LRP-RP)/Δt.
            # III - I will run isot. synthetic generator with this geometry.
            #    a) signal will appear at correct times.
            #    b) Source mechanism remains unaltered in all trayectory.
            #    c) X, Y and Z amplitud relations (normalized) are dominated by :
            #             Rᵢⱼₖ tensor  (dominated by dir of incidence, see "I -").
            #             Mᵢⱼ  tensor. (recall, signal is normalized, M₀= 1.0).

            # cpp_c: same as cpp but for currect rec_sou signal. Row dimention varies
            # depending on i-source and j-receiver i.e. number of intersections

            cpp_c     = cpp[1:n_ij,:] ;
            csh_c     = csh[1:n_ij,:] ;
            csv_c     = csv[1:n_ij,:] ;
            LRP_pp       = cpp_c[end,:]  ; # last refraction point
            LRP_sh       = csh_c[end,:]  ; # last refraction point
            LRP_sv       = csv_c[end,:]  ; # last refraction point
            RP           = receiv        ; # receiver point
            LRP_RP_pp    = RP-LRP_pp         ; # vector from LRP_pp to RP.
            LRP_RP_sh    = RP-LRP_sh         ; # vector from LRP_sh to RP.
            LRP_RP_sv    = RP-LRP_sv         ; # vector from LRP_sv to RP.
            Δ_LRP_RP_pp  = norm(LRP_RP_pp)   ; # distance from LRP_pp to RP.
            Δ_LRP_RP_sh  = norm(LRP_RP_sh)   ; # distance from LRP_sh to RP.
            Δ_LRP_RP_sv  = norm(LRP_RP_sv)   ; # distance from LRP_sv to RP.
            t_pp =  tp_out[i,j]              ; # p-phase travel time
            t_sh = tsh_out[i,j]              ; # sh-phase travel time
            t_sv = tsv_out[i,j]              ; # sv-phase travel time
            vpp_i = Δ_LRP_RP_pp/t_pp    ; # equivalent vp isot. vel
            vsh_i = Δ_LRP_RP_sh/t_sh    ; # equivalent vsh isot. vel
            vsv_i = Δ_LRP_RP_sv/t_sv    ; # equivalent vsv isot. vel

            # figure 4.20 aki-richards. Figure 4.10 aki-richards
            α_pp = take_off_angle(pp_ang)   ; # t-off ang. Where the ray comes from
            α_sh = take_off_angle(sh_ang)   ; # t-off ang.
            α_sv = take_off_angle(sv_ang)   ; # t-off ang.

            β_pp = take_in_angle(pp_ang)   ; # t-off ang. Where the ray comes from
            β_sh = take_in_angle(sh_ang)   ; # t-off ang.
            β_sv = take_in_angle(sv_ang)   ; # t-off ang.

            ϕ_pp = take_off_azimuth(LRP_RP_pp) ; # t_off azimuth
            ϕ_sh = take_off_azimuth(LRP_RP_sh) ; # t_off azimuth
            ϕ_sv = take_off_azimuth(LRP_RP_sv) ; # t_off azimuth
            # eq.4.88 Aki-Richards: Projections for P- SH- and SV- rays
            pd_pp,ϕd_pp = wave_directionals(α_pp,ϕ_pp);
            pd_sh,ϕd_sh = wave_directionals(α_sh,ϕ_sv);
            pd_sv,ϕd_sv = wave_directionals(α_sh,ϕ_sv);

            idᵣ = ID[1]; # id of receiver layer
            idₛ = ID[2]; # id of source layer

            ρᵣ   = rho[idᵣ];
            ρₛ   = rho[idₛ];

            vppᵣ = an_vel(deg2rad(α_pp),vp[idᵣ],vs[idᵣ],an[idᵣ,:],"vp");
            vshᵣ = an_vel(deg2rad(α_sh),vp[idᵣ],vs[idᵣ],an[idᵣ,:],"sh");
            vsvᵣ = an_vel(deg2rad(α_sv),vp[idᵣ],vs[idᵣ],an[idᵣ,:],"sv");
            vppₛ = an_vel(deg2rad(β_pp),vp[idₛ],vs[idₛ],an[idₛ,:],"vp");
            vshₛ = an_vel(deg2rad(β_sh),vp[idₛ],vs[idₛ],an[idₛ,:],"sh");
            vsvₛ = an_vel(deg2rad(β_sv),vp[idₛ],vs[idₛ],an[idₛ,:],"sv");


            Upp_xyz = iso_synth_generator(dt,ns,n_r,RP,t0s,LRP_pp,f0,vpp_i,vppᵣ,vppₛ,ρᵣ,ρₛ,M,"pp");
            Ush_xyz = iso_synth_generator(dt,ns,n_r,RP,t0s,LRP_sh,f0,vsh_i,vshᵣ,vshₛ,ρᵣ,ρₛ,M,"sh");
            Usv_xyz = iso_synth_generator(dt,ns,n_r,RP,t0s,LRP_sv,f0,vsv_i,vsvᵣ,vsvₛ,ρᵣ,ρₛ,M,"sv");

            # normfact = normdata(Upp_xyz + Ush_xyz + Usv_xyz)[2];
            # Upp_xyz = Upp_xyz/normfact;
            # Ush_xyz = Ush_xyz/normfact;
            # Usv_xyz = Usv_xyz/normfact;

            # theory:
            # for each ray (say pp(x,y,z), sh(x,y,z) or sv(x,y,z)), I can
            # rotate to a spherical coordinate system with axis on radial direction (r)
            # on the direction of the ray where only P lives, and a perpendicular plane,
            # where only S lives, with directions (i:inc, ϕ:azimuth).
            # If I am rotating pp(x,y,z) to pp(r,ϕ,i), I can discard pp(i), pp(ϕ), as their
            # enery is zero.
            # If I am rotating sh(x,y,z) to sh(r,ϕ,i), I can discard sh(r), sh(i), as their
            # enery is zero.
            # If I am rotating sv(x,y,z) to sv(r,ϕ,i), I can discard sv(r), sv(ϕ), as their
            # enery is zero.

            # rotate from xyz to spherical acocrding to their ray direction
            # extract radial mechanism (where pp lives)
            Upp_r = rotation(Upp_xyz,α_pp,ϕ_pp,direct,"xyz2rϕi")[:,1];
            # extract azimuthal mechanism (where sh lives)
            Ush_ϕ = rotation(Ush_xyz,α_sh,ϕ_sh,direct,"xyz2rϕi")[:,2];
            # extract inclination mechanism (where sv lives)
            Usv_i = rotation(Usv_xyz,α_sv,ϕ_sv,direct,"xyz2rϕi")[:,3];

            # input their calculated amplitudes.
            Upp_amp,Ush_amp,Usv_amp = rescale(Upp_r,Ush_ϕ,Usv_i,T_ppdiramp,T_shdiramp,T_svdiramp);

            # rotate back
            Upp_xyz = rotation(Upp_amp,α_pp,ϕ_pp,direct,"rϕi2xyz");
            # extract azimuthal mechanism (where sh lives)
            Ush_xyz = rotation(Ush_amp,α_sh,ϕ_sh,direct,"rϕi2xyz");
            # extract inclination mechanism (where sv lives)
            Usv_xyz = rotation(Usv_amp,α_sv,ϕ_sv,direct,"rϕi2xyz");

            final_signal = Upp_xyz + Ush_xyz + Usv_xyz;

         end

         signal[i,j] = normdata(final_signal)[1];
         sigpp[i,j] = Upp_xyz;
         sigsh[i,j] = Ush_xyz;
         sigsv[i,j] = Usv_xyz;
         amps[i,j]  = maxamp(Upp_xyz,Ush_xyz,Usv_xyz);
      end  # j=1,n_r
   end  # i=1,n_s
   return tp_out,tsh_out,tsv_out,signal,sigpp,sigsh,sigsv,amps;
end # function



# writedlm("signals_pp/signal_x_1.txt"    , Upp_xyz[:,1]                              );
# writedlm("signals_pp/signal_y_1.txt"    , Upp_xyz[:,2]                              );
# writedlm("signals_pp/signal_z_1.txt"    , Upp_xyz[:,3]                              );
#
# writedlm("signals_ppsh/signal_x_1.txt"  , Upp_xyz[:,1] + Ush_xyz[:,1]               );
# writedlm("signals_ppsh/signal_y_1.txt"  , Upp_xyz[:,2] + Ush_xyz[:,2]               );
# writedlm("signals_ppsh/signal_z_1.txt"  , Upp_xyz[:,3] + Ush_xyz[:,3]               );
#
# writedlm("signals_ppshsv/signal_x_1.txt", Upp_xyz[:,1] + Ush_xyz[:,1] + Usv_xyz[:,1]);
# writedlm("signals_ppshsv/signal_y_1.txt", Upp_xyz[:,2] + Ush_xyz[:,2] + Usv_xyz[:,2]);
# writedlm("signals_ppshsv/signal_z_1.txt", Upp_xyz[:,3] + Ush_xyz[:,3] + Usv_xyz[:,3]);
# error()
