function gsfactor(
    pp_traj::Array{Float64,2},
    sh_traj::Array{Float64,2},
    sv_traj::Array{Float64,2})
    # SIMPLE GEOMETRICAL SPREADING CORRECTION
    # INPUT
    # traj: TRAJECTORY INFORMATION FROM RAY TRACING CODE
    # OUTPUT
    # out: GEOMETRICAL SPREADING FACTOR TO MULTIPLY FINAL AMPLITUDE
    pp_out  = 1.0              ;
    sh_out  = 1.0              ;
    sv_out  = 1.0              ;
    lent = size(pp_traj,1)     ;

    if lent == 2 # if sou and rec in same layer
        return pp_out,sh_out,sv_out; # No geom spreading, as it is calculated
    else                             #  in synth generator routine
        for k = 1:lent-2      # last one (lent-1) is calculated in synth generator
            pp_vec = pp_traj[k,:] - pp_traj[k + 1,:]; # trajectory vectors
            sh_vec = sh_traj[k,:] - sh_traj[k + 1,:]; # trajectory vectors
            sv_vec = sv_traj[k,:] - sv_traj[k + 1,:]; # trajectory vectors
            pp_r   = norm(pp_vec)     ; # longitud of trajectory
            sh_r   = norm(sh_vec)     ;
            sv_r   = norm(sv_vec)     ;
            pp_out = pp_out*(1.0/pp_r);
            sh_out = sh_out*(1.0/sh_r);
            sv_out = sv_out*(1.0/sv_r);
        end
        return pp_out,sh_out,sv_out;
    end
end
