function dir_amplitudes(
    ang::Array{Float64,1},
    rec::Array{Float64,1},
    sou::Array{Float64,1},
    z::Array{Float64,1},
    vp::Array{Float64,1},
    vs::Array{Float64,1},
    rho::Array{Float64,1},
    opt::Int64)
    # DETERMINE TRANSMISSION COEFFICIENTS FOR DIRECT WAVE TRAVEL PATH
    # opt VARIABLE VALUES
    # 1 P TO P TRANSMISSION COEFFICIENT
    # 2 Sv TO Sv TRANSMISSION COEFFICIENT
    # 7 Sh TO Sh TRANSMISSION COEFFICIENT
    # INPUT UNITS IN SI SYSTEM WITH DENSITY IN kg/m3
    top  = hcat([0.0],z');
    nlay = length(top)   ;
    vp   = hcat(vp',[vp[end]]);
    vs   = hcat(vs',[vs[end]]);
    rpos = [norm(rec[1:2]),rec[end]];
    spos = [norm(sou[1:2]),sou[end]];

    # DIRECTION OF INITIAL RAY PROPAGATION IN THE VERTICAL DIRECTION
    dirV = sign(rpos[2] - spos[2]);

    # Id OF THE MODEL LAYER IMMEDIATELY ABOVE THE RECEIVER LOCATION
    rposId = Idmodel_above_rec(nlay,rpos,top)
    # Id OF THE MODEL LAYER IMMEDIATELY ABOVE THE SOURCE LOCATION
    rposId,
    sposId,
    dirV,
    vptraj,vstraj,rhotraj = Idmodel_above_sou(
                                            nlay,rposId,
                                            dirV,spos,
                                            top,vp,vs,rho)

    # VELOCITY VECTOR
    count = 1;
    if dirV < 0
        aux_vec     = sposId + dirV:dirV:rposId;
        aux_len     = length(aux_vec);
        vptraj_aux  = zeros(Float64,(aux_len + 1,1));
        vstraj_aux  = zeros(Float64,(aux_len + 1,1));
        rhotraj_aux = zeros(Float64,(aux_len + 1,1));
        vptraj_aux[1,:]  = vptraj;
        vstraj_aux[1,:]  = vstraj;
        rhotraj_aux[1,:] = rhotraj;
        for k in floor.(Int64,aux_vec);
            count = count + 1;
            vptraj_aux[count,1]  =  vp[k];
            vstraj_aux[count,1]  =  vs[k];
            rhotraj_aux[count,1] = rho[k];
        end
        vptraj  = [vptraj_aux  ;  vp[rposId - 1]];
        vstraj  = [vstraj_aux  ;  vs[rposId - 1]];
        rhotraj = [rhotraj_aux ; rho[rposId - 1]];

    elseif dirV > 0
        aux_vec     = sposId + dirV:dirV:rposId - 1;
        aux_len     = length(aux_vec);
        vptraj_aux  = zeros(Float64,(aux_len + 1,1));
        vstraj_aux  = zeros(Float64,(aux_len + 1,1));
        rhotraj_aux = zeros(Float64,(aux_len + 1,1));
        vptraj_aux[1,:]  = vptraj;
        vstraj_aux[1,:]  = vstraj;
        rhotraj_aux[1,:] = rhotraj;
        for k in floor.(Int64,aux_vec);
            count = count + 1;
            vptraj_aux[count,1]  = vp[k];
            vstraj_aux[count,1]  = vs[k];
            rhotraj_aux[count,1] = rho[k];
        end
        vptraj  = [vptraj_aux  ;  vp[rposId]];
        vstraj  = [vstraj_aux  ;  vs[rposId]];
        rhotraj = [rhotraj_aux ; rho[rposId]];
    end
    out_traj = hcat(vptraj,vstraj,rhotraj);

    # RAY PARAMETER
    if opt == 1
        p = sind(ang[1])/vptraj[1];
    elseif (opt == 2 || opt == 7)
        p = sind(ang[1])/vstraj[1];
    end

    ampout,phaseout = amp_at_sorce_loc(vptraj,vstraj,rhotraj,p,opt);

    ID = [rposId,sposId];
    return ampout,phaseout,out_traj,ID;
end



function Idmodel_above_rec(
    nlay::Int64,
    rpos::Array{Float64,1},
    top::Array{Float64,2})

    rposId = Array{Int64}(undef, 1, 1);
    for k = 1:nlay + 1
        if rpos[2] <= top[k]
            rposId[1] = k - 1;
            break;
        end
    end
    return rposId[1];
end

function Idmodel_above_sou(
    nlay::Int64,
    rposId::Int64,
    dirV::Float64,
    spos::Array{Float64,1},
    top::Array{Float64,2},
    vp::Array{Float64,2},
    vs::Array{Float64,2},
    rho::Array{Float64,1})

    vptraj = Array{Float64}(undef,1,1);
    vstraj = Array{Float64}(undef,1,1);
    rhotraj = Array{Float64}(undef,1,1);
    sposId = Array{Int64}(undef, 1, 1);
    aux = rposId[1];
    for k = 1:nlay + 1
        if spos[2] <= top[k]
            sposId[1] = k - 1;
            if sposId[1] == aux
                dirV = 0;
            end
            if dirV < 0
                aux = aux + 1;
            end
            vptraj[1,1]  = vp[sposId[1]];
            vstraj[1,1]  = vs[sposId[1]];
            rhotraj[1,1] = rho[sposId[1]];
            break
        end
    end
    return aux,sposId[1],dirV,vptraj,vstraj,rhotraj;
end

function  amp_at_sorce_loc(
    vptraj::Array{Float64,2},
    vstraj::Array{Float64,2},
    rhotraj::Array{Float64,2},
    p::Float64,
    opt::Int64)

    ampout   = zeros(Float64,(length(vptraj),1));
    phaseout = zeros(Float64,(length(vptraj),1));
    if length(vptraj)==1
        ampout   = zeros(Float64,(2,1));
        phaseout = zeros(Float64,(2,1));
    end
    # AMPLITUDE AT SOURCE LOCATION
    ampout[1,1]   = 1.0;
    phaseout[1,1] = 0.0;

    if length(vptraj) > 2
        for k = 2:length(vptraj)
            temp = RefTransCoef(vptraj[k-1],vstraj[k-1],rhotraj[k-1],vptraj[k],vstraj[k],rhotraj[k],p,opt);
            if isreal(temp)
                ampout[k,1]   = temp;
                phaseout[k,1] = 0   ;
            else
                ampout[k,1]   = abs(temp)  ;
                phaseout[k,1] = angle(temp);
            end
        end
        return ampout,phaseout;
    else
        ampout[2,1]   = 1;
        phaseout[2,1] = 0;
        return ampout[1:2,1],phaseout[1:2,1];
    end
end
