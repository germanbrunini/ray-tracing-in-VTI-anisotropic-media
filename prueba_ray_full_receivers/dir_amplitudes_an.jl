function dir_amplitudes_an(
    Cfs::Array{Float64,3},
    θpp::Array{Float64,1},
    θsh::Array{Float64,1},
    θsv::Array{Float64,1},
    rec::Array{Float64,1},
    sou::Array{Float64,1},
    z::Array{Float64,1},
    rho::Array{Float64,1},
    opt::String,
    dir::String)

    # DETERMINE TRANSMISSION COEFFICIENTS FOR DIRECT WAVE TRAVEL PATH
    # opt VARIABLE VALUES
    # 1 P TO P TRANSMISSION COEFFICIENT
    # 2 Sv TO Sv TRANSMISSION COEFFICIENT
    # 7 Sh TO Sh TRANSMISSION COEFFICIENT
    # INPUT UNITS IN SI SYSTEM WITH DENSITY IN kg/m3
    Zr = rec[3];
    Zs = sou[3];
    # Id OF THE MODEL LAYER IMMEDIATELY ABOVE THE RECEIVER LOCATION
    rposId = ID_layer(z,Zr);
    sposId = ID_layer(z,Zs);
    rangId = ID_range(sposId,rposId);
    ρ      = rho[rangId];
    C      = Cfs[:,:,rangId];

    ampout,phaseout = amp_at_sorce_loc_an(C,θpp,θsh,θsv,ρ,opt,dir);
    ID = [rposId,sposId];
    return ampout,phaseout,ID;
end

function amp_at_sorce_loc_an(C,
    θpp::Array{Float64,1},
    θsh::Array{Float64,1},
    θsv::Array{Float64,1},
    ρ::Array{Float64,1},
    opt::String,
    dir::String)

    lenρ = length(ρ);
    ampout   = zeros(Float64,(lenρ,1));
    phaseout = zeros(Float64,(lenρ,1));
    if lenρ==1
        ampout   = zeros(Float64,(2,1));
        phaseout = zeros(Float64,(2,1));
    end
    # AMPLITUDE AT SOURCE LOCATION
    # Firts is always 1.0. It means that starting amplitud is taken equal to 1.
    ampout[1,1]   = 1.0;
    phaseout[1,1] = 0.0; #

    if lenρ >= 2
        # If there is one or more interfaces (lenρ>=2), then ampout=[1,c1,c2,...].         
        for k = 2:lenρ
                # first medium properties: incident and reflected
                l1_pp = s1s3_PP(C[:,:,k-1],θpp[k-1],ρ[k-1],dir);
                l1_sh = s1s3_SH(C[:,:,k-1],θsh[k-1],ρ[k-1],dir);
                l1_sv = s1s3_SV(C[:,:,k-1],θsv[k-1],ρ[k-1],dir);
                # second medium properties: refracted
                l2_pp = s1s3_PP(C[:,:,k]  ,θpp[k]  ,ρ[k],dir);
                l2_sh = s1s3_SH(C[:,:,k]  ,θsh[k]  ,ρ[k],dir);
                l2_sv = s1s3_SV(C[:,:,k]  ,θsv[k]  ,ρ[k],dir);

            if opt == "pp"  # PP
                R1,R2,temp,T2 = Coef_PP(l1_pp,l2_pp,l1_sv,l2_sv);
            end
            if opt == "sv" # SV
                R1,R2,T1,temp = Coef_SV(l1_pp,l2_pp,l1_sv,l2_sv);
            end
            if opt == "sh"  # SH
                R1,R2,T1,temp = Coef_SH(l1_sh,l2_sh);
            end
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
        # If there is no interfaces (lenρ==1), then ampout=[1,1].
        ampout[2,1]   = 1.0;
        phaseout[2,1] = 0.0;
        return ampout[1:2,1],phaseout[1:2,1];
    end
end

# function Idmodel_above_rec_an(
#     nlay::Int64,
#     rpos::Array{Float64,1},
#     top::Array{Float64,2})
#
#     rposId = Array{Int64}(undef, 1, 1);
#     for k = 1:nlay + 1
#         if rpos[2] <= top[k]
#             rposId[1] = k - 1;
#             break;
#         end
#     end
#     return rposId[1];
# end
#
# function Idmodel_above_sou_an(
#     nlay::Int64,
#     rposId::Int64,
#     dirV::Float64,
#     spos::Array{Float64,1},
#     top::Array{Float64,2},
#     rho::Array{Float64,1})
#
#     rhotraj = Array{Float64}(undef,1,1);
#     sposId = Array{Int64}(undef, 1, 1);
#     aux = rposId[1];
#     for k = 1:nlay + 1
#         if spos[2] <= top[k]
#             sposId[1] = k - 1;
#             if sposId[1] == aux
#                 dirV = 0;
#             end
#             if dirV < 0
#                 aux = aux + 1;
#             end
#             rhotraj[1,1] = rho[sposId[1]];
#             break
#         end
#     end
#     return aux,sposId[1],dirV,rhotraj;
# end
