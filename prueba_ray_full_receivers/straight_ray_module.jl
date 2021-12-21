
function straight_ray(
    vp::Float64,
    vs::Float64,
    an::Array{Float64,1},
    s::Array{Float64,1},
    r::Array{Float64,1})
    # =======================================================================
    #
    # Esta rutina calcula los tiempos de arrivo de las ondas p,ssv y ssh
    # para el modelo de weak anisotropy de thompsen.
    #
    # ======================================================================
    dist  = norm(r-s);
    theta = atan(norm(r[1:2]-s[1:2])/abs(r[3]-s[3]));
    vpvs2 = (vp/vs)*(vp/vs);
    sthe = sin(theta);
    cthe = cos(theta);
    # tiempo de viaje de los rayos
    tp  = dist/(vp*(1.0 + an[2] * (sthe*cthe  )^2 + an[1] * sthe^4) );
    tsv = dist/(vs*(1.0 + vpvs2 * (an[1]-an[2])   * (sthe * cthe)^2));
    tsh = dist/(vs*(1.0 + an[3] * sthe^2))                           ;

return tp,tsv,tsh;
end #straight_ray
