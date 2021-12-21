function isotropic_ray(e::Array{Float64,1},
    v::Array{Float64,1},
    n_d::Int64,
    s::Array{Float64,1},
    r::Array{Float64,1},
    ftol::Float64,
    d_out::Array{Float64,1})

    # #====================================================================
    # #
    # # Esta subrutina calcula el rayo en un medio isotropo, utilizando
    # # la técnica propuesta en
    # #
    # # "A rapid and accurate two-point ray tracing method in
    # # horizontally layered velocity model" de Tian y Chen  #
    # # IN: e     --> espesores de las capas entre la fuente y el receptor.
    # #     v     --> velocidad de la onda, puede ser s o p
    # #     n_int --> cantida de disc. entre s y r
    # #     s     --> coordenadas de la fuente
    # #     r     --> coordenadas del receptor
    # #
    # # OUT: d_out --> coordenadas x de las intersecciones entre
    # #                el rayo y las disc. del modelo
    # #
    # #====================================================================
    #
    # #averiguo cuantas y cuales son las capas con v=v_max, y calculo el
    # #espesor e_max equivalente

    v_max,v_max_ind = findmax(v)       ; # Max vel from mod and index
    e_max           = sum(e[v.==v_max]); # Add all H with v_max. Equivalent H
    delta = norm(s[1:2]-r[1:2])        ; # Distance between source and receiver

    # calculo el q inicial
    a = 0.0;
    b = 0.0;
    q = 0.0;
    for i = 1:n_d
       a = a + (v[i]/v_max)*e[i]
       if (v[i] != v_max)
          b = b + ((v[i]/v_max)*e[i])/sqrt(1.0-(v[i]/v_max)^2)
       end
    end

    a     = a/e_max     ;
    aba_1 = (a*b)/(a-1) ;
    if (delta < aba_1);
       q = delta/a
    elseif (delta > aba_1)
       q = delta - b
    end
    # hago newton-rhapson para averiguar el parámetro q
    delta_f = F(e,v,v_max,e_max,q) #delta_f inicial
    while (abs(1.0-delta_f/delta) > ftol)
        num     = F(e,v,v_max,e_max,q)-delta;
        den     = F_p(e,v,v_max,e_max,q);
        q       = q - num/den;
        delta_f = F(e,v,v_max,e_max,q)
    end
    # calculo el tiempo de viaje
    t_out = 0.0;
    for i = 1:n_d
        raiz = sqrt(1.0-((q/(v_max*sqrt(e_max*e_max + q*q)))*v[i])^2);
       t_out = t_out + e[i]/(v[i]*raiz);
    end
    # calculo las distancias parciales hasta el receptor
    # que tambien sirven para calcular las coordenadas de la intersección
    # del rayo con las discontinuidades
    for i=1:n_d-1
       d_out[i] = F(e[1:i],v[1:i],v_max,e_max,q);
    end
    return d_out,t_out
end

function F(
    e::Array{Float64,1},
    v::Array{Float64,1},
    v_max::Float64,
    e_max::Float64,
    q::Float64)

    # !=====================================================================
    # ! Esta funcion calcula la distancia delta en función del modelo de
    # ! velocidades y del parámetro q
    # ! la función devuelve un arreglo con los valores de la distancia
    # ! horizontal recorrida por cada segmento del rayo
    # !=====================================================================
    F = 0.0;
    for i = 1:length(e)
       vimax = v[i]/v_max;
       F = F + (vimax)*e[i]*q/sqrt(e_max*e_max + (1.0-(vimax*vimax))*q^2);
    end
    return F;
  end

function F_p(
    e::Array{Float64,1},
    v::Array{Float64,1},
    v_max::Float64,
    e_max::Float64,
    q::Float64)
    # !======================================================================
    # !Esta funcion es la derivada de la anterior, mas o menos...
    # !======================================================================
    F_p = 0.0
    for i = 1:length(e)
       vimax = v[i]/v_max;
       F_p = F_p + (vimax)*e[i]/(e_max*e_max + (1.0-(vimax*vimax))*q^2)^(1.5);
    end
    F_p = F_p * e_max * e_max;
    return F_p;
end
