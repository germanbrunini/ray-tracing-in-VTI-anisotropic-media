function s1s3_PP(
    C   :: Array{Float64,2},
    θ   :: Float64,
    ρ   :: Float64,
    dir :: String)

    C₁₁   = C[1,1] ;
    C₁₃   = C[1,3] ;
    C₅₅   = C[5,5] ;
    C₃₃   = C[3,3] ;
    C₃₃⁻¹ = 1.0/C₃₃;
    C₅₅⁻¹ = 1.0/C₅₅;

    # eq.6.126 - Auxiliars
    cos²    = cosd(θ)*cosd(θ)         ;
    sin²    = sind(θ)*sind(θ)         ;
    sin²_2θ = sind(2.0*θ)*sind(2.0*θ) ;
    term1   = ((C₃₃ - C₅₅) * cos² - (C₁₁ - C₅₅) * sin²)^2.0 ;
    term2   =  (C₁₃ + C₅₅) * (C₁₃ + C₅₅) * sin²_2θ          ;
    # eq.6.126 -  Equation
    ℂ       = sqrt(term1 + term2)                             ; # eq.6.126
    # eq. 6.125 - Auxiliars (for SV with "-C")
    aux = (C₅₅ + C₁₁*sin² + C₃₃*cos² + ℂ);
    # eq. 6.125 - Equation
    v_θ   = sqrt(aux/(2.0*ρ));
    s₁    = sind(θ)/v_θ;           # eq. 6.124(1)

    s₁² = s₁*s₁;

    K₁ = ρ*(C₅₅⁻¹ + C₃₃⁻¹) + C₅₅⁻¹ * (C₁₃*C₃₃⁻¹*(C₁₃+2.0*C₅₅)-C₁₁) * s₁²; # eq.6.114(1)
    K₂ = C₃₃⁻¹*(C₁₁*s₁²-ρ)                                              ; # eq.6.114(2)
    K₃ = s₁²-ρ*C₅₅⁻¹                                                    ; # eq.6.114(3)

    con1 = 1.0/sqrt(2.0);
    if dir == "down"
        #s₃ = cosd(θ)/v_θ;       # eq. 6.124(2)
        s₃ = con1*sqrt(K₁ - vpsqrt(K₁^2-4.0*K₂*K₃))
    elseif dir == "up"
        #s₃ = -abs(cosd(θ)/v_θ); # eq. 6.124(2)
        s₃ = -con1*sqrt(K₁ - vpsqrt(K₁^2-4.0*K₂*K₃))
    end

    s₃² = s₃*s₃;
    num1 = C₅₅ * s₁² + C₃₃ * s₃² - ρ;
    num2 = C₁₁ * s₁² + C₅₅ * s₃² - ρ;
    den  = C₁₁ * s₁² + C₃₃ * s₃² + C₅₅*(s₁² + s₃²) - 2.0*ρ;
    β    =  vpsqrt(num1/den); # eq 6.116
    ζ    =  vpsqrt(num2/den); # eq 6.117
    return (C,s₁,s₃,β,ζ);
end

function s1s3_SH(
    C   :: Array{Float64,2},
    θ   :: Float64,
    ρ   :: Float64,
    dir :: String)

    cos²    = cosd(θ)*cosd(θ)         ;
    sin²    = sind(θ)*sind(θ)         ;
    sin_2θ  = sind(2.0*θ)              ;
    ρ⁻¹      = 1.0/ρ                  ;

    int_sqrt = C[4,4]*cos² + C[6,6]*sin² + C[4,6]*sin_2θ      ; # eq. 1.273
    vp_θ     = sqrt(int_sqrt * ρ⁻¹)                           ; # eq. 1.273

    s₁ = sind(θ)/vp_θ                                         ; # eq. 1.272(1)
    if dir == "down"
        s₃ = cosd(θ)/vp_θ;
    elseif dir == "up"
        s₃ = -abs(cosd(θ)/vp_θ);
    end
    β = 0.0; # irrelevant, just for output completeness
    ζ = 0.0; # irrelevant, just for output completeness
    return (C,s₁,s₃,β,ζ);
end

function s1s3_SV(
    C   :: Array{Float64,2},
    θ   :: Float64         ,
    ρ   :: Float64         ,
    dir :: String           )

    C₁₁   = C[1,1] ;
    C₁₃   = C[1,3] ;
    C₅₅   = C[5,5] ;
    C₃₃   = C[3,3] ;
    C₃₃⁻¹ = 1.0/C₃₃;
    C₅₅⁻¹ = 1.0/C₅₅;

    # eq.6.126 - Auxiliars
    cos²    = cosd(θ)*cosd(θ)         ;
    sin²    = sind(θ)*sind(θ)         ;
    sin²_2θ = sind(2.0*θ)*sind(2.0*θ) ;
    term1   = ((C₃₃ - C₅₅) * cos² - (C₁₁ - C₅₅) * sin²)^2.0 ;
    term2   =  (C₁₃ + C₅₅) * (C₁₃ + C₅₅) * sin²_2θ          ;
    # eq.6.126 -  Equation
    ℂ       = sqrt(term1 + term2)                             ; # eq.6.126

    # eq. 6.125 - Auxiliars (for SV with "-C")
    aux = (C₅₅ + C₁₁*sin² + C₃₃*cos² - ℂ);
    # eq. 6.125 - Equation
    v_θ   = sqrt(aux/(2.0*ρ));
    s₁    = sind(θ)/v_θ;           # eq. 6.124(1)

    s₁² = s₁*s₁;

    K₁ = ρ*(C₅₅⁻¹ + C₃₃⁻¹) + C₅₅⁻¹ * (C₁₃*C₃₃⁻¹*(C₁₃+2.0*C₅₅)-C₁₁) * s₁²; # eq.6.114(1)
    K₂ = C₃₃⁻¹*(C₁₁*s₁²-ρ)                                              ; # eq.6.114(2)
    K₃ = s₁²-ρ*C₅₅⁻¹                                                    ; # eq.6.114(3)

    con1 = 1.0/sqrt(2.0)
    if dir == "down"
        #s₃ = cosd(θ)/v_θ;       # eq. 6.124(2)
        s₃ = con1*sqrt(K₁ + vpsqrt(K₁^2-4.0*K₂*K₃))
    elseif dir == "up"
        #s₃ = -abs(cosd(θ)/v_θ); # eq. 6.124(2)
        s₃ = -con1*sqrt(K₁ + vpsqrt(K₁^2-4.0*K₂*K₃))
    end

    s₃² = s₃*s₃;
    num1 = C₅₅ * s₁² + C₃₃ * s₃² - ρ;
    num2 = C₁₁ * s₁² + C₅₅ * s₃² - ρ;
    den  = C₁₁ * s₁² + C₃₃ * s₃² + C₅₅*(s₁² + s₃²) - 2.0*ρ;
    β    =   vpsqrt(num1/den); # eq 6.116
    ζ    =  -vpsqrt(num2/den); # eq 6.117

    return (C,s₁,s₃,β,ζ);
end

function Coef_PP(
    l₁ᴾ :: Tuple{Array{Float64,2},Float64,Float64,Float64,Float64},
    l₂ᴾ :: Tuple{Array{Float64,2},Float64,Float64,Float64,Float64},
    l₁ˢ :: Tuple{Array{Float64,2},Float64,Float64,Float64,Float64},
    l₂ˢ :: Tuple{Array{Float64,2},Float64,Float64,Float64,Float64} )

    C₁ᴾ  = l₁ᴾ[1]; C₁ˢ  = l₁ˢ[1];
    s₁₁ᴾ = l₁ᴾ[2]; s₁₁ˢ = l₁ˢ[2];
    s₃₁ᴾ = l₁ᴾ[3]; s₃₁ˢ = l₁ˢ[3];
    β₁ᴾ  = l₁ᴾ[4]; β₁ˢ  = l₁ˢ[4];
    ζ₁ᴾ  = l₁ᴾ[5]; ζ₁ˢ  = l₁ˢ[5];

    C₂ᴾ  = l₂ᴾ[1]; C₂ˢ  = l₂ˢ[1];
    s₁₂ᴾ = l₂ᴾ[2]; s₁₂ˢ = l₂ˢ[2];
    s₃₂ᴾ = l₂ᴾ[3]; s₃₂ˢ = l₂ˢ[3];
    β₂ᴾ  = l₂ᴾ[4]; β₂ˢ  = l₂ˢ[4];
    ζ₂ᴾ  = l₂ᴾ[5]; ζ₂ˢ  = l₂ˢ[5];


    W₁ᴾ,X₁ᴾ,Z₁ᴾ = WXZ(C₁ᴾ,s₁₁ᴾ,s₃₁ᴾ,β₁ᴾ,ζ₁ᴾ); # eq. 6.121
    W₂ᴾ,X₂ᴾ,Z₂ᴾ = WXZ(C₂ᴾ,s₁₂ᴾ,s₃₂ᴾ,β₂ᴾ,ζ₂ᴾ); # eq. 6.121

    W₁ˢ,X₁ˢ,Z₁ˢ = WXZ(C₁ˢ,s₁₁ˢ,s₃₁ˢ,β₁ˢ,ζ₁ˢ); # eq. 6.121
    W₂ˢ,X₂ˢ,Z₂ˢ = WXZ(C₂ˢ,s₁₂ˢ,s₃₂ˢ,β₂ˢ,ζ₂ˢ); # eq. 6.121

    mat  = [β₁ᴾ  β₁ˢ -β₂ᴾ -β₂ˢ ;   # scattering mtrix eq. 6.139
    ζ₁ᴾ  ζ₁ˢ  ζ₂ᴾ  ζ₂ˢ ;
    Z₁ᴾ  Z₁ˢ -Z₂ᴾ -Z₂ˢ ;
    W₁ᴾ  W₁ˢ  W₂ᴾ  W₂ˢ];
    RHS  =  [-β₁ᴾ;ζ₁ᴾ;-Z₁ᴾ;W₁ᴾ] ;    # eq. 6.139

    coefs = inv(mat)*RHS;
    Rₚₚ = coefs[1];
    Rₚₛ = coefs[2];
    Tₚₚ = coefs[3];
    Tₚₛ = coefs[4];
    return Rₚₚ,Rₚₛ,Tₚₚ,Tₚₛ;
end

function Coef_SH(l₁ˢ,l₂ˢ)
    C₁ˢ  = l₁ˢ[1];
    s₁₁ˢ = l₁ˢ[2];
    s₃₁ˢ = l₁ˢ[3];

    C₂ˢ  = l₂ˢ[1];
    s₁₂ˢ = l₂ˢ[2];
    s₃₂ˢ = l₂ˢ[3];

    Zᴵ = C₁ˢ[4,6]*s₁₁ˢ + C₁ˢ[4,4]*s₃₁ˢ;  # eq.1.269
    Zᵀ = C₂ˢ[4,6]*s₁₂ˢ + C₂ˢ[4,4]*s₃₂ˢ;  # eq.1.269
    Zᴿ = - Zᴵ                         ;  # eq.1.281

    Rₛₛ = (Zᴵ- Zᵀ)/(Zᵀ-Zᴿ); # eq. 1.279(1)
    Tₛₛ = 1.0 + Rₛₛ       ; # eq. 1.277
    Rₛₚ = 0.0             ; # eq. irrelevant, NO SH reflection to P
    Tₛₚ = 0.0             ; # eq. irrelevant, NO SH refraction to P
    return Rₛₚ,Rₛₛ,Tₛₚ,Tₛₛ;
end

function Coef_SV(
    l₁ᴾ :: Tuple{Array{Float64,2},Float64,Float64,Float64,Float64},
    l₂ᴾ :: Tuple{Array{Float64,2},Float64,Float64,Float64,Float64},
    l₁ˢ :: Tuple{Array{Float64,2},Float64,Float64,Float64,Float64},
    l₂ˢ :: Tuple{Array{Float64,2},Float64,Float64,Float64,Float64})

    C₁ᴾ  = l₁ᴾ[1]; C₁ˢ  = l₁ˢ[1];
    s₁₁ᴾ = l₁ᴾ[2]; s₁₁ˢ = l₁ˢ[2];
    s₃₁ᴾ = l₁ᴾ[3]; s₃₁ˢ = l₁ˢ[3];
    β₁ᴾ  = l₁ᴾ[4]; β₁ˢ  = l₁ˢ[4];
    ζ₁ᴾ  = l₁ᴾ[5]; ζ₁ˢ  = l₁ˢ[5];

    C₂ᴾ  = l₂ᴾ[1]; C₂ˢ  = l₂ˢ[1];
    s₁₂ᴾ = l₂ᴾ[2]; s₁₂ˢ = l₂ˢ[2];
    s₃₂ᴾ = l₂ᴾ[3]; s₃₂ˢ = l₂ˢ[3];
    β₂ᴾ  = l₂ᴾ[4]; β₂ˢ  = l₂ˢ[4];
    ζ₂ᴾ  = l₂ᴾ[5]; ζ₂ˢ  = l₂ˢ[5];


    W₁ᴾ,X₁ᴾ,Z₁ᴾ = WXZ(C₁ᴾ,s₁₁ᴾ,s₃₁ᴾ,β₁ᴾ,ζ₁ᴾ); # eq. 6.121
    W₂ᴾ,X₂ᴾ,Z₂ᴾ = WXZ(C₂ᴾ,s₁₂ᴾ,s₃₂ᴾ,β₂ᴾ,ζ₂ᴾ); # eq. 6.121

    W₁ˢ,X₁ˢ,Z₁ˢ = WXZ(C₁ˢ,s₁₁ˢ,s₃₁ˢ,β₁ˢ,ζ₁ˢ); # eq. 6.121
    W₂ˢ,X₂ˢ,Z₂ˢ = WXZ(C₂ˢ,s₁₂ˢ,s₃₂ˢ,β₂ˢ,ζ₂ˢ); # eq. 6.121

    mat  = [β₁ᴾ  β₁ˢ -β₂ᴾ -β₂ˢ ;   # scattering mtrix eq. 6.139
    ζ₁ᴾ  ζ₁ˢ  ζ₂ᴾ  ζ₂ˢ ;
    Z₁ᴾ  Z₁ˢ -Z₂ᴾ -Z₂ˢ ;
    W₁ᴾ  W₁ˢ  W₂ᴾ  W₂ˢ];
    RHS  =  [-β₁ˢ;ζ₁ˢ;-Z₁ˢ;W₁ˢ] ;    # eq. 6.140

    coefs = inv(mat)*RHS;
    Rₛₚ = coefs[1];
    Rₛₛ = coefs[2];
    Tₛₚ = coefs[3];
    Tₛₛ = coefs[4];
    return Rₛₚ,Rₛₛ,Tₛₚ,Tₛₛ;
end

function WXZ(
    C  :: Array{Float64,2},
    s₁ :: Float64,
    s₃ :: Float64,
    β  :: Float64,
    ζ  :: Float64 )

    W = ζ*C[5,5]*s₁ + β*C[5,5]*s₃; # eq. 6.121
    X = β*C[1,1]*s₁ + ζ*C[1,3]*s₃; # eq. 6.122
    Z = β*C[1,3]*s₁ + ζ*C[3,3]*s₃; # eq. 6.123
    return W,X,Z;
end


function vpsqrt(arg)
    #principal value of sqrt
    if arg < 0.0
        return sqrt(Complex(arg));
    else
        return sqrt(arg);
    end
end


function Thomsen_to_Cij(
    an  :: Array{Float64,2},
    vₚ₀ :: Array{Float64,1},
    vₛ₀ :: Array{Float64,1},
    ρ   :: Array{Float64,1})
    # equations extracted from:
    # "Seismic signatures and Analysis of Reflection Data in Anisotropic Media"
    # Third Edition, by Ilya Tsvankin

    ϵ = an[:,1];
    δ = an[:,2];
    γ = an[:,3];

    len = length(ρ);
    C   = zeros(Float64,(6,6,len));
    # See Page Page 12, Eq 1.24
    for i = 1:len;
        # See Page 18, Eq, q.44 - 1.48
        vₚ² = vₚ₀[i]*vₚ₀[i]    ;
        vₛ² = vₛ₀[i]*vₛ₀[i]    ;
        A   = ρ[i]*(vₚ² - vₛ²) ;
        A²  = A*A              ;

        # anisotropic                                            # isotropic
        C[1,1,i] = ρ[i]*(2.0*ϵ[i] + 1.0)*vₚ²                   ; # λ+2μ
        C[6,6,i] = ρ[i]*(2.0*γ[i] + 1.0)*vₛ²                   ; # μ
        C[3,3,i] = ρ[i]*vₚ²                                    ; # λ+2μ
        C[5,5,i] = ρ[i]*vₛ²                                    ; # μ
        C[1,3,i] = sqrt(2.0*δ[i]*C[3,3,i]*A  + A²) - C[5,5,i]  ; # λ

        C[1,2,i] = C[1,1,i] - 2.0*C[6,6,i]                     ; # See eq. 1.24
        C[3,1,i] = C[1,3,i]                                    ; # See eq. 1.24
        C[2,1,i] = C[1,2,i]                                    ; # See eq. 1.24
        C[2,2,i] = C[1,1,i]                                    ; # See eq. 1.24
        C[3,2,i] = C[1,3,i]                                    ; # See eq. 1.24
        C[2,3,i] = C[1,3,i]                                    ; # See eq. 1.24
        C[4,4,i] = C[5,5,i]                                    ; # See eq. 1.24
    end
    # For more information, also see:
    # "Wave fields in Real Media: wave propagation in Anisotropic,
    # Anelastic, Porous and Electromagnetic Media", by José Carcione. Page 6.
    return C;
end


function coef_normalization(
    ρ    :: Array{Float64,1},
    traj :: Array{Float64,2},
    ang  :: Array{Float64,1},
    vₚ   :: Array{Float64,1},
    vₛ   :: Array{Float64,1},
    an   :: Array{Float64,2},
    ID   :: Array{Int64,1}  ,
    ray  :: String           )
    # Compute normalization from eq. 4.109 "Gretcka Microseismic Book".

    if ray == "pp" ray = "vp" end; #redifine for an_vel input

    idᵣ = ID[1]; # id of receiver layer
    idₛ = ID[2]; # id of source layer
    rangId = ID_range(idₛ,idᵣ) ;


    len    = length(rangId)    ;
    coef = zeros(Float64,(len,1));

    coef[1,1] = 1.0 # First its always 1.0
    if len==1 # No interfaces
        return [1.0,1.0]
    else
        ρ   =  ρ[rangId]        ;
        vp₀ = vₚ[rangId]        ;
        vs₀ = vₛ[rangId]        ;
        an₀ = an[rangId,:]      ;
        v  = zeros(Float64,(len,1));
        gᵥ = zeros(Float64,(len,3));
        for i = 2:len
            # anisotropic velocity
            v[i-1] = an_vel(deg2rad(ang[i-1]),vp₀[i-1],vs₀[i-1],an[i-1,:],ray);
            v[i]   = an_vel(deg2rad(ang[i])  ,vp₀[i]  ,vs₀[i]  ,an[i,:]  ,ray);
            # group velocity vector
            gᵥ[i-1,:] = v[i-1] * (traj[i-1,:]/norm(traj[i-1,:])) ;
            gᵥ[i,:]   = v[i]   * (traj[i,:]  /norm(traj[i  ,:])) ;
            num       = ρ[i]   * abs(dot(gᵥ[i,:]  ,[0,0,1])); # eq. 4.109 Gretcka Book.
            den       = ρ[i-1] * abs(dot(gᵥ[i-1,:],[0,0,1])); # eq. 4.109 Gretcka Book.
            coef[i]   = sqrt(num/den);
    end
    return coef;
    end
end
