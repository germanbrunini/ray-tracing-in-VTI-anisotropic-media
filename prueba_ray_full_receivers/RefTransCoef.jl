function RefTransCoef(
    vp1::Float64,
    vs1::Float64,
    rho1::Float64,
    vp2::Float64,
    vs2::Float64,
    rho2::Float64,
    p::Float64,
    opt::Int64)

# REFLECTION AND TRANSMISSION COEFFICIENTS IN ISOTROPIC MEDIA USING
# KNOTT-ZOEPPRITZ EQUATIONS (EXACT SOLUTION), vp1, vs1 AND rho1  ARE THE
# PROPERTIES OF THE MEDIUM WHERE THE WAVE IS INITIALLY TRAVELING
# opt VARIABLE
# 1 P TO P TRANSMISSION COEFFICIENT
# 2 Sv TO Sv TRANSMISSION COEFFICIENT
# 3 P TO p REFLECTION COEFFICIENT
# 4 P TO sv REFLECTION COEFFICIENT
# 5 Sv TO p REFLECTION COEFFICIENT
# 6 Sv TO sv REFLECTION COEFFICIENT
# 7 Sh TO Sh TRANSMISSION COEFFICIENT
# 8 Sh TO sh REFLECTION COEFFICIENT
#################################
# rho INPUT UNITS MUST BE kg/m3
# vp AND vs IN m/s
#################################

# SHEAR MODULII
mu1 = rho1*vs1^2;
mu2 = rho2*vs2^2;

Pang1 = asind(p*vp1);
Pang2 = asind(p*vp2);
Sang1 = asind(p*vs1);
Sang2 = asind(p*vs2);

cosP1 = cosd(Pang1);
cosP2 = cosd(Pang2);
sinP1 = sind(Pang1);
sinP2 = sind(Pang2);

cosS1 = cosd(Sang1);
cosS2 = cosd(Sang2);
sinS1 = sind(Sang1);
sinS2 = sind(Sang2);

eta1 = cosS1/vs1;
eta2 = cosS2/vs2;

# SEE page 144. Aki Richards For next formulas
# eq.5.38/39(1) Aki-Richards
a = rho2*(1.0 - 2.0*sinS2^2) - rho1*(1.0 - 2.0*sinS1^2);
b = rho2*(1.0 - 2.0*sinS2^2) + 2.0*rho1*sinS1^2;
c = rho1*(1.0 - 2.0*sinS1^2) + 2.0*rho2*sinS2^2;
d = 2.0*(rho2*vs2*vs2 - rho1*vs1*vs1);

E = b*cosP1/vp1 + c*cosP2/vp2;     # eq.5.39(1) Aki-Richards
F = b*cosS1/vs1 + c*cosS2/vs2;     # eq.5.39(2) Aki-Richards
G = a - d*(cosP1*cosS2)/(vp1*vs2); # eq.5.39(3) Aki-Richards
H = a - d*(cosP2*cosS1)/(vp2*vs1); # eq.5.39(4) Aki-Richards
D = E*F + G*H*p^2                ; # eq.5.39(5) Aki-Richards

if opt == 1 # P to P
    out = 2*rho1*cosP1*F/(vp2*D);
elseif opt == 2
    out = 2*rho1*cosS1*E/(vs2*D);
# elseif opt == 3
#     out = ( (b*cosP1/vp1-c*cosP2/vp2)*F-(a+d*cosP1*cosS2/(vp1*vs2))*H*p^2 )/D;
# elseif opt == 4
#     out = ( -2*(cosP1/vp1)*(a*b + c*d*cosP2*cosS2/(vp2*vs2))*p*vp1 )/(vs1*D);
# elseif opt == 5
#     out = ( -2*(cosS1/vs1)*(a*b + c*d*cosP2*cosS2/(vp2*vs2))*p*vs1 )/(vp1*D);
# elseif opt == 6
#     out = - ((b*cosS1/vs1-c*cosS2/vs2)*E-(a+d*cosP2*cosS1/(vp2*vs1))*G*p^2 )/D;
elseif opt == 7
    out = 2.0*mu1*eta1/(mu1*eta1 + mu2*eta2);
# elseif opt == 8
#     out = (mu1*eta1 - mu2*eta2)/(mu1*eta1 + mu2*eta2);
else
    println("non-valid option")
end

return out;
end
