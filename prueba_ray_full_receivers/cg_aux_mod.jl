# ==============================================================================
function func(
  x_in::Array{Float64,1},
  e::Array{Float64,1},
  s::Array{Float64,1},
  r::Array{Float64,1},
  vp::Array{Float64,1},
  vs::Array{Float64,1},
  an::Array{Float64,2},
  ray::String) # Tiempo de viaje

  len_e = size(e,1)              ;
  xin1  = size(x_in,1)           ;
  x     = zeros(Float64,xin1 + 2);

  x[1         ] = 0.0                ;
  x[2:xin1 + 1] = x_in               ;
  x[xin1   + 2] = norm(r[1:2]-s[1:2]);

  func    = 0.0                      ;
  mem_fun = zeros(Float64,(len_e,1)) ;
  for i = 2:len_e + 1
    dx   = abs(x[i-1]-x[i]);
    phi  = atan(dx/e[i-1]) ;
    denominador  = an_vel(phi,vp[i-1],vs[i-1],an[i-1,:],ray)*sin(phi);
    func         = func + dx / denominador; # sumes up all Δt
    mem_fun[i-1] = dx / denominador       ; # stores each Δt
  end
  return func,mem_fun;
end #function func

# ==============================================================================

function dfunc(
  x_in::Array{Float64,1},
  e::Array{Float64,1},
  s::Array{Float64,1},
  r::Array{Float64,1},
  vp::Array{Float64,1},
  vs::Array{Float64,1},
  an::Array{Float64,2},
  ray::String)

  # !Para calcular la derivada de un punto x[i] se necesita también el punto
  # x(i+1) y el punto x[i-1], por lo tanto se van a calcular las derivadas
  # desde i=2 a i=N-1, siendo x(N) el receptor; x[1] es la fuente.

  xin1  = size(x_in,1)           ;
  x     = zeros(Float64,xin1 + 2);
  dfunc = zeros(Float64,xin1    );

  x[1         ] = 0.0      ;
  x[2:xin1 + 1] = x_in     ;
  x[xin1   + 2] = norm(r[1:2]-s[1:2]);

  for i = 2:size(x,1)-1
    dx1   = x[i]   - x[i-1]              ;
    dx2   = x[i+1] - x[i]                ;
    phi1  = atan(dx1 / e[i-1])           ;
    phi2  = atan(dx2 / e[i])             ;
    dphi1 = e[i-1]   / (e[i-1]*e[i-1] + dx1*dx1);
    dphi2 = -e[i]    / (e[i]  *e[i]   + dx2*dx2);
    v1    = an_vel(phi1, vp[i-1], vs[i-1], an[i-1,:], ray);
    v2    = an_vel(phi2, vp[i]  , vs[i]  , an[i,:]  , ray);
    sp1   = sin(phi1);
    sp2   = sin(phi2);
    dv1  = an_dvel(phi1, vp[i-1], vs[i-1], an[i-1,:], ray);
    dv2  = an_dvel(phi2, vp[i]  , vs[i]  , an[i,:]  , ray);

    aux1  = cos(phi1)*v1 + sp1*dv1;
    aux2  = cos(phi2)*v2 + sp2*dv2;

    den1 = (sp1*v1)^2;
    den2 = (sp2*v2)^2;

    num1 = sp1*v1 - dx1*dphi1*aux1;
    num2 = sp2*v2 + dx2*dphi2*aux2;
    dfunc[i-1] = (num1/den1)-(num2/den2);

  end
  return dfunc;
end #function dfunc

# ==============================================================================


function an_vel(
  phi::Float64,
  vp::Float64,
  vs::Float64,
  an::Array{Float64,1},
  ray::String)

  an_vel = 0.0;
  if (ray=="vp")
    an_vel = vp * (1.0 + an[2]*(sin(phi)*cos(phi))^2 + an[1]*sin(phi)^4)
  elseif (ray=="sv")
    vpvs2  = (vp/vs)^2;
    an_vel = vs * (1.0 + vpvs2*(an[1]-an[2])*(sin(phi)*cos(phi))^2)
  elseif (ray=="sh")
    an_vel = vs * (1.0 + an[3] * sin(phi)^2)
  else
    @warn("input wrong key word argument in an_vel")
  end
  return an_vel;

end # function an_vel

# ==============================================================================

function an_dvel(
  phi::Float64,
  vp::Float64,
  vs::Float64,
  an::Array{Float64,1},
  ray::String)

  an_dvel=0.0;
  if (ray=="vp")
    vp2 = 2.0*vp;
    an_dvel = vp2*sin(phi)*cos(phi)*(an[2]*cos(phi)^2+(2.0*an[1]-an[2])*sin(phi)^2);
  elseif (ray=="sv")
    vs2   = 2.0*vs   ;
    vpvs2 = (vp/vs)^2;
    an_dvel = vs2*(an[1]-an[2])*(vpvs2)*sin(phi)*cos(phi)*(cos(phi)^2-sin(phi)^2);
  elseif (ray=="sh")
    vs2      = 2.0*vs;
    an_dvel = vs2*an[3]*sin(phi)*cos(phi);
  end
  return an_dvel;
end # function an_dvel
