function ID_layer(lay,z)
   Z_pos = findall(x->(x > z),lay)[1]
   return Z_pos;
end

function ID_range(s,r)
    if(s > r) # ray comming upwards
        return collect(range(s,r, step = -1));
    else
        return collect(range(s,r, step =  1));
    end
end

function direction(
   zs::Float64,
   zr::Float64)::String
   if zs > zr
      return "up";
   else
      return "down";
   end
end

function rescale(
   Upp::Array{Float64,1},
   Ush::Array{Float64,1},
   Usv::Array{Float64,1},
   App::Float64,
   Ash::Float64,
   Asv::Float64)

   # Upp_amp = (normdata(Upp)[1]) * App;
   # Ush_amp = (normdata(Ush)[1]) * Ash;
   # Usv_amp = (normdata(Usv)[1]) * Asv;

   Upp_amp = Upp * App;
   Ush_amp = Ush * Ash;
   Usv_amp = Usv * Asv;

   cero = zeros(Float64,(size(Upp,1),1));

   Upp_amp = hcat(Upp_amp,cero,cero);
   Ush_amp = hcat(cero,Ush_amp,cero);
   Usv_amp = hcat(cero,cero,Usv_amp);

   return Upp_amp,Ush_amp,Usv_amp;
   #return Upp,Ush,Usv;
end

   
function maxamp(A::Array{Float64,2},B::Array{Float64,2},C::Array{Float64,2})
   AA = findmax(A,dims=1)[1];
   BB = findmax(B,dims=1)[1];
   CC = findmax(C,dims=1)[1];
   return [AA[1],BB[1],CC[1],AA[2],BB[2],CC[2],AA[3],BB[3],CC[3]];
end

function decouple(A::Array{Array{Float64,2},2},ns::Int64,nr::Int64)

   AAx = zeros(Float64,(ns,nr));
   AAy = zeros(Float64,(ns,nr));
   AAz = zeros(Float64,(ns,nr));
   for j = 1:nr
      AAx[:,j] = A[1,j][:,1];
      AAy[:,j] = A[1,j][:,2];
      AAz[:,j] = A[1,j][:,3];
   end
   return AAx,AAy,AAz;
end

function rotation(
   U::Array{Float64,2},
   i::Float64,
   ϕ::Float64,
   direction::String,
   option::String)

   ns = size(U,1);
   if direction == "up" # means source below receiver.
      i = 180.0 - i;         # ray traveling upwards toward surf.
   end
   # page 55. Pita Ruiz.
   if option == "xyz2rϕi"
      # changing coordinates from (x,y,z), to (r,i,ϕ)
      x2sph = [sind(i)*cosd(ϕ), -sind(ϕ),  cosd(i)*cosd(ϕ)];
      y2sph = [sind(i)*sind(ϕ),  cosd(ϕ),  cosd(i)*sind(ϕ)];
      z2sph = [cosd(i)        ,   0.0   , -sind(i)        ];
      U_x2shp = zeros(Float64,(ns,3))
      U_y2shp = zeros(Float64,(ns,3))
      U_z2shp = zeros(Float64,(ns,3))
      for j = 1:ns
         U_x2shp[j,:] = U[j,1] * x2sph;
         U_y2shp[j,:] = U[j,2] * y2sph;
         U_z2shp[j,:] = U[j,3] * z2sph;
      end
      U_shp = U_x2shp + U_y2shp + U_z2shp;
      return U_shp;
   end
   if option == "rϕi2xyz"
      # changing coordinates from (r,i,ϕ) to (x,y,z)
      r2xyz = [sind(i)*cosd(ϕ), sind(i)*sind(ϕ),   cosd(i)];
      ϕ2xyz = [       -sind(ϕ),         cosd(ϕ),    0.0   ];
      i2xyz = [cosd(i)*cosd(ϕ), cosd(i)*sind(ϕ),  -sind(i)];
      U_r2xyz = zeros(Float64,(ns,3))
      U_ϕ2xyz = zeros(Float64,(ns,3))
      U_i2xyz = zeros(Float64,(ns,3))
      for j = 1:ns
         U_r2xyz[j,:] = U[j,1] * r2xyz;
         U_ϕ2xyz[j,:] = U[j,2] * ϕ2xyz;
         U_i2xyz[j,:] = U[j,3] * i2xyz;
      end
      U_xyz = U_r2xyz + U_ϕ2xyz + U_i2xyz;
      return U_xyz;
   end
end


function wave_directionals(i::Float64,ϕ::Float64)
   x = [1,0,0]; # versor x
   y = [0,1,0]; # versor y
   z = [0,0,1]; # versor z
   # eq.4.88 Aki-Richards: Projections for SH and SV
   p =  cosd(i)*cosd(ϕ)*x + cosd(i)*sind(ϕ)*y - sind(i)*z;
   ϕ = -        sind(ϕ)*x +         cosd(ϕ)*y +   0.0  *z;
   p_norm = norm(p);
   ϕ_norm = norm(ϕ);
   return p/p_norm,ϕ/ϕ_norm;
end

function take_off_angle(angles::Array{Float64,1})
   return angles[end]; # the one reaching the receiver
end

function take_in_angle(angles::Array{Float64,1})
   return angles[1]; # the one taking off the source
end

function take_off_azimuth(vect::Array{Float64,1})
   x = vect[1]               ;
   y = vect[2]               ;

   # azimuth [0,360) increase clockwise from north(x) to east(y)
   # open cuadrant, without axis
   cuad_I   = (x > 0.0) && (y > 0.0);
   cuad_II  = (x < 0.0) && (y > 0.0);
   cuad_III = (x < 0.0) && (y < 0.0);
   cuad_IV  = (x > 0.0) && (y < 0.0);
   # axis, just in case
   axis_x_pos = (x > 0.0) && (y == 0.0);
   axis_x_neg = (x < 0.0) && (y == 0.0);
   axis_y_pos = (x == 0.0) && (y > 0.0);
   axis_y_neg = (x == 0.0) && (y < 0.0);
   # origin
   origin     = (x == 0) && (y == 0)

   if     origin
      azimuth = 0.0    ;
      return azimuth   ;
   else
      hip   = norm(vect[1:2])    ;
      yasin = abs(asind(y/hip)) ;
      if axis_x_pos
         azimuth = 0.0    ;
         return azimuth
      elseif axis_x_neg
         aziumth  = 180.0 ;
         return azimuth
      elseif axis_y_pos
         azimuth = 90.0   ;
         return azimuth
      elseif axis_y_neg
         azimuth = 270.0  ;
         return azimuth
      elseif cuad_I
         azimuth = yasin  ;
         return azimuth
      elseif cuad_II
         azimuth = 180.0 - yasin;
         return azimuth
      elseif cuad_III
         azimuth = 180.0 + yasin;
         return azimuth
      elseif cuad_IV
         azimuth = 360.0 - yasin;
         return azimuth
      end
   end
end

function trayectories(
   s::Array{Float64,1},
   r::Array{Float64,1},
   cpp_in::Array{Float64,2},
   csh_in::Array{Float64,2},
   csv_in::Array{Float64,2},
   nij::Int64)
   cpp = vcat(s',cpp_in[1:nij,:],r');
   csh = vcat(s',csh_in[1:nij,:],r');
   csv = vcat(s',csv_in[1:nij,:],r');
   return cpp,csh,csv;
end
function trayectories(
   s::Array{Float64,1},
   r::Array{Float64,1})
   pp_ray = vcat(s',r');
   sh_ray = vcat(s',r');
   sv_ray = vcat(s',r');
   return pp_ray,sh_ray,sv_ray
end

function angles(
   source::Array{Float64,1},
   receiv::Array{Float64,1},
   cpp_in,csh_in,csv_in,nij)
   nij2    = nij + 2; # number of coord. (intersections + Srce. + Rec.)

   xs,ys,zs = source[1],source[2],source[3];
   xr,yr,zr = receiv[1],receiv[2],receiv[3];

   cpp = zeros(Float64,nij2,3);
   csh = zeros(Float64,nij2,3);
   csv = zeros(Float64,nij2,3);
   cpp[1,:]       = source;
   csh[1,:]       = source;
   csv[1,:]       = source;
   cpp[2:end-1,:] = cpp_in[1:nij,:];
   csh[2:end-1,:] = csh_in[1:nij,:];
   csv[2:end-1,:] = csv_in[1:nij,:];
   cpp[end,:]     = receiv;
   csh[end,:]     = receiv;
   csv[end,:]     = receiv;

   if zs > zr   # source below receiver. Ray traveling upwards
      zz = [0.0,0.0,1.0];   # normal pointing up (z grows with depth)
   else
      zz = [0.0,0.0,-1.0] ; # invert normal
   end

   grad = 180.0/pi;
   # the amount of angles is # intersec. + t-off angle: (nij + 1)
   pp_ang = zeros(Float64,nij + 1) #int*2);
   sh_ang = zeros(Float64,nij + 1) #int*2);
   sv_ang = zeros(Float64,nij + 1) #int*2);
   for i = 1: nij+1
      pp_vec = cpp[i,:] - cpp[i+1,:];
      sh_vec = csh[i,:] - csh[i+1,:];
      sv_vec = csv[i,:] - csv[i+1,:];
      d_pp_n = dot(pp_vec,zz)/norm(pp_vec);
      d_sh_n = dot(sh_vec,zz)/norm(sh_vec);
      d_sv_n = dot(sv_vec,zz)/norm(sv_vec);
      pp_ang[i] = acos(d_pp_n)*grad;
      sh_ang[i] = acos(d_sh_n)*grad;
      sv_ang[i] = acos(d_sv_n)*grad;
   end
   return pp_ang,sh_ang,sv_ang;
end

function angles(source::Array{Float64,1},receiv::Array{Float64,1})

   xs,ys,zs = source[1],source[2],source[3];
   xr,yr,zr = receiv[1],receiv[2],receiv[3];

   if zs > zr   # source below receiver. Ray traveling upwards
      # recall: z growth with depth. So normal points down.
      zz = [0.0,0.0,1.0];    # this points up
   else
      zz = [0.0,0.0,-1.0] ;  # this points down
   end
   grad = 180.0/pi;
   # Source at same interface as Receiver
   pp_ang = zeros(Float64,1);
   sh_ang = zeros(Float64,1);
   sv_ang = zeros(Float64,1);
   # 1-I  Source same depth as Receiver. Ray traveling horizontal.
   if zs == zr
      pp_ang[1] = 90.0;
      sh_ang[1] = 90.0;
      sv_ang[1] = 90.0;

      # 1-II Source below or above Receiver. Ray traveling vertical.
   elseif receiv[1:2] == source[1:2]
      pp_ang[1] = 0.0;
      sh_ang[1] = 0.0;
      sv_ang[1] = 0.0;

      # 1-III Source different depth as Receiver. Ray traveling updwds or downds.
   else
      vec =  source - receiv;
      d_vec_n = dot(vec,zz)/norm(vec);
      pp_ang[1] = acos(d_vec_n)*grad;
      sh_ang[1] = pp_ang[1];
      sv_ang[1] = pp_ang[1];
   end
   return  pp_ang,sh_ang,sv_ang;
end

function normdata(D::Array{Float64,2})

   Dmax = maximum(abs.(D));
   Dnorm = D/Dmax;
   return Dnorm, Dmax;
end

function normdata(D::Array{Float64,1})

   Dmax = maximum(abs.(D));
   Dnorm = D/Dmax;
   return Dnorm, Dmax;
end

# function angles(p_ray,sh_ray,sv_ray)
#    source = p_ray[1,:]  ;
#    rec    = p_ray[end,:];
#
#    p_rows = size(p_ray,1);
#    int    = p_rows - 2   ; # number of rows of trayectories "-2"
#    # indicates number of intersections.
#    if source[2] > rec[2]   # source below receiver. Ray traveling upwards
#       zz = [0.0,-1.0];     # normal pointing downwards (z grows with depth)
#    else
#       zz = [0.0,1.0] ;     # invert normal (pointing upwards)
#    end
#    grad = 180.0/pi;
#
#    # two scenarios:
#
#    # 1 - Source at same interface as Receiver
#    if (int==0) # Source at same interface as Receiver
#       pp_ang = zeros(Float64,int+1)
#       sh_ang = zeros(Float64,int+1)
#       sv_ang = zeros(Float64,int+1)
#
#       # 1-I  Source same depth as Receiver. Ray traveling horizontal.
#       if source[2] == rec[2]
#          pp_ang[1] = 90.0;
#          sh_ang[1] = 90.0;
#          sv_ang[1] = 90.0;
#
#          # 1-II Source below or above Receiver. Ray traveling vertical.
#       elseif source[1] == rec[1]
#          pp_ang[1] = 0.0;
#          sh_ang[1] = 0.0;
#          sv_ang[1] = 0.0;
#
#          # 1-III Source different depth as Receiver. Ray traveling updwds or downds.
#       else
#          vec =  source - rec;
#          d_vec_n = dot(vec,zz)/norm(vec);
#          pp_ang[1] = acos(d_vec_n)*grad;
#          sh_ang[1] = pp_ang[1];
#          sv_ang[1] = pp_ang[1];
#       end
#       return  pp_ang,sh_ang,sv_ang;
#
#       # 2 -  Source at different interface as Receiver
#    else
#       pp_ang = zeros(Float64,p_rows-1)#int*2);
#       sh_ang = zeros(Float64,p_rows-1)#int*2);
#       sv_ang = zeros(Float64,p_rows-1)#int*2);
#       for i = 1:p_rows-1
#          pp_vec =  p_ray[i+1,:] - p_ray[i,:] ;
#          sh_vec = sh_ray[i+1,:] - sh_ray[i,:];
#          sv_vec = sv_ray[i+1,:] - sv_ray[i,:];
#          d_pp_n = dot(pp_vec,zz)/norm(pp_vec);
#          d_sh_n = dot(sh_vec,zz)/norm(sh_vec);
#          d_sv_n = dot(sv_vec,zz)/norm(sv_vec);
#          pp_ang[i] = acos(d_pp_n)*grad;
#          sh_ang[i] = acos(d_sh_n)*grad;
#          sv_ang[i] = acos(d_sv_n)*grad;
#       end
#       return pp_ang,sh_ang,sv_ang;
#    end
#
# end



# function trayectories2D(s,r,cpp,csh,csv,nij)
#
#    nij2 = nij + 2; # number of coord. (intersections + Srce. + Rec.)
#
#    pp_ray = zeros(Float64,nij2,2);
#    sh_ray = zeros(Float64,nij2,2);
#    sv_ray = zeros(Float64,nij2,2);
#
#    rad_pp = sqrt.(cpp[1:nij,1].^2 + cpp[1:nij,2].^2);
#    rad_sh = sqrt.(csh[1:nij,1].^2 + csh[1:nij,2].^2);
#    rad_sv = sqrt.(csv[1:nij,1].^2 + csv[1:nij,2].^2);
#
#    aux1 = [norm(s[1:2]),s[3]];
#    aux2 = [norm(r[1:2]),r[3]];
#
#    pp_ray[1,:]       = aux1;
#    sh_ray[1,:]       = aux1;
#    sv_ray[1,:]       = aux1;
#    pp_ray[2:end-1,:] = hcat(rad_pp,cpp[1:nij,3])
#    sh_ray[2:end-1,:] = hcat(rad_sh,csh[1:nij,3])
#    sv_ray[2:end-1,:] = hcat(rad_sv,csv[1:nij,3])
#    pp_ray[end,:]     = aux2;
#    sh_ray[end,:]     = aux2;
#    sv_ray[end,:]     = aux2;
#
#    return pp_ray,sh_ray,sv_ray
# end
