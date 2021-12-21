function refracted_ray(
   e_in::Array{Float64,1}, # arreglo de espesores
   vp_in::Array{Float64,1},# arreglo de velocidades P
   vs_in::Array{Float64,1},# arreglo de velocidades S
   an_in::Array{Float64,2},# arreglo de anisotropias
   s::Array{Float64,1},    # coordenadas de la fuente
   r::Array{Float64,1},    # coordenadas de los receptores
   dp_out::Array{Float64,1}, # INOUT
   dsv_out::Array{Float64,1},# INOUT
   dsh_out::Array{Float64,1})# INOUT

   # !====================================================================
   # !
   # ! Ray tracing para medio anisotropo, de placas planas y paralelas.
   # ! El rayo anisotropo se consigue por medio de conjugate-gradient
   # ! haciendo bending.  Como modelo inicial para el conjugate gradient
   # ! se utiliza la solución para el rayo isotropo, obtenida segun:
   # !
   # ! "A rapid and accurate two-point ray tracing method in
   # ! horizontally layered velocity model" de Tian y Chen
   # !
   # ! La idea es buscar el parámetro de rayo sísmico p tal que el rayo
   # ! sea el de menor tiempo y una la fuente y el receptor fijos
   # ! Realmente no se busca p, sino q que es algo parecido.
   # !
   # !
   # ! IN: e_in    --> espesores de las capas entre la fuente y el receptor.
   # !     vp_in   --> velocidad vertical p
   # !     vs_in   --> velocidad vertical s
   # !     ani_in  --> anisotropia: epsilon delta, gamma
   # !     s       --> coordenadas de la fuente
   # !     r       --> coordenadas del receptor
   # !
   # ! OUT: dp_out  --> intersecciones con las disc del rayo p
   # !      dsh_out --> intersecciones con las disc del rayo sh
   # !      dsv_out --> intersecciones con las disc del rayo sv
   # !      tp      --> tiempo de viaje del rayo p
   # !      tsh     --> tiempo de viaje del rayo sh
   # !      tsv     --> tiempo de viaje del rayo sv
   # !
   # !====================================================================

   # local
   n_e    = size(e_in,1);       # cantidad de espesores
   e      = zeros(Float64,n_e);
   vp     = zeros(Float64,n_e);
   vs     = zeros(Float64,n_e);
   an     = zeros(Float64,size(an_in));
   ds_out = zeros(Float64,n_e-1);

   # si la fuente está por debajo del receptor invierto el modelo  de espesores
   # y velocidades para que tengan coherencia los indices
   if (s[3]>r[3])
      e  = reverse(e_in) ;
      vp = reverse(vp_in);
      vs = reverse(vs_in);
      an[1:n_e,:] = an_in[n_e:-1:1,:]
   else
      e  = e_in;
      vp = vp_in;
      vs = vs_in;
      an = an_in;
   end

   # llamo a las subrutinas que van a calcular el rayo para el medio  isotropo,
   # para el rayo p y para el rayo s
   ftol = 1.0e-4;
   dp_out,tp_out = isotropic_ray(e,vp,n_e,s,r,ftol,dp_out);
   ds_out,ts_out = isotropic_ray(e,vs,n_e,s,r,ftol,ds_out);

   dsv_out = ds_out; # both S phases are the same.
   dsh_out = ds_out;

   # angulo entre fuente y receptor
   H_dist  = norm(r[1:2] - s[1:2]);
   V_dist  = norm(r[3]   - s[3]  );

   if (V_dist == 0.0)
      ang_rs = 90.0;
   else
      op_adj  = H_dist/V_dist  ;
      ang_rs  = atand(op_adj)  ; # angulo de "verticalidad"
   end

   ang_tol = 1.0                  ; # angulo de tolerancia
   # si no es un rayo vertical calculo el anisotropo. Si es casi vertical
   # el anisotropo es igual al isotropo y me quedo con el ya calculado
   if (ang_rs > ang_tol)
      # A partir del rayo obtenido para medio isotropo, ajusto el rayo para medio
      # anisotropo usando bending, ajustando la posición x de las intersecciones
      # del rayo con las disc. (dp_out...) usando conjugate gradient.
      fret,iter,dp_out  = frprmn(dp_out , ftol, s, r, e, vp, vs, an, "vp");
      fret,iter,dsv_out = frprmn(dsv_out, ftol, s, r, e, vp, vs, an, "sv");
      fret,iter,dsh_out = frprmn(dsh_out, ftol, s, r, e, vp, vs, an, "sh");

      # tiempos de cada rayo en el medio anisotropo
      tpp, tpp_m = func(dp_out , e, s, r, vp, vs, an, "vp");
      tsv, tsv_m = func(dsv_out, e, s, r, vp, vs, an, "sv");
      tsh, tsh_m = func(dsh_out, e, s, r, vp, vs, an, "sh");
   end

   return  dp_out,dsv_out,dsh_out,tpp,tsv,tsh,tpp_m,tsv_m,tsh_m;
end #refracted_ray
