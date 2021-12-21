using LinearAlgebra;
using Random;
using Statistics;
using Reexport;
using DelimitedFiles;
using DSP;
using Printf;
using LaTeXStrings;
using FFTW;
using SeisMain;


include("cg_aux_mod.jl")           ;
include("conjgrad_mod.jl")         ;
include("isotropic_ray.jl")        ;
include("aniso_ray_tracing.jl")    ;
include("dir_amplitudes.jl")       ;
include("dir_amplitudes_an.jl")    ;
include("refracted_ray.jl")        ;
include("straight_ray_module.jl")  ;
include("RefTransCoef.jl")         ;
include("gsfactor.jl")             ;
include("iso_synth_generator.jl")  ;
include("auxiliars.jl")            ;
include("auxiliars_anisotropy.jl") ;

function ray_tracing_main()
  # ==================================================================
  #
  #  Ray tracing para medio anisotropo, de placas planas y paralelas.
  #
  #  Ver "A rapid and accurate two-point ray tracing method in
  #  horizontally layered velocity model" de Tian y Chen
  #
  # ==================================================================

  modelo = readdlm("dato/modelo_check_4.txt");      # read entire model;
  s      = readdlm("dato/fuente_check_4.txt");     # read source position;
  r      = readdlm("dato/receptores_check_4.txt"); # read receiver positions;
  z   = modelo[:,1]        ; # extract Z form model
  vp  = modelo[:,2]        ; # extract vp form model
  vs  = modelo[:,3]        ; # extract vs form model
  rho = modelo[:,4]*1000.0 ; # extract density from modeel
  an  = modelo[:,5:7]      ; # extract ani form model

  n_s = size(s,1); # number of sources
  n_r = size(r,1); # number of receivers

  tp  = zeros(n_s,n_r); # tiempo de viaje de la onda p
  tsv = zeros(n_s,n_r); # tiempo de viaje de la onda sv
  tsh = zeros(n_s,n_r); # tiempo de viaje de la onda sh

  # Call reytracing
  dt  = 0.001
  ns  = 300
  t0s = 0.0
  f0  = 100.0
  # vp   = [3500.0]                        # P-wave velocity [m/s]
  # vs   = 2400.0                        # S-wave velocity [m/s]
  # rho  = [1.0]                           # Density of the medium [gr/m³]
  # M = [-0.42110287930033197  0.47195575530657674 -0.6386289058383445 ;
  #       0.47195575530657674 -0.3536250535230323   0.9198789025255651 ;
  #      -0.6386289058383445   0.9198789025255651   0.9923820655701969  ]

       M = [0.0  -1.0  -0.0 ;
           -1.0   0.0   0.0 ;
            0.0   0.0   0.0  ]


  geom_in = (dt,ns,t0s,f0,M);

  tpp,tsh,tsv,
  signal,
  sigpp,sigsh,sigsv,amps = aniso_ray_tracing(z,vp,vs,rho,an,s,r,geom_in);

  ppx,ppy,ppz = decouple(sigpp,ns,n_r);
  shx,shy,shz = decouple(sigsh,ns,n_r);
  svx,svy,svz = decouple(sigsv,ns,n_r);


  writedlm("signals_pp/signal_x_1.txt",ppx);
  writedlm("signals_pp/signal_y_1.txt",ppy);
  writedlm("signals_pp/signal_z_1.txt",ppz);

  writedlm("signals_ppsh/signal_x_1.txt",ppx + shx);
  writedlm("signals_ppsh/signal_y_1.txt",ppy + shy);
  writedlm("signals_ppsh/signal_z_1.txt",ppz + shz);

  writedlm("signals_ppshsv/signal_x_1.txt",ppx + shx + svx);
  writedlm("signals_ppshsv/signal_y_1.txt",ppy + shy + svy);
  writedlm("signals_ppshsv/signal_z_1.txt",ppz + shz + svz);
  return
end
