include("cg_aux_mod.jl")
function frprmn(p::Array{Float64,1},
  ftol::Float64,
  ss::Array{Float64,1},
  rr::Array{Float64,1},
  ee::Array{Float64,1},
  vp::Array{Float64,1},
  vs::Array{Float64,1},
  an::Array{Float64,2},
  ray::String)

  ITMAX = 200;
  EPS   = 1.0e-10;
  # Given  a  starting  point  p  that is  a  vector  of  length  N,
  # Fletcher-Reeves-Polak-Ribiere  minimization  is performed  on  a
  # function func,  using its  gradient as  calculated by  a routine
  # dfunc. The convergence tolerance on  the function value is input
  # as  ftol.  Returned  quantities  are  p  (the  location  of  the
  # minimum), iter  (the number of iterations  that were performed),
  # and fret (the minimum value of the function). The routine linmin
  # is called  to perform line minimizations.   Parameters: ITMAX is
  # the maximum allowed number of  iterations; EPS is a small number
  # to  rectify  the special  case  of  converging to  exactly  zero
  # function value.

  g  = zeros(Float64,size(p));
  h  = zeros(Float64,size(p));
  xi = zeros(Float64,size(p));

  fp = func(p,ee,ss,rr,vp,vs,an,ray)[1];     # Initializations.
  xi = dfunc(p,ee,ss,rr,vp,vs,an,ray)

  g  = -xi;
  h  = g;
  xi = h;

  for its = 1:ITMAX        # Loop over iterations.

    iter = its;
    p,xi,fret = linmin(p,xi,ee,ss,rr,vp,vs,an,ray)   #Next statement is the normal return:
    if (2.0*abs(fret-fp) <= ftol*(abs(fret)+abs(fp)+EPS))
      return fret,iter,p;
    end
    fp  = fret;
    xi  = dfunc(p,ee,ss,rr,vp,vs,an,ray);
    gg  = dot(g,g)    ;
    dgg = dot(xi,xi)  ;   # This statement for Fletcher-Reeves.
    dgg = dot(xi+g,xi);   # This statement for Polak-Ribiere.
    if (gg == 0.0)        # Unlikely. If gradient is exactly zero
      return fret,iter,p;
    end
    gam = dgg/gg  ;     # then we are already done.
    g  = -xi      ;
    h  = g + gam*h;
    xi = h        ;
  end
  warning("frprmn: maximum iterations exceeded")
  return fret,iter,p;
end #SUBROUTINE frprmn

function linmin(
  p::Array{Float64,1},
  xi::Array{Float64,1},
  ee::Array{Float64,1},
  ss::Array{Float64,1},
  rr::Array{Float64,1},
  vp::Array{Float64,1},
  vs::Array{Float64,1},
  an::Array{Float64,2},
  ray::String)

  # REAL(DP), INTENT(OUT) :: fret
  # REAL(DP), DIMENSION(:), INTENT(INOUT) :: p,xi
  # real(kind=wp), dimension(:), intent(IN) :: ee,vp,vs,ss,rr
  # real(kind=wp), dimension(:,:), intent(in)::an
  # character(len=2), intent(in)::ray

  # Given an N-dimensional  point p and an  N -dimensional direction
  # xi, both  vectors of length N,  moves and resets p  to where the
  # fixed-name function func takes on  a minimum along the direction
  # xi from  p, and  replaces xi by  the actual  vector displacement
  # that p was moved. Also returns as  fret the value of func at the
  # returned  location  p.  This  is actually  all  accomplished  by
  # calling  the routines  mnbrak and  brent.  Parameter:  Tolerance
  # passed to brent.
  TOL = 1.0e-4;
  ax = 0.0;                   # Initial guess for brackets.
  xx = 1.0;
  ax,xx,bx,fa,fx,fb,p,xi =  mnbrak(ax,xx,p,xi,ee,ss,rr,vp,vs,an,ray)
  fret,p,xi,xmin         = brent(ax,xx,bx,p,xi,TOL,ee,ss,rr,vp,vs,an,ray)
  xi = xmin*xi;               # Construct the vector results to return.
  p = p+xi    ;
  return p,xi,fret;
end # linmin

#-------------------------------------------------------------------

function brent(
  ax::Float64,
  bx::Float64,
  cx::Float64,
  pp::Array{Float64,1},
  xxi::Array{Float64,1},
  tol::Float64,
  ee::Array{Float64,1},
  ss::Array{Float64,1},
  rr::Array{Float64,1},
  vp::Array{Float64,1},
  vs::Array{Float64,1},
  an::Array{Float64,2},
  ray::String)

  ITMAX = 100           ;
  CGOLD = 0.3819660     ;
  ZEPS  = 1.0e-3*eps(ax);
  # Given  a  function  func,  and given  a  bracketing  triplet  of
  # abscissas ax,  bx, cx (such  that bx is  between ax and  cx, and
  # func(bx) is less than both  func(ax) and func(cx)), this routine
  # isolates  the minimum  to a  fractional precision  of about  tol
  # using Brent's method.  #The abscissa of the  minimum is returned
  # as xmin,  and the minimum  function value is returned  as brent,
  # the returned function value.  Parameters: Maximum allowed number
  # of iterations;  golden ratio; and  a small number  that protects
  # against trying to achieve fractional accuracy for a minimum that
  # happens to be exactly zero.

  brent = 0.0;
  a = min(ax,cx); # a and b must be in ascending order, though
  b = max(ax,cx); # the input abscissas need not be.
  v = bx ;         #  Initializations...
  w = v  ;
  x = v  ;
  e = 0.0;        # This will be the distance moved on the step
  fx = func(pp+x*xxi,ee,ss,rr,vp,vs,an,ray)[1]; #before last.
  fv = fx;
  fw = fx;
  for iter = 1:ITMAX # Main program loop.
    xm   = 0.5*(a+b);
    tol1 = tol*abs(x) + ZEPS;
    tol2 = 2.0*tol1;
    if (abs(x-xm) <= (tol2-0.5*(b-a)));  #Test for done here.
      xmin  = x ; # Arrive here ready to exit with best values.
      brent = fx;
      return brent,pp,xxi,xmin;
    end
    if (abs(e) > tol1)  # Construct a trial parabolic fit.
      r = (x-w)*(fx-fv)
      q = (x-v)*(fx-fw)
      p = (x-v)*q-(x-w)*r
      q = 2.0*(q-r)
      if (q > 0.0)
        p = -p
      end
      q = abs(q)
      etemp = e
      if ((abs(p) >= abs(0.5*q*etemp)) || (p <= q*(a-x))) || (p >= q*(b-x))
        # The above conditions determine the acceptability of the
        # parabolic fit. Here it is not o.k., so we take the golden
        # section step into the larger of the two segments.
        e = merge_g(a-x,b-x,x>=xm)
        d = CGOLD*e;
      else     # Take the parabolic step.
        d = p/q
        u = x+d
        if (u-a < tol2 || b-u < tol2)
          d = sign(xm-x)*abs(tol1);
        end
      end
    else        # Take the golden section step into the larger
      e = merge_g(a-x,b-x, x >= xm ) # of the two segments.
      d = CGOLD*e
    end

    u = merge_g(x+d,x+sign(d)*abs(tol1), abs(d) >= tol1 )
    # Arrive here with d computed either from parabolic fit,
    # or else from golden section.
    fu = func(pp+u*xxi,ee,ss,rr,vp,vs,an,ray)[1];
    # This is the one function evaluation per iteration.
    if (fu <= fx)       # Now we have to decide what to do with our
      if (u >= x)      # function evaluation. Housekeeping follows:
        a = x
      else
        b = x
      end
      v ,w , x = shft(w , x , u );
      fv,fw,fx = shft(fw, fx, fu);
    else
      if (u < x)
        a = u
      else
        b = u
      end
      if (fu <= fw || w == x)
        v  = w
        fv = fw
        w  = u
        fw = fu
      elseif ((fu <= fv || v == x) || v == w)
        v  = u
        fv = fu
      end
    end
  end #Done with housekeeping. Back for another iteration.
  warning("brent: exceed maximum iterations")
  return brent,pp,xxi,xmin;
end #brent

function shft(
  b::Float64,
  c::Float64,
  d::Float64)

  a = b;
  b = c;
  c = d;
  return a,b,c
end

function mnbrak(
  ax::Float64,
  bx::Float64,
  pp::Array{Float64,1},
  xxi::Array{Float64,1},
  ee::Array{Float64,1},
  ss::Array{Float64,1},
  rr::Array{Float64,1},
  vp::Array{Float64,1},
  vs::Array{Float64,1},
  an::Array{Float64,2},
  ray::String)

  GOLD   = 1.618034;
  GLIMIT = 100.0;
  TINY   = 1.0e-20;
  # Given a function func, and  given distinct initial points ax and
  # bx, this routine searches in  the downhill direction (defined by
  # the function as evaluated at the initial points) and returns new
  # points ax, bx,  cx that bracket a minimum of  the function. Also
  # returned are  the function values  at the three points,  fa, fb,
  # and  fc.   Parameters:  GOLD  is  the  default  ratio  by  which
  # successive  intervals  are  magnified;  GLIMIT  is  the  maximum
  # magnification allowed for a parabolic-fit step.

  fa = func(pp+ax*xxi,ee,ss,rr,vp,vs,an,ray)[1];
  fb = func(pp+bx*xxi,ee,ss,rr,vp,vs,an,ray)[1];
  if (fb > fa)             # Switch roles of a and b so that we
    ax,bx = swap(ax,bx)       # can go downhill in the direction
    fa,fb = swap(fa,fb)       # from a to b.
  end
  cx = bx + GOLD*(bx-ax)           # First guess for c.
  fc = func(pp+cx*xxi,ee,ss,rr,vp,vs,an,ray)[1];
  k = 0
  while true   # Do-while-loop: Keep returning here
    k = k + 1
    if (fb < fc)
      return ax,bx,cx,fa,fb,fc,pp,xxi;
    end         # until we bracket.
    # Compute u by parabolic extrapolation from a, b, c.
    # TINY is used to prevent any possible division by zero.
    r = (bx-ax)*(fb-fc)
    q = (bx-cx)*(fb-fa)
    u_aux1 = sign(q-r)              ;
    u_aux2 = abs(max(abs(q-r),TINY));
    u_aux3 = 2*u_aux1*u_aux2          ;
    u = bx-((bx-cx)*q-(bx-ax)*r)/u_aux3;
    ulim = bx + GLIMIT*(cx-bx)
    # We won't go farther than this. Test various possibilities:
    if ((bx-u)*(u-cx) > 0.0)    # Parabolic u is between b and c: try it
      fu = func(pp+u*xxi,ee,ss,rr,vp,vs,an,ray)[1]
      if (fu < fc)             # Got a minimum between b and c.
        ax = bx;
        fa = fb;
        bx = u ;
        fb = fu;
        return ax,bx,cx,fa,fb,fc,pp,xxi;
      elseif (fu > fb)        # Got a minimum between a and u.
        cx = u
        fc = fu
        return ax,bx,cx,fa,fb,fc,pp,xxi;
      end
      u  = cx + GOLD*(cx-bx)            # Parabolic fit was no use. Use default
      fu = func(pp+u*xxi,ee,ss,rr,vp,vs,an,ray)[1];              # magnification.
    elseif ((cx-u)*(u-ulim) > 0.0)  # Parabolic fit is between c and its al-
      fu = func(pp+u*xxi,ee,ss,rr,vp,vs,an,ray)[1];              # lowed limit.
      if (fu < fc)
        bx = cx;
        cx = u;
        u  = cx + GOLD*(cx-bx);
        shft_in = func(pp+u*xxi,ee,ss,rr,vp,vs,an,ray)[1];
        fb,fc,fu = shft(fc,fu,shft_in)
      end
    elseif ((u-ulim)*(ulim-cx) >= 0.0)  # Limit parabolic u to maximum allowed value.
      u = ulim;
      fu = func(pp+u*xxi,ee,ss,rr,vp,vs,an,ray)[1];
    else
      u  = cx + GOLD*(cx-bx)                     # magnification.
      fu = func(pp+u*xxi,ee,ss,rr,vp,vs,an,ray)[1];
    end
    ax,bx,cx = shft(bx,cx,u)
    fa,fb,fc = shft(fb,fc,fu)                   # Eliminate oldest point and continue.
  end
  return ax,bx,cx,fa,fb,fc,pp,xxi;
end

function swap(
  a::Float64,
  b::Float64)
  dum = a   ;
  a   = b   ;
  b   = dum ;
  return a,b;
end # swap

function merge_g(
  a::Float64,
  b::Float64,
  cond::Bool)
  if cond
    return a;
  else
    return b;
  end
end
