from scipy.signal import butter, lfilter
#from scipy.stats import norm
from numpy import linalg as LA
from matplotlib import rc, font_manager

import numpy as np
import os

import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import matplotlib.colors as colors

rc('text', usetex=True)
plt.rc('font', family='serif')

path1 = '/home/gbrunini/Desktop/Julia_raytracing/ray_tracing_jul/data'

fdr_fig_    = "gathers/"            ;  # where to save the figures

###### carga de datos de errores
datax = np.loadtxt("signal_x_1.txt");
datay = np.loadtxt("signal_y_1.txt");
dataz = np.loadtxt("signal_z_1.txt");
# datax_n = np.loadtxt("signal_x_n.txt");
# datay_n = np.loadtxt("signal_y_n.txt");
# dataz_n = np.loadtxt("signal_z_n.txt");
# datax_d = np.loadtxt("signal_x_d.txt");
# datay_d = np.loadtxt("signal_y_d.txt");
# dataz_d = np.loadtxt("signal_z_d.txt");

# Deltat = 0.001   ;  # sampling interval [s]
# rho    = 2.70    ;  # 1.0;    # density of the medium [gr/m3]
t0s    = 0.000   ;  # -0.037; # t0 of source [s]
vp     = 5000.0  ;  # 3500.0; # P-wave velocity [m/s]
vs     = 2000.0  ;  # 2400.0; # S-wave velocity [m/s]

##########################################
##########################################
[nsx,ntx] = np.shape(datax.reshape(700,1))                ; # dimension of x-data.
print()
print("dimension of x-data")
print("number of samples   : ", nsx);
print("number of traces    : ", ntx);
[nsy,nty] = np.shape(datay.reshape(700,1))                ; # dimension of y-data.
print()
print("dimension of y-data")
print("number of samples   : ", nsy);
print("number of traces    : ", nty);
[nsz,ntz] = np.shape(dataz.reshape(700,1))                ; # dimension of z-data.
print()
print("dimension of z-data")
print("number of samples   : ", nsz);
print("number of traces    : ", ntz);


# porc_norm = 0.99    # porcentage of normalize data.
#
# maximx   = np.amax(abs(x_c))             ; # extract data maximum x normalization.
# maximy   = np.amax(abs(y_c))             ; # extract data maximum x normalization.
# maximz   = np.amax(abs(z_c))             ; # extract data maximum x normalization.
# max      = np.amax([maximx,maximy,maximz]) ;
#
# datax    = x_c/max                    ; # normalize
# datay    = y_c/max                    ; # normalize
# dataz    = z_c/max                    ; # normalize
#
# datax    = datax*porc_norm                ; # take a porcentage of normalize data.
# datay    = datay*porc_norm                ; # take a porcentage of normalize data.
# dataz    = dataz*porc_norm                ; # take a porcentage of normalize data.

ts = 0.0 ; # time begins in 0.0 ms
te = 700.0#nsx  ; # time ends at m ms

# Because data is normalize to 1 (or porc of 1)
# each trace occupies 2 vertical units (maximum).

# we need at least 2*nt vertical space in ylim
first_trace = 0        ;  # first trace is positioned at zero.
last_trace  = 2*(ntx-1);  # (nt-1) : number of intervals between nt traces
# 2*     : time 2, vertical space.
minY  = (first_trace - 1) - 0.2 ;
maxY  = (last_trace  + 1) + 0.2 ;

lw     = 1.5;  # linewidth
ttw    = 0.7   # linewidth for travel time lines
ls     = '-';  # linestyle
alp    = 0.5;  # alpha ransparency for line
alp2   = 0.8;  # alpha transparency for arriva time lines
xlab_s = 15;   # x-label size
ylab_s = 15;   # y-label size
tit_s  = 15;   # title size
tic_s  = 12;   # tics size

scale         = 3.0; # scales data (1.0: 0 scaling). Warn: makes data bigger.
yticksnum     = np.linspace(first_trace,last_trace,ntx);                # y tick numbers
ytickslabels  = np.tile(((np.linspace(1,ntx,ntx)).astype(np.int64)),1); # y tick labels

# # circle dimension
# xc1 = 0.0;
# yc1 = 100.0;
# zc1 = [200.0,160.0,120.0,80.0,40.0];
# xc2 = 50.0;
# yc2 = 0.0;
# zc2 = [200.0,160.0,120.0,80.0,40.0]
# radius = 108.0;


# # axis limits for scatter plots
# xminlim = -130.0 ;
# xmaxlim =  170.0 ;
# yminlim = -120.0 ;
# ymaxlim =  220.0 ;
# zminlim =  -10.0   ;
# zmaxlim =  260.0 ;
#
# # measures for scatter legends and markers, etc..
# markersize_1     = 8.0 ; # tamano de los receptores
# markersize_2     = 16.0 ; # tamano de los receptores
# scaterlegendsize = 12.0;
# scatterpointnum  = 1   ;
# scatermarkersize = 2.0 ;
#
#
# slope = (yc2-yc1)/(xc2-xc1);
# xline = np.linspace(-130,200,100);
# yline = slope*(xline-xc1) + yc1;
#
# zvline  = np.linspace(0,250,100);
# zhline1 = np.linspace(xc1,xc1,100);
# zhline2 = np.linspace(xc2,xc2,100);

fig2,ax2  = plt.subplots(figsize=(12,8));

for i in range(ntx):
        gath_X = scale*datax[:];
        gath_Y = scale*datay[:];
        gath_Z = scale*dataz[:];

        # gath_X_n = scale*datax_n[:,i];
        # gath_Y_n = scale*datay_n[:,i];
        # gath_Z_n = scale*dataz_n[:,i];
        # gath_X_d = scale*datax_d[:,i];
        # gath_Y_d = scale*datay_d[:,i];
        # gath_Z_d = scale*dataz_d[:,i];

        ################# gather
        # wlen  = 24.0                 ; # window length for trav time
        # trec  = np.arange(0,2*ntx,2) ;
        # tp1   = np.flipud(tp[ig,:])  ;
        # tp2   = tp1 + wlen           ;
        # ts1   = tp1*(vp/vs)          ;
        # ts2   = ts1 + wlen           ;

        # x-plot
        plt.subplot(1,3,1)
        plt.plot(gath_X  + 2*i,
         linestyle = ls,
         linewidth = lw,
         color     = 'black',
         alpha     = alp)
        # plt.gca().plot(tp1 , trec, lw = ttw,linestyle='--',color="red",alpha=alp2);
        # plt.gca().plot(tp2 , trec, lw = ttw,linestyle='--',color="red",alpha=alp2);
        # plt.gca().plot(ts1 , trec, lw = ttw,linestyle='--',color="green",alpha=alp2);
        # plt.gca().plot(ts2 , trec, lw = ttw,linestyle='--',color="green",alpha=alp2);
        left,right = plt.xlim()
        plt.xlim(left  = ts)  # adjust the left leaving right unchanged
        plt.xlim(right = te)  # adjust the right leaving left unchanged
        plt.tick_params(
         axis        = "y",
         width       = 1,
         length      = 2.5,
         direction   = "in",
         color       = "black",
         pad         = 2.5,
         labelsize   = tic_s,
         labelcolor  = "black",
         colors      = "black",
         zorder      = 20,
         bottom      = "on", top      = "off", left      = "on", right      = "off",
         labelbottom = "on", labeltop = "off", labelleft = "on", labelright = "off");
        bottom, top = plt.ylim()
        plt.ylim(bottom = minY)
        plt.ylim(top    = maxY)
        plt.yticks(yticksnum,ytickslabels);
        plt.tick_params(
         axis        = "x",
         width       = 1,
         length      = 2.5,
         direction   = "in",
         color       = "black",
         pad         = 2.5,
         labelsize   = tic_s,
         labelcolor  = "black",
         colors      = "black",
         zorder      = 20,
         bottom      = "on", top      = "off", left      = "on", right      = "off",
         labelbottom = "on", labeltop = "off", labelleft = "on", labelright = "off");
        plt.ylabel(r"receiver",fontsize   = ylab_s)
        plt.title(r"x-component",fontsize = tit_s)

        ################# gather

        # y-plot
        plt.subplot(1,3,2)
        plt.plot(gath_Y  + 2*i,
         linestyle = ls,
         linewidth = lw,
         color     = 'black',
         alpha     = alp)
        # plt.gca().plot(tp1 , trec, lw = ttw,linestyle='--',color="red",alpha=alp2);
        # plt.gca().plot(tp2 , trec, lw = ttw,linestyle='--',color="red",alpha=alp2);
        # plt.gca().plot(ts1 , trec, lw = ttw,linestyle='--',color="green",alpha=alp2);
        # plt.gca().plot(ts2 , trec, lw = ttw,linestyle='--',color="green",alpha=alp2);
        left,right = plt.xlim()
        plt.xlim(left  = ts)  # adjust the left leaving right unchanged
        plt.xlim(right = te)  # adjust the right leaving left unchanged
        plt.tick_params(
         axis        = "y",
         width       = 1,
         length      = 2.5,
         direction   = "in",
         color       = "black",
         pad         = 2.5,
         labelsize   = tic_s,
         labelcolor  = "black",
         colors      = "black",
         zorder      = 20,
         bottom      = "on", top      = "off", left      = "on", right      = "off",
         labelbottom = "on", labeltop = "off", labelleft = "on", labelright = "off");
        bottom, top = plt.ylim()
        plt.ylim(bottom = minY)
        plt.ylim(top    = maxY)
        plt.yticks(yticksnum,ytickslabels);
        plt.tick_params(
         axis        = "x",
         width       = 1,
         length      = 2.5,
         direction   = "in",
         color       = "black",
         pad         = 2.5,
         labelsize   = tic_s,
         labelcolor  = "black",
         colors      = "black",
         zorder      = 20,
         bottom      = "on", top      = "off", left      = "on", right      = "off",
         labelbottom = "on", labeltop = "off", labelleft = "on", labelright = "off");
        plt.ylabel(r"receiver",fontsize   = ylab_s)
        plt.title(r"y-component",fontsize = tit_s)

        ################# gather
        # z-plot
        plt.subplot(1,3,3)
        plt.plot(gath_Z  + 2*i,
         linestyle = ls,
         linewidth = lw,
         color     = 'black',
         alpha     = alp)
        # plt.gca().plot(tp1 , trec, lw = ttw,linestyle='--',color="red",alpha=alp2);
        # plt.gca().plot(tp2 , trec, lw = ttw,linestyle='--',color="red",alpha=alp2);
        # plt.gca().plot(ts1 , trec, lw = ttw,linestyle='--',color="green",alpha=alp2);
        # plt.gca().plot(ts2 , trec, lw = ttw,linestyle='--',color="green",alpha=alp2);
        left,right = plt.xlim()
        plt.xlim(left  = ts)  # adjust the left leaving right unchanged
        plt.xlim(right = te)  # adjust the right leaving left unchanged
        plt.tick_params(
         axis        = "y",
         width       = 1,
         length      = 2.5,
         direction   = "in",
         color       = "black",
         pad         = 2.5,
         labelsize   = tic_s,
         labelcolor  = "black",
         colors      = "black",
         zorder      = 20,
         bottom      = "on", top      = "off", left      = "on", right      = "off",
         labelbottom = "on", labeltop = "off", labelleft = "on", labelright = "off");
        bottom, top = plt.ylim()
        plt.ylim(bottom = minY)
        plt.ylim(top    = maxY)
        plt.yticks(yticksnum,ytickslabels);
        plt.tick_params(
         axis        = "x",
         width       = 1,
         length      = 2.5,
         direction   = "in",
         color       = "black",
         pad         = 2.5,
         labelsize   = tic_s,
         labelcolor  = "black",
         colors      = "black",
         zorder      = 20,
         bottom      = "on", top      = "off", left      = "on", right      = "off",
         labelbottom = "on", labeltop = "off", labelleft = "on", labelright = "off");
        plt.xlabel(r"time (ms)",fontsize  = xlab_s)
        plt.ylabel(r"receiver",fontsize   = ylab_s)
        plt.title(r"z-component",fontsize = tit_s)

# ########################### FILA 2 ###################################
# ########################### FILA 2 ###################################
#         ################# gather
#         # x - plot
#         plt.subplot(2,3,4)
#         plt.plot(gath_X_n  + 2*i,
#          linestyle = ls,
#          linewidth = lw,
#          color     = 'red',
#          alpha     = alp)
#         plt.plot(gath_X_d  + 2*i,
#          linestyle = ls,
#          linewidth = lw,
#          color     = 'black',
#          alpha     = alp)
#         # plt.gca().plot(tp1 , trec, lw = ttw,linestyle='--',color="red",alpha=alp2);
#         # plt.gca().plot(tp2 , trec, lw = ttw,linestyle='--',color="red",alpha=alp2);
#         # plt.gca().plot(ts1 , trec, lw = ttw,linestyle='--',color="green",alpha=alp2);
#         # plt.gca().plot(ts2 , trec, lw = ttw,linestyle='--',color="green",alpha=alp2);
#         left,right = plt.xlim()
#         plt.xlim(left  = ts)  # adjust the left leaving right unchanged
#         plt.xlim(right = te)  # adjust the right leaving left unchanged
#         plt.tick_params(
#          axis        = "y",
#          width       = 1,
#          length      = 2.5,
#          direction   = "in",
#          color       = "black",
#          pad         = 2.5,
#          labelsize   = tic_s,
#          labelcolor  = "black",
#          colors      = "black",
#          zorder      = 20,
#          bottom      = "on", top      = "off", left      = "on", right      = "off",
#          labelbottom = "on", labeltop = "off", labelleft = "on", labelright = "off");
#         bottom, top = plt.ylim()
#         plt.ylim(bottom = minY)
#         plt.ylim(top    = maxY)
#         plt.yticks(yticksnum,ytickslabels);
#         plt.tick_params(
#          axis        = "x",
#          width       = 1,
#          length      = 2.5,
#          direction   = "in",
#          color       = "black",
#          pad         = 2.5,
#          labelsize   = tic_s,
#          labelcolor  = "black",
#          colors      = "black",
#          zorder      = 20,
#          bottom      = "on", top      = "off", left      = "on", right      = "off",
#          labelbottom = "on", labeltop = "off", labelleft = "on", labelright = "off");
#         plt.ylabel(r"receiver",fontsize   = ylab_s)
#         plt.title(r"x-component",fontsize = tit_s)
#
#         ################# gather
#         # y-plot
#         plt.subplot(2,3,5)
#         plt.plot(gath_Y_n  + 2*i,
#          linestyle = ls,
#          linewidth = lw,
#          color     = 'red',
#          alpha     = alp)
#         plt.plot(gath_Y_d  + 2*i,
#          linestyle = ls,
#          linewidth = lw,
#          color     = 'black',
#          alpha     = alp)
#         # plt.gca().plot(tp1 , trec, lw = ttw,linestyle='--',color="red",alpha=alp2);
#         # plt.gca().plot(tp2 , trec, lw = ttw,linestyle='--',color="red",alpha=alp2);
#         # plt.gca().plot(ts1 , trec, lw = ttw,linestyle='--',color="green",alpha=alp2);
#         # plt.gca().plot(ts2 , trec, lw = ttw,linestyle='--',color="green",alpha=alp2);
#         left,right = plt.xlim()
#         plt.xlim(left  = ts)  # adjust the left leaving right unchanged
#         plt.xlim(right = te)  # adjust the right leaving left unchanged
#         plt.tick_params(
#          axis        = "y",
#          width       = 1,
#          length      = 2.5,
#          direction   = "in",
#          color       = "black",
#          pad         = 2.5,
#          labelsize   = tic_s,
#          labelcolor  = "black",
#          colors      = "black",
#          zorder      = 20,
#          bottom      = "on", top      = "off", left      = "on", right      = "off",
#          labelbottom = "on", labeltop = "off", labelleft = "on", labelright = "off");
#         bottom, top = plt.ylim()
#         plt.ylim(bottom = minY)
#         plt.ylim(top    = maxY)
#         plt.yticks(yticksnum,ytickslabels);
#         plt.tick_params(
#          axis        = "x",
#          width       = 1,
#          length      = 2.5,
#          direction   = "in",
#          color       = "black",
#          pad         = 2.5,
#          labelsize   = tic_s,
#          labelcolor  = "black",
#          colors      = "black",
#          zorder      = 20,
#          bottom      = "on", top      = "off", left      = "on", right      = "off",
#          labelbottom = "on", labeltop = "off", labelleft = "on", labelright = "off");
#         plt.ylabel(r"receiver",fontsize   = ylab_s)
#         plt.title(r"y-component",fontsize = tit_s)
#
#         ################# gather
#         # z-plot
#         plt.subplot(2,3,6)
#         plt.plot(gath_Z_n  + 2*i,
#          linestyle = ls,
#          linewidth = lw,
#          color     = 'red',
#          alpha     = alp)
#         plt.plot(gath_Z_d  + 2*i,
#          linestyle = ls,
#          linewidth = lw,
#          color     = 'black',
#          alpha     = alp)
#         # plt.gca().plot(tp1 , trec, lw = ttw,linestyle='--',color="red",alpha=alp2);
#         # plt.gca().plot(tp2 , trec, lw = ttw,linestyle='--',color="red",alpha=alp2);
#         # plt.gca().plot(ts1 , trec, lw = ttw,linestyle='--',color="green",alpha=alp2);
#         # plt.gca().plot(ts2 , trec, lw = ttw,linestyle='--',color="green",alpha=alp2);
#         left,right = plt.xlim()
#         plt.xlim(left  = ts)  # adjust the left leaving right unchanged
#         plt.xlim(right = te)  # adjust the right leaving left unchanged
#         plt.tick_params(
#          axis        = "y",
#          width       = 1,
#          length      = 2.5,
#          direction   = "in",
#          color       = "black",
#          pad         = 2.5,
#          labelsize   = tic_s,
#          labelcolor  = "black",
#          colors      = "black",
#          zorder      = 20,
#          bottom      = "on", top      = "off", left      = "on", right      = "off",
#          labelbottom = "on", labeltop = "off", labelleft = "on", labelright = "off");
#         bottom, top = plt.ylim()
#         plt.ylim(bottom = minY)
#         plt.ylim(top    = maxY)
#         plt.yticks(yticksnum,ytickslabels);
#         plt.tick_params(
#          axis        = "x",
#          width       = 1,
#          length      = 2.5,
#          direction   = "in",
#          color       = "black",
#          pad         = 2.5,
#          labelsize   = tic_s,
#          labelcolor  = "black",
#          colors      = "black",
#          zorder      = 20,
#          bottom      = "on", top      = "off", left      = "on", right      = "off",
#          labelbottom = "on", labeltop = "off", labelleft = "on", labelright = "off");
#         plt.xlabel(r"time (ms)",fontsize  = xlab_s)
#         plt.ylabel(r"receiver",fontsize   = ylab_s)
#         plt.title(r"z-component",fontsize = tit_s)

#
#
# plt.subplot(4,2,2)
# plt.scatter(pos_errorX, pos_errorY,   # pos X, pos Y
#                 c     = pos_error[cond_err,3],                      # condition error
#                 norm  = colors.Normalize(vmin=0.0, vmax=2.0),
#                 alpha = 0.5,
#                 label = r" $(\Delta e \backslash e) 100 > 175{\%}$");
# plt.legend(loc = (0.45,0.85),#"upper left",
#          scatterpoints  = scatterpointnum,
#          ncol           = 1,
#          fontsize       = scaterlegendsize,
#          markerscale    = scatermarkersize,
#          markerfirst    = True,
#          frameon        = False,
#          fancybox       = False,
#          shadow         = False);
# # create the line
# plt.plot(xline,yline, color="red",linewidth=1.5,alpha=0.5);
# # cross marking the actual event point
# plt.plot(pos_errorX[ig],pos_errorY[ig],
#  marker="*",markerfacecolor='yellow',markersize=markersize_2,markeredgecolor="black");
#
# plt.plot(xc1,yc1,marker="v", markerfacecolor='black', markersize=markersize_1, markeredgecolor="black");
# plt.plot(xc2,yc2,marker="v", markerfacecolor='black', markersize=markersize_1, markeredgecolor="black");
# circle1 = plt.Circle((xc1,yc1),radius,color = 'r',lw = 1.0, edgecolor = 'r',fill = False,alpha = 0.5);
# circle2 = plt.Circle((xc2,yc2),radius,color = 'r',lw = 1.0, edgecolor = 'r',fill = False,alpha = 0.5);
# plt.gca().add_artist(circle1);
# plt.gca().add_artist(circle2);
# cbar = plt.colorbar(extend='both', shrink=0.9);
# cbar.set_label("relative error");
# plt.xlabel("x (m)")
# plt.ylabel("y (m)")
# left,right = plt.xlim()
# plt.xlim(left  = xminlim)  # adjust the left leaving right unchanged
# plt.xlim(right = xmaxlim)  # adjust the right leaving left unchanged
# bottom, top = plt.ylim()
# plt.ylim(bottom = yminlim)
# plt.ylim(top    = ymaxlim)
# plt.tick_params(
#             axis        = "y",
#             width       = 1,
#             length      = 2.5,
#             direction   = "in",
#             color       = "black",
#             pad         = 2.5,
#             labelsize   = 12,
#             labelcolor  = "black",
#             colors      = "black",
#             zorder      = 20,
#             bottom      = "on", top      = "off", left      = "on", right      = "off",
#             labelbottom = "on", labeltop = "off", labelleft = "on", labelright = "off");
# plt.tick_params(
#             axis        = "x",
#             width       = 1,
#             length      = 2.5,
#             direction   = "in",
#             color       = "black",
#             pad         = 2.5,
#             labelsize   = 12,
#             labelcolor  = "black",
#             colors      = "black",
#             zorder      = 20,
#             bottom      = "on", top      = "off", left      = "on", right      = "off",
#             labelbottom = "on", labeltop = "off", labelleft = "on", labelright = "off");
# #plt.yticks([])
#
# plt.subplot(4,2,4)
# plt.scatter(pos_errorX, pos_errorY, # pos X, pos Y
#                 c     = pos_error[cond_err,3],
#                 norm  = colors.Normalize(vmin=0.0, vmax=2.0),
#                 alpha = 0.5,
#                 label = r" $(\Delta e \backslash e) 100 > 175{\%}$");
# plt.legend(loc = (0.45,0.85),#"upper left",
#          scatterpoints  = scatterpointnum,
#          ncol           = 1,
#          fontsize       = scaterlegendsize,
#          markerscale    = scatermarkersize,
#          markerfirst    = True,
#          frameon        = False,
#          fancybox       = False,
#          shadow         = False);
# # create the line
# plt.plot(xline,yline, color="red",linewidth=1.5,alpha=0.5);
# # cross marking the actual event point
# plt.plot(pos_errorX[ig],pos_errorY[ig],
#  marker="*", markerfacecolor='yellow',markersize=markersize_2,markeredgecolor="black");
#
# plt.plot(xc1,yc1,marker="v", markerfacecolor='black', markersize=markersize_1, markeredgecolor="black");
# plt.plot(xc2,yc2,marker="v", markerfacecolor='black', markersize=markersize_1, markeredgecolor="black");
# circle1 = plt.Circle((xc1,yc1),radius,color = 'r',lw = 1.0, edgecolor = 'r',fill = False,alpha = 0.5);
# circle2 = plt.Circle((xc2,yc2),radius,color = 'r',lw = 1.0, edgecolor = 'r',fill = False,alpha = 0.5);
# plt.gca().add_artist(circle1);
# plt.gca().add_artist(circle2);
# cbar = plt.colorbar(extend='both', shrink=0.9);
# cbar.set_label("relative error");
# plt.xlabel("x (m)")
# plt.ylabel("y (m)")
# left,right = plt.xlim()
# plt.xlim(left  = xminlim)  # adjust the left leaving right unchanged
# plt.xlim(right = xmaxlim)  # adjust the right leaving left unchanged
# bottom, top = plt.ylim()
# plt.ylim(bottom = yminlim)
# plt.ylim(top    = ymaxlim)
# plt.tick_params(
#             axis        = "y",
#             width       = 1,
#             length      = 2.5,
#             direction   = "in",
#             color       = "black",
#             pad         = 2.5,
#             labelsize   = 12,
#             labelcolor  = "black",
#             colors      = "black",
#             zorder      = 20,
#             bottom      = "on", top      = "off", left      = "on", right      = "off",
#             labelbottom = "on", labeltop = "off", labelleft = "on", labelright = "off");
# plt.tick_params(
#             axis        = "x",
#             width       = 1,
#             length      = 2.5,
#             direction   = "in",
#             color       = "black",
#             pad         = 2.5,
#             labelsize   = 12,
#             labelcolor  = "black",
#             colors      = "black",
#             zorder      = 20,
#             bottom      = "on", top      = "off", left      = "on", right      = "off",
#             labelbottom = "on", labeltop = "off", labelleft = "on", labelright = "off");
# #plt.yticks([])
#
# plt.subplot(4,2,6)
# plt.scatter(pos_errorX, pos_errorZ, # pos X, pos Z
#                 c     = pos_error[cond_err,3],
#                 norm  = colors.Normalize(vmin=0.0, vmax=2.0),
#                 alpha = 0.5,
#                 label = r"$(\Delta e \backslash e) 100 > 175{\%}$");
# plt.legend(loc = (0.45,0.85),#"upper left",
#          scatterpoints  = scatterpointnum,
#          ncol           = 1,
#          fontsize       = scaterlegendsize,
#          markerscale    = scatermarkersize,
#          markerfirst    = True,
#          frameon        = False,
#          fancybox       = False,
#          shadow         = False);
# # create the line
# plt.plot(zhline1,zvline, color="red",linewidth=1.5,alpha=0.5);
# plt.plot(zhline2,zvline, color="red",linewidth=1.5,alpha=0.5);
# # cross marking the actual event point
# plt.plot(pos_errorX[ig],pos_errorZ[ig],
#  marker="*", markerfacecolor='yellow',markersize=markersize_2,markeredgecolor="black");
#
# plt.plot(xc1,zc1[0],marker="v", markerfacecolor='black', markersize=markersize_1, markeredgecolor="black");
# plt.plot(xc1,zc1[1],marker="v", markerfacecolor='black', markersize=markersize_1, markeredgecolor="black");
# plt.plot(xc1,zc1[2],marker="v", markerfacecolor='black', markersize=markersize_1, markeredgecolor="black");
# plt.plot(xc1,zc1[3],marker="v", markerfacecolor='black', markersize=markersize_1, markeredgecolor="black");
# plt.plot(xc1,zc1[4],marker="v", markerfacecolor='black', markersize=markersize_1, markeredgecolor="black");
# plt.plot(xc2,zc2[0],marker="v", markerfacecolor='black', markersize=markersize_1, markeredgecolor="black");
# plt.plot(xc2,zc2[1],marker="v", markerfacecolor='black', markersize=markersize_1, markeredgecolor="black");
# plt.plot(xc2,zc2[2],marker="v", markerfacecolor='black', markersize=markersize_1, markeredgecolor="black");
# plt.plot(xc2,zc2[3],marker="v", markerfacecolor='black', markersize=markersize_1, markeredgecolor="black");
# plt.plot(xc2,zc2[4],marker="v", markerfacecolor='black', markersize=markersize_1, markeredgecolor="black");
#
# # circle1 = plt.Circle((xc1,yc1),radius,color = 'r',lw = 1.0, edgecolor = 'r',fill = False,alpha = 0.5);
# # circle2 = plt.Circle((xc2,yc2),radius,color = 'r',lw = 1.0, edgecolor = 'r',fill = False,alpha = 0.5);
# # plt.gca().add_artist(circle1);
# # plt.gca().add_artist(circle2);
# cbar = plt.colorbar(extend='both', shrink=0.9);
# cbar.set_label("relative error");
# plt.xlabel("x (m)")
# plt.ylabel("z (m)")
# left,right = plt.xlim()
# plt.xlim(left  = xminlim)  # adjust the left leaving right unchanged
# plt.xlim(right = xmaxlim)  # adjust the right leaving left unchanged
# bottom, top = plt.ylim()
# plt.ylim(bottom = zminlim)
# plt.ylim(top    = zmaxlim)
# plt.tick_params(
#             axis        = "y",
#             width       = 1,
#             length      = 2.5,
#             direction   = "in",
#             color       = "black",
#             pad         = 2.5,
#             labelsize   = 12,
#             labelcolor  = "black",
#             colors      = "black",
#             zorder      = 20,
#             bottom      = "on", top      = "off", left      = "on", right      = "off",
#             labelbottom = "on", labeltop = "off", labelleft = "on", labelright = "off");
# plt.tick_params(
#             axis        = "x",
#             width       = 1,
#             length      = 2.5,
#             direction   = "in",
#             color       = "black",
#             pad         = 2.5,
#             labelsize   = 12,
#             labelcolor  = "black",
#             colors      = "black",
#             zorder      = 20,
#             bottom      = "on", top      = "off", left      = "on", right      = "off",
#             labelbottom = "on", labeltop = "off", labelleft = "on", labelright = "off");
# #plt.yticks([])
#
# amplitud = xdata175[ig,3:];
# ampP      = amplitud[0:30];
# ampS     = amplitud[30:60];
#
# ampPwell1 = ampP[0:15];
# ampPwell2 = ampP[15:30];
#
# ampSwell1 = ampS[0:15];
# ampSwell2 = ampS[15:30];
#
#
# plt.subplot(4,2,(7,8))
# plt.scatter(range(len(ampP)),ampP,color="red",alpha = 0.5,marker="v",s=70);
# plt.scatter(range(len(ampS)),ampS,color="blue",alpha = 0.5,marker="*",s=70);
# # plt.legend(loc = (0.45,0.85),#"upper left",
# #          scatterpoints  = scatterpointnum,
# #          ncol           = 1,
# #          fontsize       = scaterlegendsize,
# #          markerscale    = scatermarkersize,
# #          markerfirst    = True,
# #          frameon        = False,
# #          fancybox       = False,
# #          shadow         = False);
# # # create the line
# # plt.plot(zhline1,zvline, color="red",linewidth=1.5,alpha=0.5);
# # plt.plot(zhline2,zvline, color="red",linewidth=1.5,alpha=0.5);
# # # cross marking the actual event point
# # plt.plot(pos_errorX[ig],pos_errorZ[ig],
# #  marker="*", markerfacecolor='yellow',markersize=markersize_2,markeredgecolor="black");
# #
# # plt.plot(xc1,zc1[0],marker="v", markerfacecolor='black', markersize=markersize_1, markeredgecolor="black");
# # plt.plot(xc1,zc1[1],marker="v", markerfacecolor='black', markersize=markersize_1, markeredgecolor="black");
# # plt.plot(xc1,zc1[2],marker="v", markerfacecolor='black', markersize=markersize_1, markeredgecolor="black");
# # plt.plot(xc1,zc1[3],marker="v", markerfacecolor='black', markersize=markersize_1, markeredgecolor="black");
# # plt.plot(xc1,zc1[4],marker="v", markerfacecolor='black', markersize=markersize_1, markeredgecolor="black");
# # plt.plot(xc2,zc2[0],marker="v", markerfacecolor='black', markersize=markersize_1, markeredgecolor="black");
# # plt.plot(xc2,zc2[1],marker="v", markerfacecolor='black', markersize=markersize_1, markeredgecolor="black");
# # plt.plot(xc2,zc2[2],marker="v", markerfacecolor='black', markersize=markersize_1, markeredgecolor="black");
# # plt.plot(xc2,zc2[3],marker="v", markerfacecolor='black', markersize=markersize_1, markeredgecolor="black");
# # plt.plot(xc2,zc2[4],marker="v", markerfacecolor='black', markersize=markersize_1, markeredgecolor="black");
# #
# # # circle1 = plt.Circle((xc1,yc1),radius,color = 'r',lw = 1.0, edgecolor = 'r',fill = False,alpha = 0.5);
# # # circle2 = plt.Circle((xc2,yc2),radius,color = 'r',lw = 1.0, edgecolor = 'r',fill = False,alpha = 0.5);
# # # plt.gca().add_artist(circle1);
# # # plt.gca().add_artist(circle2);
# # cbar = plt.colorbar(extend='both', shrink=0.9);
# # cbar.set_label("relative error");
# # plt.xlabel("x (m)")
# plt.ylabel("Amplitude")
# left,right = plt.xlim()
# plt.xlim(left  = -1.0)  # adjust the left leaving right unchanged
# plt.xlim(right = 30)  # adjust the right leaving left unchanged
# bottom, top = plt.ylim()
# plt.ylim(bottom = -1.1)
# plt.ylim(top    = 1.1)
# plt.xticks(np.arange(0,30,1));
# plt.tick_params(
#             axis        = "y",
#             width       = 1,
#             length      = 2.5,
#             direction   = "in",
#             color       = "black",
#             pad         = 2.5,
#             labelsize   = 12,
#             labelcolor  = "black",
#             colors      = "black",
#             zorder      = 20,
#             bottom      = "on", top      = "off", left      = "on", right      = "off",
#             labelbottom = "on", labeltop = "off", labelleft = "on", labelright = "off");
# plt.tick_params(
#             axis        = "x",
#             width       = 1,
#             length      = 2.5,
#             direction   = "in",
#             color       = "black",
#             pad         = 2.5,
#             labelsize   = 12,
#             labelcolor  = "black",
#             colors      = "black",
#             zorder      = 20,
#             bottom      = "on", top      = "off", left      = "on", right      = "off",
#             labelbottom = "on", labeltop = "off", labelleft = "on", labelright = "off");
# #plt.yticks([])

plt.tight_layout(w_pad=0.0,h_pad=0.0)

fname = fdr_fig_+'gather_xyz.png'
plt.savefig(fname, bbox_inches='tight',
dpi=300,
facecolor='w',
edgecolor='w',
orientation='portrait',
papertype=None,
transparent=False,
pad_inches=0.1,
frameon=None,
metadata=None)
# fname = fdr_fig_+'gather_xyz_'+iname+'.eps'
# plt.savefig(fname)
# plt.show()
plt.close(fig2)
