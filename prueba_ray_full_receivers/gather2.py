import os
from scipy.signal import butter, lfilter
from matplotlib import rc, font_manager
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as colors

rc('text', usetex=True)
plt.rc('font', family='serif')

path = 'signals_ppshsv/';
fdr_fig_  = "gathers_ppshsv/"
filenames   = os.listdir(path)      ;  # list of all file names in folder in path
num_files   = np.shape(filenames)[0];  # length of the list: number of files
num_gathers = int(num_files/3);

##########################################
##########################################
for i in range(num_gathers):
        iname     = str(i+1);
        print(iname)
        filenamex = path+'signal_x_'+iname+'.txt';
        filenamey = path+'signal_y_'+iname+'.txt';
        filenamez = path+'signal_z_'+iname+'.txt';

        datax     = np.loadtxt(filenamex)               ; # extract x-data.
        datay     = np.loadtxt(filenamey)               ; # extract y-data.
        dataz     = np.loadtxt(filenamez)               ; # extract z-data.
        [nsx,ntx] = np.shape(datax)                ; # dimension of x-data.
        print()
        print("dimension of x-data")
        print("number of samples   : ", nsx);
        print("number of traces    : ", ntx);
        [nsy,nty] = np.shape(datay)                ; # dimension of y-data.
        print()
        print("dimension of y-data")
        print("number of samples   : ", nsy);
        print("number of traces    : ", nty);
        [nsz,ntz] = np.shape(dataz)                ; # dimension of z-data.
        print()
        print("dimension of z-data")
        print("number of samples   : ", nsz);
        print("number of traces    : ", ntz);


        porc_norm = 0.99    # porcentage of normalize data.

        maximx   = np.amax(abs(datax))             ; # extract data maximum x normalization.
        maximy   = np.amax(abs(datay))             ; # extract data maximum x normalization.
        maximz   = np.amax(abs(dataz))             ; # extract data maximum x normalization.
        max      = np.amax([maximx,maximy,maximz]) ;

        datax    = datax/max                    ; # normalize
        datay    = datay/max                    ; # normalize
        dataz    = dataz/max                    ; # normalize

        datax    = datax*porc_norm                ; # take a porcentage of normalize data.
        datay    = datay*porc_norm                ; # take a porcentage of normalize data.
        dataz    = dataz*porc_norm                ; # take a porcentage of normalize data.

        ts = 0.0 ; # time begins in 0.0 ms
        te = 300.0#nsx  ; # time ends at m ms

        # Because data is normalize to 1 (or porc of 1)
        # each trace occupies 2 vertical units (maximum).

        # we need at least 2*nt vertical space in ylim
        first_trace = 0        ;  # first trace is positioned at zero.
        last_trace  = 2*(ntx-1);  # (nt-1) : number of intervals between nt traces
        # 2*     : time 2, vertical space.
        minY  = (first_trace - 1) - 0.2 ;
        maxY  = (last_trace  + 1) + 0.2 ;

        lw     = 1.5;  # linewidth
        ls     = '-';  # linestyle
        alp    = 0.6;  # alpha ransparency for line
        xlab_s = 15;   # x-label size
        ylab_s = 15;   # y-label size
        tit_s  = 15;   # title size
        tic_s  = 12;   # tics size

        scale  = 3.0; # scales data (1.0: 0 scaling). Warn: makes data bigger.
        yticksnum     = np.linspace(first_trace,last_trace,ntx);                # y tick numbers
        ytickslabels  = np.tile(((np.linspace(1,ntx,ntx)).astype(np.int64)),1); # y tick labels

        fig2,ax2  = plt.subplots(figsize=(6,8));

        for i in range(ntx):
                # x-plot
                gath_X = scale*datax[:,i];
                gath_Y = scale*datay[:,i];
                gath_Z = scale*dataz[:,i];

                plt.subplot(3,1,1)
                plt.plot(gath_X  + 2*i,
                linestyle = ls,
                linewidth = lw,
                color     = 'black',
                alpha     = alp)
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
                plt.subplot(3,1,2)
                # y-plot
                plt.plot(gath_Y  + 2*i,
                linestyle = ls,
                linewidth = lw,
                color     = 'black',
                alpha     = alp)
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
                plt.subplot(3,1,3)
                # z-plot
                plt.plot(gath_Z  + 2*i,
                linestyle = ls,
                linewidth = lw,
                color     = 'black',
                alpha     = alp)
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
                plt.tight_layout(w_pad=0.0,h_pad=0.2)

        fname = fdr_fig_+'gather_xyz_'+iname+'.png'
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
