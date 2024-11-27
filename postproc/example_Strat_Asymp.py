import h5py
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np


plt.rcParams['text.usetex'] = True
plt.rc('text.latex', preamble=r'\usepackage{amsmath}')

fontsize_title = 20
fontsize_labels = 20
fontsize_ticks = 12
fontsize_legend = 16

cmap = plt.colormaps["jet"]


def plot_global_quantities(times, N, H, outputDir):
    plt.figure(figsize=(7,5), constrained_layout=True)

    dNdt = np.diff(N) #/ np.diff(times)

    plt.plot(times, N, label="$\mathcal{N}$")
    plt.plot(times, H, label="$E$")

    #plt.plot(times[0:len(N)-1], dNdt/N[0:len(N)-1], label="$\mathcal{N}$")

    plt.xlabel("$t$", fontsize=fontsize_labels)
    #plt.title("$\\rm{Global~quantities~vs~time}$", fontsize=fontsize_title)
    plt.legend(fontsize=fontsize_legend, loc="upper left")
    plt.xlim([0, times[-1]])
    #plt.xlim([45,48])
    #plt.xticks([0, 2.5e3, 5e3, 7.5e3, 1e4], ["$0$", "$2.5 \\times 10^3$", "$5 \\times 10^3$", "$7.5 \\times 10^3$", "$1 \\times 10^4$"], fontsize=fontsize_ticks)
    #plt.ylim([0, 1000])
    #plt.yticks([0, 2.5e2, 5e2, 7.5e2, 1e3], ["$0$", "$2.5 \\times 10^2$", "$5 \\times 10^2$", "$7.5 \\times 10^2$", "$1 \\times 10^3$"], fontsize=fontsize_ticks)
    plt.savefig(f"{outputDir}globals.png",dpi=300)
    #plt.show()

def plot_nk(kh, kz, nk, times, outputDir, it=-1):
    KH, KZ = np.meshgrid(kh, kz)
    
    plt.figure(figsize=(6,5), constrained_layout=True)
    plt.pcolormesh(KH, KZ, np.log10(nk[it, :, :] * KH / KZ ), cmap=cmap)
    plt.xlabel("$k_h$", fontsize=fontsize_labels)
    plt.ylabel("$k_z$", fontsize=fontsize_labels)
    plt.title("$t = {:.1f}$".format(times[it],1), fontsize=fontsize_title)
    
    # Forcing parameter
    khf = 0.007/np.sqrt(2)
    kzf = 0.007/np.sqrt(2)
    of = khf/kzf
    invf = of/kzf


    # ES
    ##plt.plot([min(kh),khf],[2*kzf,2*kzf],"m--")
    #plt.text(-1.75,-1.25,"ES", color="m", fontsize=16)

    # PSI + 
    ##plt.plot([khf,max(kh)],[khf/(2*of),max(kh)/(2*of)],"k--")
    #plt.text(-1,-0.5,"PSI", color="k", fontsize=16)

    # PSI - 
    khmax = max(kz) * of / 2 
    ##plt.plot([khf,khmax],[2*khf/of,2*khmax/of],"k--")
    #x = np.array([2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 128.0, 256.0])
    #plt.plot(x*khf,x*2*khf/of,"ko")
    #plt.text(-1,-0.5,"PSI", color="k", fontsize=16)

    # ID
    invf = of/2
    #plt.plot([min(kh),max(kh)],[np.sqrt(min(kh)/invf),np.sqrt(max(kh)/invf)],"w--")
    #plt.text(-1,-0.35,"ID", color="w", fontsize=16)
    
    

    plt.clim([-8,2])
    plt.xlim([min(kh), max(kh)])
    plt.xticks(fontsize=fontsize_ticks)
    plt.xscale('log')
    plt.ylim([min(kz), max(kz)])
    plt.yticks(fontsize=fontsize_ticks)
    plt.yscale('log')
    plt.axis("square")
    cbar = plt.colorbar()
    #cbar.set_label("$e(k_h,k_z,t)$", fontsize=fontsize_labels) 
    cbar.ax.yaxis.set_ticks([-5, -2, 1, 4, 7], ["$10^{-5}$", "$10^{-2}$", "$10^1$", "$10^4$", "$10^7$"], fontsize=fontsize_ticks)
    

    plt.savefig(f"{outputDir}ek_it{it:03d}.png",dpi=300)
    #plt.show()
    plt.close()


def plot_nk_slices(kh, kz, nk, times, outputDir, it=-1):
    nkt = nk[it, :, :]
    Mh = len(kh)
    Mz = len(kz)

    # Slice in kz
    plt.figure(figsize=(7,5), constrained_layout=True)

    plt.plot(kz, nkt[int(np.floor(Mh/4)),:], c=mpl.cm.Blues(0.25), label="$k_h = {:.2e}$".format(kh[int(np.floor(Mh/4))]))
    plt.plot(kz, nkt[int(np.floor(2*Mh/4)),:], c=mpl.cm.Blues(0.5), label="$k_h = {:.2e}$".format(kh[int(np.floor(2*Mh/4))]))
    plt.plot(kz, nkt[int(np.floor(3*Mh/4)),:], c=mpl.cm.Blues(0.75), label="$k_h = {:.2e}$".format(kh[int(np.floor(3*Mh/4))]))
    plt.xlabel("$k_z$", fontsize=fontsize_labels)
    plt.ylabel("$n(k_h,k_z)$", fontsize=fontsize_labels)

    plt.xlim([min(kh), max(kh)])
    plt.xticks(fontsize=fontsize_ticks)
    plt.xscale('log')
    plt.ylim([1e-14, 1e9])
    plt.yticks(fontsize=fontsize_ticks)
    plt.yscale('log')
    plt.legend(fontsize=fontsize_legend, loc="upper right")
    plt.title("$t = {:.1f}$".format(times[it],1), fontsize=fontsize_title)


    plt.savefig(f"{outputDir}nk_slices_kz_it{it:03d}.png",dpi=100)
    #plt.show()
    plt.close()

    # Slice in kh
    plt.figure(figsize=(7,5), constrained_layout=True)

    plt.plot(kh, nkt[:,int(np.floor(Mh/4))], c=mpl.cm.Reds(0.25), label="$k_z = {:.2e}$".format(kz[int(np.floor(Mz/4))]))
    plt.plot(kh, nkt[:,int(np.floor(2*Mh/4))], c=mpl.cm.Reds(0.5), label="$k_z = {:.2e}$".format(kz[int(np.floor(2*Mz/4))]))
    plt.plot(kh, nkt[:,int(np.floor(3*Mh/4))], c=mpl.cm.Reds(0.75), label="$k_z = {:.2e}$".format(kz[int(np.floor(3*Mz/4))]))
    plt.xlabel("$k_h$", fontsize=fontsize_labels)
    plt.ylabel("$n(k_h,k_z)$", fontsize=fontsize_labels)

    plt.xlim([min(kz), max(kz)])
    plt.xticks(fontsize=fontsize_ticks)
    plt.xscale('log')
    plt.ylim([1e-14, 1e9])
    plt.yticks(fontsize=fontsize_ticks)
    plt.yscale('log')
    plt.legend(fontsize=fontsize_legend, loc="upper right")
    plt.title("$t = {:.1f}$".format(times[it],1), fontsize=fontsize_title)
    

    plt.savefig(f"{outputDir}nk_slices_kh_it{it:03d}.png",dpi=100)
    #plt.show()
    plt.close()


def plot_nk_vs_time(kh, kz, nk, times, outputDir):
    for it in range(int(np.floor(len(times)))): #int(np.floor(len(times)/10))
        plot_nk(kh, kz, nk, times, outputDir,it)  

def plot_nk_slices_vs_time(kh, kz, nk, times, outputDir):
    for it in range(int(np.floor(len(times)/10))):
        plot_nk_slices(kh, kz, nk, times, outputDir, 10*it+1)
        



def main():
    #outputDir = f"../"
    outputDir = f"/home/vincent/Data/Wavkins/WKE_Stratified_Asymp/localized_M80_kmin0.005_kmax1.0_kdinf0.005_kdsup1.0_kfh0.3_kfz0.3_lapPower4/"
    data_file = h5py.File(f"{outputDir}WKE_Stratified_Asymp_data.nc", "r")
    #print(list(data_file.keys()))
    times = data_file["Global/Times"][:]
    timesSP = data_file['Spectral/TimesSP'][:]
    #print(timesSP)
    N = data_file["Global/N"][:]
    H = data_file["Global/H"][:]
    nk = data_file["Spectral/nk"][:, :, :]
    kh = data_file["Spectral/kh"][:]
    kz = data_file["Spectral/kz"][:]
    Mh = len(kh)
    Mz = len(kz)
    
    plot_nk(kh, kz, nk, timesSP, outputDir)
    plot_global_quantities(times, N, H, outputDir)
    plot_nk_vs_time(kh, kz, nk, timesSP, outputDir)
    #plot_nk_slices(kh, kz, nk, times, outputDir)
    #plot_nk_slices_vs_time(kh, kz, nk, times, outputDir)
    #print(times)


if __name__ == "__main__":
    main()
