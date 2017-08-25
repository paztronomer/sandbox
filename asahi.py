import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#Load data and get relations
df = pd.read_table("DECam_filters_transmission.txt",sep="\s+",engine="python")
print len(df.index)
u,g,r,i,z,Y = 0,0,0,0,0,0
for idx,row in df.iterrows():
    if row["u"] > 1E-4: u += row["u"]
    if row["g"] > 1E-4: g += row["g"]
    if row["r"] > 1E-4: r += row["r"]
    if row["i"] > 1E-4: i += row["i"]
    if row["z"] > 1E-4: z += row["z"]
    if row["Y"] > 1E-4: Y += row["Y"]
integ = np.array([u,g,r,i,z,Y])
pcent = []
for i in integ:
    tmp = i/max(integ)
    pcent.append("{0:.3}".format(tmp))

#Plot
colors = ["blue","forestgreen","gold","red","magenta","darkgray"]
col = ["u","g","r","i","z","Y"]
fig,ax = plt.subplots(1,figsize=(10,6))
for idx,nm in enumerate(col):
    lab = "{0}, rel. ratio={1}".format(nm,pcent[idx])
    im = ax.plot(df["nm"],df[nm],"-",color=colors[idx],lw=1,alpha=0.8,
                label=lab)

#Setup plot
ax.set_xlabel(r"$\lambda$ [nm]",fontsize=14)
ax.set_ylabel("flux",fontsize=14)
ax.set_yscale("log")
ax.set_title("Asahi bandpass\n"+r"http://www.ctio.noao.edu/noao/sites/default/files/DECam/DECam_filters_transmission.txt")
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles,labels)

#Grid
ax.grid(color="lightgray",linestyle="dotted",dash_capstyle="round",alpha=0.7)

#xticks
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
majorLocator_x = MultipleLocator(100)
majorFormatter_x = FormatStrFormatter('%d')
minorLocator_x = MultipleLocator(10)
ax.xaxis.set_minor_locator(minorLocator_x)
ax.xaxis.set_major_locator(majorLocator_x)
ax.xaxis.set_major_formatter(majorFormatter_x)

#finally
plt.subplots_adjust(left=0.07,bottom=0.1,right=0.99,top=0.91)
plt.savefig("asahi.pdf",dpi=300,format="pdf")
plt.show()
