import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from os import path
import wavepal as wv

here = path.abspath(path.dirname(__file__))+"/"
path_to_figure_folder=here+"figures/"

data=np.genfromtxt(here+"data/ODP1148-BF-18O.txt")
myt=data[:,0]
mydata=data[:,1]

x=wv.Wavepal(myt, mydata, "Age", "$\delta{}^{18}O$", t_units="ka", mydata_units="permil")
x.check_data()
plot_timestep=x.plot_timestep()
plot_timestep.savefig(path_to_figure_folder+"timestep.pdf")
plot_timestep.close()
plot_trend=x.plot_trend(pol_degree=7)
plot_trend.savefig(path_to_figure_folder+"trend.pdf")
plot_trend.close()
x.choose_trend_degree(7)
x.trend_vectors()
x.carma_params(make_carma_fig=True,nbins=20,dpi=400,path_to_figure_folder=path_to_figure_folder)
percentile=np.zeros(2)
percentile[0]=95.
percentile[1]=99.9
x.freq_analysis(freqstep=0.0001,D=600.,percentile=percentile,n_moments=12)
plot_periodogram=x.plot_periodogram(fontsize_legend=10)
plot_periodogram.savefig(path_to_figure_folder+"periodogram.pdf")
plot_periodogram.close()
x.timefreq_analysis(permin=10.,percentile=percentile)
time_string=[0., 500., 1000., 1500., 2000., 2500., 3000., 3500., 4000., 4500., 5000., 5500., 6000.]
period_string=[10., 21., 41., 100., 200., 400., 800., 1500.]
dashed_periods=[21., 41., 100.]
plot_scalogram=x.plot_scalogram(color_cl_anal=['m','g'],fontsize_title=40,fontsize_ticks=20,fontsize_axes=30,time_string=time_string,period_string=period_string,dashed_periods=dashed_periods)
fig = plt.gcf()
fig.set_size_inches(52, 26)
plot_scalogram.savefig(path_to_figure_folder+"scalogram_w0_5,5.pdf")
plot_scalogram.close()
print "INSTALLATION COMPLETED SUCCESSFULLY"
