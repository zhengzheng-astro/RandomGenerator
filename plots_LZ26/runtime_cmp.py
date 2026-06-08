import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.patches as mpatches

fig, ax = plt.subplots()

ax.xaxis.set_tick_params(which='both', direction='in', top=True)
ax.yaxis.set_tick_params(which='both', direction='in', right=True)

ax.set_xlim(0, 20)
ax.set_ylim(5e-8, 5e-6)
ax.set_yscale('log')
ax.set_xlabel(r'$x$', fontsize=15)
ax.set_ylabel('time per draw (s)', fontsize=15)

x, t0=np.loadtxt('data/runtime0_T2.dat',unpack=True)
x, t1=np.loadtxt('data/runtime1_T2.dat',unpack=True)
x, t2=np.loadtxt('data/runtime2_T2.dat',unpack=True)
x3,t3=np.loadtxt('data/runtime3_T2.dat',unpack=True)
x, t4=np.loadtxt('data/runtime4_T2.dat',unpack=True)

ax.plot(x,t0,color='m',linewidth=0.5,alpha=1.0)
ax.plot(x,t1,color='k',linewidth=0.5,alpha=1.0)

ax.plot(x[x<=10],t2[x<=10],color='r',linewidth=0.5,alpha=1.0)
ax.plot(x[x> 10],t2[x> 10],color='r',linewidth=0.5,alpha=1.0,linestyle='dotted')

ax.plot(x3[x3<=8], t3[x3<=8],color='g',linewidth=0.5,alpha=1.0)
ax.plot(x3[x3> 8], t3[x3> 8],color='g',linewidth=0.5,alpha=1.0,linestyle='dotted')

ax.plot(x[x< 9], t4[x< 9], color='b',linewidth=0.5,alpha=1.0)
ax.plot(x[x>=9], t4[x>=9], color='b',linewidth=0.5,alpha=1.0,linestyle='dotted')

x, t0=np.loadtxt('data/runtime0_T4.dat',unpack=True)
x, t1=np.loadtxt('data/runtime1_T4.dat',unpack=True)
x, t2=np.loadtxt('data/runtime2_T4.dat',unpack=True)
x3,t3=np.loadtxt('data/runtime3_T4.dat',unpack=True)
x, t4=np.loadtxt('data/runtime4_T4.dat',unpack=True)

ax.plot(x,t0,color='m',linewidth=1,alpha=0.75)
ax.plot(x,t1,color='k',linewidth=1,alpha=0.75)

ax.plot(x[x<=10],t2[x<=10],color='r',linewidth=1,alpha=0.75)
ax.plot(x[x> 10],t2[x> 10],color='r',linewidth=1,alpha=0.75,linestyle='dotted')

ax.plot(x3[x3<=8], t3[x3<=8],color='g',linewidth=1,alpha=0.75)
ax.plot(x3[x3> 8], t3[x3> 8],color='g',linewidth=1,alpha=0.75,linestyle='dotted')

ax.plot(x[x< 9], t4[x<9],  color='b',linewidth=1,alpha=0.75)
ax.plot(x[x>=9], t4[x>=9], color='b',linewidth=1,alpha=0.75,linestyle='dotted')

x, t0=np.loadtxt('data/runtime0_T6.dat',unpack=True)
x, t1=np.loadtxt('data/runtime1_T6.dat',unpack=True)
x, t2=np.loadtxt('data/runtime2_T6.dat',unpack=True)
x3,t3=np.loadtxt('data/runtime3_T6.dat',unpack=True)
x, t4=np.loadtxt('data/runtime4_T6.dat',unpack=True)

ax.plot(x,t0,color='m',linewidth=2,alpha=0.5)
ax.plot(x,t1,color='k',linewidth=2,alpha=0.5)

ax.plot(x[x<=10],t2[x<=10],color='r',linewidth=2,alpha=0.5)
ax.plot(x[x> 10],t2[x> 10],color='r',linewidth=2,alpha=0.5,linestyle='dotted')

ax.plot(x3[x3<=8], t3[x3<=8],color='g',linewidth=2,alpha=0.5)
ax.plot(x3[x3> 8], t3[x3> 8],color='g',linewidth=2,alpha=0.5,linestyle='dotted')

ax.plot(x[x< 9], t4[x<9],  color='b',linewidth=2,alpha=0.5)
ax.plot(x[x>=9], t4[x>=9], color='b',linewidth=2,alpha=0.5,linestyle='dotted')


line0 = Line2D([0], [0], label='ZM2002',       color='m')
line2 = Line2D([0], [0], label='This Work',    color='k')
line3 = Line2D([0], [0], label='Semelin+2007', color='r')
line4 = Line2D([0], [0], label='RASCAS',       color='g')
line5 = Line2D([0], [0], label='Smith+2015',   color='b')

lineT2 = Line2D([0], [0], label=r'$T=10^2{\rm K}\, (a=4.7\times10^{-3})$',color='k',linewidth=0.5,alpha=1)
lineT4 = Line2D([0], [0], label=r'$T=10^4{\rm K}\, (a=4.7\times10^{-4})$',color='k',linewidth=1,alpha=0.75)
lineT6 = Line2D([0], [0], label=r'$T=10^6{\rm K}\, (a=4.7\times10^{-5})$',color='k',linewidth=2,alpha=0.5)

noline = Line2D([0], [0], label='', alpha=0)

legend1=ax.legend(handles=[line0,line3,line5,line4,line2],framealpha=1,loc=[0.53,0.10])

legend2=ax.legend(handles=[lineT2,lineT4,lineT6],framealpha=1,loc=[0.53,0.75])

ax.add_artist(legend1)
ax.add_artist(legend2)

#plt.savefig('runtime_cmp.pdf')
plt.savefig('runtime_cmp.png')

plt.show()

