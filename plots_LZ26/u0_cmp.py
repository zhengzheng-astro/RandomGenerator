import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.patches as mpatches

fig, ax = plt.subplots()

ax.xaxis.set_tick_params(which='both', direction='in', top=True)
ax.yaxis.set_tick_params(which='both', direction='in', right=True)

ax.set_xlim(0, 20)
ax.set_ylim(-0.1,5)
ax.set_xlabel(r'$x$', fontsize=15)
ax.set_ylabel(r'$u_0$', fontsize=15)

x=np.linspace(0,20,1000)

# T=10^2K
a=4.7e-3
lga=np.log10(a)
lna=np.log(a)

u0_0=x/(1.01+x/210.0+x*x/105.0)

u0_3=np.zeros(len(x))
for i in range(len(x)):
    if x[i]<3:
       u0_3[i]=0
    else:
       u0_3[i]=1.85-lna/6.73+np.log(np.log(x[i]))
#       u0_3[i]=1.85-lga/6.73+np.log(np.log(x[i]))

z=lga
z2=z*z
u0_4 =(      2.648963 + 2.014446*z + 0.351479*z2
     + x*(-4.058673 - 3.675859*z - 0.640003*z2
     + x*( 3.017395 + 2.117133*z + 0.370294*z2
     + x*(-0.869789 - 0.565886*z - 0.096312*z2
     + x*( 0.110987 + 0.070103*z + 0.011557*z2
     + x*(-0.005200 - 0.003240*z - 0.000519*z2)))))
     )

xcw=6.9184721+81.766279/(lna-14.651253);
u0_5=xcw-1/xcw+0.15*(x-xcw)
for i in range(len(x)):
    if 1<x[i] and x[i]<xcw: u0_5[i]=x[i]-1/(x[i]+np.exp(1-x[i]*x[i])/a)
    if x[i]<=1:             u0_5[i]=0

ax.plot(x,u0_0,color='m')

ax.plot(x[x<10], u0_3[x<10], color='r',linewidth=0.5)
ax.plot(x[x>=10],u0_3[x>=10],color='r',linewidth=0.5,linestyle='dotted')

ax.plot(x[x<8], u0_4[x<8], color='g',linewidth=0.5)
ax.plot(x[x>=8],u0_4[x>=8],color='g',linewidth=0.5,linestyle='dotted')

ax.plot(x[x<9], u0_5[x<9], color='b',linewidth=0.5)
ax.plot(x[x>=9],u0_5[x>=9],color='b',linewidth=0.5,linestyle='dotted')


# T=10^4K
a=4.7e-4
lga=np.log10(a)
lna=np.log(a)

for i in range(len(x)):
    if x[i]<3:
       u0_3[i]=0
    else:
       u0_3[i]=1.85-lna/6.73+np.log(np.log(x[i]))
#       u0_3[i]=1.85-lga/6.73+np.log(np.log(x[i]))

z=lga
z2=z*z
u0_4 =(      2.648963 + 2.014446*z + 0.351479*z2
     + x*(-4.058673 - 3.675859*z - 0.640003*z2
     + x*( 3.017395 + 2.117133*z + 0.370294*z2
     + x*(-0.869789 - 0.565886*z - 0.096312*z2
     + x*( 0.110987 + 0.070103*z + 0.011557*z2
     + x*(-0.005200 - 0.003240*z - 0.000519*z2)))))
     )

xcw=6.9184721+81.766279/(lna-14.651253);
u0_5=xcw-1/xcw+0.15*(x-xcw)
for i in range(len(x)):
    if 1<x[i] and x[i]<xcw: u0_5[i]=x[i]-1/(x[i]+np.exp(1-x[i]*x[i])/a)
    if x[i]<=1:             u0_5[i]=0


ax.plot(x[x<10], u0_3[x<10], color='r',linewidth=1,alpha=0.75)
ax.plot(x[x>=10],u0_3[x>=10],color='r',linewidth=1,alpha=0.75,linestyle='dotted')

ax.plot(x[x<8], u0_4[x<8], color='g',linewidth=1,alpha=0.75)
ax.plot(x[x>=8],u0_4[x>=8],color='g',linewidth=1,alpha=0.75,linestyle='dotted')

ax.plot(x[x<9], u0_5[x<9], color='b',linewidth=1,alpha=0.75)
ax.plot(x[x>=9],u0_5[x>=9],color='b',linewidth=1,alpha=0.75,linestyle='dotted')


# T=10^6K
a=4.7e-5
lga=np.log10(a)
lna=np.log(a)

for i in range(len(x)):
    if x[i]<3:
       u0_3[i]=0
    else:
       u0_3[i]=1.85-lna/6.73+np.log(np.log(x[i]))
#       u0_3[i]=1.85-lga/6.73+np.log(np.log(x[i]))

z=lga
z2=z*z
u0_4 =(      2.648963 + 2.014446*z + 0.351479*z2
     + x*(-4.058673 - 3.675859*z - 0.640003*z2
     + x*( 3.017395 + 2.117133*z + 0.370294*z2
     + x*(-0.869789 - 0.565886*z - 0.096312*z2
     + x*( 0.110987 + 0.070103*z + 0.011557*z2
     + x*(-0.005200 - 0.003240*z - 0.000519*z2)))))
     )

xcw=6.9184721+81.766279/(lna-14.651253);
u0_5=xcw-1/xcw+0.15*(x-xcw)
for i in range(len(x)):
    if 1<x[i] and x[i]<xcw: u0_5[i]=x[i]-1/(x[i]+np.exp(1-x[i]*x[i])/a)
    if x[i]<=1:             u0_5[i]=0


ax.plot(x[x<10], u0_3[x<10], color='r',linewidth=2,alpha=0.5)
ax.plot(x[x>=10],u0_3[x>=10],color='r',linewidth=2,alpha=0.5,linestyle='dotted')

ax.plot(x[x<8], u0_4[x<8], color='g',linewidth=2,alpha=0.5)
ax.plot(x[x>=8],u0_4[x>=8],color='g',linewidth=2,alpha=0.5,linestyle='dotted')

ax.plot(x[x<9], u0_5[x<9], color='b',linewidth=2,alpha=0.5)
ax.plot(x[x>=9],u0_5[x>=9],color='b',linewidth=2,alpha=0.5,linestyle='dotted')


line0 = Line2D([0], [0], label='ZM2002',       color='m')
line3 = Line2D([0], [0], label='Semelin+2007', color='r')
line4 = Line2D([0], [0], label='RASCAS',       color='g')
line5 = Line2D([0], [0], label='Smith+2015',   color='b')

lineT2 = Line2D([0], [0], label=r'$T=10^2{\rm K}\, (a=4.7\times10^{-3})$',color='k',linewidth=0.5,alpha=1)
lineT4 = Line2D([0], [0], label=r'$T=10^4{\rm K}\, (a=4.7\times10^{-4})$',color='k',linewidth=1,alpha=0.75)
lineT6 = Line2D([0], [0], label=r'$T=10^6{\rm K}\, (a=4.7\times10^{-5})$',color='k',linewidth=2,alpha=0.5)

noline = Line2D([0], [0], label='', alpha=0)

legend1=ax.legend(handles=[line0,line3,line5,line4],framealpha=1,loc=[0.52,0.3])
legend2=ax.legend(handles=[lineT2,lineT4,lineT6],framealpha=1,loc=[0.52,0.05])

ax.add_artist(legend1)
ax.add_artist(legend2)

#plt.savefig('u0_cmp.pdf')
plt.savefig('u0_cmp.png')

plt.show()

