import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from radar_simulate.rayleigh_functions import scatSpheroid

# electric field interaction from neighbor dipole 2 (h-pol=y, v-pol=z)
def kernal(rdip1, rdip2, ef2, alph_dip):
    x1 = rdip1[0]
    y1 = rdip1[1]
    z1 = rdip1[2]
    x2 = rdip2[0]
    y2 = rdip2[1]
    z2 = rdip2[2]

    rd = np.sqrt((x2-x1)**2.+(y2-y1)**2.+(z2-z1)**2.)
    enf_hpol = alph_dip*ef2*(3.*((y2-y1)/rd)**2.-1.)/(4.*np.pi*rd**3.)
    enf_vpol = alph_dip*ef2*(3.*((z2-z1)/rd)**2.-1.)/(4.*np.pi*rd**3.)
    return np.real(enf_hpol), np.real(enf_vpol)

# angle between two dipoles
def angle(rdip1, rdip2):
    x1 = rdip1[0]
    y1 = rdip1[1]
    z1 = rdip1[2]
    x2 = rdip2[0]
    y2 = rdip2[1]
    z2 = rdip2[2]

    rd = np.sqrt((x2-x1)**2.+(y2-y1)**2.+(z2-z1)**2.)
    ang_hpol = (3.*((y2-y1)/rd)**2.-1.)/(4.*np.pi*rd**3.)
    ang_vpol = (3.*((z2-z1)/rd)**2.-1.)/(4.*np.pi*rd**3.)
    return ang_hpol, ang_vpol

# open dda input file
data = np.genfromtxt('crystal3.0.txt', skip_header=3)
print data.shape
xint = data[:,0]
yint = data[:,1]
zint = data[:,2]
numdip = len(xint)

# estimate rayleigh gans scattering
dmax = 3.
wavl = 32.1
k = np.pi*2./wavl
eps_ice = complex(3.16835, 0.0089)
nx = np.max(xint)-np.min(xint)
nz = np.max(zint)-np.min(zint)
asp = float(nx)/float(nz)
thick = dmax/asp

print nx
dip_len = dmax/float(nx)
dip_vol = dip_len**3.

dipmom_sphere = 4.*np.pi*dip_vol*(eps_ice-1)/(eps_ice+2.)
dipmom_cube = dipmom_sphere*3./(4.*np.pi)
dipmom_tot = dipmom_cube*numdip
s = k**2.*dipmom_tot/(4.*np.pi)
print s

# check compared to sphere of same volume
vol = dip_vol*numdip
rad = (vol*3./(4.*np.pi))**(1./3.)
shh, svv = scatSpheroid(eps_ice, 2.*rad, 2.*rad, wavl)
print shh[0]

# create horizontal slice of crystal
xf = xint*dip_len
yf = yint*dip_len
zf = zint*dip_len
xf = xf-np.mean(xf)
yf = yf-np.mean(yf)
zf = zf-np.mean(zf)
print np.min(xf), np.max(xf)
print np.min(yf), np.max(yf)
print np.min(zf), np.max(zf)
zmed = np.median(zf)

x2d = xf[zf==zmed]
y2d = yf[zf==zmed]
z2d = zf[zf==zmed]
num2d = len(x2d)

'''
# apply kernal to one dipole
enf_h = np.empty([num2d])
enf_v = np.empty([num2d])
ang_h = np.empty([num2d])
ang_v = np.empty([num2d])
dipind = np.random.randint(0, num2d, size=1)
x1 = x2d[dipind]
y1 = y2d[dipind]
z1 = z2d[dipind]
r1 = np.array([x1, y1, z1])
for i in range(num2d):
    r2 = np.array([x2d[i], y2d[i], z2d[i]])
    enf_h[i], enf_v[i] = kernal(r1, r2, 1., dipmom_cube)
'''

# apply kernal to all dipoles
enf_h_tot = np.empty([num2d])
enf_v_tot = np.empty([num2d])
for i in range(num2d):
    r1 = np.array([x2d[i], y2d[i], z2d[i]])
    for j in range(num2d):
        r2 = np.array([x2d[j], y2d[j], z2d[j]])
        if (i!=j):
            enf_h, enf_v = kernal(r1, r2, 1., dipmom_cube)
            enf_h_tot[i] = enf_h_tot[i]+enf_h
            enf_v_tot[i] = enf_v_tot[i]+enf_v

# plot dipole locations
plt.scatter(x2d, y2d, s=8., c=(1.+enf_h)*100., cmap='gist_rainbow', marker='.')
#plt.scatter(x2d, y2d, s=6., c=ang_h, cmap='gist_rainbow')
plt.colorbar()

ax = plt.gca()
ax.set_aspect(1.)
plt.savefig('diploc.png')
