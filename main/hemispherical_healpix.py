import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import healpy as hp
import math

############ cone parameters
cone_angle=np.pi/2.0
radius_of_separation=np.sin(cone_angle)/np.sin((np.pi-cone_angle)/2.0) 

############ healpix parameters
nside=8


RA,DEC,longt,lat,t,et,v, bt, e_bt, btc, mabs=np.genfromtxt("data/t_and_mabs_btc_lt_15.txt", unpack=True)

# converting to radian  0<theta<pi  0<phi<2*pi
theta=(90-lat)*math.pi/180.0
phi=longt*math.pi/180.0
#***************************************

axes_num=hp.nside2npix(nside)/2

ksd=np.zeros(axes_num)
p_value=np.zeros(axes_num)

f=open("KS_hemispherical.txt","w")
f.write("#longt	  lat	KSD	p-value	  gal_count_in	 gal_count_out"+"\n")
#for i in range(hp.nside2npix(nside)):
for i in range(axes_num):
	print i
	theta_cone, phi_cone=hp.pix2ang(nside,i)
	cone_vec_x,cone_vec_y,cone_vec_z=hp.pix2vec(nside,i)
	gal_count=0
	for j in range(len(t)):
		x_gal, y_gal, z_gal=hp.ang2vec(theta[j], phi[j])
		gal_separation=np.sqrt((cone_vec_x-x_gal)**2+(cone_vec_y-y_gal)**2+(cone_vec_z-z_gal)**2) #calculating the separation of the gal from the r_cone
		if gal_separation < radius_of_separation:
			gal_count=gal_count+1
	t_in=np.zeros(gal_count)
	gal_count_out=len(t)-gal_count
	t_out=np.zeros(gal_count_out)

	gal_count=0
	gal_count_out=0
	for j in range(len(t)):
		x_gal, y_gal, z_gal=hp.ang2vec(theta[j], phi[j])
		gal_separation=np.sqrt((cone_vec_x-x_gal)**2+(cone_vec_y-y_gal)**2+(cone_vec_z-z_gal)**2) #calculating the separation of the gal from the r_cone
		if gal_separation < radius_of_separation:
			t_in[gal_count]=t[j]
			gal_count=gal_count+1
		else:
			t_out[gal_count_out]=t[j]
			gal_count_out=gal_count_out+1
	ksd[i],p_value[i]=stats.ks_2samp(t_in, t_out)
	l=phi_cone*180.0/math.pi
	b=90-(theta_cone*180.0/math.pi)
	f.write(str(l)+"\t"+str(b)+"\t"+str(ksd[i])+"\t"+str(p_value[i])+"\t"+str(gal_count)+"\t"+str(gal_count_out)+"\n")

f.close()

max_ksd_index=np.unravel_index(ksd.argmax(), ksd.shape)
theta_max, phi_max=hp.pix2ang(nside,max_ksd_index[0])
l_max=phi_max*180.0/math.pi
b_max=90-(theta_max*180.0/math.pi)

f2=open("max_KS_hemisphere.txt","w")
f2.write("direction of the largest KS"+"\n")
f2.write(str(l_max)+"\t"+str(b_max)+"\t"+str(max(ksd))+"\t"+str(p_value[max_ksd_index[0]])+"\n")
f2.close()


