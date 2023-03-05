########################################################################################################################################################################
#           Autor: Lukas Scherne, Januar 2023
########################################################################################################################################################################
import time
start_time = time.time()


import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

#eq. 4.6 - caustic
def k(R, s):
    y = (R/(2*np.sqrt(1-s**2)))*(1+s**2-2*s**4)
    z = R*s**3
    return y,z

#eq. 4.5 - reflected edge ray
def RefStrahl(R, d, t):
    p = d/R
    y = R/(2*np.sqrt(1-p**2)) + t*(1-2*p**2)
    z = t * 2*p*(np.sqrt(1-p**2))
    return y, z

def calc_zs(d,r): # ideal radius of the camera
    p = d/r
    p2 = p*p
    return R /(2*np.sqrt(1-p2)) + R/(2*np.sqrt(1-p2))*((p2*(1-np.sqrt(1-p2))*(2*p2-1)) / (1-np.sqrt(1-p2)*pow((4-3*p2),3/2)-p2*(2*p2-1)))

def calc_ys(d,r): #ideal Spotradius at given Mirrorradius
    p = d/r
    p2 = p*p
    p3 = p*p*p
    return R *(p3*(1-np.sqrt(1-p2))/((1-np.sqrt(1-p2))*pow((4-3*p2),3/2)-p2*(2*p2-1)))



#
# N : number of points, rs: mirror radius in m , ds: aperture in m
def NumericalSolution(N, rs, ds):

    s = np.linspace(0, 1, N, endpoint=False)
    yk, zk = k(rs, s)

    #rays
    t = np.linspace(-rs/2, rs/2, N)
    yrs, zrs = RefStrahl(rs, -ds/1., t)

    #calculate intersection of rays and caustic
    dis = np.sqrt((yk[:, np.newaxis] - yrs)**2 + (zk[:, np.newaxis] - zrs)**2)
    min_dis_idx = np.argmin(dis)
    ynum, znum = yk[min_dis_idx//N], zk[min_dis_idx//N]

    return ynum, znum, yrs, zrs, yk, zk


findbestspot =[]
spotindeg = []
#---------------------------- Pierre Auger Werte ------------------------------------------------
rs_pa = 3.4 #m
rk_pa = 1.743 #m
spotsize_pa = 0.5 #deg
pixelsize_pa = 1.5 #deg
aeff_pa = 3 #m^2
mirrorsize_pa =  12 #m^2
kamerasize_pa =  0.8 #m^2
theta_pa = 30 #deg
phi_pa = 30 #deg

sn_pa = np.sqrt(aeff_pa/pixelsize_pa**2)
#-----------------------------------------------------------------------------------------------
N = 3000  #number of points generated

rs = 3.4 #Radius of the mirror

theta = 10 #deg
phi = 30 #degin deg

step_size = 0.05 #m

max_app = rs
lastap = max_app

count = 0

#----------------------------------------------------------------------------------------------


i = 0.1 # m , first aperutre to look at

theta_rad = theta * np.pi/180.
phi_rad = phi * np.pi/180.


#bow string
Dxs_string = 2*np.sin(theta_rad/2)*rs
Dys_string = 2*np.sin(phi_rad/2)*rs

#bow lenght
Dxs = (theta/360)*2*np.pi*rs
Dys = (phi/360)*2*np.pi*rs

ynum = 0
znum = 0
yrs = 0
zrs = 0
yk = 0
zk = 0


Aeff_list=[]
sn_list = []
Appertur_list = []
spots_list = []
ASpiegel_list =[]
AKamera_list = []
RKamera_list = []
Ratio_toPA = []
pixels_needed_list = []
app = []
appindeg = []

#suchschleife
while i < max_app:
    count += 1
    print("----------------------------------------------------------------------")
    print("Radius of tested apertur, in m: ", np.round(i,3))


    #calc numerical NumericalSolution
    ynum,znum,yrs,zrs,yk,zk = NumericalSolution(N, rs, i)
    x = (np.round(znum,10)*1000)

    m = (yrs[N-10]-yrs[0])/(zrs[N-10]-zrs[0])

    gamma = 90 + np.arctan(m) * 180/np.pi #angle between optical axis and edge ray

    if (i==0):
        gamma = 0

    #mirror extension in deg
    delta_alpha_y = np.arcsin(i/(2*rs))*180/np.pi
    delta_alpha_x = np.arcsin(i/(2*rs))*180/np.pi



    if(gamma >= 0 and gamma <= 60):

        #spot in mm and deg
        findbestspot.append(x)
        min_pixel = np.round(2*np.arcsin((znum)/(2*ynum))*(180./np.pi),3)
        Lichtakzeptanz_in_deg = 2*np.arcsin(i/(2*rs)) * 180./np.pi

        rk = ynum
        #rk = rs/2 #fixed nominal radius

        #string lenght
        DKxs_string = 2*np.sin(theta_rad/2)*rk + x/1000
        DKys_string = 2*np.sin(phi_rad/2)*rk + x/1000

        #bow lenght
        DKxs= (theta/360)*2*np.pi*rk
        DKys= (phi/360)*2*np.pi*rk

        #FOV
        phi_kamera = (360* DKys) / (rk * 2*np.pi)
        theta_kamera = (360* DKxs) / (rk * 2*np.pi)

        #Radius and area of camera
        RKamera_list.append(np.round(rk, 4))
        AKamera = DKxs*DKys
        AKamera_list.append(AKamera)

        #area of the mirror
        ASpiegel = (Dxs+2*i)*(Dys+2*i)
        ASpiegelOhneErweiterung = Dxs*Dys
        ASpiegel_list.append(ASpiegel)

        #effektive area
        AeffGrob = np.pi*i**2 - (DKxs)*(DKys)

        if(AeffGrob < 0):
            AeffGrob = 0

        #spotsize
        spotGrob = np.round(((x/1000)*360)/(2*np.pi*rk),4)

        #signal to noise ratio
        if(spotGrob != 0):
            sn = np.sqrt(AeffGrob/((2*spotGrob)**2))
        else:
            sn = 0

        Appertur_list.append(Lichtakzeptanz_in_deg)
        Aeff_list.append(AeffGrob)

        sn_list.append(sn)

        spots_list.append(2*spotGrob)

        if ( x > 0):
            pixels_needed = AKamera / ((2*(x/1000))**2)
        if (x == 0):
            pixels_needed = 0

        pixels_needed_list.append(pixels_needed)

        #sensitivity ratio to auger
        Ratio_toPA.append(sn / sn_pa)

        print("Found optical parameters, configuration", count)
        print("  ")
        print("Ideal apertur in m: ", np.round(i,3), "light acceptence: ", 2*np.round(Lichtakzeptanz_in_deg,3) ,"deg per pixel!")
        print("spotsize in mm (radius): ", x ,"spotsize in deg (radius): ", spotGrob)
        print("camera and mirror parameters: ")
        print(" ")
        print("Mirror parameters: ")
        print("Mirrorradius, in m: ", rs)
        print("Mirror in x, in m: ", np.round(Dxs_string+2*i,4)) #i ist der Radius der Apertur
        print("Mirror in y, in m : ", np.round(Dys_string+2*i,4))
        print("Mirror in x, in deg:", np.round(phi+delta_alpha_x, 4))
        print("Mirror in y, in deg:", np.round(theta+delta_alpha_y, 4))


        print("Mirrorarea without extension, in m²: ", np.round(ASpiegelOhneErweiterung,4))
        print("Mirrorarea, in m²: ", np.round(ASpiegel,4))
        print("  ")
        print("Camera parameters:" )
        print("Cameraradius, in m: ", np.round(rk,4))
        print("Camera in x, in m: ", DKxs_string)
        print("Camera in y, in m: ", DKys_string)
        print("FOV:", theta_kamera, "x", phi_kamera )
        print("Camera area, in m²: ", (DKxs)*(DKys))
        if(pixels_needed != 0):
            print("Pixels needed:", pixels_needed)
            print("Squared Pixel with side lenght in mm :", 2*(x) )
        else:
            print("error: no pixel guess")
        print("  ")
        print("effective imaging area, light inc. 0 deg: ", np.round(AeffGrob, 4))
        print("Sqrt(Aeff/Spotsize**2), light inc. 0 deg: ", np.round(sn,4))
        print("sensitivity-ratio to PA: ", sn / sn_pa)

        if(i!=max_app):
            print("----------------------------------------------------------------------")
            print("still searching..")
        print("  ")
        i += step_size

        #Apertur in m / in deg
        app_v = i
        app.append(i)

        appindeg.append((app_v*360)/(2*np.pi*rs))


    else:
        lastap = i
        i = max_app
        print("max. apertur reached")


#save all data to .txt files
file = np.array([Aeff_list,sn_list,Appertur_list,spots_list,ASpiegel_list,AKamera_list,RKamera_list,Ratio_toPA,pixels_needed_list]).T
#np.savetxt("Kontent_r5_locked.txt",file)
print("Runtime: %s sec." % (time.time() - start_time))
