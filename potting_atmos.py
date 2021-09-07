import matplotlib.pyplot as plt
import numpy as np
import std_atmos_SI
import std_atmos_Eng

Engl_altitude = np.linspace(0,175000,1001)
altitude = np.linspace(0,70000,501)
print (altitude)

#SI
""" geop_alt = std_atmos_SI.interp(altitude)[0]
temp = std_atmos_SI.interp(altitude)[1]
pres = std_atmos_SI.interp(altitude)[2]
den = std_atmos_SI.interp(altitude)[3]
speed_of_sound = std_atmos_SI.interp(altitude)[4]

plt.close('all')
plt.figure(1)
plt.plot(altitude,geop_alt)
plt.xlabel('altitude (m)')
plt.ylabel('Geopotential Altitude (m)')
plt.show() """

#English
geop_alt = std_atmos_Eng.interp(Engl_altitude)[0]
temp = std_atmos_Eng.interp(Engl_altitude)[1]
pres = std_atmos_Eng.interp(Engl_altitude)[2]
den = std_atmos_Eng.interp(Engl_altitude)[3]
speed_of_sound = std_atmos_Eng.interp(Engl_altitude)[4]

plt.close('all')
plt.figure(2)
plt.plot(Engl_altitude,speed_of_sound)
plt.xlabel('altitude (ft)')
plt.ylabel('Speed of Sound (ft/s)')
plt.show()
