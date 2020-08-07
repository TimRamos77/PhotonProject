import numpy as np
import math
import random

class Simulation:

    n1 = 1.5
    n2 = 1.0
    theta_critical = (math.asin(n2 / n1))

    detected_photon = 0
    #for airgap, air_gap = True
    #for no airgap, air_gap = False

    def photon(self, v, Ro, l, w, h, detector, air_gap, rec = 0):

        if rec > 900:
            return False

        for i in range(1):
            if v[i] > 0.:
                x = l/2
                t = (x - Ro[i]) / v[i]
                y = Ro[i+1] + v[i+1] * t
                z = Ro[i+2] + v[i+2] * t
                R = np.array([x, y, z])  #finds the coordinates where the photon hits the x = l/2 plane
            if v[i] < 0.:
                x = -l/2
                t = (x - Ro[i]) / v[i]
                y = Ro[i+1] + v[i+1] * t
                z = Ro[i+2] + v[i+2] * t
                R = np.array([x, y, z])  #finds the coordinates where the photon hits the x = -l/2 plane
            if (-l/2 <= R[i] <= l/2) and (-w/2 <= R[i+1] <= w/2) and (-h/2 <= R[i+2] <= h/2): #checks to see if any of those two points are within the boundaries of the box
                theta_i = (math.acos(abs(v[i]) / math.sqrt(v[i] ** 2 + v[i+1] ** 2 + v[i+2] ** 2)))
                if theta_i > theta_critical:
                    if air_gap:
                        v[i] *= -1
                        Ro = R
                        return photon(v, Ro, l, w, h, detector, air_gap, rec + 1)
                    if not(air_gap):
                        if (detector == 1 and x == l / 2) or (detector == 2 and x == -l / 2):
                            return True
                        else:
                            v[i] *= -1
                            Ro = R
                            return photon(v, Ro, l, w, h, detector, air_gap, rec + 1)
                else:
                    theta_t = math.asin(n1 * math.sin(theta_i))
                    r_perp = (n1 * math.cos(theta_i) - n2 * math.cos(theta_t)) / (n1 * math.cos(theta_i) + n2 * math.cos(theta_t))
                    r_para = (n1 * math.cos(theta_t) - n2 * math.cos(theta_i)) / (n1 * math.cos(theta_t) + n2 * math.cos(theta_i))
                    r_ave = (abs(r_perp) + abs(r_para)) / 2
                    Reflectance = (r_ave ** 2) * 100
                    Transmittance = 100 - Reflectance
                    select_path = random.uniform(0, 101)
                    Min = min(Transmittance, Reflectance)
                    if air_gap:
                        if 0 <= select_path <= Min:
                            if Min == Reflectance:
                                v[i] *= -1
                                Ro = R
                                return photon(v, Ro, l, w, h, detector, air_gap, rec + 1)
                            if Min == Transmittance:
                                if (detector == 1 and x == l / 2) or (detector == 2 and x == -l / 2):
                                    return True
                                else:
                                    return False
                        else:
                            if Min == Reflectance:
                                if (detector == 1 and x == l / 2) or (detector == 2 and x == -l / 2):
                                    return True
                                else:
                                    return False
                            if Min == Transmittance:
                                v[i] *= -1
                                Ro = R
                                return photon(v, Ro, l, w, h, detector, air_gap, rec + 1)
                    if not air_gap:
                        if (detector == 1 and x == l / 2) or (detector == 2 and x == -l / 2):
                            return True
                        else:
                            if 0 <= select_path <= Min:
                                if Min == Reflectance:
                                    v[i] *= -1
                                    Ro = R
                                    return photon(v, Ro, l, w, h, detector, air_gap, rec + 1)
                                if Min == Transmittance:
                                    return False
                            else:
                                if Min == Reflectance:
                                    return False
                                if Min == Transmittance:
                                    v[i] *= -1
                                    Ro = R
                                    return photon(v, Ro, l, w, h, detector, air_gap, rec + 1)


        for i in range(1,2):
            if v[i] > 0.:
                y = w/2
                t = (y - Ro[i]) / v[i]
                x = Ro[i-1] + v[i-1] * t
                z = Ro[i+1] + v[i+1] * t
                R = np.array([x, y, z])
            if v[i] < 0.:
                y = -w/2
                t = (y - Ro[i]) / v[i]
                x = Ro[i-1] + v[i-1] * t
                z = Ro[i+1] + v[i+1] * t
                R = np.array([x, y, z])
            if (-l/2 <= R[i-1] <= l/2) and (-w/2 <= R[i] <= w/2) and (-h/2 <= R[i+1] <= h/2):
                theta_i = (math.acos(abs(v[i]) / math.sqrt(v[i-1] ** 2 + v[i] ** 2 + v[i+1] ** 2)))
                if theta_i > theta_critical:
                    if air_gap:
                        v[i] *= -1
                        Ro = R
                        return photon(v, Ro, l, w, h, detector, air_gap, rec + 1)
                    if not(air_gap):
                        if (detector == 3 and y == w / 2) or (detector == 4 and y == -w / 2):
                            return True
                        else:
                            v[i] *= -1
                            Ro = R
                            return photon(v, Ro, l, w, h, detector, air_gap, rec + 1)
                else:
                    theta_t = math.asin(n1 * math.sin(theta_i))
                    r_perp = (n1 * math.cos(theta_i) - n2 * math.cos(theta_t)) / (n1 * math.cos(theta_i) + n2 * math.cos(theta_t))
                    r_para = (n1 * math.cos(theta_t) - n2 * math.cos(theta_i)) / (n1 * math.cos(theta_t) + n2 * math.cos(theta_i))
                    r_ave = (abs(r_perp) + abs(r_para)) / 2
                    Reflectance = (r_ave ** 2) * 100
                    Transmittance = 100 - Reflectance
                    select_path = random.uniform(0, 100)
                    Min = min(Transmittance, Reflectance)
                    if air_gap:
                        if 0 <= select_path <= Min:
                            if Min == Reflectance:
                                v[i] *= -1
                                Ro = R
                                return photon(v, Ro, l, w, h, detector, air_gap, rec + 1)
                            if Min == Transmittance:
                                if (detector == 3 and y == w / 2) or (detector == 4 and y == -w / 2):
                                    return True
                                else:
                                    return False
                        else:
                            if Min == Reflectance:
                                if (detector == 3 and y == w / 2 ) or (detector == 4 and y == -w / 2):
                                    return True
                                else:
                                    return False
                            if Min == Transmittance:
                                v[i] *= -1
                                Ro = R
                                return photon(v, Ro, l, w, h, detector, air_gap, rec + 1)
                    if not air_gap:
                        if (detector == 3 and y == w / 2) or (detector == 4 and y == -w / 2):
                            return True
                        else:
                            if 0 <= select_path <= Min:
                                if Min == Reflectance:
                                    v[i] *= -1
                                    Ro = R
                                    return photon(v, Ro, l, w, h, detector, air_gap, rec + 1)
                                if Min == Transmittance:
                                    return False
                            else:
                                if Min == Reflectance:
                                    return False
                                if Min == Transmittance:
                                    v[i] *= -1
                                    Ro = R
                                    return photon(v, Ro, l, w, h, detector, air_gap, rec + 1)
        for i in range(2,3):
            if v[i] > 0.:
                z = h/2
                t = (z - Ro[i]) / v[i]
                x = Ro[i-2] + v[i-2] * t
                y = Ro[i-1] + v[i-1] * t
                R = np.array([x, y, z])
            if v[i] < 0.:
                z = -h/2
                t = (z - Ro[i]) / v[i]
                x = Ro[i-2] + v[i-2] * t
                y = Ro[i-1] + v[i-1] * t
                R = np.array([x, y, z])
            if (-l/2 <= R[i-2] <= l/2) and (-w/2 <= R[i-1] <= w/2) and (-h/2 <= R[i] <= h/2):
                theta_i = (math.acos(abs(v[i]) / math.sqrt(v[i-2] ** 2 + v[i-1] ** 2 + v[i] ** 2)))
                if theta_i > theta_critical:
                    if air_gap:
                        v[i] *= -1
                        Ro = R
                        return photon(v, Ro, l, w, h, detector, air_gap, rec + 1)
                    if not(air_gap):
                        if (detector == 5 and z == h / 2) or (detector == 6 and z == -h / 2):
                            return True
                        else:
                            v[i] *= -1
                            Ro = R
                            return photon(v, Ro, l, w, h, detector, air_gap, rec + 1)
                else:
                    theta_t = math.asin(n1 * math.sin(theta_i))
                    r_perp = (n1 * math.cos(theta_i) - n2 * math.cos(theta_t)) / (n1 * math.cos(theta_i) + n2 * math.cos(theta_t))
                    r_para = (n1 * math.cos(theta_t) - n2 * math.cos(theta_i)) / (n1 * math.cos(theta_t) + n2 * math.cos(theta_i))
                    r_ave = (abs(r_perp) + abs(r_para)) / 2
                    Reflectance = (r_ave ** 2) * 100
                    Transmittance = 100 - Reflectance
                    select_path = random.uniform(0, 101)
                    Min = min(Transmittance, Reflectance)
                    if air_gap:
                        if 0 <= select_path <= Min:
                            if Min == Reflectance:
                                v[i] *= -1
                                Ro = R
                                return photon(v, Ro, l, w, h, detector, air_gap, rec + 1)
                            if Min == Transmittance:
                                if (detector == 5 and z == h / 2) or (detector == 6 and z == -h / 2):
                                    return True
                                else:
                                    return False
                        else:
                            if Min == Reflectance:
                                if (detector == 5 and z == h / 2) or (detector == 6 and z == -h / 2):
                                    return True
                                else:
                                    return False
                            if Min == Transmittance:
                                v[i] *= -1
                                Ro = R
                                return photon(v, Ro, l, w, h, detector, air_gap, rec + 1)
                    if not air_gap:
                        if (detector == 5 and z == h / 2) or (detector == 6 and z == -h / 2):
                            return True
                        else:
                            if 0 <= select_path <= Min:
                                if Min == Reflectance:
                                    v[i] *= -1
                                    Ro = R
                                    return photon(v, Ro, l, w, h, detector, air_gap, rec + 1)
                                if Min == Transmittance:
                                    return False
                            else:
                                if Min == Reflectance:
                                    return False
                                if Min == Transmittance:
                                    v[i] *= -1
                                    Ro = R
                                    return photon(v, Ro, l, w, h, detector, air_gap, rec + 1)



    t = None
    length = None

    Xoy = 0
    Xoz = 0
    phi_line = random.uniform(-math.pi/2, math.pi/2)   #-pi/2 <= phi_line <= pi/2
    theta_line = math.acos(random.uniform(-1., 1.))  #0 <= theta_line <= pi
    l = 20
    w = 10
    h = 4

    tmin_x = -l/2
    tmax_x = l/2
    tmin_y = ((-w/2)/math.tan(phi_line)) + Xoy
    tmax_y = ((w/2)/math.tan(phi_line)) + Xoy
    tmin_z = ((-h/2)/math.tan((math.pi)/2 - theta_line)) + Xoz
    tmax_z = ((h/2)/math.tan((math.pi)/2 - theta_line)) + Xoz

    ts = [tmin_x, tmax_x, tmin_y, tmax_y, tmin_z, tmax_z]
    ts.sort()

    for i in range(2,4):
        if ((-l/2 <= ts[i] <= l/2) and (-w/2 <= (ts[i]-Xoy)*math.tan(phi_line) <= w/2) and
                (-h/2 <= (ts[i]-Xoz)*math.tan((math.pi/2)-theta_line) <= h/2)):
            t = random.uniform(ts[2], ts[3])
            length = (ts[3] - ts[2]) * math.sqrt(1 + (math.tan(phi_line)) ** 2 + (math.tan((math.pi / 2) - theta_line)) ** 2)
        else:
            print("t is not defined")
            break


    for n in range(1000):
        if t is None:
            break

        phi_initial = random.uniform(0., 2*math.pi)
        theta_initial = math.acos(random.uniform(-1., 1.))
        Vx = math.sin(theta_initial) * math.cos(phi_initial)
        Vy = math.sin(theta_initial) * math.sin(phi_initial)
        Vz = math.cos(theta_initial)

        t = random.uniform(ts[2], ts[3])
        Rox = t #x(t)
        Roy = (t-Xoy) * math.tan(phi_line) #y(t)
        Roz = (t-Xoz) * math.tan((math.pi/2)-theta_line) #z(t)

        detected_photon += int(photon(np.array([Vx, Vy, Vz]), np.array([Rox, Roy, Roz]), l, w, h, 1, False))


    print("Hits the detector " + str(detected_photon) + " times")
    print("Does not hit the detector " + str(1000 - detected_photon) + " times")
    print("Light efficiency = " + str(detected_photon/10) + "%")
    print(length)
    print(phi_line * (180/math.pi))
    print(theta_line * (180/math.pi))
    print(Xoy)
    print(Xoz)



    #if detector is placed at the positive x plane, detector = 1
    #if detector is placed at the negative x plane, detector = 2
    #if detector is placed at the positive y plane, detector = 3
    #if detector is placed at the negative y plane, detector = 4
    #if detector is placed at the positive z plane, detector = 5
    #if detector is placed at the negative z plane, detector = 6




