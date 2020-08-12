import math
import numpy as np
import random

#T, phi_line, theta_line, Xoy, Xoz, V, Ro



class Simulation:

    def __init__(self, l, w, h, n1, n2, detector, air_gap, rec=0, phi_line = math.pi/4,
                 theta_line = math.pi/4, Xoy = 0, Xoz = 0):
        self.l = l
        self.w = w
        self.h = h
        self.n1 = n1
        self.n2 = n2
        self.detector = detector
        self.air_gap = air_gap
        self.rec = rec
        self.phi_line = phi_line
        self.theta_line = theta_line
        self.Xoy = Xoy
        self.Xoz = Xoz
        self.theta_critical = (math.asin(n2 / n1))


    def random_line(self):

        Tmin_x = -self.l / 2
        Tmax_x = self.l / 2
        Tmin_y = ((-self.w / 2) / math.tan(self.phi_line)) + self.Xoy
        Tmax_y = ((self.w / 2) / math.tan(self.phi_line)) + self.Xoy
        Tmin_z = ((-self.h / 2) / math.tan((math.pi) / 2 - self.theta_line)) + self.Xoz
        Tmax_z = ((self.h / 2) / math.tan((math.pi) / 2 - self.theta_line)) + self.Xoz
        # checks the T boundaries of each x, y, and z based on the dimensions of the box, l, w, and h

        self.Ts = [Tmin_x, Tmax_x, Tmin_y, Tmax_y, Tmin_z, Tmax_z]
        self.Ts.sort()
        # sorts out all possible T values; the two middle values after sorting are the Tmin and Tmax that would satisfy
        # the range of T's of the generated random line

        for i in range(2, 4):
            if ((-self.l / 2 <= self.Ts[i] <= self.l / 2) and (
                    -self.w / 2 <= (self.Ts[i] - self.Xoy) * math.tan(self.phi_line) <= self.w / 2) and
                    (-self.h / 2 <= (self.Ts[i] - self.Xoz) * math.tan((math.pi / 2) - self.theta_line) <= self.h / 2)):
                self.T = random.uniform(self.Ts[2], self.Ts[3])
                self.length = (self.Ts[3] - self.Ts[2]) * math.sqrt(
                    1 + (math.tan(self.phi_line)) ** 2 + (math.tan((math.pi / 2) - self.theta_line)) ** 2)
            else:
                print("T is not defined")
                return False


    def photon(self, V, Ro):

        self.rec += 1

        if self.rec > 900:
            return False

        for i in range(1):
            if V[i] > 0.:
                x = self.l/2
                t = (x - Ro[i]) / V[i]
                y = Ro[i+1] + V[i+1] * t
                z = Ro[i+2] + V[i+2] * t
                R = np.array([x, y, z])  #finds the coordinates where the photon hits the x = l/2 plane
            if V[i] < 0.:
                x = -self.l/2
                t = (x - Ro[i]) / V[i]
                y = Ro[i+1] + V[i+1] * t
                z = Ro[i+2] + V[i+2] * t
                R = np.array([x, y, z])  #finds the coordinates where the photon hits the x = -l/2 plane
            if (-self.l/2 <= R[i] <= self.l/2) and (-self.w/2 <= R[i+1] <= self.w/2) and (-self.h/2 <= R[i+2] <= self.h/2): #checks to see if any of those two points are within the boundaries of the box
                theta_i = (math.acos(abs(V[i]) / math.sqrt(V[i] ** 2 + V[i+1] ** 2 + V[i+2] ** 2)))
                if theta_i > self.theta_critical:
                    if self.air_gap:
                        V[i] *= -1
                        Ro = R
                        return self.photon(V, Ro)
                    if not(self.air_gap):
                        if (self.detector == 1 and x == self.l / 2) or (self.detector == 2 and x == -self.l / 2):
                            return True
                        else:
                            V[i] *= -1
                            Ro = R
                            return self.photon(V, Ro)
                else:
                    theta_t = math.asin(self.n1 * math.sin(theta_i))
                    r_perp = (self.n1 * math.cos(theta_i) - self.n2 * math.cos(theta_t)) / (self.n1 * math.cos(theta_i) + self.n2 * math.cos(theta_t))
                    r_para = (self.n1 * math.cos(theta_t) - self.n2 * math.cos(theta_i)) / (self.n1 * math.cos(theta_t) + self.n2 * math.cos(theta_i))
                    r_ave = (abs(r_perp) + abs(r_para)) / 2
                    Reflectance = (r_ave ** 2) * 100
                    Transmittance = 100 - Reflectance
                    select_path = random.uniform(0, 101)
                    Min = min(Transmittance, Reflectance)
                    if self.air_gap:
                        if 0 <= select_path <= Min:
                            if Min == Reflectance:
                                V[i] *= -1
                                Ro = R
                                return self.photon(V, Ro)
                            if Min == Transmittance:
                                if (self.detector == 1 and x == self.l / 2) or (self.detector == 2 and x == -self.l / 2):
                                    return True
                                else:
                                    return False
                        else:
                            if Min == Reflectance:
                                if (self.detector == 1 and x == self.l / 2) or (self.detector == 2 and x == -self.l / 2):
                                    return True
                                else:
                                    return False
                            if Min == Transmittance:
                                V[i] *= -1
                                Ro = R
                                return self.photon(V, Ro)
                    if not self.air_gap:
                        if (self.detector == 1 and x == self.l / 2) or (self.detector == 2 and x == -self.l / 2):
                            return True
                        else:
                            if 0 <= select_path <= Min:
                                if Min == Reflectance:
                                    V[i] *= -1
                                    Ro = R
                                    return self.photon(V, Ro)
                                if Min == Transmittance:
                                    return False
                            else:
                                if Min == Reflectance:
                                    return False
                                if Min == Transmittance:
                                    V[i] *= -1
                                    Ro = R
                                    return self.photon(V, Ro)


        for i in range(1,2):
            if V[i] > 0.:
                y = self.w/2
                t = (y - Ro[i]) / V[i]
                x = Ro[i-1] + V[i-1] * t
                z = Ro[i+1] + V[i+1] * t
                R = np.array([x, y, z])
            if V[i] < 0.:
                y = -self.w/2
                t = (y - Ro[i]) / V[i]
                x = Ro[i-1] + V[i-1] * t
                z = Ro[i+1] + V[i+1] * t
                R = np.array([x, y, z])
            if (-self.l/2 <= R[i-1] <= self.l/2) and (-self.w/2 <= R[i] <= self.w/2) and (-self.h/2 <= R[i+1] <= self.h/2):
                theta_i = (math.acos(abs(V[i]) / math.sqrt(V[i-1] ** 2 + V[i] ** 2 + V[i+1] ** 2)))
                if theta_i > self.theta_critical:
                    if self.air_gap:
                        V[i] *= -1
                        Ro = R
                        return self.photon(V, Ro)
                    if not(self.air_gap):
                        if (self.detector == 3 and y == self.w / 2) or (self.detector == 4 and y == -self.w / 2):
                            return True
                        else:
                            V[i] *= -1
                            Ro = R
                            return self.photon(V, Ro)
                else:
                    theta_t = math.asin(self.n1 * math.sin(theta_i))
                    r_perp = (self.n1 * math.cos(theta_i) - self.n2 * math.cos(theta_t)) / (self.n1 * math.cos(theta_i) + self.n2 * math.cos(theta_t))
                    r_para = (self.n1 * math.cos(theta_t) - self.n2 * math.cos(theta_i)) / (self.n1 * math.cos(theta_t) + self.n2 * math.cos(theta_i))
                    r_ave = (abs(r_perp) + abs(r_para)) / 2
                    Reflectance = (r_ave ** 2) * 100
                    Transmittance = 100 - Reflectance
                    select_path = random.uniform(0, 100)
                    Min = min(Transmittance, Reflectance)
                    if self.air_gap:
                        if 0 <= select_path <= Min:
                            if Min == Reflectance:
                                V[i] *= -1
                                Ro = R
                                return self.photon(V, Ro)
                            if Min == Transmittance:
                                if (self.detector == 3 and y == self.w / 2) or (self.detector == 4 and y == -self.w / 2):
                                    return True
                                else:
                                    return False
                        else:
                            if Min == Reflectance:
                                if (self.detector == 3 and y == self.w / 2 ) or (self.detector == 4 and y == -self.w / 2):
                                    return True
                                else:
                                    return False
                            if Min == Transmittance:
                                V[i] *= -1
                                Ro = R
                                return self.photon(V, Ro)
                    if not self.air_gap:
                        if (self.detector == 3 and y == self.w / 2) or (self.detector == 4 and y == -self.w / 2):
                            return True
                        else:
                            if 0 <= select_path <= Min:
                                if Min == Reflectance:
                                    V[i] *= -1
                                    Ro = R
                                    return self.photon(V, Ro)
                                if Min == Transmittance:
                                    return False
                            else:
                                if Min == Reflectance:
                                    return False
                                if Min == Transmittance:
                                    V[i] *= -1
                                    Ro = R
                                    return self.photon(V, Ro)

        for i in range(2,3):
            if V[i] > 0.:
                z = self.h/2
                t = (z - Ro[i]) / V[i]
                x = Ro[i-2] + V[i-2] * t
                y = Ro[i-1] + V[i-1] * t
                R = np.array([x, y, z])
            if V[i] < 0.:
                z = -self.h/2
                t = (z - Ro[i]) / V[i]
                x = Ro[i-2] + V[i-2] * t
                y = Ro[i-1] + V[i-1] * t
                R = np.array([x, y, z])
            if (-self.l/2 <= R[i-2] <= self.l/2) and (-self.w/2 <= R[i-1] <= self.w/2) and (-self.h/2 <= R[i] <= self.h/2):
                theta_i = (math.acos(abs(V[i]) / math.sqrt(V[i-2] ** 2 + V[i-1] ** 2 + V[i] ** 2)))
                if theta_i > self.theta_critical:
                    if self.air_gap:
                        V[i] *= -1
                        Ro = R
                        return self.photon(V, Ro)
                    if not(self.air_gap):
                        if (self.detector == 5 and z == self.h / 2) or (self.detector == 6 and z == -self.h / 2):
                            return True
                        else:
                            V[i] *= -1
                            Ro = R
                            return self.photon(V, Ro)
                else:
                    theta_t = math.asin(self.n1 * math.sin(theta_i))
                    r_perp = (self.n1 * math.cos(theta_i) - self.n2 * math.cos(theta_t)) / (self.n1 * math.cos(theta_i) + self.n2 * math.cos(theta_t))
                    r_para = (self.n1 * math.cos(theta_t) - self.n2 * math.cos(theta_i)) / (self.n1 * math.cos(theta_t) + self.n2 * math.cos(theta_i))
                    r_ave = (abs(r_perp) + abs(r_para)) / 2
                    Reflectance = (r_ave ** 2) * 100
                    Transmittance = 100 - Reflectance
                    select_path = random.uniform(0, 101)
                    Min = min(Transmittance, Reflectance)
                    if self.air_gap:
                        if 0 <= select_path <= Min:
                            if Min == Reflectance:
                                V[i] *= -1
                                Ro = R
                                return self.photon(V, Ro)
                            if Min == Transmittance:
                                if (self.detector == 5 and z == self.h / 2) or (self.detector == 6 and z == -self.h / 2):
                                    return True
                                else:
                                    return False
                        else:
                            if Min == Reflectance:
                                if (self.detector == 5 and z == self.h / 2) or (self.detector == 6 and z == -self.h / 2):
                                    return True
                                else:
                                    return False
                            if Min == Transmittance:
                                V[i] *= -1
                                Ro = R
                                return self.photon(V, Ro)
                    if not self.air_gap:
                        if (self.detector == 5 and z == self.h / 2) or (self.detector == 6 and z == -self.h / 2):
                            return True
                        else:
                            if 0 <= select_path <= Min:
                                if Min == Reflectance:
                                    V[i] *= -1
                                    Ro = R
                                    return self.photon(V, Ro)
                                if Min == Transmittance:
                                    return False
                            else:
                                if Min == Reflectance:
                                    return False
                                if Min == Transmittance:
                                    V[i] *= -1
                                    Ro = R
                                    return self.photon(V, Ro)


    def run(self, detected_photon=0):
        self.random_line()
        for n in range(1000):
            if self.random_line() == False:
                break

            phi_initial = random.uniform(0., 2 * math.pi)
            theta_initial = math.acos(random.uniform(-1., 1.))
            Vx = math.sin(theta_initial) * math.cos(phi_initial)
            Vy = math.sin(theta_initial) * math.sin(phi_initial)
            Vz = math.cos(theta_initial)
            V = np.array([Vx, Vy, Vz])
            # generates random direction (V) based on random angles of phi and theta

            self.T = random.uniform(self.Ts[2], self.Ts[3])
            Rox = self.T  # x(t)
            Roy = (self.T - self.Xoy) * math.tan(self.phi_line)  # y(t)
            Roz = (self.T - self.Xoz) * math.tan((math.pi / 2) - self.theta_line)  # z(t)
            Ro = np.array([Rox, Roy, Roz])
            # generates random point (Ro) on a line based on the angle phi and theta of that line

            detected_photon += int(self.photon(V, Ro))

        print("Hits the detector " + str(detected_photon) + " times")
        print("Does not hit the detector " + str(1000 - detected_photon) + " times")
        print("Light efficiency = " + str(detected_photon / 10) + "%")

sim1 = Simulation(20.0, 10.0, 4.0, 1.5, 1.0, 1, False)

sim1.run()

