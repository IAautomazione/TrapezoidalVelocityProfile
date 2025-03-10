"""
Author: Archetti Ivan
Date: 01/12/2024

Class to calculate different motion profile
"""
#=========================================================================================================================

import numpy as np
from numpy import sin, cos, pi


#=========================================================================================================================

class MotionProfile():
    def __init__(self):
        pass
        

#=========================================================================================================================

    def trapezoidal_MP(self, Ds=1, s0=0, Dt=1, ti=0, n_points=100, shape=[0.2, 0.6, 0.2]) -> tuple:
        """
        Calculate a trapezoidal law of motion (with constant acceleration) of a straight path.
        Sum of shape values must be always 1
        
        :param Ds: total ammount of space to travel
               Dt: travel time
               s0: initial space
               ti: initial time
               n_points: number of points for calculation
               shape: parameters that determinate the shape of the law of motion
        :return: a tuple of 4 arrays (time, space, speed, acceleration) during the LOM
        """
        
        try:
            # check for trapezoidal coefficients
            trapezoid = np.round((sum(shape)), 3) 
            if trapezoid != 1:
                raise ValueError(f"\n\n Shape's parameters are not correct ({trapezoid})!\n\n") 
        except Exception as err:
            print("-"*50)
            print("Invalid values! Shape is set to [0.2, 0.6, 0.2]")
            print(f"{err}")
            print("-"*50)
            shape = [0.2, 0.6, 0.2]
            
        # definition of time slots
        rnd = 4
        dt = Dt/(n_points)
        t1 = int(n_points*np.round(shape[0], rnd))
        t2 = int(n_points*np.round(shape[1], rnd))
        t3 = int(n_points*np.round(shape[2], rnd))

        t = np.array([]) # array of time points
        s = np.array([]) # array of space points 
        v = np.array([]) # array of speed points 
        a = np.array([]) # array of acceleration points 

        # initialization of the magnitudes
        space = 0
        v_max = 2*Ds/((shape[0]+2*shape[1]+shape[2])*Dt)

        # first section
        if shape[0] != 0:
            v_tmp = 0
            a_max_pos = v_max/(shape[0]*Dt)
            for i in range(t1):
                t = np.append(t, ti + dt*i)
                s = np.append(s, space + s0)
                v = np.append(v, v_tmp) 
                a = np.append(a, a_max_pos)
                space = space + v_tmp*dt + 1/2*a_max_pos*dt**2
                v_tmp = v_tmp + a_max_pos * dt
        
        # second section
        if shape[1] != 0:
            v_tmp = v_max
            for i in range(t1, t1+t2):
                t = np.append(t, ti + dt*i)
                s = np.append(s, space + s0)
                v = np.append(v, v_tmp)
                a = np.append(a, 0)
                space = space + v_tmp*dt
            
        # third section
        if shape[2] != 0:
            v_tmp = v_max
            a_max_neg = -v_max/(shape[2]*Dt)
            for i in range(t1+t2, t1+t2+t3):
                t = np.append(t, ti + dt*i)
                s = np.append(s, space + s0)
                v = np.append(v, v_tmp)
                a = np.append(a, a_max_neg)
                space = space + v_tmp*dt + 1/2*a_max_neg*dt**2
                v_tmp = v_tmp + a_max_neg * dt
                
        return t, s, v, a

#=========================================================================================================================

    def get_line_angle(self, P0=(0,0), P1=(1,1)) -> float:
        """
        Calculate the angle of a line between two points
        
        :param P0: tuple of starting point
               P1: tuple of end point
        :return: angle of the line (float)
        """

        return np.arctan((P1[1]-P0[1])/(P1[0]-P0[0]))

#------------------------------------------------------------------------------------------------------------------------

    def get_line_position(self, s=np.array([0]), angle=pi/2) -> tuple:
        """
        Calculate the coordinates (x,y) of points from an array s
        
        :param s: array of positions
               angle: angle of the line (float)
        :return: a tuple of 2 arrays (position x and y) 
        """

        sx = np.array([])
        sy = np.array([])

        for i in s:
            sx = np.append(sx, i*cos(angle))
            sy = np.append(sy, i*sin(angle))
        
        return sx, sy
        
#=========================================================================================================================

    def get_line_speed(self, v=np.array([0]), angle=pi/2) -> tuple:
        """
        Calculate speed in x and y of a vector v
        
        :param s: vector of speed
               angle: angle of the line (float)
        :return: a tuple of 2 arrays (speed x and y)  
        """

        vx = np.array([])
        vy = np.array([])

        for i in v:
            vx = np.append(vx, i*cos(angle))
            vy = np.append(vy, i*sin(angle))
        
        return vx, vy
    
#------------------------------------------------------------------------------------------------------------------------

    def get_line_acceleration(self, a=np.array([0]), angle=pi/2) -> tuple:
        """
        Calculate acceleration in x and y of a vector a
        
        :param s: vector of acceleration
               angle: angle of the line (float)
        :return: a tuple of 2 arrays (acceleration x and y)    
        """

        ax = np.array([])
        ay = np.array([])

        for i in a:
            ax = np.append(ax, i*cos(angle))
            ay = np.append(ay, i*sin(angle))
        
        return ax, ay
    
#=========================================================================================================================

    def arc_motion(self, kin_param, radius=1, start_angle=0) -> tuple:
        """
        Calculate a law of motion (by kin_param) for an arc.
        
        :param kin_param: a tuple of 4 array (time, space, speed, acceleration) calculateted
                          from an equivalent straight path with same length
               radius: radius of the arc
               start_angle: starting angle of the arc
        :return: a tuple of 6 arrays (time, space, speed, acceleration, tangential and centripetal acceleration) during the LOM
        """

        # check radius value
        if radius == 0:
            print(f"\n\n Radius must be not 0!\n The value will be forced to 1 \n") 
            radius = 1

        # calculation of the kinematic magnitudes
        t, s_tmp, v_tmp, a_tmp = kin_param

        # initialization
        s = np.array([s_tmp[0]])
        v = np.array([])
        ac = np.array([])
        at = np.array([])
        a = np.array([])

        speed_angle = pi/2 - start_angle
        dv_x0 = -v_tmp[0]*cos(speed_angle)
        dv_y0 = v_tmp[0]*sin(speed_angle)
        v = np.append(v, (dv_x0**2 + dv_y0**2)**(1/2))

        dat_x0 = -a_tmp[0]*cos(speed_angle)
        dat_y0 = a_tmp[0]*sin(speed_angle)
        at = np.append(at, (dat_x0**2 + dat_y0**2)**(1/2))
        
        dac_x0 = -v_tmp[0]**2/radius*sin(speed_angle)
        dac_y0 = -v_tmp[0]**2/radius*cos(speed_angle)
        ac = np.append(ac, (dac_x0**2 + dac_y0**2)**(1/2))
        
        a = np.append(a, (at[-1]**2 + ac[-1]**2)**(1/2))

        d_theta_prev = 0

        # initialization of x and y components
        self.arc_sx = np.array([radius*cos(start_angle)])
        self.arc_sy = np.array([radius*sin(start_angle)])
        self.arc_vx = np.array([dv_x0])
        self.arc_vy = np.array([dv_y0])
        self.arc_atx = np.array([dat_x0])
        self.arc_aty = np.array([dat_y0])
        self.arc_acx = np.array([dac_x0])
        self.arc_acy = np.array([dac_y0])
        self.arc_ax = np.array([(dat_x0**2 + dac_x0**2)**(1/2)])
        self.arc_ay = np.array([(dat_y0**2 + dac_y0**2)**(1/2)])

        for si, sf, vi, ai in zip(s_tmp[:-1], s_tmp[1:], v_tmp[1:], a_tmp[1:]):

            # angle fraction
            d_theta = (sf-si)/radius

            # angle associated to the chord
            pos_angle = d_theta_prev + d_theta/2 
        
            # chord length
            ds = 2*radius*sin(d_theta/2)

            # horizontal and vertical displacement increment
            dx = ds*sin(start_angle + pos_angle)
            dy = ds*cos(start_angle + pos_angle)

            # calculation of the horizontal and vertical speed
            speed_angle = pi/2 - (d_theta + d_theta_prev)

            # horizontal and vertical speed increment
            dv_x = -vi*cos(speed_angle)
            dv_y = vi*sin(speed_angle)
            
            # horizontal and vertical tangential acceleration increment
            dat_x = -ai*cos(speed_angle)
            dat_y = ai*sin(speed_angle)

            # horizontal and vertical centripetal acceleration increment
            dac_x = -vi**2/radius*sin(speed_angle)
            dac_y = -vi**2/radius*cos(speed_angle)
            
            # update vectors and angle increment
            d_theta_prev = d_theta + d_theta_prev

            s = np.append(s, s[-1]+(dx**2 + dy**2)**(1/2))
            v = np.append(v, (dv_x**2 + dv_y**2)**(1/2))
            at = np.append(at, np.sign(ai)*(dat_x**2 + dat_y**2)**(1/2))
            ac = np.append(ac, (dac_x**2 + dac_y**2)**(1/2))

            a = np.append(a, (at[-1]**2 + ac[-1]**2)**(1/2)) #  total acceleration


            # update x any kinematic components
            self.arc_sx = np.append(self.arc_sx, self.arc_sx[-1] - dx)
            self.arc_sy = np.append(self.arc_sy, self.arc_sy[-1] + dy)
            self.arc_vx = np.append(self.arc_vx, dv_x)
            self.arc_vy = np.append(self.arc_vy, dv_y)
            self.arc_atx = np.append(self.arc_atx, dat_x)
            self.arc_aty = np.append(self.arc_aty, dat_y)
            self.arc_acx = np.append(self.arc_acx, dac_x)
            self.arc_acy = np.append(self.arc_acy, dac_y)
            self.arc_ax = np.append(self.arc_ax, (dat_x**2 + dac_x**2)**(1/2))
            self.arc_ay = np.append(self.arc_ay, (dat_y**2 + dac_y**2)**(1/2))
                      
        return t, s, v, a, at, ac

#------------------------------------------------------------------------------------------------------------------------

    def get_arc_position(self) -> tuple:
        """
        Get position in x and y
        :return: a tuple of 2 arrays (position x and y) 
        """

        try:
            self.arc_sx
            self.arc_sy
        except AttributeError:
            print("Vector doesn't exist!")
            self.arc_sx = np.array([])
            self.arc_sy = np.array([])
        finally:
            sx = self.arc_sx
            sy = self.arc_sy
       
        return sx, sy
    
    #=========================================================================================================================

    def get_arc_speed(self) -> tuple:
        """
        Get speed in x and y
        :return: a tuple of 2 arrays (speed x and y) 
        """

        try:
            self.arc_vx
            self.arc_vy
        except AttributeError:
            print("Vector doesn't exist!")
            self.arc_vx = np.array([])
            self.arc_vy = np.array([])
        finally:
            vx = self.arc_vx
            vy = self.arc_vy
       
        return vx, vy
    
#------------------------------------------------------------------------------------------------------------------------

    def get_arc_tan_acceleration(self) -> tuple:
        """
        Get tangential acceleration in x and y
        :return: a tuple of 2 arrays (tangential acceleration x and y) 
        """

        try:
            self.arc_atx
            self.arc_aty
        except AttributeError:
            print("Vector doesn't exist!")
            self.arc_atx = np.array([])
            self.arc_aty = np.array([])
        finally:
            atx = self.arc_atx
            aty = self.arc_aty
       
        return atx, aty
    
#------------------------------------------------------------------------------------------------------------------------

    def get_arc_centr_acceleration(self) -> tuple:
        """
        Get centripetal acceleration in x and y
        :return: a tuple of 2 arrays (centripetal acceleration x and y) 
        """

        try:
            self.arc_acx
            self.arc_acy
        except AttributeError:
            print("Vector doesn't exist!")
            self.arc_acx = np.array([])
            self.arc_acy = np.array([])
        finally:
            acx = self.arc_acx
            acy = self.arc_acy
       
        return acx, acy
    
#------------------------------------------------------------------------------------------------------------------------

    def get_arc_total_acceleration(self) -> tuple:
        """
        Get total acceleration in x and y
        :return: a tuple of 2 arrays (total acceleration x and y) 
        """

        try:
            self.arc_ax
            self.arc_ay
        except AttributeError:
            print("Vector doesn't exist!")
            self.arc_ax = np.array([])
            self.arc_ay = np.array([])
        finally:
            ax = self.arc_ax
            ay = self.arc_ay
       
        return ax, ay
    
#=========================================================================================================================
    # ******************** creare procedura per ottenere la legge di moto da un percorso generico ***********
    def polygon_path_motion(self, kin_param, fillet_radius=1, n_points=1000, shapes=[[0.2, 0.6, 0.2]]) -> tuple:
        """
        Calculate a law of motion for a polygon that has straight and curved parts alternated 
        with same amount of points.
        
        :param kin_param: a tuple of 4 array (time, space, speed, acceleration) calculateted
                          from an equivalent straight path with same length
               fillet_radius: radius of the arc
               
        :return: 
        """
        pass
        """for i in shapes:

            

            self.trapezoidal_LOM(self, Ds=, s0=, Dt=, ti=, n_points=n_points, shape=i)

            self.arc_motion(self, kin_param, radius=1, start_angle=0)"""