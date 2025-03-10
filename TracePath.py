"""
Author: Archetti Ivan
Date: 08/02/2025

Class to draw path of a body on x, y plan
"""
#=========================================================================================================================

import numpy as np
from numpy import sin, cos, pi


#=========================================================================================================================

class TracePath():
    def __init__(self):
        pass
        
    
        
#=========================================================================================================================

    def trace_line(self, P_start=(0, 0), P_end=(1, 1), n_points=100) -> tuple:
        """
        From P_start to P_end create two vectors of n_points for x and y direction

        :param P_start: point of start
               P_end: point of end
               n_points: number of points for calculation
        :return 2 arrays of n_points
        """
        # controll of number of points
        if n_points <= 0:
            print("n_points must at least 1!")
            n_points = 2

        x = np.linspace(P_start[0], P_end[0], n_points) 
        y = np.linspace(P_start[1], P_end[1], n_points)
        
        return x, y
    
#=========================================================================================================================

    def trace_arc(self, center=(0, 0), radius=1, start_angle=0, angle=pi/2, n_points=100) -> tuple:
        """
        From start_angle till angle create a vector of n_points that describe an arc in x y directions

        :param center: a tuple of the arc center
               radius: radius of the arc
               start_angle: start angle of the arc
               angle: magnitude of the arc angle
               n_points: number of points for calculation
        :return 2 arrays of n_points
        """
        # initialisation
        R = radius

        # controll of number of points
        if n_points <= 0:
            print("n_points must at least 1!")
            n_points = 3

        d_theta = angle/n_points
        ds = 2*R*sin(d_theta/2)

        x = np.array([R*cos(start_angle) + center[0]])
        y = np.array([R*sin(start_angle) + center[1]])
     
        for i in range(0, n_points-1):
            pos_angle = start_angle + i*d_theta + d_theta/2
            
            dx = ds*sin(pos_angle)
            dy = ds*cos(pos_angle)

            x = np.append(x, x[-1] - dx)
            y = np.append(y, y[-1] + dy)

        return x, y

#=========================================================================================================================

    def trace_rounded_rectangle(self, center=(0, 0), rectangle_angle=0, fillet_radius=1, heigth=2, width=1, n_points=100) -> tuple:
        """
        Create a rectangle with rounded corners

        :param center: a tuple of the rectangle center
               rectangle_angle: inclination agle of the rectangle
               radius: radius of the rounded corners
               length: length of the rectangle
               width: width of the rectangle
               n_points: number of points for calculation
        :return 2 arrays of n_points
        """
        # initialisation
        R = fillet_radius
        w = width
        h = heigth
        r_angle = rectangle_angle
        x = np.array([])
        y = np.array([])
        
        # controll of the radius
        if R < 0 or R > w/2 or R > h/2:
            print("Radius must be greater than 0 but smaller than half legth or half width.")
            R = 0

        wp = w - 2*R                        # width without fillet
        hp = h - 2*R                        # height without fillet
        s = ((w/2-R)**2+(h/2-R)**2)**(1/2)  # distance between rectangle center and fillet center
        beta = np.arctan((h/2-R)/(w/2-R))   # angle between line of rectangle center/fillet center and a horizontal line

        dx = h/2*sin(r_angle)               # horizontal displacement of the side due to rectangle angle 
        dy = h/2*(1-cos(r_angle))           # vertical displacement of the side due to rectangle angle

        P_start = (-wp/2*cos(r_angle) + dx + center[0], 
                   -h/2-wp/2*sin(r_angle) + dy + center[1])

        for i in range(4):
            # check if is a horizontal or vertical side of the rectangle
            if i % 2 == 0:
                P_end = (P_start[0] + wp*cos(r_angle + i*pi/2), P_start[1] + wp*sin(r_angle + i*pi/2))
            else:
                P_end = (P_start[0] + hp*cos(r_angle + i*pi/2), P_start[1] + hp*sin(r_angle + i*pi/2))

            # check if is right side or left side of rectangle
            if i < 2:
                s_angle = 0
            else:
                s_angle = pi
            
            pos_x, pos_y = self.trace_line(P_start=P_start, P_end=P_end, n_points=n_points) 

            x = np.append(x, pos_x)
            y = np.append(y, pos_y)

            r_center = (s*cos(r_angle + (-1)**(i+1)*beta + s_angle) + center[0], 
                        s*sin(r_angle + (-1)**(i+1)*beta + s_angle) + center[1])
            pos_x, pos_y = self.trace_arc(radius=R, center=r_center, start_angle=-pi/2 + r_angle + i*pi/2, angle=pi/2, n_points=n_points)
 
            x = np.append(x, pos_x)
            y = np.append(y, pos_y)

            P_start = (x[-1], y[-1])
        
        return x, y
    
#=========================================================================================================================

    def trace_milling_path(self, P_start=(0, 0), n_step=1, height=2, width=1, n_points=100) -> tuple:
        """
        Create a tool-like homogeneous path on a surface

        :param start: a tuple of the starting point
               radius: radius of the rounded corners
               length: length of the rectangle
               width: width of the rectangle
               n_points: number of points for calculation
        :return 2 arrays of n_points
        """
        # initialisation
        x = np.array([])
        y = np.array([])

        w = width
        h = height

        P_end = (0,0)

        # controll of n_step value
        if n_step<1:
            print("Step must be greater 1 or at least 1!")
            n_step = 1
        if not isinstance(n_step, int):
            print("n_step must be an integer!")
            n_step = int(np.round(n_step, 0))


        step = h/n_step           # calculation of the step lenght
        R = step/2                # calculation of the radius fillet

        c = np.sign(w)*np.sign(h) # verify if the path running clockwise or counterclockwise
        
        for i in range(n_step+1):
            # controll the parameters of the first path
            if i == 0:
                # first line
                k = 0
                j = 1
            elif i > 0 and i < n_step:
                # between the firsr and last line
                k = 2
                j = 2
            else:
                # last line
                j = 1

            # line
            P_start = (P_end[0], P_end[1] + k*R)
            P_end = (P_start[0] + (-1)**(i)*(w-c*j*R), P_start[1])
            pos_x, pos_y = self.trace_line(P_start=P_start, P_end=P_end, n_points=n_points) 
            x = np.append(x, pos_x)
            y = np.append(y, pos_y)
            
            if i < n_step:
                # arc
                center = (P_end[0], P_end[1] + R)
                angle = (-1)**(i)*pi*c
                pos_x, pos_y = self.trace_arc(center=center, radius=R, start_angle=-pi/2, angle=angle, n_points=n_points)
                x = np.append(x, pos_x)
                y = np.append(y, pos_y)

        return x, y
    
    #=========================================================================================================================

    def trace_regular_rounded_polygon(self, center=(0, 0), polygon_angle=0, radius=1, n_sides=3, fillet_radius=0.4, n_points=100) -> tuple:
        """
        Create a rectangle with rounded corners

        :param center: a tuple of the rectangle center
               polygon_angle: inclination angle of the polygon
               radius: radius of the circle that circumscribes the polygon
               n_sides: number of polygon sides
               fillet_radius: radius of the rounded corner
               n_points: number of points for calculation
        :return 2 arrays of n_points
        """
        # initialisation
        R = radius
        N = n_sides
        r = fillet_radius
        p_angle = polygon_angle
        x = np.array([])
        y = np.array([])

        # controll number of sides
        if N < 3:
            print("Number of sides must be at least 3!")
            N = 3

        # controll of the fillet radius
        if r < 0 or r > R:
            print("Fillet radius must be greater than 0 and smaller than radius.")
            r = 0

        # calculation of the geometric parameters of the polygon
        theta = (N-2)/N*pi                    # internal angle of polygon
        beta = 2*pi/N                         # angle of fillet_radius
        alpha = 0                             # slope of each side respect to the horizontal line

        h = R*np.sqrt(1-(sin(pi/N))**2)       # distance from circle's center and poligon center side
        self.l_polygon = 2*np.sqrt(R**2-h**2) # lenght of poligon side
        hp = h - r                            # vertical ddistance from the polygon center to the fillet_radius center
        Lp = hp*np.tan(beta/2)                # half side of the polygon
        s = hp/cos(beta/2)                    # distance from the polygon center to the fillet_radius center

        dx = h*sin(p_angle)                   # horizontal displacement of the side due to polygon angle
        dy = h*(1-cos(p_angle))               # verical displacement of the side due to polygon angle

        P_start = (-Lp*cos(alpha+p_angle) + center[0] + dx, 
                   -h-Lp*sin(alpha+p_angle) + center[1] + dy)

        for i in range(N):
            P_end = (P_start[0] + 2*Lp*cos(alpha+p_angle), P_start[1] + 2*Lp*sin(alpha+p_angle))
            pos_x, pos_y = self.trace_line(P_start=P_start, P_end=P_end, n_points=n_points)
            alpha = alpha + (pi - theta)
            
            x = np.append(x, pos_x)
            y = np.append(y, pos_y)

            fillet_center = (s*sin(beta/2+i*beta + p_angle) + center[0] , -s*cos(beta/2+i*beta + p_angle) + center[1])
            start_angle = -pi/2 + p_angle + i*beta
            pos_x, pos_y = self.trace_arc(center=fillet_center, radius=r, start_angle=start_angle, angle=beta, n_points=n_points)
            
            x = np.append(x, pos_x)
            y = np.append(y, pos_y)

            P_start = (x[-1], y[-1])
        return x, y

    #------------------------------------------------------------------------------------------------------------------------

    def get_polygon_lenght_side(self) -> float:
        """
        Get lenght of polygon side without fillet

        :return one float of the length
        """
        try:
            self.l_polygon
        except AttributeError:
            print("Polygon not present.")
            self.l_polygon = 0
        finally:
            l_polygon = self.l_polygon
       
        return l_polygon
    
    #=========================================================================================================================