"""
Author: Archetti Ivan
Date: 25/01/2025

Class to plot law of motion
"""
#=========================================================================================================================

import numpy as np
from matplotlib import pyplot as plt


#=========================================================================================================================

class PlotMotionProfile():
    def __init__(self, um_axes=("time", "space", "vel", "acc")):
        self.time = um_axes[0]
        self.space = um_axes[1]
        self.vel = um_axes[2]
        self.acc = um_axes[3]

    #=========================================================================================================================

    def plot_kinematic_value(self, title="", t=np.array([]), x=np.array([]), um=[], color="black", label=[]) -> None:
        """
        Show a single array of a kinematic value
        
        :param title: title of each graph and plot page
               t: array of time
               x: array of a "generic quantity"
               um: unit of measurement 
        """

        plt.plot(t, x, label=label, color=color)  
        plt.grid(True)
        plt.xlabel(f"{label[0]} {um[0]}", {'size': 12})  
        plt.ylabel(f"{label[1]} {um[1]}", {'size': 12})  
        plt.title(title, {'size': 16, 'color':'darkred'})  

        plt.show()

#=========================================================================================================================

    def plot_motion_profile(self, title=["", "", "", ""], t=np.array([]), s=np.array([]), v=np.array([]), a=np.array([]), amax=(0,0), vmax=(0,0), um=[]) -> None:
        """
        Show space, speed and acceleration of a law of motion 
        
        :param title: title of each graph and plot page
               t: array of time
               s: array of space
               v: array of speed
               a: array of acceleration
               amax: tuple that defines the lines of minimum a maximum of acceleration
               vmax: tuple that defines the lines of minimum a maximum of speed
               um: unit of measurement of each magnitude
        """
        
        grafic_area = plt.figure(figsize=(12, 8))
        grafic_area.suptitle(title[0], fontsize=24)
        
        axes_space = grafic_area.add_subplot(3, 1, 1)
        axes_space.plot(t, s, color="black", label=self.space)
        axes_space.grid(True)
        axes_space.set_title(title[1], {'size': 16, 'color':'darkred'})
        axes_space.set_xlabel(f"{self.time} {um[0]}", {'size': 12})
        axes_space.set_ylabel(f"{self.space} {um[1]}", {'size': 12})      

        # Add a lateral summary
        pos_s_max = np.argmax(s)
        t_max = t[pos_s_max]
        s_max = s[pos_s_max]
        pos_s_min = np.argmin(s)
        t_min = t[pos_s_min]
        s_min = s[pos_s_min]
        text = f"Max: ({np.round(t_max, 2)}, {np.round(s_max, 2)})\nMin: ({np.round(t_min, 2)}, {np.round(s_min, 2)})\nStop: ({np.round(t[-1], 2)}, {np.round(s[-1], 2)})"
        axes_space.text(1.03, 0.8, text, transform=axes_space.transAxes, fontsize=12, verticalalignment='top', 
                        bbox=dict(facecolor='lightyellow', edgecolor='black', boxstyle='round,pad=0.5'))
        
        # add vertical lines
        axes_space.vlines(x=[t[0]], ymin=0, ymax=s[0], colors='k', linestyles='--')
        axes_space.vlines(x=[t[-1]], ymin=0, ymax=s[-1], colors='k', linestyles='--')
        axes_space.legend()
        
        axes_vel = grafic_area.add_subplot(3, 1, 2)
        axes_vel.plot(t, v, color="green", label=self.vel)
        axes_vel.grid(True)
        axes_vel.set_title(title[2], {'size': 16, 'color':'darkred'})
        axes_vel.set_xlabel(f"{self.time} {um[0]}", {'size': 12})
        axes_vel.set_ylabel(f"{self.vel} {um[2]}", {'size': 12})

        # Add a lateral summary
        pos_v_max = np.argmax(v)
        t_max = t[pos_v_max]
        v_max = v[pos_v_max]
        pos_v_min = np.argmin(v)
        t_min = t[pos_v_min]
        v_min = v[pos_v_min]
        text = f"Max: ({np.round(t_max, 2)}, {np.round(v_max, 2)})\nMin: ({np.round(t_min, 2)}, {np.round(v_min, 2)})\nStop: ({np.round(t[-1], 2)}, {np.round(v[-1], 2)})"
        axes_vel.text(1.03, 0.8, text, transform=axes_vel.transAxes, fontsize=12, verticalalignment='top', 
                        bbox=dict(facecolor='lightyellow', edgecolor='black', boxstyle='round,pad=0.5'))
        
        # add acceleration limits
        if vmax[0] != 0:
            axes_vel.hlines(y=[vmax[0]], xmin=t[0], xmax=t[-1], colors='#658437', linestyles='--', label="max set")
        if vmax[1] != 0:
            axes_vel.hlines(y=[vmax[1]], xmin=t[0], xmax=t[-1], colors='#a7c66c', linestyles='--', label="min set")
        
        # add vertical lines
        axes_vel.vlines(x=[t[0]], ymin=0, ymax=v[0], colors='k', linestyles='--')
        axes_vel.vlines(x=[t[-1]], ymin=0, ymax=v[-1], colors='k', linestyles='--')
        axes_vel.legend()

        axes_acc = grafic_area.add_subplot(3, 1, 3)
        axes_acc.plot(t, a, color="red", label=self.acc)
        axes_acc.grid(True)
        axes_acc.set_title(title[3], {'size': 16, 'color':'darkred'})
        axes_acc.set_xlabel(f"{self.time} {um[0]}", {'size': 12})
        axes_acc.set_ylabel(f"{self.acc} {um[3]}", {'size': 12})

        # add acceleration limits
        if amax[0] != 0:
            axes_acc.hlines(y=[amax[0]], xmin=t[0], xmax=t[-1], colors='orange', linestyles='--', label="max set")
        if amax[1] != 0:
            axes_acc.hlines(y=[amax[1]], xmin=t[0], xmax=t[-1], colors='brown', linestyles='--', label="min set")

        # Add a lateral summary
        pos_a_max = np.argmax(a)
        t_max = t[pos_a_max]
        a_max = a[pos_a_max]
        pos_a_min = np.argmin(a)
        t_min = t[pos_a_min]
        a_min = a[pos_a_min]
        text = f"Max: ({np.round(t_max, 2)}, {np.round(a_max, 2)})\nMin: ({np.round(t_min, 2)}, {np.round(a_min, 2)})\nStop: ({np.round(t[-1], 2)}, {np.round(a[-1], 2)})"
        axes_acc.text(1.03, 0.8, text, transform=axes_acc.transAxes, fontsize=12, verticalalignment='top', 
                        bbox=dict(facecolor='lightyellow', edgecolor='black', boxstyle='round,pad=0.5'))
        
        # add verticle lines
        axes_acc.vlines(x=[t[0]], ymin=0, ymax=a[0], colors='k', linestyles='--')
        axes_acc.vlines(x=[t[-1]], ymin=0, ymax=a[-1], colors='k', linestyles='--')
        axes_acc.legend()

        grafic_area.tight_layout()
        
        plt.show()


#=========================================================================================================================

    def plot_path(self, title="", x=np.array([0]), y=np.array([0]), um=[]) -> None:
        """
        Show the track of the point

        :param title: title of plot page
               x: array of x space
               y: array of y space
               um: unit of measurement of each magnitude

        """

        plt.plot(x, y, color="black", label=self.space)
        plt.grid(True)
        plt.title(title, {'size': 16, 'color':'darkred'})
        plt.xlabel(f"x {um[0]}", {'size': 12})
        plt.ylabel(f"y {um[1]}", {'size': 12})
        plt.legend()
      
        plt.show()

        
