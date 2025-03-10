"""
Autore: Archetti Ivan
Data: 10/01/2025

Messa in opera delle leggi di moto e composizione di moti
"""


from MotionProfile import MotionProfile
from PlotMotionProfile import PlotMotionProfile
from TracePath import TracePath
import numpy as np
from numpy import pi

# =======================================================================================================================

def traccia_linea_generica():
    # definizione punti
    punto_partenza = (1,1)
    punto_arrivo = (2,3)

    n_punti = 1000

    um1 = ("Tempo", "Spazio", "Velocità", "Accelerazione")    
    LDM = MotionProfile()
    PlotLDM = PlotMotionProfile(um_axes=um1)
    TracciaPercorso = TracePath()

    # discretizzo il percorso
    pos_x, pos_y = TracciaPercorso.trace_line(punto_partenza, punto_arrivo, n_points=n_punti)
    
    # mostro il percorso
    um = ("[m]", "[m]")
    PlotLDM.plot_path(x=pos_x, y=pos_y, um=um)

    # calcolo la legge di moto e mostro i grafici
    Ds_x = punto_arrivo[0] - punto_partenza[0]
    Ds_y = punto_arrivo[1] - punto_partenza[1]
    Ds = (Ds_x**2 + Ds_y**2)**(1/2)
    Dt = 1
    t_i = 0
    s0 = 0
    shape = [0.2, 0.6, 0.2]

    # spostamento in modulo
    t, s, v, a  = LDM.trapezoidal_MP(Ds=Ds, Dt=Dt, ti=t_i, s0=s0, n_points=n_punti, shape=shape)

    um_x = ("[s]", "[m]", "$[\\frac{m}{s}]$", "$[\\frac{m}{s^{2}}]$")
    titoli_x = ["          Moto", "Spazio", "Velocità", "Accelerazione"] 
    PlotLDM.plot_motion_profile(title=titoli_x, t=t, s=s, v=v, a=a, amax=(10,-10), vmax=(2,1), um=um_x)

    # spostamento, velocità ed accelerazione in x e y
    angle = LDM.get_line_angle(P0=punto_partenza, P1=punto_arrivo)
    sx, sy = LDM.get_line_position(s=s, angle=angle)
    vx, vy = LDM.get_line_speed(v=v, angle=angle)
    ax, ay = LDM.get_line_acceleration(a=a, angle=angle)

    um_x = ("[s]", "[m]", "$[\\frac{m}{s}]$", "$[\\frac{m}{s^{2}}]$")
    titoli_x = ["          Moto in x", "Spazio in x", "Velocità in x", "Accelerazione in x"] 
    PlotLDM.plot_motion_profile(title=titoli_x, t=t, s=sx, v=vx, a=ax, amax=(0,0), vmax=(0,0), um=um_x)
    
    um_y = ("[s]", "[m]", "$[\\frac{m}{s}]$", "$[\\frac{m}{s^{2}}]$")
    titoli_y = ["          Moto in y", "Spazio in y", "Velocità in y", "Accelerazione in y"] 
    PlotLDM.plot_motion_profile(title=titoli_y, t=t, s=sy, v=vy, a=ay, amax=(0,0), vmax=(0,0), um=um_y)

    


# =======================================================================================================================

def traccia_arco():
    # parametri della curva
    R = 1
    centro = (0, 0)
    angolo_partenza = 0
    angolo = pi/2
    
    n_punti = 1000

    um1 = ("Tempo", "Spazio", "Velocità", "Accelerazione")    
    LDM = MotionProfile()
    PlotLDM = PlotMotionProfile(um_axes=um1)
    TracciaPercorso = TracePath()

    # discretizzo il percorso
    pos_x, pos_y = TracciaPercorso.trace_arc(radius=R, center=centro, start_angle=angolo_partenza, angle=angolo, n_points=n_punti)
    
    # mostro il percorso
    um = ("[m]", "[m]")
    #PlotLDM.plot_path(x=pos_x, y=pos_y, um=um)

    # calcolo la legge di moto
    Dt = 1
    ti = 0
    shape = [0.2, 0.6, 0.2]
    s0 = angolo_partenza*R
    Ds = angolo*R
    par_cinematici = list(LDM.trapezoidal_MP(Ds=Ds, Dt=Dt, ti=ti, s0=s0, n_points=n_punti, shape=shape))
    t, s, v, a, at, ac = LDM.arc_motion(par_cinematici, radius=R, start_angle=angolo_partenza)

    um1 = ["[s]", "[m]", "$[\\frac{m}{s}]$", "$[\\frac{m}{s^{2}}]$"]
    titoli = ["          Moto curvilineo", "Spazio", "Velocità", "Accelerazione totale"]
    PlotLDM.plot_motion_profile(title=titoli, t=t, s=s, v=v, a=a, um=um1)

    # spostamento, velocità ed accelerazione in x e y
    sx, sy = LDM.get_arc_position()
    vx, vy = LDM.get_arc_speed()
    atx, aty = LDM.get_arc_tan_acceleration()
    acx, acy = LDM.get_arc_centr_acceleration()
    ax, ay = LDM.get_arc_total_acceleration()

    """PlotMotionProfile.plot_kinematic_value(title="Accelerazione centripeta", t=t, x=ac, um=["s","$[\\frac{m}{s^{2}}]$"], 
                                 label=["tempo", "Accelerazione centripeta"], color="RED")
    
    PlotMotionProfile.plot_kinematic_value(title="Accelerazione tangenziale", t=t, x=at, um=["s","$[\\frac{m}{s^{2}}]$"], 
                                 label=["tempo", "Accelerazione tangenziale"], color="RED")"""

    """um_y = ("[s]", "[m]", "$[\\frac{m}{s}]$", "$[\\frac{m}{s^{2}}]$")
    titoli_y = ["          Moto in y", "Spazio in y", "Velocità in y", "Accelerazione in y"] 
    PlotMotionProfile.plot_motion_profile(title=titoli_y, t=t, s=sy, v=vy, a=ay, amax=(0,0), vmax=(0,0), um=um_y)"""

    """PlotMotionProfile.plot_kinematic_value(title="Accelerazione totale", t=t, x=a, um=["s","$[\\frac{m}{s^{2}}]$"], 
                                 label=["tempo", "Accelerazione totale"], color="RED")"""
    
    um_y = ("[s]", "[m]", "$[\\frac{m}{s}]$", "$[\\frac{m}{s^{2}}]$")
    titoli_y = ["          Moto in y", "Spazio in y", "Velocità in y", "Accelerazione totale in y"] 
    PlotMotionProfile.plot_motion_profile(title=titoli_y, t=t, s=sy, v=vy, a=ay, amax=(0,0), vmax=(0,0), um=um_y)

    um_x = ("[s]", "[m]", "$[\\frac{m}{s}]$", "$[\\frac{m}{s^{2}}]$")
    titoli_y = ["          Moto in x", "Spazio in x", "Velocità in x", "Accelerazione totale in x"] 
    PlotMotionProfile.plot_motion_profile(title=titoli_y, t=t, s=sx, v=vx, a=ax, amax=(0,0), vmax=(0,0), um=um_x)


# =======================================================================================================================

def traccia_percorso():
    # prove legge di moto
    um1 = ("Tempo", "Spazio", "Velocità", "Accelerazione")    
    LDM = MotionProfile()
    PlotLDM = PlotMotionProfile(um_axes=um1)
    TracciaPercorso = TracePath()

    n_punti = 1000

    # geometrie del percorso
    h = 1
    w = 2
    R = 0.4

    pos_x, pos_y = TracciaPercorso.trace_rounded_rectangle(center=(0,0), rectangle_angle=0, fillet_radius=R, heigth=h, width=w, n_points=n_punti)

    # mostro il percorso
    um2 = ("[m]", "[m]")
    PlotMotionProfile.plot_path(x=pos_x, y=pos_y, um=um2)

    # -------------------------------------------------------------------------------------------------------
    # calcolo la legge di moto

    # inizializzo le grandezze cinematiche
    t = np.array([])
    s = np.array([])
    v = np.array([])
    a = np.array([])

    # tracciatura circuito come ciclo for
    angolo_curve = pi/2
    segmenti = (w, angolo_curve*R, h, angolo_curve*R, w, angolo_curve*R, h, angolo_curve*R)
    forma_acc = [0.2, 0.8, 0] 
    forma_cost = [0, 1, 0]     
    forma_dec = [0, 0.8, 0.2] 
    v_max = 1
    

    # inizializzazione
    t0 = 0
    s0 = 0
    
    for i, Ds in enumerate(segmenti):
        if i == 0:
            forma = forma_acc
        elif 0 < i < len(segmenti) - 1:
            forma = forma_cost
        elif i == len(segmenti) - 1:
            forma = forma_dec
        
        Dt = Ds/(v_max*(forma[0]/2 + forma[1] + forma[2]/2))

        ti, si, vi, ai = LDM.trapezoidal_MP(Ds=Ds, Dt=Dt, ti=t0, s0=s0, n_points=n_punti, shape=forma)
        
        if i % 2 == 1:
            par_cinematici = ti, si, vi, ai
            ti, si, vi, ai, ait, aic = LDM.arc_motion(par_cinematici, radius=R, start_angle=0)

        s0 = si[-1]
        t0 = ti[-1]

        t = np.concatenate([t, ti])
        s = np.concatenate([s, si])
        v = np.concatenate([v, vi])
        a = np.concatenate([a, ai])

    # traccio la legge di moto complessiva
    um1 = ["[s]", "[m]", "$[\\frac{m}{s}]$", "$[\\frac{m}{s^{2}}]$"]
    titoli = ["          Moto rettilineo", "Spazio percorso", "Velocità percorso", "Accelerazione percorso"]
    PlotMotionProfile.plot_motion_profile(title=titoli, t=t, s=s, v=v, a=a, um=um1)

# =======================================================================================================================

def traccia_percorso_utensile():
    # prove legge di moto
    um1 = ("Tempo", "Spazio", "Velocità", "Accelerazione")    
    LDM = MotionProfile()
    PlotLDM = PlotMotionProfile(um_axes=um1)
    TracciaPercorso = TracePath()

    n_punti = 1000

    # geometrie del percorso
    w = 10
    h = 20

    pos_x, pos_y = TracciaPercorso.trace_milling_path(P_start=(0,0), n_step=10, height=h, width=w, n_points=n_punti)

    # mostro il percorso
    um2 = ("[m]", "[m]")
    PlotLDM.plot_path(x=pos_x, y=pos_y, um=um2)

# =======================================================================================================================

def traccia_poligono():
    um1 = ("Tempo", "Spazio", "Velocità", "Accelerazione")
    LDM = MotionProfile()
    PlotLDM = PlotMotionProfile(um_axes=um1)
    TracciaPercorso = TracePath()

    n_punti = 1000

    pos_x, pos_y = TracciaPercorso.trace_regular_rounded_polygon(center=(0,0), n_sides=6, fillet_radius=0.4, radius=1, polygon_angle=0, n_points=n_punti)
    um2 = ("[m]", "[m]")
    PlotLDM.plot_path(x=pos_x, y=pos_y, um=um2)

    

# =======================================================================================================================



if __name__ == "__main__":
    scelta = 0
    
    if scelta == 0:
        traccia_linea_generica()
    elif scelta == 1:
        traccia_arco()
    elif scelta == 2:
        traccia_percorso()
    elif scelta == 3:
        traccia_percorso_utensile()
    elif scelta == 4:
        traccia_poligono()

