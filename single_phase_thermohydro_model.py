from dolfin import *
# from mshr import *
from datetime import datetime
import numpy as np
import csv
import os
import sys

# Setting up dolfin parameters
if has_linear_algebra_backend("Epetra"):
    parameters["linear_algebra_backend"] = "Epetra"

# Setting up solver parameters
parameters["form_compiler"]["cpp_optimize"] = True
ffc_options = {"quadrature_degree": 6, "optimize": True}

'''
Verbose:
If true, multiple information will be printed in the terminal, regarding the
mesh, mesh quality, time evolution of quantites of interest, and so on.
'''
VERBOSE = True

'''
Heating Method:
HUC - Heating Up Curve - Drying curve of Gong et al., 1995.
ISO - IS O834:2014 Fire Curve - Fire simulation with the ISO 834 fire curve.
'''
# HEATING_CASE = 'ISO'
HEATING_CASE = 'HUC'


'''
Sorption Isotherm Method:
numerical - Numerical Derivative of non-smoothed sorption isotherm.
analytical - Analytical Derivative of non-smoothed sorption isotherm.
smooth - Analytical Derivative of smoothed sorption isotherm.
smooth_num - Numerical Derivative of smoothed sorption isotherm.
'''
SORPTION_METHOD = 'smooth'

'''
Batch Mode Options:
single - Single simulation mode, no extra arguments needed tfor running.
convergence - If set, the script needs to be called following:
'python single_phase_thermohydro_model.py dt nx SORPTION_METHOD'.
sensitivity - If set, the script needs to be called following:
'python single_phase_thermohydro_model.py lambda_c K_0 porosity w_c'.
For both convergence and sensitivity, it is recommended to run through the
batch_convergece.py or batch_sensitivity.py files.
'''
BATCH_MODE = 'single'
RADIATION = 0

'''
Domain:
1D case used on convergence and sensitivity analysis or 2D fibers case.
NOTE: For running the 2D fibers case, first run the 'generate_2D_meshes.py' file,
to create the necessary files.
'''
N_DIM = 2


if BATCH_MODE == 'convergence':
    dt = eval(sys.argv[1])
    nx = eval(sys.argv[2])
    SORPTION_METHOD = sys.argv[3]
    print(dt, nx, SORPTION_METHOD)
    set_log_level(50)
    VERBOSE = False
    base_dir = './convergence_analysis/'
    if not os.path.exists(base_dir):
        os.mkdir(base_dir)
    saving_dir = (base_dir + '/case_' + HEATING_CASE + '_' + SORPTION_METHOD
                  + '_1D_dt_' + str(dt) + '_nx_' + str(nx))
    if not os.path.exists(saving_dir):
        os.mkdir(saving_dir)

elif BATCH_MODE == 'sensitivity':
    lambda_c = eval(sys.argv[1])
    K_0 = eval(sys.argv[2])
    porosity = eval(sys.argv[3])
    wt_c = eval(sys.argv[4])
    print(lambda_c, K_0, porosity, wt_c)
    set_log_level(20)
    VERBOSE = False
    base_dir = './sensitivity_analysis/'
    if not os.path.exists(base_dir):
        os.mkdir(base_dir)
    saving_dir = (base_dir + '/case_' + HEATING_CASE + '_lambda_c_'
                  + str(lambda_c) + '_K_0_' + str(K_0) + '_porosity_'
                  + str(porosity) + '_w_c_' + str(wt_c))
    if not os.path.exists(saving_dir):
        os.mkdir(saving_dir)


elif BATCH_MODE == 'single':
    if N_DIM == 1:
        '''
        Define the proper name for the desired analysis by substituting the
        '_my_case' string with proper identifiers.
        '''
        base_dir = './1D_case/'
        if not os.path.exists(base_dir):
            os.mkdir(base_dir)
        saving_dir = (base_dir + 'case_norad_' + HEATING_CASE + '_my_case')
        if not os.path.exists(saving_dir):
            os.mkdir(saving_dir)

    elif N_DIM == 2:
        '''
        For the 2D case, both the size of the inclusions, INC_SIZE and the
        formulation strategy, FORMULATION, need to be defined.
          FORMULATION options:
          primal - primal formulation defined by Equations (23) and (24).
          mixed - mixed formulationn defined by Equations () and ().

          INC_SIZE options:
          short - shorter fibers, for details see 'generate_2D_meshes.py'.
          long - longer fibers, for details see 'generate_2D_meshes.py'.

        '''
        FORMULATION = 'mixed'
        INC_SIZE = 'short'

        base_dir = './2D_case/'
        if not os.path.exists(base_dir):
            os.mkdir(base_dir)
        saving_dir = (base_dir + 'case_' + HEATING_CASE + '_formulation_'
                      + FORMULATION + '_inc_size_' + INC_SIZE)
        if not os.path.exists(saving_dir):
            os.mkdir(saving_dir)

else:
    print('Invalid BATCH_MODE. Choose properly, options are: single, \
           convergence and sensitivity')

# Set the log level of dolfin, default is 20, for no information, use 50.
set_log_level(20)
if VERBOSE:
    print('    |-Running with the following parameters set:')
    print('    |--BATCH_MODE=' + BATCH_MODE)
    print('    |--HEATING_CASE=' + HEATING_CASE)
    print('    |--SORPTION_METHOD=' + SORPTION_METHOD)
    print('    |--N_DIM=' + str(N_DIM))
    print('\n')
    set_log_level(50)  # otherwise dolfin info will populate the terminal.

# Materials parameters and constants
DeltaH_d = Constant(0.0)
g = Constant(9.81)
alpha_T = Constant(5.0)
beta = Constant(1e-6)
p_inf = Constant(2850.0)
T_inf = Constant(298.15)
if BATCH_MODE != 'sensitivity':
    lambda_c = Constant(1.67)
    K_0 = Constant(1e-12)
    w_c = Constant(300.0)
    w_0 = Constant(100)
rho_c = Constant(2200)
C_pc = Constant(1100)
C_pw = Constant(4100.0)

if BATCH_MODE == 'sensitivity':
    set_log_level(50)
    rho_l = 1000
    rho_s = 2200
    w_0 = porosity * rho_l
    wt_c = wt_c / 100
    w_c = rho_s * wt_c
    print(w_0, w_c)

# Defining the coefficients of time discretization
if BATCH_MODE != 'convergence':
    if HEATING_CASE == 'HUC':
        dt = 60
    elif HEATING_CASE == 'ISO':
        dt = 5
    else:
        print('Invalid HEATING_CASE. Choose properly, options: HUC and ISO')

theta = Constant(0.)
alpha = Constant(1.0)
t = 0
if HEATING_CASE == 'HUC':
    T_total = 30 * 3600

elif HEATING_CASE == 'ISO':
    T_total = 1 * 3600

# Mesh definition
if BATCH_MODE != 'convergence':
    if N_DIM == 1:
        nx = 200
        lx = 0.2
        mesh = IntervalMesh(nx, 0.0, lx)
        # Boundaries definition
        boundaries = MeshFunction('size_t', mesh, mesh.topology().dim() - 1)
        boundaries.set_all(0)

        class left(SubDomain):
            def inside(self, x, on_boundary):
                return abs(x[0]) < DOLFIN_EPS and on_boundary

        class right(SubDomain):
            def inside(self, x, on_boundary):
                return abs(x[0] - lx) < DOLFIN_EPS and on_boundary

        left = left()
        right = right()
        left.mark(boundaries, 1)
        right.mark(boundaries, 2)
        ds = Measure("ds", domain=mesh, subdomain_data=boundaries)
        dx = Measure("dx")

        if VERBOSE:
            print("    |-Mesh done:")
            print("    |--Number of vertices = " + str(mesh.num_vertices()))
            print("    |--Number of cells = " + str(mesh.num_cells()))
            print("    |--Cell size hmax,hmin = %.3g %.3g" % (mesh.hmax(),
                                                              mesh.hmin()))
    if N_DIM == 2:
        lx = ly = 0.2
        if INC_SIZE == 'long':
            mesh = Mesh('./meshes/long_test.xml')
            sub_dom_xml = MeshFunction('int', mesh,
                                       './meshes/long_test_gmsh-physical.xml')
            sub_dom = MeshFunction("size_t", mesh, mesh.topology().dim())
            sub_dom.array()[:] = sub_dom_xml.array()[:]

        elif INC_SIZE == 'short':
            mesh = Mesh('./meshes/short_test.xml')
            sub_dom_xml = MeshFunction('int', mesh,
                                       './meshes/short_test_gmsh-physical.xml')
            sub_dom = MeshFunction("size_t", mesh, mesh.topology().dim())
            sub_dom.array()[:] = sub_dom_xml.array()[:]

        if VERBOSE:
            rat = MeshQuality.radius_ratio_min_max(mesh)
            print("    |-Mesh done")
            print("    |--Number of vertices = " + str(mesh.num_vertices()))
            print("    |--Number of cells = " + str(mesh.num_cells()))
            print("    |--Cell size hmax,hmin = %.3g %.3g" % (mesh.hmax(),
                                                              mesh.hmin()))
            print("    |--Radius ratio RRmax,RRmin = %.3g %.3g" % (max(rat),
                                                                   min(rat)))
        # Boundaries definition
        boundaries = MeshFunction('size_t', mesh, mesh.topology().dim() - 1)
        boundaries.set_all(0)

        class left(SubDomain):
            def inside(self, x, on_boundary):
                return abs(x[0]) < 1e-5 and on_boundary

        class right(SubDomain):
            def inside(self, x, on_boundary):
                return abs(x[0] - lx) < 1e-5 and on_boundary

        class top(SubDomain):
            def inside(self, x, on_boundary):
                return abs(x[1] - ly) < 1e-5 and on_boundary

        class down(SubDomain):
            def inside(self, x, on_boundary):
                return abs(x[1]) < 1e-5 and on_boundary

        left = left()
        right = right()
        top = top()
        down = down()
        left.mark(boundaries, 1)
        right.mark(boundaries, 2)
        top.mark(boundaries, 3)
        down.mark(boundaries, 4)
        ds = Measure("ds", domain=mesh, subdomain_data=boundaries)
        dx = Measure("dx")(subdomain_data=sub_dom)

        File(saving_dir + '/inc_size_' + INC_SIZE + '_form'
             + FORMULATION + '_mesh.pvd') << mesh
        File(saving_dir + '/inc_size_' + INC_SIZE + '_form'
             + FORMULATION + '_sub_dom.pvd') << sub_dom
        File(saving_dir + '/inc_size_' + INC_SIZE + '_form'
             + FORMULATION + '_bcs.pvd') << boundaries

elif BATCH_MODE == 'convergence':
    set_log_level(50)
    lx = 0.2
    mesh = IntervalMesh(nx, 0.0, lx)
    # Boundaries definition
    boundaries = MeshFunction('size_t', mesh, mesh.topology().dim() - 1)
    boundaries.set_all(0)

    class left(SubDomain):
        def inside(self, x, on_boundary):
            return abs(x[0]) < DOLFIN_EPS and on_boundary

    class right(SubDomain):
        def inside(self, x, on_boundary):
            return abs(x[0] - lx) < DOLFIN_EPS and on_boundary

    left = left()
    right = right()
    left.mark(boundaries, 1)
    right.mark(boundaries, 2)
    ds = Measure("ds", domain=mesh, subdomain_data=boundaries)
    if VERBOSE:
        print("    |-Mesh done:")
        print("    |--Number of vertices = " + str(mesh.num_vertices()))
        print("    |--Number of cells = " + str(mesh.num_cells()))
        print("    |--Cell size hmax,hmin = %.3g %.3g" % (mesh.hmax(),
                                                          mesh.hmin()))

# Heating case and Dirichlet BC definitions
if HEATING_CASE == 'HUC':
    rate_ = 30
    T_HUC = Expression('(t-t_0)*rate/3600.0 + T_0', degree=1, t=t, t_0=0,
                       rate=rate_, T_0=298.15)

elif HEATING_CASE == 'ISO':
    T_HUC = Expression('298.15 + 345.0*log10(8.0*(t/60)+1.0)', degree=1, t=t)


# Materials Parameters
# Antoine pressure
def p_sat(T):
    x = T - 273.15
    A1 = 133.322365
    B1 = 8.07131
    C1 = 1730.63
    D1 = 233.426
    A2 = 133.322365
    B2 = 8.14019
    C2 = 1810.94
    D2 = 244.485
    return conditional(le(T, 373.15), A1 * 10.**(B1 - (C1 / (D1 + x))),
                       A2 * 10.**(B2 - (C2 / (D2 + x))))


# Relative humidity
def phi(p, T):
    return p / p_sat(T)


# Properties of concrete
# Dehydration water
def w_d(T):
    x = T - 273.15
    A1 = 18.494166017699435
    A2 = -0.5666796585938574
    A3 = 0.007323190020922426
    x0 = 267.84667550724026
    dx = 17.338882746686192
    value = A1 + (A2 - A1) / (1 + exp((x - x0) / dx)) + A3 * x
    return (value + abs(value)) / 2


# Sorption isotherm
def m(T):
    return 1.04 - ((T - 263.15)**2) / ((T - 263.15)**2
                    + 22.34 * (298.15 - 263.15)**2)


def w_1(p, T):
    return w_c * ((w_0 / w_c) * (phi(p, T)))**(1 / m(T))


def w_2(p, T):
    return w_c * (0.037 * (phi(p, T) - 1.04) +
                  0.3335 * (1.0 - ((T - 273.15)**2.0) / (3.6e5)))


if SORPTION_METHOD == 'analytical':
    # Analytical derivative of sorption isotherm
    def w(p, T):
        w_096 = w_1(0.96 * p_sat(T), T)
        w_104 = w_2(1.04 * p_sat(T), T)
        w_int = w_096 + ((w_104 - w_096) * (phi(p, T) - 0.96) / 0.08)
        return conditional(le(T, 647.15),
                           conditional(le(phi(p, T), 0.96), w_1(p, T),
                                       conditional(le(phi(p, T), 1.04), w_int,
                                                   w_2(p, T))), 0.0)

    def dP_sdT(p, T):
        x = T - 273.15
        A1 = 531279.337303094
        B1 = 8.07131
        C1 = 1730.63
        D1 = 233.426
        A2 = 555933.390207996
        B2 = 8.14019
        C2 = 1810.94
        D2 = 244.485
        return conditional(le(T, 373.15),
                           A1 * 10.0**(B1 - C1 / (x + D1)) / (x + D1)**2,
                           A2 * 10.0**(B2 - C2 / (x + D2)) / (x + D2)**2)

    def dmdT(p, T):
        return -(((-2 * T + 526.3) * (T - 263.15)**2 /
                ((T - 263.15)**2 + 22.34 * (-263.15 + 298.15)**2)**2) +
                (2 * T - 526.3) / ((T - 263.15)**2 + 22.34 * (-263.15 + 298.15)**2))


    def dphidp(p, T):
        return 1 / p_sat(T)


    def dphidT(p, T):
        return -p * dP_sdT(p, T) / p_sat(T)**2


    def dw_1dP(p, T):
        return w_c * (w_0 * phi(p, T) / w_c)**(1 / m(T)) * \
        dphidp(p, T) / (phi(p, T) * m(T))


    def dw_intdP(p, T):
        return 12.5 * (- w_c * (p * w_0 / (w_c * p_sat(T)))**(1 / m(T)) + \
                    w_c * (0.037 * p / p_sat(T) - 9.26388888888889e-7 * \
                    (T - 273.15)**2 + 0.29502)) * dphidp(p, T) + \
                    (0.037 * w_c / p_sat(T) - w_c * \
                     (p * w_0 / (w_c * p_sat(T)))**(1 / m(T)) / \
                     (p * m(T))) * (12.5 * phi(p, T) - 12.0) + \
                      w_c * (p * w_0 / (w_c * p_sat(T)))**(1 / m(T)) / (p * m(T))

    def dw_2dP(p, T):
        return 0.037 * w_c * dphidp(p, T)


    def dwdp(p, T):
        return conditional(le(T, 647.15),
                conditional(le(phi(p, T), 0.96), dw_1dP(p, T),
                            conditional(le(phi(p, T), 1.04),
                                          dw_intdP(p, T),
                                          dw_2dP(p, T))), 0.0)


    def dw_1dT(p, T):
        return w_c * (w_0 * phi(p, T) / w_c)**(1 / m(T)) * \
        (-ln(w_0 * phi(p, T) / w_c) * dmdT(p, T) / m(T)**2 +
         dphidT(p, T) / (phi(p, T) * m(T)))


    def dw_intdT(p, T):
        return w_c * (p * w_0 / (w_c * p_sat(T)))**(1 / m(T)) * \
           (-ln(p * w_0 / (w_c * p_sat(T))) * dmdT(p, T) / m(T)**2 -
            dP_sdT(p, T) / (p_sat(T) * m(T))) + 12.5 * \
            (-w_c * (p * w_0 / (w_c * p_sat(T)))**(1 / m(T)) + \
             w_c * (0.037 * p / p_sat(T) - 9.26388888888889e-7 * \
                    (T - 273.15)**2 + 0.29502)) * dphidT(p, T) + \
            (-w_c * (p * w_0 / (w_c * p_sat(T)))**(1 / m(T)) * \
             (-ln(p * w_0 / (w_c * p_sat(T))) * dmdT(p, T) / m(T)**2 \
              - dP_sdT(p, T) / (p_sat(T) * m(T))) + w_c * \
             (-0.037 * p * dP_sdT(p, T) / p_sat(T)**2 - \
              1.85277777777778e-6 * T + 0.00050608625)) * (12.5 * phi(p, T) - 12.0)



    def dw_2dT(p, T):
        return w_c * (-1.8527777778e-6 * T + 0.037 * dphidT(p, T) + 0.00050608625)


    def dwdT(p, T):
        return conditional(le(T, 647.15),
                conditional(le(phi(p, T), 0.96), dw_1dT(p, T),
                            conditional(le(phi(p, T), 1.04),
                                          dw_intdT(p, T),
                                          dw_2dT(p, T))), 0.0)

    # Central differentiation formula for sorption isotherm
    def dwdt(p, T, p_n, T_n, p_n_theta=0, T_n_theta=0):
        if p_n_theta == 0:
            p_n_theta = p_n
            T_n_theta = T_n
        dPdt = (p - p_n) / dt
        dTdt = (T - T_n) / dt
        return dwdp(p_n_theta, T_n_theta) * dPdt + dwdT(p_n_theta, T_n_theta) * dTdt


elif SORPTION_METHOD == 'numerical':
    def w(p, T):
        w_096 = w_1(0.96 * p_sat(T), T)
        w_104 = w_2(1.04 * p_sat(T), T)
        w_int = w_096 + ((w_104 - w_096) * (phi(p, T) - 0.96) / 0.08)
        return conditional(le(T, 647.15),
                       conditional(le(phi(p, T), 0.96), w_1(p, T),
                                   conditional(le(phi(p, T), 1.04), w_int,
                                               w_2(p, T))), 0.0)

    def dwdt(p, T, p_n, T_n, p_n_theta=0, T_n_theta=0):
        if p_n_theta == 0:
            p_n_theta = p_n
            T_n_theta = T_n
        delta = 1e-8
        dP = delta * p_n_theta
        dT = delta * T_n_theta
        dwdP = ((w(p_n_theta + dP / 2, T_n_theta) - w(p_n_theta - dP / 2, T_n_theta)) / (dP))
        dwdT = ((w(p_n_theta, T_n_theta + dT / 2) - w(p_n_theta, T_n_theta - dT / 2)) / (dT))
        dPdt = (p - p_n) / dt
        dTdt = (T - T_n) / dt
        return dwdP * dPdt + dwdT * dTdt


elif SORPTION_METHOD == 'smooth_num':
    def w(p, T):
        w_int = -3887.0 * w_c * (0.96 * w_0 / w_c)**(1 / m(T)) - \
                169.0 * w_c * (0.96 * w_0 / w_c)**(1 / m(T)) / m(T) + \
                3888.0 * w_c * ((-(T - 273.15)**2 / 360000 + 1.0) * 0.3335) \
                - 149.76 * w_c * 0.037 + \
                phi(p, T)**3 * (3906.25 * w_c * (0.96 * w_0 / w_c)**(1 / m(T)) + \
                162.760416666667 * w_c * (0.96 * w_0 / w_c)**(1 / m(T)) / m(T) - \
                3906.25 * w_c * ((-(T - 273.15)**2 / 360000 + 1.0) * 0.3335) + \
                156.25 * w_c * 0.037) + \
                phi(p, T)**2 * (-11718.75 * w_c * (0.96 * w_0 / w_c)**(1 / m(T))\
                - 494.791666666667 * w_c * (0.96 * w_0 / w_c)**(1 / m(T)) / m(T) +\
                11718.75 * w_c * ((-(T - 273.15)**2 / 360000 + 1.0) * 0.3335) - \
                462.5 * w_c * 0.037) + \
                phi(p, T) * (11700.0 * w_c * (0.96 * w_0 / w_c)**(1/m(T)) + \
                501.041666666667 * w_c * (0.96 * w_0 / w_c)**(1 / m(T)) / m(T) - \
                11700.0 * w_c * ((-(T - 273.15)**2 / 360000 + 1.0) * 0.3335) + \
                456.0 * w_c * 0.037)
    
        return conditional(le(T, 647.15),
                           conditional(le(phi(p, T), 0.96), w_1(p, T),
                                       conditional(le(phi(p, T), 1.04), w_int,
                                                   w_2(p, T))), 0.0)

    def dwdt(p, T, p_n, T_n, p_n_theta=0, T_n_theta=0):
        if p_n_theta == 0:
            p_n_theta = p_n
            T_n_theta = T_n        
        delta = 1e-8
        dP = delta * p_n_theta
        dT = delta * T_n_theta
        dwdP = ((w(p_n_theta + dP / 2, T_n_theta) - w(p_n_theta - dP / 2, T_n_theta)) / (dP))
        dwdT = ((w(p_n_theta, T_n_theta + dT / 2) - w(p_n_theta, T_n_theta - dT / 2)) / (dT))
        dPdt = (p - p_n) / dt
        dTdt = (T - T_n) / dt
        return dwdP * dPdt + dwdT * dTdt


elif SORPTION_METHOD == 'smooth':
    def w(p, T):
        w_int = -3887.0 * w_c * (0.96 * w_0 / w_c)**(1 / m(T)) - \
                169.0 * w_c * (0.96 * w_0 / w_c)**(1 / m(T)) / m(T) + \
                3888.0 * w_c * ((-(T - 273.15)**2 / 360000 + 1.0) * 0.3335) \
                - 149.76 * w_c * 0.037 + \
                phi(p, T)**3 * (3906.25 * w_c * (0.96 * w_0 / w_c)**(1 / m(T)) + \
                162.760416666667 * w_c * (0.96 * w_0 / w_c)**(1 / m(T)) / m(T) - \
                3906.25 * w_c * ((-(T - 273.15)**2 / 360000 + 1.0) * 0.3335) + \
                156.25 * w_c * 0.037) + \
                phi(p, T)**2 * (-11718.75 * w_c * (0.96 * w_0 / w_c)**(1 / m(T))\
                - 494.791666666667 * w_c * (0.96 * w_0 / w_c)**(1 / m(T)) / m(T) +\
                11718.75 * w_c * ((-(T - 273.15)**2 / 360000 + 1.0) * 0.3335) - \
                462.5 * w_c * 0.037) + \
                phi(p, T) * (11700.0 * w_c * (0.96 * w_0 / w_c)**(1/m(T)) + \
                501.041666666667 * w_c * (0.96 * w_0 / w_c)**(1 / m(T)) / m(T) - \
                11700.0 * w_c * ((-(T - 273.15)**2 / 360000 + 1.0) * 0.3335) + \
                456.0 * w_c * 0.037)
    
        return conditional(le(T, 647.15),
                           conditional(le(phi(p, T), 0.96), w_1(p, T),
                                       conditional(le(phi(p, T), 1.04), w_int,
                                                   w_2(p, T))), 0.0)
    
    
    # Analytical derivative of sorption isotherm
    def dP_sdT(p, T):
        x = T - 273.15
        A1 = 531279.337303094
        B1 = 8.07131
        C1 = 1730.63
        D1 = 233.426
        A2 = 555933.390207996
        B2 = 8.14019
        C2 = 1810.94
        D2 = 244.485
        return conditional(le(T, 373.15),
                           A1 * 10.0**(B1 - C1 / (x + D1)) / (x + D1)**2,
                           A2 * 10.0**(B2 - C2 / (x + D2)) / (x + D2)**2)
    
    
    def dmdT(p, T):
        return -(((-2 * T + 526.3) * (T - 263.15)**2 /
                  ((T - 263.15)**2 + 22.34 * (-263.15 + 298.15)**2)**2) +
                 (2 * T - 526.3) / ((T - 263.15)**2 + 22.34 * (-263.15 + 298.15)**2))
    
    
    def dphidp(p, T):
        return 1 / p_sat(T)
    
    
    def dphidT(p, T):
        return -p * dP_sdT(p, T) / p_sat(T)**2
    
    
    def dw_1dP(p, T):
        return w_c * (w_0 * phi(p, T) / w_c)**(1 / m(T)) * \
            dphidp(p, T) / (phi(p, T) * m(T))
    
    
    def dw_intdP(p, T):
        return 3 * phi(p, T)**2 * dphidp(p, T) * \
        (3906.25 * w_c * (0.96 * w_0 / w_c)**(1 / m(T)) +\
        162.760416666667 * w_c * (0.96 * w_0 / w_c)**(1 / m(T)) / m(T) - \
        3906.25 * w_c * ((-(T - 273.15)**2 / 360000 + 1.0) * 0.3335) + \
        156.25 * w_c * 0.037) + 2 * phi(p, T) * dphidp(p, T) * \
        (-11718.75 * w_c * (0.96 * w_0 / w_c)**(1 / m(T)) - \
        494.791666666667 * w_c * (0.96 * w_0 / w_c)**(1 / m(T)) / m(T) + \
        11718.75 * w_c * ((-(T - 273.15)**2 / 360000 + 1.0) * 0.3335) - \
        462.5 * w_c * 0.037) + dphidp(p, T) * (11700.0 * w_c * \
        (0.96 * w_0 / w_c)**(1 / m(T)) + 501.041666666667 * w_c * \
        (0.96 * w_0 / w_c)**(1 / m(T)) / m(T) - 11700.0 * w_c * \
        ((-(T - 273.15)**2 / 360000 + 1.0) * 0.3335) + 456.0 * w_c * 0.037)
    
    def dw_2dP(p, T):
        return 0.037 * w_c * dphidp(p, T)
    
    
    def dwdp(p, T):
        return conditional(le(T, 647.15),
                    conditional(le(phi(p, T), 0.96), dw_1dP(p, T),
                                conditional(le(phi(p, T), 1.04),
                                              dw_intdP(p, T),
                                              dw_2dP(p, T))), 0.0)
    
    
    def dw_1dT(p, T):
        return w_c * (w_0 * phi(p, T) / w_c)**(1 / m(T)) * \
            (-ln(w_0 * phi(p, T) / w_c) * dmdT(p, T) / m(T)**2 +
             dphidT(p, T) / (phi(p, T) * m(T)))
    
    
    def dw_intdT(p, T):
        return 3887.0 * w_c * (0.96 * w_0 / w_c)**(1 / m(T)) * \
                ln(0.96 * w_0 / w_c) * dmdT(p, T) / m(T)**2 + \
                169.0 * w_c * (0.96 * w_0 / w_c)**(1/m(T))*dmdT(p, T) / m(T)**2\
                + 169.0 * w_c * (0.96 * w_0 / w_c)**(1 / m(T)) * \
                ln(0.96 * w_0 / w_c) * dmdT(p, T) / m(T)**3 + \
                3888.0 * w_c * (-T / 180000 + 0.0015175) * 0.3335 + \
                phi(p, T)**3 * (-3906.25 * w_c * (0.96 * w_0 / w_c)**(1 / m(T))*\
                ln(0.96 * w_0 / w_c) * dmdT(p, T) / m(T)**2 - \
                162.760416666667 * w_c * (0.96 * w_0 / w_c)**(1 / m(T)) * \
                dmdT(p, T) / m(T)**2 - 162.760416666667 * w_c * \
                (0.96 * w_0 / w_c)**(1 / m(T)) * ln(0.96 * w_0 / w_c) * \
                dmdT(p, T) / m(T)**3 - 3906.25 * w_c * (-T / 180000 + 0.0015175) * 0.3335)\
                + 3 * phi(p, T)**2 * dphidT(p, T) * \
                (3906.25 * w_c * (0.96 * w_0 / w_c)**(1 / m(T)) + \
                162.760416666667 * w_c * (0.96 * w_0 / w_c)**(1 / m(T)) / m(T) - \
                3906.25 * w_c*((-(T - 273.15)**2 / 360000 + 1.0) * 0.3335) + \
                156.25 * w_c * 0.037) + phi(p, T)**2 * (11718.75 * w_c * \
                (0.96 * w_0 / w_c)**(1 / m(T)) * ln(0.96 * w_0/w_c) * \
                dmdT(p, T) / m(T)**2 + 494.791666666667 * w_c * \
                (0.96 * w_0 / w_c)**(1 / m(T)) * dmdT(p, T) / m(T)**2 + \
                494.791666666667 * w_c * (0.96 * w_0 / w_c)**(1 / m(T)) * \
                ln(0.96 * w_0 / w_c) * dmdT(p, T) / m(T)**3 + \
                11718.75 * w_c * (-T / 180000 + 0.0015175) * 0.3335) + \
                2*phi(p, T) * dphidT(p, T) * (-11718.75 * w_c * \
                (0.96 * w_0 / w_c)**(1 / m(T)) - 494.791666666667 * w_c * \
                (0.96 * w_0 / w_c)**(1 / m(T)) / m(T) + 11718.75 * w_c * \
                ((-(T - 273.15)**2 / 360000 + 1.0) * 0.3335) - 462.5 * w_c * 0.037)\
                + phi(p, T) * (-11700.0 * w_c * (0.96 * w_0 / w_c)**(1 / m(T)) * \
                ln(0.96 * w_0 / w_c) * dmdT(p, T) / m(T)**2 - 501.041666666667 *\
                w_c * (0.96 * w_0 / w_c)**(1 / m(T)) * dmdT(p, T) / m(T)**2 - \
                501.041666666667 * w_c * (0.96 * w_0 / w_c)**(1 / m(T)) * \
                ln(0.96 * w_0 / w_c) * dmdT(p, T)/m(T)**3 - 11700.0 * w_c * \
                (-T / 180000 + 0.0015175) * 0.3335) + dphidT(p, T) * \
                (11700.0 * w_c * (0.96 * w_0 / w_c)**(1 / m(T)) + \
                501.041666666667 * w_c * (0.96 * w_0 / w_c)**(1 / m(T)) / m(T) - \
                11700.0 * w_c * ((-(T - 273.15)**2 / 360000 + 1.0) * 0.3335) + \
                456.0*w_c*0.037)
    
    
    
    def dw_2dT(p, T):
        return w_c * (-1.8527777778e-6 * T + 0.037 * dphidT(p, T) + 0.00050608625)

    def dwdT(p, T):
        return conditional(le(T, 647.15),
                           conditional(le(phi(p, T), 0.96), dw_1dT(p, T),
                                       conditional(le(phi(p, T), 1.04),
                                                   dw_intdT(p, T),
                                                   dw_2dT(p, T))), 0.0)

    # Central differentiation formula for sorption isotherm
    def dwdt(p, T, p_n, T_n, p_n_theta=0, T_n_theta=0):
        if p_n_theta == 0:
            p_n_theta = p_n
            T_n_theta = T_n
        dPdt = (p - p_n) / dt
        dTdt = (T - T_n) / dt
        return (dwdp(p_n_theta, T_n_theta) * dPdt
                + dwdT(p_n_theta, T_n_theta) * dTdt)

# Latent heat of evaporation of water
def DeltaH_e(T):
    return conditional(le(T, 647.3), 3.5e5 * ((647.3 - T)**(1 / 3)), 0.0)


# Permeability
def K_t(T):
    return (T - 273.15) * 0.95 / 70. - 0.28928571


def f1(T, phi):
    return conditional(lt(phi, 1),
                       (K_t(T) + (1. - K_t(T)) / (1. + (4. * (1. - phi))**4.)),
                       1)


def f2(T):
    return exp(2700 * (1. / (273.15 + 25.) - 1. / (T)))


def f3(T):
    return exp(((T - 273.15) - 95.)
               / (0.881 + 0.214 * ((T - 273.15) - 95.)))


def K(p, T):
    return conditional(le(T, 368.15), K_0 * f1(T, phi(p, T)) * f2(T),
                       K_0 * 5.6 * f3(T))


if N_DIM == 2:
    K_ratio = 1e6

    def lambda_inclusion(p, T):
        return 0.0262

    def rho_inclusion(p, T):
        return 1.2754

    def C_pinclusion(p, T):
        return 1006

    def K_inclusion(p, T):
        return Constant(1e-12 * 1e6)

    def K_matrix(p, T):
        return conditional(le(T, 368.15), K_0 * f1(T, phi(p, T)) * f2(T),
                           K_0 * 5.6 * f3(T))

# Variational formulation in residual form
if N_DIM == 1:
    orderp = 1
    orderT = 1
    FE_p = FiniteElement("P", mesh.ufl_cell(), orderp)
    FE_T = FiniteElement("P", mesh.ufl_cell(), orderT)

    element = MixedElement([FE_p, FE_T])
    X = FunctionSpace(mesh, element)
    Tst = TestFunction(X)
    w_h, v_h = split(Tst)
    P_0 = p_inf
    p_n = interpolate(P_0, X.sub(0).collapse())
    T_0 = Constant(298.15)
    T_n = interpolate(T_0, X.sub(1).collapse())
    u = Function(X)
    p, T = split(u)
    T_n_theta = theta * T + (1 - theta) * T_n
    p_n_theta = theta * p + (1 - theta) * p_n

    T_n_alpha = alpha * T + (1 - alpha) * T_n
    p_n_alpha = alpha * p + (1 - alpha) * p_n

    if not RADIATION:
        bc1 = DirichletBC(X.sub(1), T_HUC, left)
        bcs = [bc1]
    else:
        bcs=[]
    # Mass balance equation equation
    Resp = dwdt(p, T, p_n, T_n, p_n_theta, T_n_theta) * w_h * dx
    Resp += (K(p_n_theta, T_n_theta) / g) * inner(nabla_grad(p_n_alpha),
                                                  nabla_grad(w_h)) * dx
    Resp += - ((w_d(T) - w_d(T_n)) / dt) * w_h * dx
    Resp += beta * (p - p_inf) * w_h * (ds(1) + ds(2))

    # Energy balance equation
    ResT = rho_c * C_pc * ((T - T_n) / dt) * v_h * dx
    ResT += lambda_c * inner(nabla_grad(T_n_alpha), nabla_grad(v_h)) * dx
    ResT += DeltaH_d * ((w_d(T) - w_d(T_n)) / dt) * v_h * dx
    ResT += (- DeltaH_e(T_n) * dwdt(p, T, p_n, T_n, p_n_theta, T_n_theta)
             * v_h * dx)
    ResT += (C_pw * (K(p_n_theta, T_n_theta) / g)
             * inner(nabla_grad(p_n_theta), nabla_grad(T_n_alpha)) * v_h * dx)

    if RADIATION:
        epsilon = 0.8
        sigma_sb = 5.67e-8
        ResT += (epsilon * sigma_sb * (T_n**3 * T - T_HUC**4) + alpha_T * (T - T_HUC)) * v_h * (ds(1))
    # Boundary conditions for energy balance equation
    ResT += alpha_T * (T - T_inf) * v_h * (ds(2))

    # Total residual
    Res = ResT + Resp

    # Jacobian
    Jac = derivative(Res, u)

elif N_DIM == 2:
    if FORMULATION == 'primal':
        # Finite element and space functions definition
        orderp = 1
        orderT = 1
        FE_p = FiniteElement("P", mesh.ufl_cell(), orderp)
        FE_T = FiniteElement("P", mesh.ufl_cell(), orderT)

        element = MixedElement([FE_p, FE_T])
        X = FunctionSpace(mesh, element)

        Tst = TestFunction(X)
        w_h, v_h = split(Tst)

        u = Function(X)
        p, T = split(u)
        u_n = Function(X)
        p_n, T_n = split(u_n)
        P_0 = p_inf
        p_n = interpolate(P_0, X.sub(0).collapse())
        T_0 = Constant(298.15)
        T_n = interpolate(T_0, X.sub(1).collapse())

        bc1 = DirichletBC(X.sub(1), T_HUC, left)
        bcs = [bc1]

        # Mass balance equation
        Resp = dwdt(p, T, p_n, T_n) * w_h * (dx(1) + dx(2))
        Resp += ((K_matrix(p_n, T_n) / g) * inner(nabla_grad(p),
                                                  nabla_grad(w_h)) * dx(1)
                + (K_inclusion(p_n, T_n) / g) * inner(nabla_grad(p),
                                                      nabla_grad(w_h)) * dx(2))
        Resp += - ((w_d(T) - w_d(T_n)) / dt) * w_h * (dx(1) + dx(2))
        Resp += beta * (p - p_inf) * w_h * (ds(1) + ds(2) + ds(3) + ds(4))

        # Energy balance equation
        ResT = (rho_c * C_pc * ((T - T_n) / dt) * v_h * dx(1)
                + rho_inclusion(p_n, T_n) * C_pinclusion(p_n, T_n)
                * ((T - T_n) / dt) * v_h * dx(2))
        ResT += (lambda_c * inner(nabla_grad(T), nabla_grad(v_h)) * dx(1)
                 + lambda_inclusion(p_n, T_n) * inner(nabla_grad(T),
                                                 nabla_grad(v_h)) * dx(2))
        ResT += DeltaH_d * ((w_d(T) - w_d(T_n)) / dt) * v_h * (dx(1) + dx(2))
        ResT += - DeltaH_e(T_n) * dwdt(p, T, p_n, T_n) * v_h * (dx(1) + dx(2))
        ResT += (C_pw * (K_matrix(p_n, T_n) / g) * inner(nabla_grad(p_n),
                                                        nabla_grad(T_n))
                 * v_h * dx(1) + C_pw * (K_inclusion(p_n, T_n) / g) 
                 * inner(nabla_grad(p_n), nabla_grad(T_n)) * v_h * dx(2))
        ResT = ResT + (alpha_T * (T - T_inf)) * v_h * (ds(2))
        # Total residual
        Res = ResT + Resp

        # Jacobian
        Jac = derivative(Res, u)

    if FORMULATION == 'mixed':
        # Mixed element and space functions definition
        orderT = 1
        orderJ = 1
        orderp = 0
        P1 = FiniteElement("P", mesh.ufl_cell(), orderT)
        EM = FiniteElement("RT", mesh.ufl_cell(), orderJ)
        DG = FiniteElement("DG", mesh.ufl_cell(), orderp)

        element = MixedElement([EM, DG, P1])
        X = FunctionSpace(mesh, element)

        Tst = TestFunction(X)
        (tau, w_h, v_h) = split(Tst)

        u = Function(X)
        (Jota, p, T) = split(u)

        u_n = Function(X)
        (Jota_n, p_n, T_n) = split(u_n)

        # Initial conditions
        P_0 = p_inf
        p_n = interpolate(P_0, X.sub(1).collapse())
        T_0 = Constant(298.15)
        T_n = interpolate(T_0, X.sub(2).collapse())

        bc1 = DirichletBC(X.sub(2), T_HUC, left)
        bcs = [bc1]

        # Variational formulation in residual form
        # Mass balance equation equation
        Resp = (dwdt(p, T, p_n, T_n) * w_h * dx(1)
                + dwdt(p, T, p_n, T_n) * w_h * dx(2))
        Resp += div(Jota) * w_h * dx(1) + div(Jota) * w_h * dx(2)
        Resp += (- ((w_d(T) - w_d(T_n)) / dt) * w_h * dx(1)
                 - ((w_d(T) - w_d(T_n)) / dt) * w_h * dx(2))

        # Exterior normal to each facet
        nn = FacetNormal(mesh)

        ResJ = ((g / K_matrix(p_n, T_n)) * inner(Jota, tau) * dx(1)
                + (g / K_inclusion(p_n, T_n)) * inner(Jota, tau) * dx(2))
        ResJ += - div(tau) * p * dx(1) - div(tau) * p * dx(2)
        ResJ += (dot(tau,nn) * ((1 / beta) * dot(Jota,nn) + p_inf)
                 * (ds(1) + ds(2) + ds(3) + ds(4)))

        # Energy balance equation
        ResT = (rho_c * C_pc * ((T - T_n) / dt) * v_h * dx(1)
                + rho_inclusion(p_n, T_n) * C_pinclusion(p_n, T_n)
                * ((T - T_n) / dt) * v_h * dx(2))
        ResT += (lambda_c * inner(nabla_grad(T), nabla_grad(v_h)) * dx(1)
                 + lambda_inclusion(p_n, T_n) * inner(nabla_grad(T),
                                                      nabla_grad(v_h)) * dx(2))
        ResT += (DeltaH_d * ((w_d(T) - w_d(T_n)) / dt) * v_h * dx(1)
                 + DeltaH_d * ((w_d(T) - w_d(T_n)) / dt) * v_h * dx(2))
        ResT += (- DeltaH_e(T_n) * dwdt(p, T, p_n, T_n) * v_h * dx(1)
                 - DeltaH_e(T_n) * dwdt(p, T, p_n, T_n) * v_h * dx(2))
        ResT += (- C_pw * inner(Jota, nabla_grad(T_n)) * v_h * dx(1)
                 - C_pw * inner(Jota, nabla_grad(T_n)) * v_h * dx(2))

        # Boundary conditions for energy balance equation
        ResT += (alpha_T * (T - T_inf)) * v_h * (ds(2) + ds(3) + ds(4))

        # Total residual
        Res = ResT + Resp + ResJ

        # Jacobian
        Jac = derivative(Res, u)

t = 0

# IO Files - Postprocessing
time = []
convergence = []

filex = XDMFFile(saving_dir + '/results.xdmf')
filex.parameters['functions_share_mesh'] = True
filex.parameters['rewrite_function_mesh'] = False
filex.parameters["flush_output"] = True

# Nonlinear problem and solver parameters
problem = NonlinearVariationalProblem(Res, u, bcs, Jac, ffc_options)
solver = NonlinearVariationalSolver(problem)
prm = solver.parameters
# prm["newton_solver"]["relaxation_parameter"] = 0.8
prm["newton_solver"]["absolute_tolerance"] = 1E-8
prm["newton_solver"]["relative_tolerance"] = 1E-15
prm["newton_solver"]["maximum_iterations"] = 7
prm['newton_solver']['error_on_nonconvergence'] = True

if N_DIM == 1:
    # IO Files - Postprocessing
    q_st_domain = []
    q_ev_domain = []
    q_conv_domain = []
    q_dehyd_domain = []
    w_domain = []
    time = []
    convergence = []

    # Projections
    P1 = FiniteElement('P', mesh.ufl_cell(), 1)
    Q = FunctionSpace(mesh, P1)
    q = TrialFunction(Q)
    r = TestFunction(Q)

    a_r_q = q * r * dx
    L_w = w(p_n, T_n) * r * dx
    L_h = phi(p_n, T_n) * r * dx
    L_a = (K(p_n, T_n)) * r * dx
    w_proj = Function(Q, name='Water=w(p, T)', label='w')
    h_proj = Function(Q, name='phi=p/p_sat(T)', label='phi')
    a_proj = Function(Q, name='Permeability=a(p, T)', label='a')
    nt = 0

    # Calculate time of simulation
    f = open(saving_dir + '/time_evolution.csv', 'w')
    writer = csv.writer(f, delimiter='\t')
    startTime = datetime.now()
    # Loop over time steps
    if HEATING_CASE == 'HUC':
        freq_out = 100
    elif HEATING_CASE == 'ISO':
        freq_out = 5
    nt = 0
    print(f'    Percent Complete\tTime (h)\tT_min (K) \
               \tT_max (K)\tp_min (Pa)\tp_max: (Pa)')
    while t <= T_total:
        # Solve non-linear problem
        T_HUC.t = t
        if t < 5.833 * 3600:
            T_HUC.rate = 30
            T_HUC.T_0 = 298.15
        elif (t > 5.833 * 3600) & (t < 15.833 * 3600):
            T_HUC.rate = 0
            T_HUC.T_0 = 473.15
        elif (t > 15.833 * 3600):
            T_HUC.rate = 30
            T_HUC.t_0 = 15.833 * 3600
            T_HUC.T_0 = 473.15
        n, conv = solver.solve()

        # Integrals: calculating energy  and water

        # Integration over domain of rhos*Cs*T
        q_st = assemble(rho_c * C_pc * ((T - T_n) / dt) * dx,
                        form_compiler_parameters=ffc_options)
        q_st_domain.append(q_st)

        # Integration over domain of div(gradp, gradT)
        q_conv = assemble(C_pw * (K(p_n, T_n) / g)
                          * inner(nabla_grad(p_n), nabla_grad(T_n)) * dx,
                          form_compiler_parameters=ffc_options)
        q_conv_domain.append(q_conv)

        # Integration over domain of C_a*dwdt
        q_ev = assemble(- DeltaH_e(T_n)
                        * dwdt(p, T, p_n, T_n, p_n_theta, T_n_theta)
                        * dx, form_compiler_parameters=ffc_options)
        q_ev_domain.append(q_ev)

        # Integration over domain of h_d*dwddt
        q_dehyd = assemble(- DeltaH_d * ((w_d(T) - w_d(T_n)) / dt) * dx,
                           form_compiler_parameters=ffc_options)
        q_dehyd_domain.append(q_dehyd)

        # integration over domain of water quantity
        wat = assemble((w(p_n, T_n)) * dx,
                       form_compiler_parameters=ffc_options)
        w_domain.append(wat)
        convergence.append(n)
        time.append(t)

        # Solve projections
        solve(a_r_q == L_w, w_proj, [], form_compiler_parameters=ffc_options)
        solve(a_r_q == L_h, h_proj, [], form_compiler_parameters=ffc_options)
        solve(a_r_q == L_a, a_proj, [], form_compiler_parameters=ffc_options)

        # Update solution with last computed value
        (_P, _T) = u.split(True)
        p_n.vector()[:] = _P.vector()
        T_n.vector()[:] = _T.vector()
        T_n_theta = theta * T + (1 - theta) * T_n
        p_n_theta = theta * p + (1 - theta) * p_n
        T_n_alpha = alpha * T + (1 - alpha) * T_n
        p_n_alpha = alpha * p + (1 - alpha) * p_n
        p_max = max(p_n.vector()[:]) / 1e6
        p_min = min(p_n.vector()[:]) / 1e6
        T_max = max(T_n.vector()[:])
        T_min = min(T_n.vector()[:])

        writer.writerow([t, q_st, q_conv, q_ev, q_dehyd,
                         wat, n, T_max, p_max])

        if (nt % freq_out == 0):
            _P.rename("Pressure [Pa]", "p")
            _T.rename("Temperature [K]", "T")
            filex.write(_P, t)
            filex.write(_T, t)
            filex.write(w_proj, t)
            filex.write(h_proj, t)
            filex.write(a_proj, t)
        print(f'    {t / T_total * 100:3.2f}%\t\t{t/3600:.2f}\t\t{T_min:.2f} \
               \t\t{T_max:.2f}\t\t{p_min:.4f}\t\t{p_max:.4f}', end='\r')
        nt += 1
        t += dt
    # End loop over time steps
    f.close()
    time_delta = datetime.now() - startTime
    print('Simulation time: ', str(time_delta))

elif N_DIM == 2:
    # Projections to save the results
    P1 = FiniteElement("P", mesh.ufl_cell(), 1)
    Q = FunctionSpace(mesh, P1)
    q = TrialFunction(Q)
    r = TestFunction(Q)

    a_r_q = q * r * (dx(1) + dx(2))
    L_w = w(p_n, T_n) * r * (dx(1) + dx(2))
    L_h = phi(p_n, T_n) * r * (dx(1) + dx(2))
    L_a = ((K_matrix(p_n, T_n) * r * dx(1)))
    w_proj = Function(Q, name='Water=w(p, T)', label='w')
    h_proj = Function(Q, name='phi=p/p_sat(T)', label='phi')
    a_proj = Function(Q, name='Permeability=a(p, T)', label='a')

    # Loop over time steps
    freq_out = 5
    nt = 0
    # Calculate time of simulation
    f = open(saving_dir + '/time_evolution.csv', 'w')
    writer = csv.writer(f, delimiter='\t')
    startTime = datetime.now()
    while t <= T_total:
        print('Progress: ' + str(round(t / T_total * 100, 2)) + '%', end='\r')
        # Solve non-linear problem
        T_HUC.t = t
        if HEATING_CASE == 'HUC':
            if t < 10 * 3600:
                T_HUC.rate = rate_
                T_HUC.T_0 = 298.15

        n, conv = solver.solve()

        # integration over domain of water quantity
        convergence.append(n)
        time.append(t)

        # Solve projections
        solve(a_r_q == L_w, w_proj, [], form_compiler_parameters=ffc_options)
        solve(a_r_q == L_h, h_proj, [], form_compiler_parameters=ffc_options)
        solve(a_r_q == L_a, a_proj, [], form_compiler_parameters=ffc_options)

        # Update solution with last computed value
        if FORMULATION == 'primal':
            _P, _T = u.split(True)

        if FORMULATION == 'mixed':
            _Jota, _P, _T = u.split(True)

        p_n.vector()[:] = _P.vector()
        T_n.vector()[:] = _T.vector()
        p_max = max(p_n.vector()[:]) / 1e6
        p_min = min(p_n.vector()[:]) / 1e6
        T_max = max(T_n.vector()[:])
        T_min = min(T_n.vector()[:])
        if VERBOSE:
            print('  T_min:', T_min, 'K', '\n  T_max:', T_max, 'K')
            print('  p_min:', p_min, 'MPa', '\n p_max:', p_max, 'MPa')

        writer.writerow([t, n, T_max, p_max])

        if (nt % freq_out == 0):
            if FORMULATION == 'primal':
                _P.rename("Pressure [Pa]", "p")
                _T.rename("Temperature [K]", "T")
                w_proj.rename("Sorption [Kg/m3]", "w")
                h_proj.rename("Rel. Humidity [-]", "phi")
                a_proj.rename("Peremeability [m/s]", "a")
                filex.write(_P, t)
                filex.write(_T, t)
                filex.write(w_proj, t)
                filex.write(h_proj, t)
                filex.write(a_proj, t)

            if FORMULATION == 'mixed':
                _Jota.rename("Mass Flux [kg/(m2s)]", "J")
                _P.rename("Pressure [Pa]", "p")
                _T.rename("Temperature [K]", "T")
                w_proj.rename("Sorption [Kg/m3]", "w")
                h_proj.rename("Rel. Humidity [-]", "phi")
                a_proj.rename("Peremeability [m/s]", "a")
                filex.write(_Jota, t)
                filex.write(_P, t)
                filex.write(_T, t)
                filex.write(w_proj, t)
                filex.write(h_proj, t)
                filex.write(a_proj, t)

        nt += 1
        t += dt
    # End loop over time steps
    f.close()
    time_delta = datetime.now() - startTime
    print('Simulation time: ', str(time_delta))
