import matplotlib.pyplot as plt
from scipy.spatial.distance import euclidean
import numpy as np
import os
import pygmsh
import meshio

'''
generate_2D_meshes.py is a helper script aiming to create 2D meshes with
elliptic inclusions. The user can define the limits of sizes and aspect ratio
of the ellpises, and test the effect of increasing porisity
(area of inclusions divided by the total area). The results of the companion
paper, Moreira et al. "Finite element implementation of a thermohygro model for
porous materials at high temperatures using the FEniCS platform", 2020 were
obtained with a constant porosity of roughly 6%. The meshes can be user 
defined, but default values are LONG (4.5cm) and SHORT(2cm).
'''
SIZE = 'LONG'
# SIZE = 'SHORT'

base_dir = './meshes/'
if not os.path.exists(base_dir):
    os.mkdir(base_dir)

class eli_particle:
    def __init__(self, position, a, b, theta):
        self.pos = position
        self.a = a
        self.b = b
        self.radius = max([a, b])
        self.theta = theta


def inv_porosity(n, lx, ly):
    return lx * ly * n


def ellipse_area(a, b):
    return np.pi * a * b


def paramteric_ellipse(pos, a, b, theta, t):
    x = pos[0] + a * np.sin(t) * np.cos(theta) - b * np.cos(t) * np.sin(theta)
    y = pos[1] + b * np.cos(t) * np.cos(theta) + a * np.sin(t) * np.sin(theta)
    return (x, y)


def newpos_eli(lx, ly, max_radius, max_ar):
    x_lim = np.array([0, lx])
    y_lim = np.array([0, ly])
    x = np.random.uniform(x_lim[0] + 0.2 * max_radius,
                          x_lim[1] - 0.2 * max_radius)
    y = np.random.uniform(y_lim[0] + 0.2 * max_radius,
                          y_lim[1] - 0.2 * max_radius)
    return [x, y]


def eli_porosity(n, radii, aspect_ratio, thetas,
                 lx, ly, sep_factor=1, lim_iters=150, verbose=False):
    pore_area = inv_porosity(n, lx, ly)
    theta = np.random.choice(thetas)
    max_radius = radii.max()
    max_ar = max(aspect_ratio)
    pos = newpos_eli(lx, ly, max_radius, max_ar)
    a = np.random.choice(radii)
    b = a * np.random.choice(aspect_ratio)
    particles = [eli_particle(pos, a, b, theta)]
    inc_area = ellipse_area(particles[0].a, particles[0].b)
    j = 0
    while inc_area < pore_area:
        # print(inc_area)
        a = np.random.choice(radii)
        b = a * np.random.choice(aspect_ratio)
        radius = max([a, b])
        pos = newpos_eli(lx, ly, max_radius, max_ar)
        theta = np.random.choice(thetas)
        try:
            sort_dist = sorted(particles, key=lambda p: euclidean(
                p.pos, pos) - (p.radius + radius))
            mindist = euclidean(
                sort_dist[0].pos, pos) - (sort_dist[0].radius + radius)
            mindist_radius = sort_dist[0].radius
            i = 0
            while (mindist <= sep_factor * (mindist_radius + radius)):
                # print('oi')
                pos = newpos_eli(lx, ly, max_radius, max_ar)
                sort_dist = sorted(particles, key=lambda p: euclidean(
                    p.pos, pos) - (p.radius + radius))
                mindist = euclidean(
                    sort_dist[0].pos, pos) - (sort_dist[0].radius + radius)
                mindist_radius = sort_dist[0].radius
                i += 1
                if i > lim_iters:
                    if verbose:
                        print(f'Non Overlapping Failed in Iteration {j}')
                    mindist = 2 * sep_factor * (mindist_radius + radius)
        except:
            if verbose:
                print('Problem Found in Rearrangement Algorithm')
        particles.append(eli_particle(pos, a, b, theta))
        inc_area += ellipse_area(a, b)
        j += 1
        if j > lim_iters:
            old_area = inc_area
            inc_area = 1 * lx * ly
        # print('This circle', circle_area(particle(pos, radius).radius)
    if verbose:
        if j < lim_iters:
            print(
                f'Number of Iterations: {j} \nDefined \
                \\ Obtained Porisity: {round(n*100, 2)}% \
                \\ {round(inc_area/(lx*ly)*100, 2)}%')
        if j > lim_iters:
            print(
                f'Number of Iterations: {j} \nDefined \\ \
                 Obtained Porisity: {round(n*100, 2)}% \\ \
                 {round(old_area/(lx*ly)*100, 2)}%')
    return particles


if SIZE == 'LONG':
    np.random.seed(2303)
    mu_radii = 0.045
    mu_ab_ratio = 0.025  # 0.050

elif SIZE == 'SHORT':
    np.random.seed(2303)
    mu_radii = 0.02
    mu_ab_ratio = 0.045  # 0.095

sigma_radii = mu_radii * 0.2
normal_radii = np.random.normal(mu_radii, sigma_radii, 50)

sigma_ab_ratio = mu_ab_ratio * 0.05
normal_ab_ratio = np.random.normal(mu_ab_ratio, sigma_ab_ratio, 50)

particles = eli_porosity(0.06, normal_radii, normal_ab_ratio,
                         np.linspace(0, 2 * np.pi),
                         0.2, 0.2, lim_iters=50, sep_factor=0.001)
t = np.linspace(0, 2 * np.pi, 20)
lx = 0.2
ly = 0.2
plt.figure(figsize=(6, 6))
for i, ptc in enumerate(particles):
    plt.plot(*paramteric_ellipse(ptc.pos, ptc.a, ptc.b,
                                 ptc.theta, t), label=i, c='k', ls='-')
plt.xlim(-0.0005, 0.2005)
plt.ylim(-0.0005, 0.2005)
square = np.transpose([[0, 0], [lx, 0], [lx, ly], [0, ly], [0, 0]])
plt.plot(*square, c='k')
plt.xlim(-0.0015, 0.2005)
plt.ylim(-0.0005, 0.2015)
plt.axis('off')
plt.show()

geom = pygmsh.opencascade.Geometry()
# geom.add_raw_code('Mesh.Algorithm = 6;')  # Specify the Frontal Algorithm

square = geom.add_rectangle((0, 0, 0), lx, ly, 0)
p_mid = geom.add_point((0, 0.1, 0), lcar=0.09)

ellips = []
for particle in particles[:]:
    if particle.a > particle.b:
        ellips.append(geom.add_disk((*particle.pos, 0),
                                    radius0=abs(particle.a), radius1=abs(particle.b)))
    else:
        ellips.append(geom.add_disk((*particle.pos, 0),
                                    radius0=abs(particle.b), radius1=abs(particle.a)))

for ellip, particle in zip(ellips, particles[:]):
    geom.rotate(ellip, (*particle.pos, 0), particle.theta, (0, 0, 1))

all_ellips = geom.boolean_union(ellips, delete_other=True)
ellips_inside = geom.boolean_intersection(
    [all_ellips, square], delete_other=False)
square_mask = geom.boolean_difference(
    [square], [ellips_inside], delete_first=True, delete_other=False)

geom.add_physical(ellips_inside, label=2)
geom.add_physical(square_mask, label=1)
# pts_ellips = [i.replace('[][]', '[]')
#               for i in ellips_inside.char_length_code(l_car_ellips)][0][:9]
# pts_square = [i.replace('[][]', '[]')
#               for i in square_mask.char_length_code(l_car_square)][0][:9]

# Use Field to control mesh size
geom.add_raw_code('Field[1] = Distance;')
geom.add_raw_code('lns_bo2[] = Abs(Boundary{Surface{bo2[]}; });')
geom.add_raw_code('Field[1].EdgesList = {lns_bo2[]};')
geom.add_raw_code('Field[1].NNodesByEdge = 1000;')
geom.add_raw_code('Field[2] = Threshold;')
geom.add_raw_code('Field[2].IField = 1;')
geom.add_raw_code('Field[2].LcMin = 0.0015;')
# geom.add_raw_code('Field[2].LcMax = 0.0225;')
geom.add_raw_code('Field[2].LcMax = 0.0150;')
geom.add_raw_code('Field[2].DistMin = 2*0.0015;')
geom.add_raw_code('Field[2].DistMax = 10*0.0025;')
# geom.add_raw_code('Field[2].LcMin = 0.0025;')
# geom.add_raw_code('Field[2].LcMax = 0.0225;')
# geom.add_raw_code('Field[2].DistMin = 2*0.0025;')
# geom.add_raw_code('Field[2].DistMax = 10*0.0025;')
geom.add_raw_code('Field[3] = Min;')
geom.add_raw_code('Field[3].FieldsList = {2};')
geom.add_raw_code('Background Field = 3;')
# geom.add_raw_code('Mesh.Algorithm = 6;')
geom.add_raw_code('Mesh.CharacteristicLengthFactor = 1;')
geom.add_raw_code('Mesh.CharacteristicLengthExtendFromBoundary = 0;')
geom.add_raw_code('Mesh.Smoothing = 40;')

if SIZE == 'LONG':
    mesh = pygmsh.generate_mesh(
        geom, geo_filename="./meshes/long_inclusions.geo", prune_z_0=True)
    meshio.write("./meshes/long_test.vtk", mesh)
    meshio.write("./meshes/long_test.xml", mesh)
    os.rename('./meshes/long_test_gmshgeometrical.xml',
              r'meshes/long_test_gmsh-geometrical.xml')
    os.rename('./meshes/long_test_gmshphysical.xml',
              r'meshes/long_test_gmsh-physical.xml')
elif SIZE == 'SHORT':
    mesh = pygmsh.generate_mesh(
        geom, geo_filename="./meshes/short_inclusions.geo", prune_z_0=True)
    meshio.write("./meshes/short_test.vtk", mesh)
    meshio.write("./meshes/short_test.xml", mesh)
    os.rename('./meshes/short_test_gmshgeometrical.xml',
              './meshes/short_test_gmsh-geometrical.xml')
    os.rename('./meshes/short_test_gmshphysical.xml',
              './meshes/short_test_gmsh-physical.xml')
