import imp

# various try statements to deal with lumerical API paths being different for different computers #
try:
    imp.load_source("lumapi", "C:/Program Files/Lumerical/v211/api/python/lumapi.py")
try:
    imp.load_source("lumapi", "C:/Program Files/Lumerical/v241/api/python/lumapi.py")
except FileNotFoundError:
    print('Check the Lumerical API import path!')

import lumapi

import numpy as np
import pickle
import matplotlib.pyplot as plt
import sys
import time
import json
import os

# You may ask, why Hope are you using an accursed format such as pckl??? Answer: to keep track of all complex numbers in the arrays.
# Complex numbers are not json serializable
# I'm totally leaving this as a note to myself in case I forget this for the third time.

##################################################################################
fdtd = lumapi.FDTD()
fdtd.setresource("FDTD", 1, "processes", "1")
fdtd.setresource("FDTD", 1, "threads", 16)


def sim_vol(tsim, coord, mesh_res, xspan, yspan, zspan, mesh_accuracy):
    x = coord[0]
    y = coord[1]
    z = coord[2]
    fdtd.addfdtd()
    fdtd.set("simulation time", tsim)
    if mesh_res == 0:
        fdtd.set("mesh type", 1)
        fdtd.set("mesh accuracy", mesh_accuracy)
    else:
        fdtd.set("mesh type", 3)
        fdtd.set("define x mesh by", 2)
        fdtd.set("define y mesh by", 2)
        fdtd.set("define z mesh by", 2)
        fdtd.set("dx", mesh_res)
        fdtd.set("dy", mesh_res)
        fdtd.set("dz", mesh_res)
    fdtd.set("dimension", "3D")
    fdtd.set("x", x)
    fdtd.set("y", y)
    fdtd.set("z", z)
    fdtd.set("x span", xspan)
    fdtd.set("y span", yspan)
    fdtd.set("z span", zspan)
    fdtd.set("x min bc", "PML")
    fdtd.set("x max bc", "PML")
    fdtd.set("y min bc", "PML")
    fdtd.set("y max bc", "PML")
    fdtd.set("z min bc", "PML")
    fdtd.set("z max bc", "PML")
    print("Simulation Volume Added")
    return ()


def sim_surface_xy(tsim, coord, mesh_res, xspan, yspan, mesh_accuracy):
    x = coord[0]
    y = coord[1]
    z = coord[2]
    fdtd.addfdtd()
    fdtd.set("simulation time", tsim)
    fdtd.set("dimension", "2D")
    if mesh_res == 0:
        fdtd.set("mesh type", 1)
        fdtd.set("mesh accuracy", mesh_accuracy)
    else:
        fdtd.set("mesh type", 3)
        fdtd.set("define x mesh by", 2)
        fdtd.set("define y mesh by", 2)
        fdtd.set("dx", mesh_res)
        fdtd.set("dy", mesh_res)
    fdtd.set("x", x)
    fdtd.set("y", y)
    fdtd.set("z", z)
    fdtd.set("x span", xspan)
    fdtd.set("y span", yspan)
    fdtd.set("x min bc", "PML")
    fdtd.set("x max bc", "PML")
    fdtd.set("y min bc", "PML")
    fdtd.set("y max bc", "PML")
    print("Simulation Surface (xy) Added")
    return ()


def sim_symmetrical_vol(tsim, coord, mesh_res, xspan, yspan, zspan, mesh_accracy):
    x = coord[0]
    y = coord[1]
    z = coord[2]
    fdtd.addfdtd()
    fdtd.set("simulation time", tsim)
    if mesh_res == 0:
        fdtd.set("mesh type", 1)
        fdtd.set("mesh accuracy", mesh_accracy)
    else:
        fdtd.set("define x mesh by", 1)
        fdtd.set("define y mesh by", 1)
        fdtd.set("define z mesh by", 1)
        fdtd.set("mesh cells per wavelength", mesh_res)
    fdtd.set("dimension", "3D")
    fdtd.set("x", x)
    fdtd.set("y", y)
    fdtd.set("z", z)
    print(xspan)
    fdtd.set("x span", xspan)
    fdtd.set("y span", yspan)
    fdtd.set("z span", zspan)
    fdtd.set("x min bc", "Symmetric")
    fdtd.set("x max bc", "PML")
    fdtd.set("y min bc", "Anti-Symmetric")
    fdtd.set("y max bc", "PML")
    fdtd.set("z min bc", "PML")
    fdtd.set("z max bc", "PML")
    print("Simulation Volume Added")
    return ()


def add_mesh(coord, xspan, yspan, zspan, dx, dy, dz):
    x = coord[0]
    y = coord[1]
    z = coord[2]
    fdtd.addmesh()
    fdtd.set("x", x)
    fdtd.set("y", y)
    fdtd.set("z", z)
    fdtd.set("x span", xspan)
    fdtd.set("y span", yspan)
    fdtd.set("z span", zspan)
    fdtd.set("override x mesh", 1)
    fdtd.set("override y mesh", 1)
    fdtd.set("override z mesh", 1)
    fdtd.set("dx", dx)
    fdtd.set("dy", dy)
    fdtd.set("dz", dz)
    print(f"Manual mesh region added with step size: {dx},{dy},{dz}")
    return ()


def create_diamond():
    fdtd.eval("temp = addmaterial('Dielectric');")
    fdtd.eval("setmaterial(temp,'name','diamond');")
    fdtd.eval("setmaterial('diamond','Permittivity', 2.4114 * 2.4114);")
    fdtd.eval("setmaterial('diamond','color',[0;0.6;0.8;1]);")
    print("Created Diamond Material Entry")
    return "diamond"


def create_sio2():
    print("Assigned SiO2 Material Entry")
    return "SiO2 (Glass) - Palik"


def create_si():
    print("Assigned Si Material Entry")
    return "Si (Silicon) - Palik"


def create_hbn():
    fdtd.eval("temp = addmaterial('Dielectric');")
    fdtd.eval("setmaterial(temp,'name','hbn');")
    permitivities = [1.69**2, 1.69**2, 2.17**2]
    fdtd.eval("setmaterial('hbn', 'Anisotropy', 1);")  # diagonal
    fdtd.eval(f"setmaterial('hbn','Permittivity', {permitivities});")
    fdtd.eval("setmaterial('hbn','color',[0;0.6;1;0.5]);")
    print("Created hBN Material Entry")
    return "hbn"


def create_sic():
    fdtd.eval("temp = addmaterial('Dielectric');")
    fdtd.eval("setmaterial(temp,'name','sic');")
    fdtd.eval("setmaterial('sic','Permittivity', 2.64 * 2.64);")
    fdtd.eval("setmaterial('sic','color',[0;0.6;1;0.5]);")
    print("Created SiC Material Entry")
    return "sic"


def create_al2o3():
    fdtd.eval("temp = addmaterial('Dielectric');")
    fdtd.eval("setmaterial(temp,'name','al2o3');")
    fdtd.eval("setmaterial('al2o3','Permittivity', 1.75 * 1.75);")
    fdtd.eval("setmaterial('al2o3','color',[0;0.6;1;0.5]);")
    print("Created Al2O3 Material Entry")
    return "al2o3"


def create_air():
    fdtd.eval("temp = addmaterial('Dielectric');")
    fdtd.eval("setmaterial(temp,'name','air');")
    fdtd.eval("setmaterial('air','Permittivity', 1 * 1);")
    fdtd.eval("setmaterial('air','color',[0;0.1;0.1;0.1]);")
    print("Created Air Material Entry")
    return "air"


def create_sin():
    print("Assigned SiN Material Entry")
    return "Si3N4 (Silicon Nitride) - Phillip"


def create_sin_2p2():
    fdtd.eval("temp = addmaterial('Dielectric');")
    fdtd.eval("setmaterial(temp,'name','sin2p2');")
    fdtd.eval("setmaterial('sin2p2','Permittivity', 2.2 * 2.2);")
    fdtd.eval("setmaterial('sin2p2','color',[0.85;0;0.7;0.2]);")
    print("Created SiN n=2.2 Material Entry")
    return "sin2p2"


def create_etch():
    print("Assigned Etch Material Entry")
    return "etch"


def create_ln_wvg():
    fdtd.eval("temp = addmaterial('Dielectric');")
    fdtd.eval("setmaterial(temp,'name','ln_wvg');")
    fdtd.eval("setmaterial('ln_wvg','Permittivity', 2.2204 * 2.2204);")
    print("Created LN Wvg Material Entry")
    return "ln_wvg"


def create_ln_substrate():
    fdtd.eval("temp = addmaterial('Dielectric');")
    fdtd.eval("setmaterial(temp,'name','ln_substrate');")
    fdtd.eval("setmaterial('ln_substrate','Permittivity', 2.29 * 2.29);")
    print("Created LN Substrate Material Entry")
    return "ln_substrate"


def substrate(name, coord, material, thick, xspan, yspan):
    x = coord[0]
    y = coord[1]
    z = coord[2]
    fdtd.addrect()
    fdtd.set("name", name)
    fdtd.set("material", material)
    fdtd.set("x", x)
    fdtd.set("y", y)
    fdtd.set("z", z)
    fdtd.set("x span", xspan)
    fdtd.set("y span", yspan)
    fdtd.set("z span", thick)
    print("Created Substrate with Material {}".format(material))
    return ()


def disk(name, coord, material, thick, diameter):
    x = coord[0]
    y = coord[1]
    z = coord[2]
    fdtd.addcircle()
    fdtd.set("name", name)
    fdtd.set("material", material)
    fdtd.set("x", x)
    fdtd.set("y", y)
    fdtd.set("z", z)
    fdtd.set("radius", diameter / 2)
    fdtd.set("z span", thick)
    print("Created Disk Resonator with Material {}".format(material))
    return ()


def ring(name, coord, material, width, thick, diameter):
    x = coord[0]
    y = coord[1]
    z = coord[2]
    fdtd.addcircle()
    fdtd.set("name", name)
    fdtd.set("material", material)
    fdtd.set("x", x)
    fdtd.set("y", y)
    fdtd.set("z", z)
    fdtd.set("radius", diameter / 2)
    fdtd.set("z span", thick)

    fdtd.addcircle()
    fdtd.set("name", name)
    fdtd.set("material", "etch")
    fdtd.set("x", x)
    fdtd.set("y", y)
    fdtd.set("z", z)
    fdtd.set("radius", (diameter - 2 * width) / 2)
    fdtd.set("z span", thick)

    print("Created Ring Resonator with Material {}".format(material))
    return ()


def rectWaveguide(name, coord, material, thick, xspan, yspan):
    x = coord[0]
    y = coord[1]
    z = coord[2]
    fdtd.addrect()
    fdtd.set("name", name)
    fdtd.set("material", material)
    fdtd.set("x", x)
    fdtd.set("y", y)
    fdtd.set("z", z)
    fdtd.set("x span", xspan)
    fdtd.set("y span", yspan)
    fdtd.set("z span", thick)
    print("Created Rectangular Waveguide with Material {}".format(material))
    return ()


def trapWaveguide(name, coord, material, thick, xspan, yspan, angle):
    x = coord[0]
    y = coord[1]
    z = coord[2]
    off = thick * np.tan(angle * np.pi / 180)
    lx_top = yspan - 2 * off
    fdtd.addobject("isos_trpzd_extpoly")
    fdtd.set("name", name)
    fdtd.set("material", material)
    fdtd.set("lx top", lx_top)
    fdtd.set("lx base", yspan)
    fdtd.set("y span", thick)
    fdtd.set("z span", xspan)
    fdtd.set("first axis", "x")
    fdtd.set("rotation 1", 90)
    fdtd.set("second axis", "z")
    fdtd.set("rotation 2", 90)
    fdtd.set("x", x)
    fdtd.set("y", y)
    fdtd.set("z", z)
    print("Created Trapezoidal Waveguide with Material {}".format(material))
    return ()


def taper(name, coord, material, thick, maxWidth, minWidth, xspan, flip):
    x = coord[0]
    y = coord[1]
    z = coord[2]
    fdtd.addobject("isos_trpzd_extpoly")
    fdtd.set("name", name)
    fdtd.set("material", material)
    fdtd.set("lx top", minWidth)
    fdtd.set("lx base", maxWidth)
    fdtd.set("y span", xspan)
    fdtd.set("z span", thick)
    if flip == 1:
        fdtd.set("first axis", "z")
        fdtd.set("rotation 1", 90)
    else:
        fdtd.set("first axis", "z")
        fdtd.set("rotation 1", 270)
    fdtd.set("x", x)
    fdtd.set("y", y)
    fdtd.set("z", z)
    print("Created Adiabatic Taper {}".format(material))
    return ()


def trap_taper(name, coord, material, thick, maxBot, minBot, xspan, angle, flip):
    x = coord[0]
    y = coord[1]
    z = coord[2]

    fdtd.addstructuregroup()
    fdtd.adduserprop("thickness", 2, thick)
    fdtd.adduserprop("angle_side", 0, angle)
    fdtd.adduserprop("width_bot_l", 2, maxBot)
    fdtd.adduserprop("width_bot_r", 2, minBot)
    fdtd.adduserprop("len", 2, xspan)
    fdtd.adduserprop("material", 5, material)

    ######################################################################
    solid_script = "deleteall; \n"
    solid_script += "delta_w = thickness*tan(angle_side*pi/180); \n"
    solid_script += "width_top_l = width_bot_l - 2*delta_w; \n"
    solid_script += "width_top_r = width_bot_r - 2*delta_w; \n"
    solid_script += "zmin = -thickness/2; \n"
    solid_script += "zmax = thickness/2; \n"
    solid_script += "xmin = -len/2; \n"  # 'xmin = 20e-6-len; \n'
    solid_script += "xmax = len/2; \n"  # 'xmax = 20e-6; \n'
    solid_script += "ymin_bot_l = -width_bot_l/2; \n"
    solid_script += "ymax_bot_l = width_bot_l/2; \n"
    solid_script += "ymin_bot_r = -width_bot_r/2; \n"
    solid_script += "ymax_bot_r = width_bot_r/2; \n"
    solid_script += "ymin_top_l = -width_top_l/2; \n"
    solid_script += "ymax_top_l = width_top_l/2; \n"
    solid_script += "ymin_top_r = -width_top_r/2; \n"
    solid_script += "ymax_top_r = width_top_r/2; \n"
    solid_script += "vtx=    [xmin,ymin_bot_l,zmin;\n        xmax,ymin_bot_r,zmin;\n        xmax,ymax_bot_r,zmin;\n        xmin,ymax_bot_l,zmin;\n        xmin,ymin_top_l,zmax;\n        xmax,ymin_top_r,zmax;\n        xmax,ymax_top_r,zmax;\n        xmin,ymax_top_l,zmax]; \n"
    solid_script += "a = cell(6); \n"
    solid_script += "for(i = 1:6){\n    a{i} = cell(1);\n} \n"
    solid_script += "a{1}{1} = [1,4,3,2]; \n"
    solid_script += "a{2}{1} = [1,2,6,5]; \n"
    solid_script += "a{3}{1} = [2,3,7,6]; \n"
    solid_script += "a{4}{1} = [3,4,8,7]; \n"
    solid_script += "a{5}{1} = [1,5,8,4]; \n"
    solid_script += "a{6}{1} = [5,6,7,8]; \n"
    solid_script += "addplanarsolid(vtx,a); \n"
    solid_script += 'if (material=="<Object defined dielectric>"){\n    setnamed("solid", "index",index);} \n'
    solid_script += 'else{\n    setnamed("solid", "material",material);\n} \n'
    ######################################################################

    fdtd.set("x", x)
    fdtd.set("y", y)
    fdtd.set("z", z)
    fdtd.set("script", solid_script)
    fdtd.set("name", name)
    if flip == 1:
        fdtd.set("first axis", "z")
        fdtd.set("rotation 1", 180)
    print("Created trapezoidal taper waveguide")
    return ()


def fiber_taper(name, coord, material, taper_theta, taper_length, rotation=None):
    x = coord[0]
    y = coord[1]
    z = coord[2]

    fdtd.addobject("cone")
    fdtd.set("name", name)
    fdtd.set("material", material)
    fdtd.set("x", x)
    fdtd.set("y", y)
    fdtd.set("z", z)
    fdtd.set("z span", taper_length)
    fdtd.set("theta", taper_theta)
    if rotation is not None:
        fdtd.set("first axis", "x")
        fdtd.set("second axis", "y")
        fdtd.set("third axis", "z")
        fdtd.set("rotation 1", rotation[0])
        fdtd.set("rotation 2", rotation[1])
        fdtd.set("rotation 3", rotation[2])
    print("Created Conical Fiber Taper")
    return ()


def mirror_cavity_holes(
    name, coord, thick, mirror_dx, dia_1, dia_2, num_holes, start_pos
):
    x = coord[0]
    y = coord[1]
    z = coord[2]

    fdtd.addstructuregroup()
    fdtd.set("name", name)  # group all holes in one structure group
    fdtd.set("x", x)
    fdtd.set("y", y)
    fdtd.set("z", z)

    for i in range(num_holes):
        # adding holes symmetrically
        fdtd.addcircle()
        fdtd.set("material", "etch")
        fdtd.set("make ellipsoid", 1)
        fdtd.set("name", "ellipse")
        fdtd.set("radius", dia_1 / 2)
        fdtd.set("radius 2", dia_2 / 2)
        fdtd.set("z span", thick)
        fdtd.set("x", mirror_dx * (i + 1) + start_pos)
        fdtd.addtogroup(name)

        fdtd.addcircle()
        fdtd.set("material", "etch")
        fdtd.set("make ellipsoid", 1)
        fdtd.set("name", "ellipse")
        fdtd.set("radius", dia_1 / 2)
        fdtd.set("radius 2", dia_2 / 2)
        fdtd.set("z span", thick)
        fdtd.set("x", -(mirror_dx * (i + 1) + start_pos))
        fdtd.addtogroup(name)

    extents = (num_holes + 1) * mirror_dx + start_pos

    print("Mirror holes created.")
    print(f"Array extent {extents} is returned.")
    return extents


def gaussian_cavity_holes(
    name, coord, thick, mirror_dx, dia_1, dia_2, gauss_std, num_holes, preview
):
    x = coord[0]
    y = coord[1]
    z = coord[2]

    fdtd.addstructuregroup()
    fdtd.set("name", name)  # group all holes in one structure group
    fdtd.set("x", x)
    fdtd.set("y", y)
    fdtd.set("z", z)

    if num_holes % 2 == 0:
        point_list = [0.5 + i for i in range(num_holes // 2)]

        gauss_list = [np.exp(-(i**2) / 2 / gauss_std**2) for i in point_list]
        gauss_norm_list = [i - min(gauss_list) for i in gauss_list]

        perturb_list = [mirror_dx * (1 - i) for i in gauss_norm_list]

        if preview == 1:
            position_list = [-i for i in point_list[::-1]] + point_list
            holes_list = gauss_norm_list[::-1] + gauss_norm_list
            preview_list = [mirror_dx * (1 - i) for i in holes_list]

            plt.plot(position_list, preview_list, marker="o")
            plt.show()

    else:
        point_list = [i for i in range(num_holes // 2 + 1)]

        gauss_list = [np.exp(-(i**2) / 2 / gauss_std**2) for i in point_list]
        gauss_norm_list = [i - min(gauss_list) for i in gauss_list]

        perturb_list = [mirror_dx * (1 - i) for i in gauss_norm_list]

        if preview == 1:
            gauss_reversed = gauss_norm_list[::-1]
            position_list = [-i for i in point_list[::-1]] + point_list[1::]
            holes_list = gauss_reversed[0:-1] + gauss_norm_list
            preview_list = [mirror_dx * (1 - i) for i in holes_list]

            plt.plot(position_list, preview_list, marker="o")
            plt.show()

    if num_holes % 2 == 1:
        fdtd.addcircle()
        fdtd.set("material", "etch")
        fdtd.set("make ellipsoid", 1)
        fdtd.set("name", "ellipse")
        fdtd.set("radius", dia_1 / 2)
        fdtd.set("radius 2", dia_2 / 2)
        fdtd.set("z span", thick)
        fdtd.set("x", 0)
        fdtd.addtogroup(name)

        perturb_list = perturb_list[1::]

        # calculating extent of the cavity hole array, in one direction
        extents = sum(perturb_list)
        # print(extents)

        indexing = 0

    else:
        indexing = -perturb_list[0] / 2

        # calculating extent of the cavity hole array, in one direction
        extents = sum(perturb_list) - perturb_list[0] / 2
        # print(extents)

    for i in perturb_list:
        # adding holes symmetrically
        fdtd.addcircle()
        fdtd.set("material", "etch")
        fdtd.set("make ellipsoid", 1)
        fdtd.set("name", "ellipse")
        fdtd.set("radius", dia_1 / 2)
        fdtd.set("radius 2", dia_2 / 2)
        fdtd.set("z span", thick)
        fdtd.set("x", indexing + i)
        fdtd.addtogroup(name)

        fdtd.addcircle()
        fdtd.set("material", "etch")
        fdtd.set("make ellipsoid", 1)
        fdtd.set("name", "ellipse")
        fdtd.set("radius", dia_1 / 2)
        fdtd.set("radius 2", dia_2 / 2)
        fdtd.set("z span", thick)
        fdtd.set("x", -(indexing + i))
        fdtd.addtogroup(name)

        indexing += i

    print("Gaussian Tapered Cavity Holes Created.")
    print(f"Array extent {extents} is returned.")
    return extents


def quadratic_cavity_holes(
    name, coord, thick, mirror_dx, dia_1, dia_2, quad_constant, num_holes, preview
):
    x = coord[0]
    y = coord[1]
    z = coord[2]

    fdtd.addstructuregroup()
    fdtd.set("name", name)  # group all holes in one structure group
    fdtd.set("x", x)
    fdtd.set("y", y)
    fdtd.set("z", z)

    if num_holes % 2 == 0:
        point_list = [0.5 + i for i in range(num_holes // 2)]

        quad_list = [quad_constant * i**2 for i in point_list]
        quad_norm_list = [i - max(quad_list) for i in quad_list]

        if max(quad_list) > 1:
            print("Error, Choose smaller quadratic constant")
            sys.exit()

        perturb_list = [mirror_dx * (1 + i) for i in quad_norm_list]

        if preview == 1:
            position_list = [-i for i in point_list[::-1]] + point_list
            holes_list = quad_norm_list[::-1] + quad_norm_list
            preview_list = [mirror_dx * (1 + i) for i in holes_list]

            plt.plot(position_list, preview_list, marker="o")
            plt.show()

    else:
        point_list = [i for i in range(num_holes // 2 + 1)]

        quad_list = [quad_constant * i**2 for i in point_list]
        quad_norm_list = [i - max(quad_list) for i in quad_list]

        if max(quad_list) > 1:
            print("Error, Choose smaller quadratic constant")
            sys.exit()

        perturb_list = [mirror_dx * (1 + i) for i in quad_norm_list]

        if preview == 1:
            gauss_reversed = quad_norm_list[::-1]
            position_list = [-i for i in point_list[::-1]] + point_list[1::]
            holes_list = gauss_reversed[0:-1] + quad_norm_list
            preview_list = [mirror_dx * (1 + i) for i in holes_list]

            plt.plot(position_list, preview_list, marker="o")
            plt.show()

    if num_holes % 2 == 1:
        fdtd.addcircle()
        fdtd.set("material", "etch")
        fdtd.set("make ellipsoid", 1)
        fdtd.set("name", "ellipse")
        fdtd.set("radius", dia_1 / 2)
        fdtd.set("radius 2", dia_2 / 2)
        fdtd.set("z span", thick)
        fdtd.set("x", 0)
        fdtd.addtogroup(name)

        perturb_list = perturb_list[1::]

        # calculating extent of the cavity hole array, in one direction
        extents = sum(perturb_list)
        # print(extents)

        indexing = 0

    else:
        indexing = -perturb_list[0] / 2

        # calculating extent of the cavity hole array, in one direction
        extents = sum(perturb_list) - perturb_list[0] / 2
        # print(extents)

    for i in perturb_list:
        # adding holes symmetrically
        fdtd.addcircle()
        fdtd.set("material", "etch")
        fdtd.set("make ellipsoid", 1)
        fdtd.set("name", "ellipse")
        fdtd.set("radius", dia_1 / 2)
        fdtd.set("radius 2", dia_2 / 2)
        fdtd.set("z span", thick)
        fdtd.set("x", indexing + i)
        fdtd.addtogroup(name)

        fdtd.addcircle()
        fdtd.set("material", "etch")
        fdtd.set("make ellipsoid", 1)
        fdtd.set("name", "ellipse")
        fdtd.set("radius", dia_1 / 2)
        fdtd.set("radius 2", dia_2 / 2)
        fdtd.set("z span", thick)
        fdtd.set("x", -(indexing + i))
        fdtd.addtogroup(name)

        indexing += i

    print("Quadratic Tapered Cavity Holes Created.")
    print(f"Array extent {extents} is returned.")
    return extents


def importgds(name, coord, material, thick, filename, cellname, layer):
    ### Origin will correspond to that defined in the GDS ###
    fdtd.gdsimport(filename, cellname, layer, material, -thick / 2, thick / 2)
    fdtd.set("name", name)
    fdtd.set("x", coord[0])
    fdtd.set("y", coord[1])
    fdtd.set("z", coord[2])
    print("Created GDS Imported Object with Material {}".format(material))


def dipole(coord, wave_center, wave_width, theta, phi, name=""):
    x = coord[0]
    y = coord[1]
    z = coord[2]
    fdtd.adddipole()
    if len(name) > 0:
        fdtd.set("name", name)
    fdtd.set("x", x)
    fdtd.set("y", y)
    fdtd.set("z", z)
    fdtd.set("wavelength start", wave_center - wave_width / 2)
    fdtd.set("wavelength stop", wave_center + wave_width / 2)
    fdtd.set("theta", theta)
    fdtd.set("phi", phi)
    # fdtd.set("dipole type","electric dipole")
    print("Created Dipole Emitter at Coordinates {}".format(coord))
    return ()


def gaussian_x(
    name,
    coord,
    height,
    width,
    mode,
    direction,
    wave_center,
    wave_width,
    theta=0,
    bend_radius=None,
    bend_orientation=None,
    bend_location=[],
    mode_number=1,
):
    x = coord[0]
    y = coord[1]
    z = coord[2]
    fdtd.addmode()
    fdtd.set("name", name)
    fdtd.set("injection axis", "x")
    fdtd.set("direction", direction)
    fdtd.set("x", x)
    fdtd.set("y", y)
    fdtd.set("z", z)
    fdtd.set("y span", width)
    fdtd.set("z span", height)
    if bend_radius != None:
        fdtd.set("bent waveguide", 1)
        fdtd.set("bend radius", bend_radius)
        if bend_orientation != None:
            fdtd.set("bend orientation", bend_orientation)
        if bend_location != []:
            fdtd.set("bend location", "user specified")
            fdtd.set("bend location x", bend_location[0])
            fdtd.set("bend location y", bend_location[1])
            fdtd.set("bend location z", bend_location[2])
    fdtd.set("mode selection", mode)
    if mode != "user select":
        fdtd.set("selected mode number", str(mode_number))
    fdtd.set("wavelength start", wave_center - wave_width / 2)
    fdtd.set("wavelength stop", wave_center + wave_width / 2)
    fdtd.set("theta", theta)
    print("Created Mode Source (x) at Coordinates {}".format(coord))
    return ()


def gaussian_y(
    name,
    coord,
    height,
    width,
    mode,
    direction,
    wave_center,
    wave_width,
    theta=0,
    bend_radius=None,
    bend_orientation=None,
    bend_location=[],
    mode_number=1,
):
    x = coord[0]
    y = coord[1]
    z = coord[2]
    fdtd.addmode()
    fdtd.set("name", name)
    fdtd.set("injection axis", "y")
    fdtd.set("direction", direction)
    fdtd.set("x", x)
    fdtd.set("y", y)
    fdtd.set("z", z)
    fdtd.set("x span", width)
    fdtd.set("z span", height)
    if bend_radius != None:
        fdtd.set("bent waveguide", 1)
        fdtd.set("bend radius", bend_radius)
        if bend_orientation != None:
            fdtd.set("bend orientation", bend_orientation)
        if bend_location != []:
            fdtd.set("bend location", "user specified")
            fdtd.set("bend location x", bend_location[0])
            fdtd.set("bend location y", bend_location[1])
            fdtd.set("bend location z", bend_location[2])
    fdtd.set("mode selection", mode)
    if mode != "user select":
        fdtd.set("selected mode number", str(mode_number))
    fdtd.set("wavelength start", wave_center - wave_width / 2)
    fdtd.set("wavelength stop", wave_center + wave_width / 2)
    fdtd.set("theta", theta)
    print("Created Mode Source (y) at Coordinates {}".format(coord))
    return ()


def gaussian_z(
    name,
    coord,
    height,
    width,
    mode,
    direction,
    wave_center,
    wave_width,
    theta=0,
    bend_radius=None,
    bend_orientation=None,
    bend_location=[],
    mode_number=1,
):
    x = coord[0]
    y = coord[1]
    z = coord[2]
    fdtd.addmode()
    fdtd.set("name", name)
    fdtd.set("direction", direction)
    fdtd.set("injection axis", "z")
    fdtd.set("x", x)
    fdtd.set("y", y)
    fdtd.set("z", z)
    fdtd.set("x span", width)
    fdtd.set("y span", height)
    if bend_radius != None:
        fdtd.set("bent waveguide", 1)
        fdtd.set("bend radius", bend_radius)
        if bend_orientation != None:
            fdtd.set("bend orientation", bend_orientation)
        if bend_location != []:
            fdtd.set("bend location", "user specified")
            fdtd.set("bend location x", bend_location[0])
            fdtd.set("bend location y", bend_location[1])
            fdtd.set("bend location z", bend_location[2])
    fdtd.set("mode selection", mode)
    if mode != "user select":
        fdtd.set("selected mode number", str(mode_number))
    fdtd.set("wavelength start", wave_center - wave_width / 2)
    fdtd.set("wavelength stop", wave_center + wave_width / 2)
    fdtd.set("theta", theta)
    print("Created Mode Source (z) at Coordinates {}".format(coord))
    return ()


def power_xy(coord, height, width, name):
    x = coord[0]
    y = coord[1]
    z = coord[2]
    fdtd.addpower()
    fdtd.set("monitor type", 7)
    fdtd.set("name", name)
    fdtd.set("x", x)
    fdtd.set("y", y)
    fdtd.set("z", z)
    fdtd.set("x span", width)
    fdtd.set("y span", height)
    return ()


def power_xz(coord, height, width, name):
    x = coord[0]
    y = coord[1]
    z = coord[2]
    fdtd.addpower()
    fdtd.set("monitor type", 6)
    fdtd.set("name", name)
    fdtd.set("x", x)
    fdtd.set("y", y)
    fdtd.set("z", z)
    fdtd.set("x span", width)
    fdtd.set("z span", height)
    return ()


def power_yz(coord, height, width, name):
    x = coord[0]
    y = coord[1]
    z = coord[2]
    fdtd.addpower()
    fdtd.set("monitor type", 5)
    fdtd.set("name", name)
    fdtd.set("x", x)
    fdtd.set("y", y)
    fdtd.set("z", z)
    fdtd.set("y span", width)
    fdtd.set("z span", height)
    return ()


def profile_xy(
    coord, height, width, bool, freq_point=None, wave_center=None, wave_width=None
):
    x = coord[0]
    y = coord[1]
    z = coord[2]
    fdtd.addprofile()
    fdtd.set("monitor type", 7)
    fdtd.set("x", x)
    fdtd.set("y", y)
    fdtd.set("z", z)
    fdtd.set("x span", width)
    fdtd.set("y span", height)
    fdtd.set("override global monitor settings", bool)
    if bool == 1 and freq_point != None:
        fdtd.set("frequency points", freq_point)
    if bool == 1 and wave_center != None:
        fdtd.set("use source limits", 0)
        fdtd.set("wavelength center", wave_center * 1e6)
        fdtd.set("wavelength span", wave_width * 1e6)
    print("Created Profile Monitor (xy) at Coordinates {}".format(coord))
    return ()


def profile_xz(
    coord, height, width, bool, freq_point=None, wave_center=None, wave_width=None
):
    x = coord[0]
    y = coord[1]
    z = coord[2]
    fdtd.addprofile()
    fdtd.set("monitor type", 6)
    fdtd.set("x", x)
    fdtd.set("y", y)
    fdtd.set("z", z)
    fdtd.set("x span", width)
    fdtd.set("z span", height)
    fdtd.set("override global monitor settings", bool)
    if bool == 1 and freq_point != None:
        fdtd.set("frequency points", freq_point)
    if bool == 1 and wave_center != None:
        fdtd.set("use source limits", 0)
        fdtd.set("wavelength center", wave_center * 1e6)
        fdtd.set("wavelength span", wave_width * 1e6)
    print("Created Profile Monitor (xz) at Coordinates {}".format(coord))
    return ()


def profile_yz(
    coord, height, width, bool, freq_point=None, wave_center=None, wave_width=None
):
    x = coord[0]
    y = coord[1]
    z = coord[2]
    fdtd.addprofile()
    fdtd.set("monitor type", 5)
    fdtd.set("x", x)
    fdtd.set("y", y)
    fdtd.set("z", z)
    fdtd.set("y span", width)
    fdtd.set("z span", height)
    fdtd.set("override global monitor settings", bool)
    if bool == 1 and freq_point != None:
        fdtd.set("frequency points", freq_point)
    if bool == 1 and wave_center != None:
        fdtd.set("use source limits", 0)
        fdtd.set("wavelength center", wave_center * 1e6)
        fdtd.set("wavelength span", wave_width * 1e6)
    print("Created Profile Monitor (yz) at Coordinates {}".format(coord))


def line_x(coord, width, bool, freq_point=None, wave_center=None, wave_width=None):
    x = coord[0]
    y = coord[1]
    z = coord[2]
    fdtd.addprofile()
    fdtd.set("monitor type", 2)
    fdtd.set("x", x)
    fdtd.set("y", y)
    fdtd.set("z", z)
    fdtd.set("x span", width)
    fdtd.set("override global monitor settings", bool)
    if bool == 1 and freq_point != None:
        fdtd.set("frequency points", freq_point)
    if bool == 1 and wave_center != None:
        fdtd.set("use source limits", 0)
        fdtd.set("wavelength center", wave_center * 1e6)
        fdtd.set("wavelength span", wave_width * 1e6)
    print("Created Line Monitor (x) at Coordinates {}".format(coord))


def point_monitor(coord, bool, freq_point=None, wave_center=None, wave_width=None):
    x = coord[0]
    y = coord[1]
    z = coord[2]
    fdtd.addprofile()
    fdtd.set("monitor type", 1)
    fdtd.set("x", x)
    fdtd.set("y", y)
    fdtd.set("z", z)
    fdtd.set("override global monitor settings", bool)
    if bool == 1 and freq_point != None:
        fdtd.set("frequency points", freq_point)
    if bool == 1 and wave_center != None:
        fdtd.set("use source limits", 0)
        fdtd.set("wavelength center", wave_center * 1e6)
        fdtd.set("wavelength span", wave_width * 1e6)
    print("Created Point Monitor at Coordinates {}".format(coord))


def volume_monitor(
    coord, xspan, yspan, zspan, bool, freq_point=None, wave_center=None, wave_width=None
):
    x = coord[0]
    y = coord[1]
    z = coord[2]
    fdtd.addprofile()
    fdtd.set("monitor type", 8)
    fdtd.set("x", x)
    fdtd.set("y", y)
    fdtd.set("z", z)
    fdtd.set("x span", xspan)
    fdtd.set("y span", yspan)
    fdtd.set("z span", zspan)
    fdtd.set("override global monitor settings", bool)
    if bool == 1 and freq_point != None:
        fdtd.set("frequency points", freq_point)
    if bool == 1 and wave_center != None:
        fdtd.set("use source limits", 0)
        fdtd.set("wavelength center", wave_center * 1e6)
        fdtd.set("wavelength span", wave_width * 1e6)
    print("Created Volume Monitor at Coordinates {}".format(coord))


def movie_xy(coord, height, width):
    x = coord[0]
    y = coord[1]
    z = coord[2]
    fdtd.addmovie()
    fdtd.set("monitor type", 3)
    fdtd.set("x", x)
    fdtd.set("y", y)
    fdtd.set("z", z)
    fdtd.set("x span", width)
    fdtd.set("y span", height)
    print("Created Movie Monitor (xy) at Coordinates {}".format(coord))
    return ()


def movie_xz(coord, height, width):
    x = coord[0]
    y = coord[1]
    z = coord[2]
    fdtd.addmovie()
    fdtd.set("monitor type", 2)
    fdtd.set("x", x)
    fdtd.set("y", y)
    fdtd.set("z", z)
    fdtd.set("x span", width)
    fdtd.set("z span", height)
    print("Created Movie Monitor (xz) at Coordinates {}".format(coord))
    return ()


def movie_yz(coord, height, width):
    x = coord[0]
    y = coord[1]
    z = coord[2]
    fdtd.addmovie()
    fdtd.set("monitor type", 1)
    fdtd.set("x", x)
    fdtd.set("y", y)
    fdtd.set("z", z)
    fdtd.set("y span", width)
    fdtd.set("z span", height)
    print("Created Movie Monitor (yz) at Coordinates {}".format(coord))
    return ()


def time_monitor(coord):
    x = coord[0]
    y = coord[1]
    z = coord[2]
    fdtd.addtime()
    fdtd.set("monitor type", 1)
    fdtd.set("x", x)
    fdtd.set("y", y)
    fdtd.set("z", z)
    print("Created Time (point) Monitor at Coordiantes {}".format(coord))
    return ()


def spectrum_monitor(coord, freq_points):
    x = coord[0]
    y = coord[1]
    z = coord[2]
    fdtd.addpower()
    fdtd.set("monitor type", 1)
    fdtd.set("x", x)
    fdtd.set("y", y)
    fdtd.set("z", z)
    fdtd.set("override global monitor settings", 1)
    fdtd.set("frequency points", freq_points)
    print("Created Power (point) Monitor at Coordiantes {}".format(coord))
    return ()


def modeexpansion_xy_monitor(
    coord,
    xspan,
    yspan,
    mode,
    input,
    bend_radius=None,
    bend_orientation=None,
    bend_location=[],
    mode_number=1,
    name="",
):
    x = coord[0]
    y = coord[1]
    z = coord[2]
    fdtd.addmodeexpansion()
    fdtd.set("monitor type", 3)
    if len(name) > 0:
        fdtd.set("name", name)
    fdtd.set("x", x)
    fdtd.set("y", y)
    fdtd.set("z", z)
    fdtd.set("x span", xspan)
    fdtd.set("y span", yspan)
    fdtd.setexpansion("input", input)
    if bend_radius != None:
        fdtd.set("bent waveguide", 1)
        fdtd.set("bend radius", bend_radius)
        if bend_orientation != None:
            fdtd.set("bend orientation", bend_orientation)
        if bend_location != []:
            fdtd.set("bend location", "user specified")
            fdtd.set("bend location x", bend_location[0])
            fdtd.set("bend location y", bend_location[1])
            fdtd.set("bend location z", bend_location[2])
    fdtd.set("mode selection", mode)
    if mode == "user select":
        fdtd.set("selected mode number", str(mode_number))
    print(f"Mode Expansion xy Set for {mode}, with input {input}")


def modeexpansion_xz_monitor(
    coord,
    xspan,
    zspan,
    mode,
    input,
    bend_radius=None,
    bend_orientation=None,
    bend_location=[],
    mode_number=1,
    name="",
):
    x = coord[0]
    y = coord[1]
    z = coord[2]
    fdtd.addmodeexpansion()
    fdtd.set("monitor type", 2)
    if len(name) > 0:
        fdtd.set("name", name)
    fdtd.set("x", x)
    fdtd.set("y", y)
    fdtd.set("z", z)
    fdtd.set("x span", xspan)
    fdtd.set("z span", zspan)
    fdtd.setexpansion("input", input)
    if bend_radius != None:
        fdtd.set("bent waveguide", 1)
        fdtd.set("bend radius", bend_radius)
        if bend_orientation != None:
            fdtd.set("bend orientation", bend_orientation)
        if bend_location != []:
            fdtd.set("bend location", "user specified")
            fdtd.set("bend location x", bend_location[0])
            fdtd.set("bend location y", bend_location[1])
            fdtd.set("bend location z", bend_location[2])
    fdtd.set("mode selection", mode)
    if mode == "user select":
        fdtd.set("selected mode number", str(mode_number))
    print(f"Mode Expansion xz Set for {mode}, with input {input}")


def modeexpansion_yz_monitor(
    coord,
    xspan,
    zspan,
    mode,
    input,
    bend_radius=None,
    bend_orientation=None,
    bend_location=[],
    mode_number=1,
    name="",
):
    x = coord[0]
    y = coord[1]
    z = coord[2]
    fdtd.addmodeexpansion()
    fdtd.set("monitor type", 1)
    if len(name) > 0:
        fdtd.set("name", name)
    fdtd.set("x", x)
    fdtd.set("y", y)
    fdtd.set("z", z)
    fdtd.set("y span", xspan)
    fdtd.set("z span", zspan)
    fdtd.setexpansion("input", input)
    if bend_radius != None:
        fdtd.set("bent waveguide", 1)
        fdtd.set("bend radius", bend_radius)
        if bend_orientation != None:
            fdtd.set("bend orientation", bend_orientation)
        if bend_location != []:
            fdtd.set("bend location", "user specified")
            fdtd.set("bend location x", bend_location[0])
            fdtd.set("bend location y", bend_location[1])
            fdtd.set("bend location z", bend_location[2])
    fdtd.set("mode selection", mode)
    if mode == "user select":
        fdtd.set("selected mode number", str(mode_number))
    print(f"Mode Expansion yz Set for {mode}, with input {input}")


def qanalysis(coord, height, width, wave_min, wave_max):
    x = coord[0]
    y = coord[1]
    z = coord[2]

    c = 299792458
    fmin = c / wave_max
    fmax = c / wave_min

    fdtd.addobject("Qanalysis")
    fdtd.set("x", x)
    fdtd.set("y", y)

    fdtd.set("x span", width)
    fdtd.set("y span", height)
    fdtd.set("f min", fmin)
    fdtd.set("f max", fmax)

    print("Created QAnalysis Group Object")
    return ()


def Lowqanalysis(coord, height, width):
    x = coord[0]
    y = coord[1]
    z = coord[2]

    fdtd.addobject("lowQanalysis")
    fdtd.set("x", x)
    fdtd.set("y", y)

    fdtd.set("x span", width)
    fdtd.set("y span", height)

    print("Created LowQAnalysis Group Object")
    return ()


def profile_save(filename, objname, resultname):
    """For monitors and mode sources"""
    data_file = filename + f"_{objname}_{resultname}.pkl"
    keys = filename + f"_{objname}_{resultname}_keys.txt"

    data = fdtd.getresult(objname, resultname)

    for k, v in data.items():
        if isinstance(v, np.ndarray):
            data[k] = v.tolist()
        else:
            # print(type(v))
            pass

    data_keys = list(data.keys())
    g = open(keys, "w")
    g.write(f"Keys for file: {data_file}\n")
    for i in data_keys:
        g.write(f"{i}\n")
    g.close()

    f = open(data_file, "wb")
    pickle.dump(data, f)
    f.close()


def profile_save_single(filename, objname, resultname, slice_num):
    data_file = filename + f"_{objname}_{resultname}_{slice_num}.pkl"
    keys = filename + f"_{objname}_{resultname}_keys.txt"

    data = fdtd.getresult(objname, resultname)

    # val = np.transpose(np.squeeze(val), (2,3,1,0))

    for k, v in data.items():
        if isinstance(v, np.ndarray):
            if k == resultname:
                v = np.transpose(np.squeeze(v), (2, 3, 1, 0))
                data[k] = v[slice_num].tolist()
            elif k == "lambda":
                data[k] = v[slice_num].tolist()
            else:
                data[k] = v.tolist()
        else:
            # print(type(v))
            pass

    data_keys = list(data.keys())
    g = open(keys, "w")
    g.write(f"Keys for file: {data_file}\n")
    for i in data_keys:
        g.write(f"{i}\n")
    g.close()

    f = open(data_file, "wb")
    pickle.dump(data, f)
    f.close()


def profile_rawData_save(filename, objname):
    data_file = filename + f"_{objname}_rawData.pkl"
    keys = filename + f"_{objname}_rawData_keys.txt"

    x = fdtd.getdata(objname, "x")
    y = fdtd.getdata(objname, "y")
    z = fdtd.getdata(objname, "z")
    f = fdtd.getdata(objname, "f")
    Ex = fdtd.getdata(objname, "Ex")
    Ey = fdtd.getdata(objname, "Ey")
    Ez = fdtd.getdata(objname, "Ez")
    Hx = fdtd.getdata(objname, "Hx")
    Hy = fdtd.getdata(objname, "Hy")
    Hz = fdtd.getdata(objname, "Hz")

    g = open(keys, "w")
    g.write(f"Keys for file: {data_file}\n")
    g.write("x\n")
    g.write("y\n")
    g.write("z\n")
    g.write("f\n")
    g.write("Ex\n")
    g.write("Ey\n")
    g.write("Ez\n")
    g.write("Hx\n")
    g.write("Hy\n")
    g.write("Hz\n")
    g.close()

    h = open(data_file, "wb")
    data_dict = {
        "x": x,
        "y": y,
        "z": z,
        "f": f,
        "Ex": Ex,
        "Ey": Ey,
        "Ez": Ez,
        "Hx": Hx,
        "Hy": Hy,
        "Hz": Hz,
    }
    pickle.dump(data_dict, h)
    h.close()


def spectrum_save(filename, objname, resultname):
    data_file = filename + f"_{objname}_{resultname}.pkl"
    keys = filename + f"_{objname}_{resultname}_keys.txt"

    data = fdtd.getresult(objname, resultname)
    for k, v in data.items():
        if isinstance(v, np.ndarray):
            data[k] = v.tolist()
        else:
            # print(type(v))
            pass

    data_keys = list(data.keys())
    g = open(keys, "w")
    g.write(f"Keys for file: {data_file}\n")
    for i in data_keys:
        g.write(f"{i}\n")
    g.close()

    f = open(data_file, "wb")
    pickle.dump(data, f)
    f.close()


def load_file(filename):
    fdtd.load(filename)
    print(f"Loaded Sim File: {filename}")


def save_file(filename):
    fdtd.save(filename)
    print(f"Saved Sim File: {filename}")


def QAnalysis_run():
    try:
        fdtd.runanalysis("Qanalysis")
        time.sleep(3)
        print("QAnalysis Run")
        ret = 1
    except Exception as e:
        print("QAnalysis Fail")
        ret = 0
    return ret


def LowQAnalysis_run():
    try:
        fdtd.runanalysis("lowQanalysis")
        time.sleep(3)
        print("Low QAnalysis Run")
        ret = 1
    except Exception as e:
        print("Low QAnalysis Fail")
        ret = 0
    return ret


def farfield_run_save(monitor_name, filename, wave_center, wave_width, freq_point=None):
    if freq_point == None:
        freq_point = 5
    wave_list = np.linspace(
        wave_center - wave_width / 2, wave_center + wave_width / 2, freq_point
    )
    ux = fdtd.farfieldux(monitor_name, 1)
    uy = fdtd.farfielduy(monitor_name, 1)
    farfield_slices = []
    for i in range(freq_point):
        fdtd.farfieldfilter(1)
        E = fdtd.farfield3d(monitor_name, i + 1)
        farfield_slices.append(E.tolist())
    results_dict = {
        "wavelength": wave_list.tolist(),
        "farfield_E": farfield_slices,
        "ux": np.squeeze(ux).tolist(),
        "uy": np.squeeze(uy).tolist(),
    }
    saveJson(results_dict, filename + f"_{monitor_name}_farfield3D")


def modeexpansion_run_save(monitor_name, filename):
    data_file = filename + f"_{monitor_name}_modeExpansion.pkl"
    keys = filename + f"_{monitor_name}_modeExpansion_keys.txt"

    data = fdtd.getresult(monitor_name, "expansion for input")
    for k, v in data.items():
        if isinstance(v, np.ndarray):
            data[k] = v.tolist()
        else:
            # print(type(v))
            pass

    data_keys = list(data.keys())
    g = open(keys, "w")
    g.write(f"Keys for file: {data_file}\n")
    for i in data_keys:
        g.write(f"{i}\n")
    g.close()

    f = open(data_file, "wb")
    pickle.dump(data, f)
    f.close()


def QAnalysis_save(run_bool, filename):
    if run_bool == 1:
        data = fdtd.getresult("Qanalysis", "Q")
        wave = np.squeeze(data["lambda"]).tolist()
        quality = data["Q"].tolist()
        results_dict = {"wavelength": wave, "Q": quality}
        saveJson(results_dict, filename + "_QAnalysis")


def LowQAnalysis_save(run_bool, filename):
    if run_bool == 1:
        data = fdtd.getresult("lowQanalysis", "Q")
        wave = np.squeeze(data["lambda"]).tolist()
        quality = data["Q"].tolist()

        results_dict = {"wavelength": wave, "Q": quality}
        saveJson(results_dict, filename + "_LowQAnalysis")


def saveJson(dictionary, save_dir):
    json_object = json.dumps(dictionary, indent=4)
    with open(save_dir + ".json", "w") as outfile:
        outfile.write(json_object)


def CheckMakeDirectory(path, dir):
    full = os.path.join(path, dir)
    if os.path.exists(full) == False:
        os.mkdir(full)
    else:
        print("Folder Exists")
    return full


def createSubFolders(path: str):
    """Create organizational subfolders"""
    simpath = CheckMakeDirectory(path, "sim/")
    savepath = CheckMakeDirectory(path, "results/")
    datapath = CheckMakeDirectory(path, "data/")
    paramspath = CheckMakeDirectory(path, "params/")

    return (simpath, savepath, datapath, paramspath)
