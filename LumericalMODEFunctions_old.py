import imp

imp.load_source("lumapi", "C:/Program Files/Lumerical/v211/api/python/lumapi.py")
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

mode = lumapi.MODE()


def sim_yNormal(coord, width, height, cell_num=[], mesh_res=[], bc_type="Metal"):
    x = coord[0]
    y = coord[1]
    z = coord[2]
    bc_dict = {"PML": 1, "Metal": 2, "Periodic": 3}
    mode.addfde()
    mode.set("solver type", "2D Y normal")
    if len(cell_num) != 0:
        mode.set("mesh cells x", cell_num[0])
        mode.set("mesh cells z", cell_num[1])
    if len(mesh_res) != 0:
        mode.set("define x mesh by", "maximum mesh step")
        mode.set("define z mesh by", "maximum mesh step")
        mode.set("dx", mesh_res[0])
        mode.set("dz", mesh_res[1])
    mode.set("x", x)
    mode.set("y", y)
    mode.set("z", z)
    mode.set("x span", width)
    mode.set("z span", height)
    mode.set("x min bc", bc_dict[bc_type])
    mode.set("x max bc", bc_dict[bc_type])
    mode.set("z min bc", bc_dict[bc_type])
    mode.set("z max bc", bc_dict[bc_type])
    print("Simulation Region Y-Normal Added")
    return ()


def sim_xNormal(coord, width, height, cell_num=[], mesh_res=[], bc_type="Metal"):
    x = coord[0]
    y = coord[1]
    z = coord[2]
    bc_dict = {"PML": 1, "Metal": 2, "Periodic": 3}
    mode.addfde()
    mode.set("solver type", "2D X normal")
    if len(cell_num) != 0:
        mode.set("mesh cells y", cell_num[0])
        mode.set("mesh cells z", cell_num[1])
    if len(mesh_res) != 0:
        mode.set("define y mesh by", "maximum mesh step")
        mode.set("define z mesh by", "maximum mesh step")
        mode.set("dy", mesh_res[0])
        mode.set("dz", mesh_res[1])
    mode.set("x", x)
    mode.set("y", y)
    mode.set("z", z)
    mode.set("y span", width)
    mode.set("z span", height)
    mode.set("y min bc", bc_dict[bc_type])
    mode.set("y max bc", bc_dict[bc_type])
    mode.set("z min bc", bc_dict[bc_type])
    mode.set("z max bc", bc_dict[bc_type])
    print("Simulation Region X-Normal Added")
    return ()


def sim_zNormal(coord, width, height, cell_num=[], mesh_res=[], bc_type="Metal"):
    x = coord[0]
    y = coord[1]
    z = coord[2]
    bc_dict = {"PML": 1, "Metal": 2, "Periodic": 3}
    mode.addfde()
    mode.set("solver type", "2D Z normal")
    if len(cell_num) != 0:
        mode.set("mesh cells x", cell_num[0])
        mode.set("mesh cells y", cell_num[1])
    if len(mesh_res) != 0:
        mode.set("define x mesh by", "maximum mesh step")
        mode.set("define y mesh by", "maximum mesh step")
        mode.set("dx", mesh_res[0])
        mode.set("dy", mesh_res[1])
    mode.set("x", x)
    mode.set("y", y)
    mode.set("z", z)
    mode.set("x span", width)
    mode.set("y span", height)
    mode.set("x min bc", bc_dict[bc_type])
    mode.set("x max bc", bc_dict[bc_type])
    mode.set("y min bc", bc_dict[bc_type])
    mode.set("y max bc", bc_dict[bc_type])
    print("Simulation Region Z-Normal Added")
    return ()


def eigenSolver(wave, num_modes):
    mode.select("FDE")
    mode.set("wavelength", wave)
    mode.set("number of trial modes", num_modes)
    print("Solver Settings Set")


def eigenSolverBentWaveguide(bend_raduis, bend_orientation, bend_location=[]):
    mode.select("FDE")
    mode.set("bent waveguide", 1)
    mode.set("bend radius", bend_raduis)
    mode.set("bend orientation", bend_orientation)
    if len(bend_location) != 0:
        mode.set("bend location", "user specified")
        mode.set("bend location (x)", bend_location[0])
        mode.set("bend location (y)", bend_location[1])
        mode.set("bend location (z)", bend_location[2])
    print("Solver Settings Set for Bent Waveguide")


def create_diamond():
    mode.eval("temp = addmaterial('Dielectric');")
    mode.eval("setmaterial(temp,'name','diamond');")
    mode.eval("setmaterial('diamond','Permittivity', 2.4114 * 2.4114);")
    mode.eval("setmaterial('diamond','color',[0;0.6;0.8;1]);")
    print("Created Diamond Material Entry")
    return "diamond"


def create_sio2():
    print("Assigned SiO2 Material Entry")
    return "SiO2 (Glass) - Palik"


def create_si():
    print("Assigned Si Material Entry")
    return "Si (Silicon) - Palik"


def create_sic():
    mode.eval("temp = addmaterial('Dielectric');")
    mode.eval("setmaterial(temp,'name','sic');")
    mode.eval("setmaterial('sic','Permittivity', 2.64 * 2.64);")
    mode.eval("setmaterial('sic','color',[0;0.6;1;0.5]);")
    print("Created SiC Material Entry")
    return "sic"


def create_sin():
    print("Assigned SiN Material Entry")
    return "Si3N4 (Silicon Nitride) - Phillip"


def create_etch():
    print("Assigned Etch Material Entry")
    return "etch"


def rectWaveguide(name, coord, material, thick, xspan, yspan):
    x = coord[0]
    y = coord[1]
    z = coord[2]
    mode.addrect()
    mode.set("name", name)
    mode.set("material", material)
    mode.set("x", x)
    mode.set("y", y)
    mode.set("z", z)
    mode.set("x span", xspan)
    mode.set("y span", yspan)
    mode.set("z span", thick)
    print("Created Rectangular Waveguide with Material {}".format(material))
    return ()


### PARAMETER FUNCTIONS START HERE ###


def sim_yNormal_params(
    params_dict, coord, width, heigth, cell_num=[], mesh_res=[], bc_type="Metal"
):
    component = "sim_yNormal"
    dict = {}
    dict[f"{component}_coord"] = coord
    dict[f"{component}_width"] = width
    dict[f"{component}_height"] = heigth
    dict[f"{component}_cell_num"] = cell_num
    dict[f"{component}_mesh_res"] = mesh_res
    dict[f"{component}_bc_type"] = bc_type

    params_dict[f"{component}"] = dict
    return params_dict


def sim_xNormal_params(
    params_dict, coord, width, heigth, cell_num=[], mesh_res=[], bc_type="Metal"
):
    component = "sim_xNormal"
    dict = {}
    dict[f"{component}_coord"] = coord
    dict[f"{component}_width"] = width
    dict[f"{component}_height"] = heigth
    dict[f"{component}_cell_num"] = cell_num
    dict[f"{component}_mesh_res"] = mesh_res
    dict[f"{component}_bc_type"] = bc_type

    params_dict[f"{component}"] = dict
    return params_dict


def sim_zNormal_params(
    params_dict, coord, width, heigth, cell_num=[], mesh_res=[], bc_type="Metal"
):
    component = "sim_zNormal"
    dict = {}
    dict[f"{component}_coord"] = coord
    dict[f"{component}_width"] = width
    dict[f"{component}_height"] = heigth
    dict[f"{component}_cell_num"] = cell_num
    dict[f"{component}_mesh_res"] = mesh_res
    dict[f"{component}_bc_type"] = bc_type

    params_dict[f"{component}"] = dict
    return params_dict


def rectWaveguide_params(params_dict, count, coord, material, thick, xspan, yspan):
    component = "rectWaveguide"
    dict = {}
    dict[f"{component}_coord_{count}"] = coord
    dict[f"{component}_material_{count}"] = material
    dict[f"{component}_thick_{count}"] = thick
    dict[f"{component}_xspan_{count}"] = xspan
    dict[f"{component}_yspan_{count}"] = yspan
    params_dict[f"{component}_{count}"] = dict
    return params_dict


### DATA EXTRACTION AND PROCESSING FUNCTIONS ###


def saveModes(filename, resultname, mode_num=[]):
    if len(mode_num) == 0:
        mode_num = range(int(mode.nummodes()))
    for i in mode_num:
        mode.select(f"FDE::data::mode{i+1}")
        data = mode.getresult(f"mode{i+1}", resultname)

        data_file = filename + f"_mode{i}_{resultname}.pkl"
        keys = filename + f"_mode{i}_{resultname}_keys.txt"

        for k, v in data.items():
            if isinstance(v, np.ndarray):
                data[k] = v.tolist()
            else:
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


def saveData(filename, mode_num=[]):
    if len(mode_num) == 0:
        mode_num = range(int(mode.nummodes()))

    index_real_array = []
    index_img_array = []
    loss_array = []
    TE_array = []
    num_list = []

    for i in mode_num:
        mode.select(f"FDE::data::mode{i+1}")

        index = mode.getdata(f"mode{i+1}", "neff")[0]
        loss = mode.getdata(f"mode{i+1}", "loss")
        TE = mode.getdata(f"mode{i+1}", "TE polarization fraction")

        index_real_array.append(index.real[0])
        index_img_array.append(index.imag[0])
        loss_array.append(loss)
        TE_array.append(TE)
        num_list.append(i)

    data_dict = {
        "mode_num": num_list,
        "effective_index_real": index_real_array,
        "effective_index,img": index_img_array,
        "loss": loss_array,
        "TE_polarization": TE_array,
    }

    saveJson(data_dict, filename + "_data")


def mode_E_plot_xz(filename, savename, modename, resultname, display):
    profile_file = filename + f"_{modename}_{resultname}.pkl"
    data_file = open(profile_file, "rb")

    data = pickle.load(data_file)
    data_keys = list(data.keys())
    print(data_keys)

    # first extact dimensions
    x = np.squeeze(data["x"])
    x_min = min(x) * 1e6
    x_max = max(x) * 1e6

    z = np.squeeze(data["z"])
    z_min = min(z) * 1e6
    z_max = max(z) * 1e6

    # extract mode profile E
    E = np.asarray(data["E"])
    E = np.transpose(np.squeeze(E))
    E_array = np.sqrt(abs(E[0]) ** 2 + abs(E[1]) ** 2 + abs(E[2]) ** 2)

    plt.imshow(E_array, origin="lower", extent=[x_min, x_max, z_min, z_max])
    plt.colorbar()
    plt.ylabel("um")
    plt.xlabel("um")
    plt.savefig(savename + f"_{modename}_E.png")
    if display == 1:
        plt.show()
    plt.clf()

    # # extract mode profile H
    # H = np.asarray(data['H'])
    # H = np.transpose(np.squeeze(H))
    # H_array = np.sqrt(abs(H[0])**2+abs(H[1])**2+abs(H[2])**2)

    # plt.imshow(H_array, origin='lower', extent=[x_min, x_max, y_min, y_max])
    # plt.colorbar()
    # plt.ylabel('um')
    # plt.xlabel('um')
    # plt.savefig(filename+f'_{modename}_H.svg')
    # if display == 1:
    #     plt.show()

    # # extract effective index profile
    # index = np.asarray(data['index'])
    # index = np.transpose(np.squeeze(index))
    # index_array = np.sqrt(index[0]**2+index[1]**2+index[2]**2)

    # plt.imshow(index_array, origin='lower', extent=[x_min, x_max, y_min, y_max])
    # plt.ylabel('um')
    # plt.xlable('um')
    # plt.colorbar()
    # plt.savefig(filename+f'_{modename}_Idx.svg')
    # if display == 1:
    #     plt.show()


#### FILE IO FUNCTIONS START HERE ###


def load_file(filename):
    mode.load(filename)
    print(f"Loaded Sim File: {filename}")


def save_file(filename):
    mode.save(filename)
    print(f"Saved Sim File: {filename}")


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
