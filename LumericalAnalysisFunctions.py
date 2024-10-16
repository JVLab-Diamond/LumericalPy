import matplotlib.pyplot as plt
from skimage.feature import peak_local_max
import numpy as np
import pickle
import json
import os

##################################################################################


def mode_neff_extract(filename, objname, resultname):
    profile_file = filename + f"_{objname}_{resultname}.pkl"
    data_file = open(profile_file, "rb")

    data = pickle.load(data_file)
    data_keys = list(data.keys())
    print(data_keys)

    neff = data["neff"]
    return neff


def mode_profile_plot_xy(filename, savename, objname, resultname, display):
    profile_file = filename + f"_{objname}_{resultname}.pkl"
    data_file = open(profile_file, "rb")

    data = pickle.load(data_file)
    data_keys = list(data.keys())
    print(data_keys)

    # first extact dimensions
    x = np.squeeze(data["x"])
    x_min = min(x) * 1e6
    x_max = max(x) * 1e6

    y = np.squeeze(data["y"])
    y_min = min(y) * 1e6
    y_max = max(y) * 1e6

    # extract mode profile E
    E = np.asarray(data["E"])
    E = np.transpose(np.squeeze(E))
    E_array = np.sqrt(abs(E[0]) ** 2 + abs(E[1]) ** 2 + abs(E[2]) ** 2)

    plt.imshow(E_array, origin="lower", extent=[x_min, x_max, y_min, y_max])
    plt.colorbar()
    plt.ylabel("um")
    plt.xlabel("um")
    plt.savefig(savename + f"_{objname}_E.svg")
    if display == 1:
        plt.show()

    # extract mode profile H
    H = np.asarray(data["H"])
    H = np.transpose(np.squeeze(H))
    H_array = np.sqrt(abs(H[0]) ** 2 + abs(H[1]) ** 2 + abs(H[2]) ** 2)

    plt.imshow(H_array, origin="lower", extent=[x_min, x_max, y_min, y_max])
    plt.colorbar()
    plt.ylabel("um")
    plt.xlabel("um")
    plt.savefig(savename + f"_{objname}_H.svg")
    if display == 1:
        plt.show()

    # extract effective index profile
    index = np.asarray(data["index"])
    index = np.transpose(np.squeeze(index))
    index_array = np.sqrt(index[0] ** 2 + index[1] ** 2 + index[2] ** 2)

    plt.imshow(index_array, origin="lower", extent=[x_min, x_max, y_min, y_max])
    plt.ylabel("um")
    plt.xlable("um")
    plt.colorbar()
    plt.savefig(savename + f"_{objname}_Idx.svg")
    if display == 1:
        plt.show()


def mode_profile_plot_xz(filename, objname, resultname, display):
    profile_file = filename + f"_{objname}_{resultname}.pkl"
    data_file = open(profile_file, "rb")

    data = pickle.load(data_file)
    data_keys = list(data.keys())
    print(data_keys)

    # first extact dimensions
    x = np.squeeze(data["x"])
    x_min = min(x) * 1e6
    x_max = max(x) * 1e6

    y = np.squeeze(data["z"])
    y_min = min(y) * 1e6
    y_max = max(y) * 1e6

    # extract mode profile E
    E = np.asarray(data["E"])
    E = np.transpose(np.squeeze(E))
    E_array = np.sqrt(abs(E[0]) ** 2 + abs(E[1]) ** 2 + abs(E[2]) ** 2)

    plt.imshow(E_array, origin="lower", extent=[x_min, x_max, y_min, y_max])
    plt.colorbar()
    plt.ylabel("um")
    plt.xlabel("um")
    plt.savefig(filename + f"_{objname}_E.svg")
    if display == 1:
        plt.show()

    # extract mode profile H
    H = np.asarray(data["H"])
    H = np.transpose(np.squeeze(H))
    H_array = np.sqrt(abs(H[0]) ** 2 + abs(H[1]) ** 2 + abs(H[2]) ** 2)

    plt.imshow(H_array, origin="lower", extent=[x_min, x_max, y_min, y_max])
    plt.colorbar()
    plt.ylabel("um")
    plt.xlabel("um")
    plt.savefig(filename + f"_{objname}_H.svg")
    if display == 1:
        plt.show()

    # extract effective index profile
    index = np.asarray(data["index"])
    index = np.transpose(np.squeeze(index))
    index_array = np.sqrt(index[0] ** 2 + index[1] ** 2 + index[2] ** 2)

    plt.imshow(index_array, origin="lower", extent=[x_min, x_max, y_min, y_max])
    plt.ylabel("um")
    plt.xlabel("um")
    plt.colorbar()
    plt.savefig(filename + f"_{objname}_Idx.svg")
    if display == 1:
        plt.show()


def mode_profile_plot_yz(filename, objname, resultname, display):
    profile_file = filename + f"_{objname}_{resultname}.pkl"
    data_file = open(profile_file, "rb")

    data = pickle.load(data_file)
    data_keys = list(data.keys())
    print(data_keys)

    # first extact dimensions
    x = np.squeeze(data["y"])
    x_min = min(x) * 1e6
    x_max = max(x) * 1e6

    y = np.squeeze(data["z"])
    y_min = min(y) * 1e6
    y_max = max(y) * 1e6

    # extract mode profile E
    E = np.asarray(data["E"])
    E = np.transpose(np.squeeze(E))
    E_array = np.sqrt(abs(E[0]) ** 2 + abs(E[1]) ** 2 + abs(E[2]) ** 2)

    plt.imshow(E_array, origin="lower", extent=[x_min, x_max, y_min, y_max])
    plt.colorbar()
    plt.ylabel("um")
    plt.xlabel("um")
    plt.title("E")
    plt.savefig(filename + f"_{objname}_E.svg")
    if display == 1:
        plt.show()

    # extract mode profile H
    H = np.asarray(data["H"])
    H = np.transpose(np.squeeze(H))
    H_array = np.sqrt(abs(H[0]) ** 2 + abs(H[1]) ** 2 + abs(H[2]) ** 2)

    plt.imshow(H_array, origin="lower", extent=[x_min, x_max, y_min, y_max])
    plt.colorbar()
    plt.ylabel("um")
    plt.xlabel("um")
    plt.title("H")
    plt.savefig(filename + f"_{objname}_H.svg")
    if display == 1:
        plt.show()

    # extract effective index profile
    index = np.asarray(data["index"])
    index = np.transpose(np.squeeze(index))
    index_array = np.sqrt(index[0] ** 2 + index[1] ** 2 + index[2] ** 2)

    plt.imshow(index_array, origin="lower", extent=[x_min, x_max, y_min, y_max])
    plt.ylabel("um")
    plt.xlabel("um")
    plt.title("Index")
    plt.colorbar()
    plt.savefig(filename + f"_{objname}_Idx.svg")
    if display == 1:
        plt.show()

    plt.clf()


def profile_unpack(filename, objname, resultname):
    profile_file = filename + f"_{objname}_{resultname}.pkl"
    data_file = open(profile_file, "rb")

    data = pickle.load(data_file)
    data_keys = list(data.keys())
    print(data_keys)

    # first extact dimensions
    x = np.squeeze(data["y"])

    y = np.squeeze(data["z"])

    # extract mode profile E
    E = np.asarray(data["E"])
    E = np.transpose(np.squeeze(E))
    E_array = np.sqrt(abs(E[0]) ** 2 + abs(E[1]) ** 2 + abs(E[2]) ** 2)

    # extract mode profile H
    H = np.asarray(data["H"])
    H = np.transpose(np.squeeze(H))
    H_array = np.sqrt(abs(H[0]) ** 2 + abs(H[1]) ** 2 + abs(H[2]) ** 2)

    # extract effective index profile
    index = np.asarray(data["index"])
    index = np.transpose(np.squeeze(index))
    index_array = np.sqrt(index[0] ** 2 + index[1] ** 2 + index[2] ** 2)

    data_dict = {
        "x": x.tolist(),
        "y": y.tolist(),
        "E": E_array.tolist(),
        "H": H_array.tolist(),
        "index": index_array.tolist(),
    }
    return data_dict


def monitor_profile_plot_xy(filename, objname, resultname, savename, display):
    profile_file = filename + f"_{objname}_{resultname}.pkl"
    data_file = open(profile_file, "rb")

    data = pickle.load(data_file)
    data_keys = list(data.keys())
    print(data_keys)

    # first extact dimensions
    x = np.squeeze(data["x"])
    x_min = min(x) * 1e6
    x_max = max(x) * 1e6

    y = np.squeeze(data["y"])
    y_min = min(y) * 1e6
    y_max = max(y) * 1e6

    wave = np.squeeze(data["lambda"])

    # extract mode profile E
    val = np.asarray(data[resultname])
    print(val.shape)
    val = np.transpose(np.squeeze(val), (2, 3, 1, 0))
    print(val.shape)

    for i, j in zip(wave, val):
        val_array = np.sqrt(abs(j[0]) ** 2 + abs(j[1]) ** 2 + abs(j[2]) ** 2)
        plt.imshow(val_array, origin="lower", extent=[x_min, x_max, y_min, y_max])
        plt.colorbar()
        plt.title(f"{resultname} {int(i*1e9)}nm")
        plt.ylabel("um")
        plt.xlabel("um")
        plt.savefig(savename + f"_{objname}_{resultname}_{int(i*1e9)}.png")
        if display == 1:
            plt.show()
        plt.close()

    return (wave, val)


def monitor_profile_plot_xy_single(
    filename, objname, resultname, savename, display, slice_num
):
    profile_file = filename + f"_{objname}_{resultname}_{slice_num}.pkl"
    data_file = open(profile_file, "rb")

    data = pickle.load(data_file)
    data_keys = list(data.keys())
    print(data_keys)

    # first extact dimensions
    x = np.squeeze(data["x"])
    x_min = min(x) * 1e6
    x_max = max(x) * 1e6

    y = np.squeeze(data["y"])
    y_min = min(y) * 1e6
    y_max = max(y) * 1e6

    wave = np.squeeze(data["lambda"])

    # extract mode profile E
    val = np.asarray(data[resultname])

    val_array = np.sqrt(abs(val[0]) ** 2 + abs(val[1]) ** 2 + abs(val[2]) ** 2)
    plt.imshow(val_array, origin="lower", extent=[x_min, x_max, y_min, y_max])
    plt.colorbar()
    plt.title(f"{resultname} {int(wave*1e9)} nm")
    plt.ylabel("um")
    plt.xlabel("um")
    plt.savefig(savename + f"_{objname}_{resultname}_{slice_num}.png")
    if display == 1:
        plt.show()
    plt.close()

    return (wave, val_array)


def profile_lobe_counting(array, min_distance, threshold_factor):
    detected_peaks = peak_local_max(
        array,
        min_distance=min_distance,
        threshold_abs=threshold_factor * np.mean(array),
    )
    peak_y, peak_x = detected_peaks.T

    # plt.imshow(array, origin='lower')
    # plt.scatter(peak_x, peak_y)
    # plt.colorbar()
    # plt.ylabel('um')
    # plt.xlabel('um')
    # plt.show()

    return (peak_x, peak_y, len(peak_x))


def monitor_profile_plot_yz(filename, objname, resultname, savename, display):
    profile_file = filename + f"_{objname}_{resultname}.pkl"
    data_file = open(profile_file, "rb")

    data = pickle.load(data_file)
    data_keys = list(data.keys())
    print(data_keys)

    # first extact dimensions
    x = np.squeeze(data["y"])
    x_min = min(x) * 1e6
    x_max = max(x) * 1e6

    y = np.squeeze(data["z"])
    y_min = min(y) * 1e6
    y_max = max(y) * 1e6

    wave = np.squeeze(data["lambda"])

    # extract mode profile E
    val = np.asarray(data[resultname])
    val = np.transpose(np.squeeze(val), (2, 3, 1, 0))

    for i, j in zip(wave, val):
        val_array = np.sqrt(abs(j[0]) ** 2 + abs(j[1]) ** 2 + abs(j[2]) ** 2)
        plt.imshow(val_array, origin="lower", extent=[x_min, x_max, y_min, y_max])
        plt.colorbar()
        plt.title(f"{resultname} {int(i*1e9)}nm")
        plt.ylabel("um")
        plt.xlabel("um")
        plt.savefig(savename + f"_{objname}_{resultname}_{int(i*1e9)}.png")
        if display == 1:
            plt.show()
        plt.close()

    return (wave, val)


def monitor_profile_plot_xz(filename, objname, resultname, savename, display):
    profile_file = filename + f"_{objname}_{resultname}.pkl"
    data_file = open(profile_file, "rb")

    data = pickle.load(data_file)
    data_keys = list(data.keys())
    print(data_keys)

    # first extact dimensions
    x = np.squeeze(data["x"])
    x_min = min(x) * 1e6
    x_max = max(x) * 1e6

    y = np.squeeze(data["z"])
    y_min = min(y) * 1e6
    y_max = max(y) * 1e6

    wave = np.squeeze(data["lambda"])

    # extract mode profile E
    val = np.asarray(data[resultname])
    val = np.transpose(np.squeeze(val), (2, 3, 1, 0))

    for i, j in zip(wave, val):
        val_array = np.sqrt(abs(j[0]) ** 2 + abs(j[1]) ** 2 + abs(j[2]) ** 2)
        plt.imshow(
            val_array,
            origin="lower",
            extent=[x_min, x_max, y_min, y_max],
            aspect=8,
            cmap="inferno",
        )
        # plt.colorbar()
        # plt.title(f'{resultname} {int(i*1e9)}nm')
        plt.ylabel("um")
        plt.xlabel("um")
        plt.savefig(
            savename + f"_{objname}_{resultname}_{int(i*1e9)}.png",
            dpi=300,
            transparent=True,
            bbox_inches="tight",
        )
        if display == 1:
            plt.show()
        plt.close()

    return (wave, val)


def monitor_plot_1D(filename, objname, resultname, savename, display):
    profile_file = filename + f"_{objname}_{resultname}.pkl"
    data_file = open(profile_file, "rb")

    data = pickle.load(data_file)
    data_keys = list(data.keys())
    print(data_keys)

    wave = np.squeeze(data["lambda"])

    # extract mode profile E
    val = np.asarray(data[resultname])
    val = np.squeeze(val).T

    val_array = np.sqrt(np.abs(val[0]) ** 2 + np.abs(val[1]) ** 2 + np.abs(val[2]) ** 2)
    plt.plot(wave, val_array)
    # plt.colorbar()
    # plt.title(f'{resultname} {int(i*1e9)}nm')
    plt.ylabel("abu")
    plt.xlabel("Wavelength")
    plt.savefig(
        savename + f"_{objname}_{resultname}.png",
        dpi=300,
        transparent=True,
        bbox_inches="tight",
    )
    if display == 1:
        plt.show()
    plt.close()

    return (wave, val_array)


def monitor_transmission_plot(filename, objname, savename, display):
    profile_file = filename + f"_{objname}_T.pkl"
    data_file = open(profile_file, "rb")

    data = pickle.load(data_file)
    data_keys = list(data.keys())
    print(data_keys)

    wave = np.squeeze(data["lambda"])
    transmission = np.squeeze(data["T"])

    plt.plot(wave * 1e9, transmission, marker="o")
    plt.ylabel("Transmission")
    plt.xlabel("Wavelength (nm)")
    plt.savefig(savename + f"_{objname}_T.png")
    if display == 1:
        plt.show()
    plt.close()

    data_dict = {"wavelength": wave.tolist(), "transmission": transmission.tolist()}

    saveJson(data_dict, savename + f"_{objname}_T")

    return (wave, transmission)


def monitor_spectrum_plot(filename, objname, savename, propname, display):
    profile_file = filename + f"_{objname}_{propname}.pkl"
    data_file = open(profile_file, "rb")

    data = pickle.load(data_file)
    data_keys = list(data.keys())
    print(data_keys)

    wave = np.squeeze(data["lambda"])
    spectrum = np.squeeze(data[propname])

    plt.plot(wave * 1e9, abs(spectrum))
    plt.ylabel(propname)
    plt.xlabel("Wavelength (nm)")
    plt.savefig(savename + f"_{objname}_{propname}.png")
    if display == 1:
        plt.show()
    plt.close()

    return (wave, spectrum)


def farfield3d_plot(filename, objname, save_name, slice_num=None):
    f = open(filename + f"_{objname}_farfield3d.json")
    data = json.load(f)
    ux = data["ux"]
    uy = data["uy"]

    # convert cartesian meshgrid to polar
    cartesian_x, cartesian_y = np.meshgrid(ux, uy)
    radius_list = np.asarray(
        [np.sqrt(i**2 + j**2) for i, j in zip(cartesian_x, cartesian_y)]
    )
    theta_list = []
    for i, j in zip(cartesian_x, cartesian_y):
        theta_row = []
        for ii, jj in zip(i, j):
            if ii == 0 and jj == 0:
                angle = 0
            else:
                angle = np.arccos(ii / np.sqrt(ii**2 + jj**2))
            if jj < 0:
                angle += np.pi
            theta_row.append(angle)

        theta_list.append(theta_row)
    theta_list = np.asarray(theta_list)

    half = len(cartesian_x) // 2
    if slice_num == None:
        for i in range(len(data["wavelength"])):
            wave = round(data["wavelength"][i] * 1e9, 3)
            E = np.asarray(data["farfield_E"][i])

            E_chunk_1 = np.flip(E[0:half, 0:half], 1)
            E_chunk_2 = np.flip(E[0:half, half:], 1)

            E_arrange = np.zeros((len(E), len(E[0])))

            E_arrange[0:half, 0:half] = E_chunk_2
            E_arrange[0:half, half:] = E_chunk_1
            E_arrange[half:, 0:half] = E[half:, 0:half]
            E_arrange[half:, half:] = E[half:, half:]

            E_plot = np.flip(np.flip(E_arrange, 0), 1)

            # E[0:half, half:] = np.zeros((half,half))

            # wrap_radius = np.asarray(radius_list[half][half:]).flatten()
            # wrap_z = np.asarray(E[half][half:]).flatten()
            # wrap_theta = np.asarray([2*np.pi for i in wrap_z]).flatten()

            # theta_plot = np.append(zipped[0], wrap_theta)
            # radius_plot = np.append(zipped[1], wrap_radius)
            # z_plot = np.append(zipped[2], wrap_z)

            fig = plt.figure()
            plt.subplot(projection="polar")
            plt.ylim(0, max(ux))
            # plt.tricontourf([i+3*np.pi/2 for i in theta_plot], radius_plot, z_plot, levels=256, cmap='jet')
            # plt.tricontourf([i+np.pi/2 for i in theta_list.flatten()], radius_list.flatten(), E.flatten())
            plt.tricontourf(
                [i + np.pi / 2 for i in theta_list.flatten()],
                radius_list.flatten(),
                E_plot.flatten(),
                levels=256,
                cmap="jet",
            )
            plt.grid(alpha=0.9)
            plt.savefig(save_name + f"_farfield3d_{int(wave)}.png")
            plt.close()


def modeExpansion_plot(filename, objname, save_name, display):
    profile_file = filename + f"_{objname}_modeExpansion.pkl"
    data_file = open(profile_file, "rb")

    data = pickle.load(data_file)
    data_keys = list(data.keys())
    print(data_keys)

    wave = np.squeeze(data["lambda"])
    T_total = np.squeeze(data["T_total"])
    T_forward = np.squeeze(data["T_forward"])
    T_backward = np.squeeze(data["T_backward"])
    T_net = np.squeeze(data["T_net"])

    plt.plot(wave * 1e9, T_total, marker="o", label="T total")
    plt.plot(wave * 1e9, T_forward, marker="o", label="T foward")
    plt.plot(wave * 1e9, T_backward, marker="o", label="T backward")
    plt.plot(wave * 1e9, T_net, marker="o", label="T net")
    plt.ylabel("Transmission")
    plt.xlabel("Wavelength (nm)")
    plt.legend()
    plt.savefig(save_name + f"_{objname}_modeExpansion.png")
    if display == 1:
        plt.show()
    plt.close()

    data_dict = {
        "wavelength": wave.tolist(),
        "T_total": T_total.tolist(),
        "T_forward": T_forward.tolist(),
        "T_backward": T_backward.tolist(),
        "T_net": T_net.tolist(),
    }

    saveJson(data_dict, save_name + f"_{objname}_modeExpansion")

    return (wave, [T_total, T_forward, T_backward, T_net])


def overlap_integral(file1, file2, slice_num1, slice_num2, display):
    f1 = open(file1)
    f2 = open(file2)

    data1 = json.load(f1)
    data2 = json.load(f2)

    E1 = np.asarray(data1["farfield_E"][slice_num1])
    E2 = np.asarray(data2["farfield_E"][slice_num2])

    if display == 1:
        ux = data1["ux"]
        uy = data1["uy"]

        # convert cartesian meshgrid to polar
        cartesian_x, cartesian_y = np.meshgrid(ux, uy)
        radius_list = np.asarray(
            [np.sqrt(i**2 + j**2) for i, j in zip(cartesian_x, cartesian_y)]
        )
        theta_list = []
        for i, j in zip(cartesian_x, cartesian_y):
            theta_row = []
            for ii, jj in zip(i, j):
                if ii == 0 and jj == 0:
                    angle = 0
                else:
                    angle = np.arccos(ii / np.sqrt(ii**2 + jj**2))
                if jj < 0:
                    angle += np.pi
                theta_row.append(angle)

            theta_list.append(theta_row)
        theta_list = np.asarray(theta_list)

        half = len(cartesian_x) // 2

        E1_chunk_1 = np.flip(E1[0:half, 0:half], 1)
        E1_chunk_2 = np.flip(E1[0:half, half:], 1)

        E1_arrange = np.zeros((len(E1), len(E1[0])))

        E1_arrange[0:half, 0:half] = E1_chunk_2
        E1_arrange[0:half, half:] = E1_chunk_1
        E1_arrange[half:, 0:half] = E1[half:, 0:half]
        E1_arrange[half:, half:] = E1[half:, half:]

        E1_plot = np.flip(np.flip(E1_arrange, 0), 1)

        E2_chunk_1 = np.flip(E2[0:half, 0:half], 1)
        E2_chunk_2 = np.flip(E2[0:half, half:], 1)

        E2_arrange = np.zeros((len(E2), len(E2[0])))

        E2_arrange[0:half, 0:half] = E2_chunk_2
        E2_arrange[0:half, half:] = E2_chunk_1
        E2_arrange[half:, 0:half] = E2[half:, 0:half]
        E2_arrange[half:, half:] = E2[half:, half:]

        E2_plot = np.flip(np.flip(E2_arrange, 0), 1)

        fig = plt.figure()
        plt.subplot(projection="polar")
        plt.ylim(0, max(ux))
        plt.tricontourf(
            [i + np.pi / 2 for i in theta_list.flatten()],
            radius_list.flatten(),
            E1_plot.flatten(),
            levels=256,
            cmap="jet",
            alpha=0.5,
        )
        plt.tricontourf(
            [i + np.pi / 2 for i in theta_list.flatten()],
            radius_list.flatten(),
            E2_plot.flatten(),
            levels=256,
            cmap="jet",
            alpha=0.3,
        )
        plt.grid(alpha=0.9)
        plt.show()

    E1_sum = np.sum(E1)
    E2_sum = np.sum(E2)

    numerator = np.sum(E1 * E2) ** 2
    denominator = np.sum(E2**2) * np.sum(E1**2)

    overlap = numerator / denominator

    return overlap


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
