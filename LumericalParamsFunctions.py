import pandas as pd


def sim_vol_params(
    params_dict, count, tsim, coord, accuracy, mesh_res, xspan, yspan, zspan
):
    component = "sim"
    dict = {}
    dict[f"{component}_t_{count}"] = tsim
    dict[f"{component}_coord_{count}"] = coord
    dict[f"{component}_sim_accuracy_{count}"] = accuracy
    dict[f"mesh_res_{count}"] = mesh_res
    dict[f"{component}_xspan_{count}"] = xspan
    dict[f"{component}_yspan_{count}"] = yspan
    dict[f"{component}_zspan_{count}"] = zspan

    params_dict[f"{component}_{count}"] = dict
    return params_dict


def sim_surface_xy_params(
    params_dict, count, tsim, coord, accuracy, mesh_res, xspan, yspan
):
    component = "sim_surface_xy"
    dict = {}
    dict[f"{component}_t_{count}"] = tsim
    dict[f"{component}_coord_{count}"] = coord
    dict[f"{component}_sim_accuracy_{count}"] = accuracy
    dict[f"mesh_res_{count}"] = mesh_res
    dict[f"{component}_xspan_{count}"] = xspan
    dict[f"{component}_yspan_{count}"] = yspan

    params_dict[f"{component}_{count}"] = dict
    return params_dict


def mesh_params(params_dict, count, coord, xspan, yspan, zspan, dx, dy, dz):
    component = "mesh"
    dict = {}
    dict[f"{component}_coord_{count}"] = coord
    dict[f"{component}_xspan_{count}"] = xspan
    dict[f"{component}_yspan_{count}"] = yspan
    dict[f"{component}_zspan_{count}"] = zspan
    dict[f"{component}_dx_{count}"] = dx
    dict[f"{component}_dy_{count}"] = dy
    dict[f"{component}_dz_{count}"] = dz

    params_dict[f"{component}_{count}"] = dict
    return params_dict


def substrate_params(params_dict, count, coord, material, thick, xspan, yspan):
    component = "substrate"
    dict = {}
    dict[f"{component}_coord_{count}"] = coord
    dict[f"{component}_material_{count}"] = material
    dict[f"{component}_thick_{count}"] = thick
    dict[f"{component}_xspan_{count}"] = xspan
    dict[f"{component}_yspan_{count}"] = yspan
    params_dict[f"{component}_{count}"] = dict
    return params_dict


def disk_params(params_dict, count, coord, material, thick, diameter):
    component = "disk"
    dict = {}
    dict[f"{component}_coord_{count}"] = coord
    dict[f"{component}_material_{count}"] = material
    dict[f"{component}_thick_{count}"] = thick
    dict[f"{component}_diameter_{count}"] = diameter
    params_dict[f"{component}_{count}"] = dict
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


def dipole_params(params_dict, count, coord, wave_center, wave_width, theta, phi):
    component = "dipole"
    dict = {}
    dict[f"{component}_coord_{count}"] = coord
    dict[f"{component}_wave_center_{count}"] = wave_center
    dict[f"{component}_wave_width_{count}"] = wave_width
    dict[f"{component}_theta_{count}"] = theta
    dict[f"{component}_phi_{count}"] = phi

    params_dict[f"{component}_{count}"] = dict
    return params_dict


def gaussian_x_params(
    params_dict, count, coord, height, width, mode, wave_center, wave_width, theta
):
    component = "gaussian_x"
    dict = {}
    dict[f"{component}_height_{count}"] = height
    dict[f"{component}_coord_{count}"] = coord
    dict[f"{component}_width_{count}"] = width
    dict[f"{component}_mode_{count}"] = mode
    dict[f"{component}_wave_center_{count}"] = wave_center
    dict[f"{component}_wave_width_{count}"] = wave_width
    dict[f"{component}_mode_theta_{count}"] = theta

    params_dict[f"{component}_{count}"] = dict
    return params_dict


def gaussian_y_params(
    params_dict, count, coord, height, width, mode, wave_center, wave_width, theta
):
    component = "gaussian_y"
    dict = {}
    dict[f"{component}_height_{count}"] = height
    dict[f"{component}_coord_{count}"] = coord
    dict[f"{component}_width_{count}"] = width
    dict[f"{component}_mode_{count}"] = mode
    dict[f"{component}_wave_center_{count}"] = wave_center
    dict[f"{component}_wave_width_{count}"] = wave_width
    dict[f"{component}_mode_theta_{count}"] = theta

    params_dict[f"{component}_{count}"] = dict
    return params_dict


def profile_xy_params(
    params_dict,
    count,
    coord,
    height,
    width,
    bool,
    freq_point=None,
    wave_center=None,
    wave_width=None,
):
    component = "profile_xy"
    dict = {}
    dict[f"{component}_coord_{count}"] = coord
    dict[f"{component}_height_{count}"] = height
    dict[f"{component}_width_{count}"] = width
    dict[f"{component}_bool_{count}"] = bool
    if bool == 1 and freq_point != None:
        dict[f"{component}_freq_point_{count}"] = freq_point
    if bool == 1 and wave_center != None:
        dict[f"{component}_wave_center_{count}"] = wave_center
        dict[f"{component}_wave_width_{count}"] = wave_width

    params_dict[f"{component}_{count}"] = dict
    return params_dict


def profile_xz_params(
    params_dict,
    count,
    coord,
    height,
    width,
    bool,
    freq_point=None,
    wave_center=None,
    wave_width=None,
):
    component = "profile_xz"
    dict = {}
    dict[f"{component}_coord_{count}"] = coord
    dict[f"{component}_height_{count}"] = height
    dict[f"{component}_width_{count}"] = width
    dict[f"{component}_bool_{count}"] = bool
    if bool == 1 and freq_point != None:
        dict[f"{component}_freq_point_{count}"] = freq_point
    if bool == 1 and wave_center != None:
        dict[f"{component}_wave_center_{count}"] = wave_center
        dict[f"{component}_wave_width_{count}"] = wave_width

    params_dict[f"{component}_{count}"] = dict
    return params_dict


def profile_yz_params(
    params_dict,
    count,
    coord,
    height,
    width,
    bool,
    freq_point=None,
    wave_center=None,
    wave_width=None,
):
    component = "profile_yz"
    dict = {}
    dict[f"{component}_coord_{count}"] = coord
    dict[f"{component}_height_{count}"] = height
    dict[f"{component}_width_{count}"] = width
    dict[f"{component}_bool_{count}"] = bool
    if bool == 1 and freq_point != None:
        dict[f"{component}_freq_point_{count}"] = freq_point
    if bool == 1 and wave_center != None:
        dict[f"{component}_wave_center_{count}"] = wave_center
        dict[f"{component}_wave_width_{count}"] = wave_width

    params_dict[f"{component}_{count}"] = dict
    return params_dict


def line_x_params(
    params_dict,
    count,
    coord,
    width,
    bool,
    freq_point=None,
    wave_center=None,
    wave_width=None,
):
    component = "line_xz"
    dict = {}
    dict[f"{component}_coord_{count}"] = coord
    dict[f"{component}_width_{count}"] = width
    dict[f"{component}_bool_{count}"] = bool
    if bool == 1 and freq_point != None:
        dict[f"{component}_freq_point_{count}"] = freq_point
    if bool == 1 and wave_center != None:
        dict[f"{component}_wave_center_{count}"] = wave_center
        dict[f"{component}_wave_width_{count}"] = wave_width

    params_dict[f"{component}_{count}"] = dict
    return params_dict


def point_monitor_params(
    params_dict, count, coord, bool, freq_point=None, wave_center=None, wave_width=None
):
    component = "point_profile"
    dict = {}
    dict[f"{component}_coord_{count}"] = coord
    dict[f"{component}_bool_{count}"] = bool
    if bool == 1 and freq_point != None:
        dict[f"{component}_freq_point_{count}"] = freq_point
    if bool == 1 and wave_center != None:
        dict[f"{component}_wave_center_{count}"] = wave_center
        dict[f"{component}_wave_width_{count}"] = wave_width

    params_dict[f"{component}_{count}"] = dict
    return params_dict


def volume_monitor_params(
    params_dict,
    count,
    coord,
    xspan,
    yspan,
    zspan,
    bool,
    freq_point=None,
    wave_center=None,
    wave_width=None,
):
    component = "volume_profile"
    dict = {}
    dict[f"{component}_coord_{count}"] = coord
    dict[f"{component}_xspan_{count}"] = xspan
    dict[f"{component}_yspan_{count}"] = yspan
    dict[f"{component}_zspan_{count}"] = zspan
    dict[f"{component}_bool_{count}"] = bool
    if bool == 1 and freq_point != None:
        dict[f"{component}_freq_point_{count}"] = freq_point
    if bool == 1 and wave_center != None:
        dict[f"{component}_wave_center_{count}"] = wave_center
        dict[f"{component}_wave_width_{count}"] = wave_width

    params_dict[f"{component}_{count}"] = dict
    return params_dict


def movie_xy_params(params_dict, count, coord, height, width):
    component = "movie_xy"
    dict = {}
    dict[f"{component}_coord_{count}"] = coord
    dict[f"{component}_height_{count}"] = height
    dict[f"{component}_width_{count}"] = width

    params_dict[f"{component}_{count}"] = dict
    return params_dict


def movie_xz_params(params_dict, count, coord, height, width):
    component = "movie_xz"
    dict = {}
    dict[f"{component}_coord_{count}"] = coord
    dict[f"{component}_height_{count}"] = height
    dict[f"{component}_width_{count}"] = width

    params_dict[f"{component}_{count}"] = dict
    return params_dict


def movie_yz_params(params_dict, count, coord, height, width):
    component = "movie_yz"
    dict = {}
    dict[f"{component}_coord_{count}"] = coord
    dict[f"{component}_height_{count}"] = height
    dict[f"{component}_width_{count}"] = width

    params_dict[f"{component}_{count}"] = dict
    return params_dict


def qanalysis_params(params_dict, count, coord, height, width, wave_min, wave_max):
    component = "qanalysis"
    dict = {}
    dict[f"{component}_coord_{count}"] = coord
    dict[f"{component}_height_{count}"] = height
    dict[f"{component}_width_{count}"] = width
    dict[f"{component}_wave_min_{count}"] = wave_min
    dict[f"{component}_wave_max_{count}"] = wave_max

    params_dict[f"{component}_{count}"] = dict
    return params_dict


def Lowqanalysis_params(params_dict, count, coord, height, width):
    component = "Lowqanalysis"
    dict = {}
    dict[f"{component}_coord_{count}"] = coord
    dict[f"{component}_height_{count}"] = height
    dict[f"{component}_width_{count}"] = width

    params_dict[f"{component}_{count}"] = dict
    return params_dict


def importgds_params(
    params_dict, count, name, coord, material, thick, filename, cellname, layer
):
    component = "importgds"
    dict = {}
    dict[f"{component}_name_{count}"] = name
    dict[f"{component}_coord_{count}"] = coord
    dict[f"{component}_material_{count}"] = material
    dict[f"{component}_thick_{count}"] = thick
    dict[f"{component}_filename_{count}"] = filename
    dict[f"{component}_cellname_{count}"] = cellname
    dict[f"{component}_layer_{count}"] = layer

    params_dict[f"{component}_{count}"] = dict
    return params_dict


def saveDataFrameCSV(
    count: int, params_dict: dict, save_path: str, dataframe: pd.DataFrame | None = None
):
    # dataframe will have to be a multiindex because of nested dictionary #
    # restructure dict to cast into dataframe #
    reform = {
        (outerKey, innerKey): values
        for outerKey, innerDict in params_dict.items()
        for innerKey, values in innerDict.items()
    }

    # create datafram from multiindex #
    params_df = pd.DataFrame.from_dict(reform, orient="index").transpose()
    params_df.columns = pd.MultiIndex.from_tuples(params_df.columns)

    # 'append' dataframe to existing #
    if dataframe is not None:
        params_df = pd.concat([dataframe, params_df], ignore_index=True)

    params_csv = params_df.to_csv(save_path + ".csv")

    return params_df
