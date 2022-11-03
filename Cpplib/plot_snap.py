#!/usr/bin/env python3

import numpy as np
import json as js
import sys
import os
sys.path.insert(1, os.path.join(sys.path[0], '../Python'))

import utils.visualizer as uv
import argparse


def get_Z(dzi, rho_i, rho_w, dzs=None, rho_s=None):
    
    with_snow = (dzs is not None) and (rho_s is not None)
    
    bottomline = (-np.dot(dzi, rho_i) - (np.dot(dzs, rho_s) if with_snow else 0))/rho_w
    
    Z_i = np.concatenate(([bottomline], dzi)).cumsum()
    Z_i = np.append(Z_i, Z_i[[-1]])
    Z_i[1:-1] -= [dz/2 for dz in dzi]
    
    if with_snow:
        Z_s = np.concatenate((Z_i[[-1]], dzs)).cumsum()
        Z_s = np.append(Z_s, Z_s[[-1]])
        Z_s[1:-1] -= [dz/2 for dz in dzs]
        return Z_i, Z_s
    
    else:
        return Z_i

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description='Plotting a snap from json files.')
    parser.add_argument('-if', '--icefiles', nargs='*',
                        help="jsons' names for ice")
    parser.add_argument('-sf', '--snowfiles', nargs='*', default=[],
                        help="jsons' names for snow")
    parser.add_argument('-ip', '--iceprefix', default='',
                        help='prefix for names of ice jsons')
    parser.add_argument('-sp', '--snowprefix', default='',
                        help='prefix for names of snow jsons')
    parser.add_argument('-o', '--output', default='output.png',
                        help='name for the output file (default format is .png)')

    args = parser.parse_args()
    
    json_pairs = {}

    for ice_path in args.icefiles:
        time = int(ice_path[ice_path.rfind('/') + len(args.iceprefix) + 1:-5])
        with open(ice_path) as json_ice:
            json_pairs.setdefault(time, {})['ice'] = js.load(json_ice)
                
    for snow_path in args.snowfiles:
        time = int(snow_path[snow_path.rfind('/') + len(args.snowprefix) + 1:-5])
        with open(snow_path) as json_snow:
            json_pairs.setdefault(time, {})['snow'] = js.load(json_snow)
    
    Z_i_list = []
    Z_s_list = []
    T_i_list = []
    T_s_list = []
    
    for time, pair in json_pairs.items():
        obj_ice = pair['ice']
        if 'snow' not in pair: 
            Z_i = get_Z(obj_ice["cells_thickness_array"], obj_ice["cells data"]["cells_density_array"], 1023.0)
            T_i = [obj_ice["single data"]["down_temperature"]] + obj_ice["cells data"]["cells_temperature_array"] \
                + [obj_ice["single data"]["up_temperature"]]
            Z_i_list.append(Z_i)
            T_i_list.append(T_i)
            Z_s_list.append(None)
            T_s_list.append(None)
            
        else:
            obj_snow = pair['snow']
            Z_i, Z_s = get_Z(obj_ice["cells_thickness_array"], obj_ice["cells data"]["cells_density_array"], 1023.0,
                         obj_snow["cells_thickness_array"], obj_snow["cells data"]["cells_density_array"])
            T_i = [obj_ice["single data"]["down_temperature"]] + obj_ice["cells data"]["cells_temperature_array"] \
                + [obj_ice["single data"]["up_temperature"]]
            T_s = [obj_snow["single data"]["down_temperature"]] + obj_snow["cells data"]["cells_temperature_array"] \
                + [obj_snow["single data"]["up_temperature"]]
            Z_i_list.append(Z_i)
            T_i_list.append(T_i)
            Z_s_list.append(Z_s)
            T_s_list.append(T_s)
            
        uv.plot_snap(Z_i_list, Z_s_list, T_i_list, T_s_list, names=["iteration %d"%time for time in json_pairs.keys()],
                     savepath=args.output)
        
