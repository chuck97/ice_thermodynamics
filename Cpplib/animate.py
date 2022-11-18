#!/usr/bin/env python3

import sys
import os
sys.path.insert(1, os.path.join(sys.path[0], '../Python'))

import utils.engine as ue
import utils.visualizer as uv
import argparse

#TODO: сделать нормальный парсер через subparser

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plotting an animation from json files.')
    parser.add_argument('-i', '--input', action='append', nargs='+',
                        help='Three values: (path, ice prefix, snow prefix(optional)). Can be used multiple times.')
#     parser.add_argument('-ip', '--iceprefix',
#                         help='prefix for names of ice jsons')
#     parser.add_argument('-sp', '--snowprefix',
#                         help='prefix for names of snow jsons')
    parser.add_argument('-o', '--output', default='output.mp4',
                        help='Name for the output file (default format is .mp4).')
    
    args = parser.parse_args()
    processes = []
    for proc_args in args.input:
        if len(proc_args) == 2:
            processes.append(ue.get_process_from_json_folder(path=proc_args[0], ice_prefix=proc_args[1]))
        elif len(proc_args) == 3:
            processes.append(ue.get_process_from_json_folder(path=proc_args[0], ice_prefix=proc_args[1],
                                                             snow_prefix=proc_args[2]))
#     process = ue.get_process_from_json_folder(path=args.input, ice_prefix=args.iceprefix, snow_prefix=args.snowprefix)
    uv.animate(processes, savepath=(args.output if '.' in args.output else args.output + '.mp4'))