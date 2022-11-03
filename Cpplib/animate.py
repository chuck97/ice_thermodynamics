#!/usr/bin/env python3

import sys
import os
sys.path.insert(1, os.path.join(sys.path[0], '../Python'))

import utils.engine as ue
import utils.visualizer as uv
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plotting an animation from json files.')
    parser.add_argument('-i', '--input', default='.',
                        help='folder in which json files are contained (current folder is default)')
    parser.add_argument('-ip', '--iceprefix',
                        help='prefix for names of ice jsons')
    parser.add_argument('-sp', '--snowprefix',
                        help='prefix for names of snow jsons')
    parser.add_argument('-o', '--output', default='output.mp4',
                        help='name for the output file (default format is .mp4)')
    
    args = parser.parse_args()
    process = ue.get_process_from_json_folder(path=args.input, ice_prefix=args.iceprefix, snow_prefix=args.snowprefix)
    uv.animate([process], savepath=(args.output if '.' in args.output else args.output + '.mp4'))