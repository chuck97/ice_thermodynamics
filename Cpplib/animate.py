#!/usr/bin/env python3

import sys
import os
sys.path.insert(1, os.path.join(sys.path[0], '../Python'))

import utils.engine as ue
import utils.visualizer as uv
import argparse


if __name__ == "__main__":
    
    class CustomAppend(argparse._AppendAction):
        def __call__(self, parser, namespace, values, option_string=None):
            if not (len(values) in (3, 5)):
                raise argparse.ArgumentError(self, "%s takes 3 or 5 values, %d given"%(option_string, len(values)))

            val_dict = {}
            val_dict['path'] = values[0]
            if values[1] not in ('ip', 'sp'):
                raise argparse.ArgumentError(self, "key shold be either 'ip' or 'sp', %s given"%values[1])
            val_dict[values[1]] = values[2]

            if len(values) == 5:
                if values[3] not in ('ip', 'sp'):
                    raise argparse.ArgumentError(self, "key shold be either 'ip' or 'sp', '%s' given"%values[3])
                val_dict[values[3]] = values[4]
            super(CustomAppend, self).__call__(parser, namespace, val_dict, option_string)
            
    class CustomFormatter(argparse.HelpFormatter):
        def _format_args(self, action, default_metavar):
            if action.__class__ is CustomAppend:
                return 'INPUT ([ip ICE_PREFIX] [sp SNOW_PREFIX])'
            else:
                return super(CustomFormatter, self)._format_args(action, default_metavar)
    
parser = argparse.ArgumentParser(description='Plotting an animation from json files.', formatter_class=CustomFormatter)
parser.add_argument('-i', '--input', action=CustomAppend, nargs='+',
                    help='Specifies input with prefixes')
parser.add_argument('-flt', '--floating', action='store_true',
                    help='Enable floating (ground level will be sea level)')
parser.add_argument('-o', '--output', default='output.mp4',
                    help='Name for the output file (default format is .mp4).')

args = parser.parse_args()
processes = []
for proc_args in args.input:
    if 'ip' in proc_args and 'sp' in proc_args:
        processes.append(ue.get_process_from_json_folder(path=proc_args['path'], ice_prefix=proc_args['ip'],
                                                         snow_prefix=proc_args['sp']))
    elif 'ip' in proc_args:
        processes.append(ue.get_process_from_json_folder(path=proc_args['path'], ice_prefix=proc_args['ip']))
    else:
        processes.append(ue.get_process_from_json_folder(path=proc_args['path'], snow_prefix=proc_args['sp']))

uv.animate(processes, floating=args.floating, savepath=(args.output if '.' in args.output else args.output + '.mp4'))