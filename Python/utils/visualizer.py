import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as anim
from matplotlib import cm as mcm
from matplotlib.ticker import FormatStrFormatter
from matplotlib.collections import LineCollection, PathCollection
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.cm import ScalarMappable
matplotlib.rcParams['animation.embed_limit'] = 5000
plt.rcParams["animation.html"] = "jshtml"
from utils.engine import rho_w as water_density, rho_s as snow_density
import itertools


def get_Z(process, rho_s, rho_w):
    
    dzi, dzs, rho_i = process.ice_dz_history, process.snow_dz_history, process.ice_density_history
    
    bottomlines = np.array([(-np.dot(dzi_t, rho_i_t) - rho_s*sum(dzs_t))/rho_w \
                            for dzi_t, dzs_t, rho_i_t in zip(dzi, dzs, rho_i)])
    
    Z_i = np.concatenate((bottomlines.reshape(-1, 1), dzi), axis=1).cumsum(axis=1)
    Z_i = np.append(Z_i, Z_i[:, [-1]], axis=1)
    Z_i[:, 1:-1] -= dzi/2
    
    Z_s = np.concatenate((Z_i[:, [-1]], dzs), axis=1).cumsum(axis=1)
    Z_s = np.append(Z_s, Z_s[:, [-1]], axis=1)
    Z_s[:, 1:-1] -= dzs/2
    
    return Z_i, Z_s

def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap

def animate(processes,
            rho_water=water_density, rho_snow=snow_density,
            clip_start=None, clip_end=None,
            t_min=None, t_max=None,
            cmap=None, names=None,
            savepath=None, dpi=None):
    
    assert len(set(process.get_length() for process in processes)) == 1, "Length of processes are not equal!"
    
    if names is not None: 
        assert len(names) == len(processes), "Number of names and number of processes should be equal!"
    
    def run(data):

        Z_ice_line_procs, Z_snow_line_procs, proc_data, frame_number = data
        
        for Z_ice_line, Z_snow_line, \
            (_, _,
             time,
             T_oi, T_ice, T_is, T_snow, T_sa,
             rho_ice_line, has_snow), \
            line_ice, line_snow, markers in zip(Z_ice_line_procs, Z_snow_line_procs, proc_data,
                                                lines_ice, lines_snow, markers_ice):
        
            T_ice_full = np.concatenate(([T_oi], T_ice, [T_is]))

            if has_snow:
                T_snow_full = np.concatenate(([T_is], T_snow, [T_sa]))
                line_snow.set_data(T_snow_full, Z_snow_line)
            else:
                line_snow.set_data([], [])

            points_ice = [[[T, Z]] for T, Z in zip(T_ice_full, Z_ice_line)]
            segments_ice = np.concatenate([points_ice[:-1], points_ice[1:]], axis=1)
            line_ice.set_segments(segments_ice)
            line_ice.set_array(T_ice_full)             
            markers.set_offsets([[T, Z] for T, Z in zip(T_ice_full, Z_ice_line)])
            markers.set_array(T_ice_full)
            ax.set_title('Time: %.2f hours'%(time/3600), size=20)
    #     ax.set_yticks(np.insert(-process.ice_dz_history.cumsum(), 0, 0))
    
        print("Rendered {:.0%} ({}/{}) frames".format(frame_number/frames_count, frame_number, frames_count), end="\r")

        return lines_ice + lines_snow + markers_ice
    
    if clip_end:
        clip_end += 1
    frames_count = (clip_end if clip_end else len(processes[0].timeline) + 1) \
                 - (clip_start if clip_start else 0)
    
    all_Z_i = []
    all_Z_s = []
    for process in processes:
        Z_i, Z_s = get_Z(process, rho_snow, rho_water)
        all_Z_i.append(Z_i)
        all_Z_s.append(Z_s)
    z_min = min([Z_i[:, 0].min() for Z_i in all_Z_i])
    z_max = max([Z_s[:, -1].max() for Z_s in all_Z_s])
    temp_pairs = [(process.get_temp('min'), process.get_temp('max')) for process in processes]
    t_min = (t_min if t_min is not None else min(temp_pair[0] for temp_pair in temp_pairs))
    t_max = (t_max if t_max is not None else max(temp_pair[1] for temp_pair in temp_pairs))
    
    cmaps = ['Blues', 'Purples', 'Greens', 'cool', 'winter'] 
        
#     if cmap is None:
        
#         cdict = {'red':  [[0.0, 0.0, 0.0],
#                           [1.0, 0.6, 0.6]],
#                  'green':  [[0.0, 0.2, 0.2],
#                             [1.0, 0.85, 0.85]],
#                  'blue':  [[0.0, 0.5, 0.5],
#                            [1.0, 1.0, 1.0]]}
#         cmap = LinearSegmentedColormap('MyCmap', segmentdata=cdict)

    fig = plt.figure(figsize=(15, 15))
    ax = fig.add_subplot(xlim=(t_min - (t_max - t_min)*0.1, t_max + (t_max - t_min)*0.1),
                         ylim=(z_min - (z_max - z_min)*0.1, z_max + (z_max - z_min)*0.1))
    lines_ice = []
    markers_ice = []
    lines_snow = []
    for i, temp_pair in enumerate(temp_pairs):
        norm = plt.Normalize(*temp_pair)
        cmap = truncate_colormap(plt.get_cmap(cmaps[i%len(cmaps)]), 0.3, 1)
        lc = LineCollection([], linewidths=3, cmap=cmap, norm=norm, animated=True)
        line_ice = ax.add_collection(lc)
        lines_ice.append(line_ice)
        markers = ax.scatter([], [], s=60, animated=True)
        markers.set_cmap(cmap)
        markers_ice.append(markers)
        line_snow, = ax.plot([], [], color='grey', lw=3, marker='o', animated=True)
        lines_snow.append(line_snow)
    ax.axhline(ls='--', lw=3, color='navy')
    ax.set_xlabel(r'T, $^o C$', size=20)
    ax.set_ylabel('Z, m', size=20)
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax.tick_params(axis='both', labelsize=15)
    ax.grid()

    pad = 0
    for i, line_ice in enumerate(lines_ice):
        cbar = fig.colorbar(line_ice, orientation='horizontal', pad=pad)
        if names is None:
            cbar.ax.set_xlabel('process %d'%(len(lines_ice)-i), size=15)
        else:
            cbar.ax.set_xlabel(names[i], size=15)
        cbar.ax.tick_params(labelsize=15)
        pad += 0.08
        
    animation = anim.FuncAnimation(fig, run, zip(zip(*[Z_i[clip_start:clip_end] for Z_i in all_Z_i]),
                                                 zip(*[Z_s[clip_start:clip_end] for Z_s in all_Z_s]),
                                                 zip(*[process[clip_start:clip_end].get_zip() for process in processes]),
                                                 range(1, frames_count + 1)
                                                ),
                                   save_count=frames_count, interval=30, blit=True)
    
    if savepath:
        animation.save(savepath, dpi=dpi)
        
    return animation
    
    
def paint_snow(arr, num_of_cells, snow_color):
    ind_snow = (arr != [0, 0, 0, 1]).nonzero()[0].max() + 1
    arr[ind_snow:(ind_snow+num_of_cells)] = snow_color
    return arr


def timeseries_img(process, rho_water, rho_snow,
                   y_points=50, x_ticks=11, y_ticks=6,
                   t_min=None, t_max=None,
                   cmap_ice='Blues', cmap_snow='Greys', color_empty='black',
                   savepath=None):

    Z_i, Z_s = get_Z(process, rho_snow, rho_water)
    z_min = min(Z_i[:, 0])
    z_max = max(Z_s[:, -1])
    
    vmin = (t_min if t_min is not None else process.get_temp('min'))
    vmax = (t_max if t_max is not None else process.get_temp('max'))

    z_mesh = np.tile(np.linspace(z_min, z_max, y_points + 1), (process.get_length(), 1))
    ice_filter = ((Z_i[:, [0]] > z_mesh) | (z_mesh >= Z_i[:, [-1]]))[:, ::-1].T
    snow_filter = ((Z_s[:, [0]] > z_mesh) | (z_mesh >= Z_s[:, [-1]]))[:, ::-1].T
    
    temp_mesh = np.array([np.interp(x=Z_new,
                                    xp=np.concatenate((Z_i_line, (Z_s_line[1:] if has_snow else []))),
                                    fp=np.concatenate(([T_oi], T_i_line, [T_is],
                                                       (T_s_line if has_snow else []),
                                                       ([T_sa] if has_snow else [])
                                                      )),
                                    left=0.0, right=0.0
                                   ) \
                          for Z_new, Z_i_line, Z_s_line, T_oi, T_i_line, T_is, T_s_line, T_sa, has_snow \
                          in zip(z_mesh, Z_i, Z_s,
                                 process.oi_temp_history, process.ice_temp_history, process.is_temp_history,
                                 process.snow_temp_history, process.sa_temp_history,
                                 process.snow_presence_history
                                )
                         ])[:, ::-1].T

    ice_masked = np.ma.masked_array(temp_mesh, ice_filter)
    snow_masked = np.ma.masked_array(temp_mesh, snow_filter)
    empty_mask = np.ma.masked_array(temp_mesh, ice_filter & snow_filter)

    fig = plt.figure(figsize=(40, 20))
    ax = fig.add_subplot()
    curr_cmap = mcm.get_cmap().copy()
    curr_cmap.set_bad(color_empty)
    img_empty = ax.imshow(empty_mask, cmap=curr_cmap, aspect='auto')
    img_ice = ax.imshow(ice_masked, cmap=cmap_ice, aspect='auto')#, interpolation='bilinear')
    img_snow = ax.imshow(snow_masked, cmap=cmap_snow, aspect='auto')#, interpolation='bilinear')
    ax.set_xlabel('t, hours', size=20)
    ax.set_ylabel('Z, m', size=20)
    ax.set_xticks(np.linspace(0, process.get_length() - 1, x_ticks, dtype=int))
    ax.set_yticks(np.linspace(0, y_points, y_ticks, dtype=int))
    ax.set_xticklabels(round(process.timeline[int(tick)]/3600, 2) for tick in ax.get_xticks())
    ax.set_yticklabels(np.round(np.linspace(z_max, z_min, y_ticks), 2))
    ax.axhline(y_points*z_max/(z_max-z_min), lw=3, ls='--', color='c')
    ax.tick_params(axis='both', labelsize=15)
    
    if process.snow_presence_history.any():
        cbar_snow = fig.colorbar(img_snow, orientation='horizontal', pad=-0.01)
        cbar_snow.ax.set_xlabel(r'Temperature of snow, $^o C$', size=20)
        cbar_snow.ax.tick_params(labelsize=15)
        
#     cbar_ice = fig.colorbar(img_ice, orientation='horizontal', pad=0.08)
#     cbar_ice.ax.set_xlabel(r'Temperature of ice, $^o C$', size=20)
#     cbar_ice.ax.tick_params(labelsize=15)
    
    if savepath is not None:
        fig.savefig(savepath, bbox_inches='tight')

        
    return temp_mesh, ice_filter, snow_filter