import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as anim
import matplotlib.cm as mcm
import matplotlib.colors as clr
from matplotlib.ticker import FormatStrFormatter
from matplotlib.collections import LineCollection
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
    new_cmap = clr.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap

def animate(processes,
            rho_snow=snow_density, rho_water=water_density,
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

    
def discretize_cmap(cmap, bounds=None, vmin=None, vmax=None, prec=None):
    
    assert bounds is not None or (vmin is not None and vmax is not None), \
           "Either bounds or limits should be set!"
    
    cmap_base = mcm.get_cmap(cmap)
    
    if bounds is None:
        if prec is None:
            prec = (vmax - vmin)*.1
        bounds = np.arange(vmin, vmax, prec)
        if bounds[-1] != vmax:
            bounds = np.concatenate((bounds, [vmax]))
            
    norm_custom = clr.BoundaryNorm(bounds, cmap_base.N, extend='both')
    
    return cmap_base, norm_custom, bounds


def timeseries_img(process, rho_snow=snow_density, rho_water=water_density,
                   y_points=100,
                   tmin_ice=None, tmax_ice=None, prec_ice=None, bounds_ice=None,
                   tmin_snow=None, tmax_snow=None, prec_snow=None, bounds_snow=None,
                   cmap_ice='Blues', cmap_snow='Greys', color_empty=[208, 245, 226],
                   savepath=None):

    Z_i, Z_s = get_Z(process, rho_snow, rho_water)
    z_min = min(Z_i[:, 0])
    z_max = max(Z_s[:, -1])

    Z_mesh = np.linspace(z_min, z_max, y_points + 1)
    temp_mesh = np.array([np.interp(x=Z_mesh,
                                    xp=np.concatenate((Z_i_line, (Z_s_line[1:] if has_snow else []))),
                                    fp=np.concatenate(([T_oi], T_i_line, [T_is],
                                                       (T_s_line if has_snow else []),
                                                       ([T_sa] if has_snow else [])
                                                      )),
                                    left=np.nan, right=np.nan
                                   ) \
                          for Z_i_line, Z_s_line, T_oi, T_i_line, T_is, T_s_line, T_sa, has_snow \
                          in zip(Z_i, Z_s,
                                 process.oi_temp_history, process.ice_temp_history, process.is_temp_history,
                                 process.snow_temp_history, process.sa_temp_history,
                                 process.snow_presence_history
                                )
                         ])

    nan_filter = np.isnan(temp_mesh)
    ice_filter = ((Z_i[:, [0]] <= Z_mesh) & (Z_mesh <= Z_i[:, [-1]]))
    snow_filter = ((Z_s[:, [0]] < Z_mesh) & (Z_mesh <= Z_s[:, [-1]]))
    
    img = np.array([[[np.nan]*4]*(y_points+1) for i in range(process.get_length())])
    img[nan_filter] = np.concatenate((np.array(color_empty)/256, [1.0]))

    if bounds_ice is None:
        if tmin_ice is None:
            tmin_ice = np.nanmin(process.ice_temp_history)
        if tmax_ice is None:
            tmax_ice = np.nanmax(process.ice_temp_history)
    cmap_ice, norm_ice, bounds_ice = discretize_cmap(cmap_ice, bounds_ice, tmin_ice, tmax_ice, prec_ice)
    img[ice_filter] = cmap_ice(norm_ice(temp_mesh[ice_filter]))

    if process.snow_temp_history.any():
        if bounds_snow is None:
            if tmin_snow is None:
                tmin_snow = np.nanmin(process.snow_temp_history)
            if tmax_snow is None:
                tmax_snow = np.nanmax(process.snow_temp_history)
        cmap_snow, norm_snow, bounds_snow = discretize_cmap(cmap_snow, bounds_snow, tmin_snow, tmax_snow, prec_snow)
        img[snow_filter] = cmap_snow(norm_snow(temp_mesh[snow_filter]))
        
    fig = plt.figure(figsize=(30, 10))
    ax = plt.axes()
    ax.imshow(img.swapaxes(0, 1), aspect='auto', origin='lower',
              extent=[process.timeline[0]/3600, process.timeline[-1]/3600, Z_mesh[0], Z_mesh[-1]])
    ax.set_xlabel('t, hours', size=20)
    ax.set_ylabel('Z, m', size=20)
    ax.tick_params(axis='both', labelsize=15)
    ax.axhline(lw=3, ls='--', color='c')

    cax_ice = fig.add_axes([ax.get_position().x0, ax.get_position().y0 - 0.15,
                            ax.get_position().width, 0.05])
    fig.colorbar(mcm.ScalarMappable(norm_ice, cmap_ice), cax=cax_ice,
                 orientation='horizontal', ticks=bounds_ice)
    cax_ice.set_xlabel('ice', size=20)
    cax_ice.tick_params(axis='x', labelsize=15)

    if process.snow_temp_history.any():
        cax_snow = fig.add_axes([ax.get_position().x0, ax.get_position().y0 - 0.3,
                                 ax.get_position().width, 0.05])
        fig.colorbar(mcm.ScalarMappable(norm_snow, cmap_snow), cax=cax_snow,
                     orientation='horizontal', ticks=bounds_snow)
        cax_snow.set_xlabel('snow', size=20)
        cax_snow.tick_params(axis='x', labelsize=15)
    
    if savepath is not None:
        fig.savefig(savepath, bbox_inches='tight')
