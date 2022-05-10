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
import itertools

def get_Z(dzi, dzs, rho_i, rho_s, rho_w):
    
    
    bottomlines = np.array([(-np.dot(dzi_t, rho_i_t) - rho_s*sum(dzs_t))/rho_w \
                            for dzi_t, dzs_t, rho_i_t in zip(dzi, dzs, rho_i)])
    
    Z_i = np.concatenate((bottomlines.reshape(-1, 1), dzi), axis=1).cumsum(axis=1)
    Z_i = np.append(Z_i, Z_i[:, [-1]], axis=1)
    Z_i[:, 1:-1] -= dzi/2
    
    Z_s = np.concatenate((Z_i[:, [-1]], dzs), axis=1).cumsum(axis=1)
    Z_s = np.append(Z_s, Z_s[:, [-1]], axis=1)
    Z_s[:, 1:-1] -= dzs/2
    
    return Z_i, Z_s

def animate(dz_ice_arr, dz_snow_arr,
            temp_oi, temp_ice, temp_is, temp_snow, temp_sa,
            time_arr,
            rho_ice_arr, rho_water, rho_snow,
            has_snow_filter,
            clip_start=0, clip_end=-2,
            t_min=None, t_max=None, cmap=None, savepath=None):
    
    dz_ice_arr = np.array(dz_ice_arr)
    dz_snow_arr = np.array(dz_snow_arr)
    temp_oi = np.array(temp_oi)
    temp_ice = np.array(temp_ice)
    temp_is = np.array(temp_is)
    temp_snow = np.array(temp_snow)
    temp_sa = np.array(temp_sa)
    time_arr = np.array(time_arr)
    rho_ice_arr = np.array(rho_ice_arr)
    
    def run(data):

        Z_ice_line, Z_snow_line,\
        T_oi, T_ice, T_is, T_snow, T_sa,\
        time,\
        rho_ice_line, rho_water,\
        has_snow = data
        
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
    #     ax.set_yticks(np.insert(-dz_ice_arr.cumsum(), 0, 0))

        return line_ice, line_snow, markers
    
    
    Z_i, Z_s = get_Z(dz_ice_arr, dz_snow_arr, rho_ice_arr, rho_snow, rho_water)
    z_min = min(Z_i[:, 0])
    z_max = max(Z_s[:, -1])
    temp_min = min(temp_arr.min() for temp_arr in [temp_oi, temp_ice, temp_is, temp_snow, temp_sa] if temp_arr.size > 0)
    temp_max = max(temp_arr.max() for temp_arr in [temp_oi, temp_ice, temp_is, temp_snow, temp_sa] if temp_arr.size > 0)
     
    t_min = (t_min if t_min is not None else temp_min)
    t_max = (t_max if t_max is not None else temp_max)
    
    if cmap is None:
        
        cdict = {'red':  [[0.0,  0.0, 0.0],
                  [1.0,  0.6, 0.6]],
         'green':  [[0.0,  0.2, 0.2],
                    [1.0,  0.85, 0.85]],
         'blue':  [[0.0,  0.5, 0.5],
                   [1.0,  1.0, 1.0]]}
        cmap = LinearSegmentedColormap('MyCmap', segmentdata=cdict)

    fig = plt.figure(figsize=(15, 10))
    ax = fig.add_subplot(xlim=(t_min - (t_max - t_min)*0.1, t_max + (t_max - t_min)*0.1),
                         ylim=(z_min - (z_max - z_min)*0.1, z_max + (z_max - z_min)*0.1))
    norm = plt.Normalize(t_max, t_min)
    lc = LineCollection([], linewidths=3, cmap=cmap, norm=norm, animated=True)
    line_ice = ax.add_collection(lc)
    line_snow, = ax.plot([], [], color='grey', lw=3, marker='o', animated=True)
    markers = ax.scatter([], [], s=60, animated=True)
    markers.set_cmap(cmap)
    ax.axhline(ls='--', lw=3, color='navy')
    ax.set_xlabel(r'T, $^o C$', size=20)
    ax.set_ylabel('Z, m', size=20)
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax.tick_params(axis='both', labelsize=15)
    ax.grid()
    fig.colorbar(line_ice)
        
    animation = anim.FuncAnimation(fig, run, zip(Z_i[clip_start:clip_end+1], Z_s[clip_start:clip_end+1],
                                                 temp_oi[clip_start:clip_end+1], temp_ice[clip_start:clip_end+1], 
                                                 temp_is[clip_start:clip_end+1], temp_snow[clip_start:clip_end+1], 
                                                 temp_sa[clip_start:clip_end+1],
                                                 time_arr[clip_start:clip_end+1],
                                                 rho_ice_arr[clip_start:clip_end+1], 
                                                 [rho_water]*dz_ice_arr[clip_start:clip_end+1].shape[0],
                                                 has_snow_filter[clip_start:clip_end+1]),
                                   save_count=dz_ice_arr[clip_start:clip_end+1].shape[0], interval=30, blit=True)
    
    if savepath:
        animation.save(savepath)
        
    return animation
    
    
def paint_snow(arr, num_of_cells, snow_color):
    ind_snow = (arr != [0, 0, 0, 1]).nonzero()[0].max() + 1
    arr[ind_snow:(ind_snow+num_of_cells)] = snow_color
    return arr


def timeseries_img(dz_ice_arr, dz_snow_arr,
                   temp_oi, temp_ice, temp_is, temp_snow, temp_sa,
                   time_arr,
                   rho_ice_arr, rho_water, rho_snow,
                   has_snow_filter,
#                    snow_color = [0.5019607843137255,
#                                  0.5019607843137255,
#                                  0.5019607843137255,
#                                  1.0],
                   y_points=50, max_aspect=0.66, x_ticks=11, y_ticks=6,
                   t_min=None, t_max=None,
                   cmap_ice='Blues', cmap_snow='Greys', color_empty='black',
                   savepath=None):

    dz_ice_arr = np.array(dz_ice_arr)
    dz_snow_arr = np.array(dz_snow_arr)
    temp_oi = np.array(temp_oi)
    temp_ice = np.array(temp_ice)
    temp_is = np.array(temp_is)
    temp_snow = np.array(temp_snow)
    temp_sa = np.array(temp_sa)
    time_arr = np.array(time_arr)
    rho_ice_arr = np.array(rho_ice_arr)

    Z_i, Z_s = get_Z(dz_ice_arr, dz_snow_arr, rho_ice_arr, rho_snow, rho_water)
    z_min = min(Z_i[:, 0])
    z_max = max(Z_s[:, -1])
    temp_min = min(temp_arr.min() for temp_arr in [temp_oi, temp_ice, temp_is, temp_snow, temp_sa] if temp_arr.size > 0)
    temp_max = max(temp_arr.max() for temp_arr in [temp_oi, temp_ice, temp_is, temp_snow, temp_sa] if temp_arr.size > 0)
    
    y_points = min(y_points, int(max_aspect*len(time_arr)))
    
    vmin = (t_min if t_min is not None else temp_min)
    vmax = (t_max if t_max is not None else temp_max)

    z_mesh = np.tile(np.linspace(z_min, z_max, y_points + 1), (len(time_arr), 1))
    ice_filter = ((Z_i[:, [0]] > z_mesh) | (z_mesh >= Z_i[:, [-1]]))[:, ::-1].T
    snow_filter = ((Z_s[:, [0]] > z_mesh) | (z_mesh >= Z_s[:, [-1]]))[:, ::-1].T
    
    temp_mesh = np.array([np.interp(x=np.linspace(z_min, z_max, y_points + 1),
                                    xp=np.concatenate((Z_i_line, (Z_s_line[1:] if has_snow else []))),
                                    fp=np.concatenate(([T_oi], T_i_line, [T_is],
                                                       (T_s_line if has_snow else []),
                                                       ([T_sa] if has_snow else []))),
                                    left=0.0, right=0.0)\
                          for Z_i_line, Z_s_line, T_oi, T_i_line, T_is, T_s_line, T_sa, has_snow\
                          in zip(Z_i, Z_s, temp_oi, temp_ice, temp_is, temp_snow, temp_sa, has_snow_filter)]
                        )[:, ::-1].T

    ice_masked = np.ma.masked_array(temp_mesh, ice_filter)
    snow_masked = np.ma.masked_array(temp_mesh, snow_filter)
    empty_mask = np.ma.masked_array(temp_mesh, ice_filter & snow_filter)

    fig = plt.figure(figsize=(15, 17))
    ax = fig.add_subplot()
    curr_cmap = mcm.get_cmap().copy()
    curr_cmap.set_bad(color_empty)
    img_empty = ax.imshow(empty_mask, cmap=curr_cmap)
    img_ice = ax.imshow(ice_masked, cmap=cmap_ice)
    img_snow = ax.imshow(snow_masked, cmap=cmap_snow)
    ax.set_xlabel('t, hours', size=20)
    ax.set_ylabel('Z, m', size=20)
    ax.set_xticks(np.linspace(0, len(time_arr) - 1, x_ticks, dtype=int))
    ax.set_yticks(np.linspace(0, y_points, y_ticks, dtype=int))
    ax.set_xticklabels(round(time_arr[int(tick)]/3600, 2) for tick in ax.get_xticks())
    ax.set_yticklabels(np.round(np.linspace(z_max, z_min, y_ticks), 2))
    ax.axhline(y_points*z_max/(z_max-z_min), lw=3, ls='--', color='c')
    ax.tick_params(axis='both', labelsize=15)
    
    if np.array(has_snow_filter).any():
        cbar_snow = fig.colorbar(img_snow, orientation='horizontal', pad=-0.01)
        cbar_snow.ax.set_xlabel(r'T_snow, $^o C$', size=20)
        cbar_snow.ax.tick_params(labelsize=15)
        
    cbar_ice = fig.colorbar(img_ice, orientation='horizontal', pad=0.08)
    cbar_ice.ax.set_xlabel(r'T_ice, $^o C$', size=20)
    cbar_ice.ax.tick_params(labelsize=15)
    
    if savepath is not None:
        fig.savefig(savepath, bbox_inches='tight')
