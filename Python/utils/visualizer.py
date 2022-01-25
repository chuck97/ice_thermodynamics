import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as anim
from matplotlib.ticker import FormatStrFormatter
from matplotlib.collections import LineCollection, PathCollection
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.cm import ScalarMappable
plt.rcParams["animation.html"] = "jshtml"
matplotlib.rcParams['animation.embed_limit'] = 2000
import itertools


def animate(dz_arr, dz_snow_arr,
            temp_arr, temp_snow,
            time_arr,
            rho_ice_arr, rho_water, rho_snow,
            t_min=None, t_max=None, cmap=None, savepath=None):
    
    dz_arr = np.array(dz_arr)
    dz_snow_arr = np.array(dz_snow_arr)
    temp_arr = np.array(temp_arr)
    temp_snow = np.array(temp_snow)
    time_arr = np.array(time_arr)
    rho_ice_arr = np.array(rho_ice_arr)
    
    def run(data):

        dz_line, dz_snow_line, T_line, T_snow, time, rho_ice_line, rho_water = data

        Z = (-dz_line).cumsum() + dz_line/2
        Z = np.insert(Z, 0, 0)
        Z = np.append(Z, -dz_line.sum())
        Z = Z[::-1] - Z[-1] - np.dot(rho_ice_line, dz_line)/rho_water - rho_snow*sum(dz_snow_line)/rho_water
        points = np.array([T_line, Z]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)
        
        Z_snow = dz_snow_line.cumsum() - dz_snow_line/2
        Z_snow = np.insert(Z_snow, 0, 0)
        Z_snow = np.append(Z_snow, dz_snow_line.sum())
        Z_snow = Z_snow + Z[-1]
        
        line.set_segments(segments)
        line.set_array(T_line)
        line_snow.set_data(np.insert(T_snow, 0, T_line[-1]), Z_snow)
        markers.set_offsets([[T, z] for T, z in zip(T_line, Z)])
        markers.set_array(T_line)
        ax.set_title('Time: %.2f hours'%(time/3600), size=20)
    #     ax.set_yticks(np.insert(-dz_arr.cumsum(), 0, 0))

        return line, line_snow, markers
    
    # потом исправить при динамической плотности (формула будет (\int_0^h \rho(x) dx)/rho_water)
#     rho_ice = rho_ice_list[0]
    
#     h_max = max([sum(proc_line[0]) for proc_line in proc_list])
#     z_min = -h_max*rho_ice/rho_water
#     z_max = z_min + h_max
    z_mins = np.array([(-np.dot(rho_ice_line, dz_line) - rho_snow*sum(dz_snow_line))/rho_water\
                      for rho_ice_line, dz_line, dz_snow_line in zip(rho_ice_arr, dz_arr, dz_snow_arr)])
    z_min = min(z_mins)
    z_max = max(z_mins + dz_arr.sum(axis=1) + dz_snow_arr.sum(axis=1))
    temp_min = min(temp_arr.min(), temp_snow.min())
    temp_max = max(temp_arr.max(), temp_snow.max())
     
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
    line = ax.add_collection(lc)
    line_snow, = ax.plot([], [], color='grey', lw=3, marker='o', animated=True)
    markers = ax.scatter([], [], s=60, animated=True)
    markers.set_cmap(cmap)
    ax.axhline(ls='--', lw=3, color='navy')
    ax.set_xlabel(r'T, $^o C$', size=20)
    ax.set_ylabel('Z, m', size=20)
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax.tick_params(axis='both', labelsize=15)
    ax.grid()
    fig.colorbar(line)
    
    print(dz_arr.shape, dz_snow_arr.shape, temp_arr.shape, temp_snow.shape, time_arr.shape, rho_ice_arr.shape)
    animation = anim.FuncAnimation(fig, run, zip(dz_arr, dz_snow_arr,
                                                 temp_arr, temp_snow,
                                                 time_arr,
                                                 rho_ice_arr, [rho_water]*dz_arr.shape[0]),
                                   interval=30, blit=True)
    
    if savepath:
        animation.save(savepath)
        
    else:
        return animation


def timeseries_img(dz_arr, dz_snow_arr,
                   temp_arr, temp_snow,
                   time_arr,
                   rho_ice_arr, rho_water, rho_snow,
                   y_points=100, x_ticks=11, y_ticks=6, t_min=None, t_max=None, savepath=None):

    dz_arr = np.array(dz_arr)
    dz_snow_arr = np.array(dz_snow_arr)
    temp_arr = np.array(temp_arr)
    temp_snow = np.array(temp_snow)
    timeline = np.array(time_arr)
    rho_ice_arr = np.array(rho_ice_arr)

    z_mins = np.array([(-np.dot(rho_ice_line, dz_line) - rho_snow*sum(dz_snow_line))/rho_water\
                      for rho_ice_line, dz_line, dz_snow_line in zip(rho_ice_arr, dz_arr, dz_snow_arr)])
    z_min = min(z_mins)
    z_max = max(z_mins + dz_arr.sum(axis=1) + dz_snow_arr.sum(axis=1))
    temp_min = min(T_arr.min(), T_snow.min())
    temp_max = max(T_arr.max(), T_snow.max())
    
    if cmap is None:
        
        cdict = {'red':  [[0.0,  0.0, 0.0],
                  [1.0,  0.6, 0.6]],
         'green':  [[0.0,  0.2, 0.2],
                    [1.0,  0.85, 0.85]],
         'blue':  [[0.0,  0.5, 0.5],
                   [1.0,  1.0, 1.0]]}
        cmap = LinearSegmentedColormap('MyCmap', segmentdata=cdict)
        cmap.set_over('black')
        cmap.set_under('black')

    
    Z = dz_arr.cumsum(axis=1) - dz_arr/2
    Z = np.insert(Z, 0, 0, axis=1)
    Z = np.append(Z, dz_arr.sum(axis=1).reshape(-1, 1), axis=1)
    Z += ((rho_ice_arr/rho_water - 1)*Z[:, -1]).reshape(-1, 1)

    proc_image = np.array([np.interp(x=np.linspace(z_min, z_max, y_points+1),
                                     xp=Z_line, fp=T_line, left=0.0, right=0.0)\
                           for Z_line, T_line in zip(Z, T_arr)]
                         )[:, ::-1].T

    fig = plt.figure(figsize=(15, 10))
    ax = fig.add_subplot()
    img = ax.imshow(proc_image, cmap=cmap,
                    vmin=(t_min if t_min is not None else temp_min),
                    vmax=(t_max if t_max is not None else temp_max))
    ax.set_xlabel('t, hours', size=20)
    ax.set_ylabel('Z, m', size=20)
    ax.set_xticks(np.linspace(0, len(timeline) - 1, x_ticks, dtype=int))
    ax.set_yticks(np.linspace(0, y_points, y_ticks, dtype=int))
    ax.set_xticklabels(round(timeline[int(tick)]/3600, 2) for tick in ax.get_xticks())
    ax.set_yticklabels(np.round(np.linspace(z_max, z_min, y_ticks + 1), 2))
    ax.axhline(y_points*z_max/dz_arr.sum(axis=1).max(), lw=3, ls='--', color='c')
    ax.tick_params(axis='both', labelsize=15)
    cbar = fig.colorbar(img, orientation='horizontal', pad=0.11)
    cbar.ax.set_xlabel(r'T, $^o C$', size=20)
    cbar.ax.tick_params(labelsize=15)
    
    if savepath is not None:
        fig.savefig(savepath, bbox_inches='tight')
        
def test():
    print(6)