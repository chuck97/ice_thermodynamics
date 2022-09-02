import numpy as np
import scipy.ndimage as spndim
import pandas as pd
import datetime as dt
import itertools
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as anim
import matplotlib.cm as mcm
import matplotlib.colors as clr
import matplotlib.ticker as mtk
import matplotlib.dates as mdt
from matplotlib.collections import LineCollection
matplotlib.rcParams['animation.embed_limit'] = 5000
plt.rcParams["animation.html"] = "jshtml"

from utils.engine import rho_w as water_density, rho_s as snow_density


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


def get_bounds(vmin, vmax, step):
    
    if step is None:
        step = (vmax - vmin)*.1
        
    bounds = np.arange(vmin, vmax, step)
    if bounds[-1] != vmax:
        bounds = np.concatenate((bounds, [vmax]))
        
    return bounds


def discretize_cmap(cmap, bounds=None, vmin=None, vmax=None, step=None):
    
    assert bounds is not None or (vmin is not None and vmax is not None), \
           "Either bounds or limits should be set!"
    
    cmap_base = mcm.get_cmap(cmap)
    
    if bounds is None:
        bounds = get_bounds(vmin, vmax, step)
            
    norm_custom = clr.BoundaryNorm(bounds, cmap_base.N, extend='both')
    
    return cmap_base, norm_custom, bounds


def gauss_filter_with_nans(arr, sigma):
    
    nan_msk = np.isnan(arr)

    loss = np.zeros(arr.shape)
    loss[nan_msk] = 1
    loss = spndim.gaussian_filter(
            loss, sigma=sigma, mode='constant', cval=1)

    gauss = arr.copy()
    gauss[nan_msk] = 0
    gauss = spndim.gaussian_filter(
            gauss, sigma=sigma, mode='constant', cval=0)
    gauss[nan_msk] = np.nan

    gauss += loss * arr

    return gauss


def animate(processes,
            rho_snow=snow_density, rho_water=water_density,
            figsize=(25, 20),
            t_min=None, t_max=None,
            cmap=None, names=None,
            interval=50, savepath=None, dpi=None):
    
    assert len(set(process.get_length() for process in processes)) == 1, "Length of processes are not equal!"
    
    if names is not None: 
        assert len(names) == len(processes), "Number of names and number of processes should be equal!"
    
    def run(data):

        frame_number, (Z_ice_line_procs, Z_snow_line_procs, proc_data) = data
        
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
            ax.set_title('Time: {}'.format(time), size=25)
    #     ax.set_yticks(np.insert(-process.ice_dz_history.cumsum(), 0, 0))
    
        print("Rendered {:.0%} ({}/{}) frames".format(frame_number/frames_count, frame_number, frames_count), end="\r")

        return lines_ice + lines_snow + markers_ice
    
    frames_count = processes[0].get_length()
    
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

    fig, axes = plt.subplots(nrows=len(processes) + 1, figsize=figsize, \
                             gridspec_kw={"height_ratios":[1] + [0.05]*len(processes)})
    ax = axes[0]
    ax.set_xlim(t_min - (t_max - t_min)*0.1, t_max + (t_max - t_min)*0.1)
    ax.set_ylim(z_min - (z_max - z_min)*0.1, z_max + (z_max - z_min)*0.1)
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
    ax.yaxis.set_major_formatter(mtk.FormatStrFormatter('%.2f'))
    ax.tick_params(labelsize=15)
    ax.grid()
    
    for i, (cax, line_ice) in enumerate(zip(axes[1:], lines_ice)):
        cbar = fig.colorbar(line_ice, cax=cax, orientation='horizontal')
        if names is None:
            cbar.ax.set_xlabel('process %d'%i, size=20)
        else:
            cbar.ax.set_xlabel(names[i], size=20)
        cbar.ax.tick_params(labelsize=15)
        
    animation = anim.FuncAnimation(fig, run, enumerate(zip(zip(*[Z_i for Z_i in all_Z_i]),
                                                           zip(*[Z_s for Z_s in all_Z_s]),
                                                           zip(*[process.get_zip() for process in processes])),
                                                       start=1),
                                   save_count=frames_count, interval=interval, blit=True)
    
    if savepath:
        animation.save(savepath, dpi=dpi)
        
    return animation


def timeseries_img(process, rho_snow=snow_density, rho_water=water_density,
                   figsize=(30, 10), y_points=100,
                   mode='hours', year=1997, x_ticks=None, y_ticks=None,
                   tmin_ice=None, tmax_ice=None, step_ice=None, bounds_ice=None,
                   tmin_snow=None, tmax_snow=None, step_snow=None, bounds_snow=None,
                   cmap_ice='Blues', cmap_snow='Greys', color_waterline='c', color_empty=[208, 245, 226],
                   savepath=None, dpi=200):

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
    cmap_ice, norm_ice, bounds_ice = discretize_cmap(cmap_ice, bounds_ice, tmin_ice, tmax_ice, step_ice)
    img[ice_filter] = cmap_ice(norm_ice(temp_mesh[ice_filter]))

    if process.snow_presence_history.any():
        if bounds_snow is None:
            if tmin_snow is None:
                tmin_snow = np.nanmin(process.snow_temp_history)
            if tmax_snow is None:
                tmax_snow = np.nanmax(process.snow_temp_history)
        cmap_snow, norm_snow, bounds_snow = discretize_cmap(cmap_snow, bounds_snow, tmin_snow, tmax_snow, step_snow)
        img[snow_filter] = cmap_snow(norm_snow(temp_mesh[snow_filter]))
        
    fig = plt.figure(figsize=figsize)
    ax = plt.axes()
    
    if mode == 'month':
        daystart = (dt.datetime(year, 1, 1) - dt.datetime(1970, 1, 1)).days
        T_axis = daystart + process.timeline
        ax.xaxis.set_major_locator(mdt.MonthLocator())
        ax.xaxis.set_major_formatter(mdt.DateFormatter('%b'))
        xlabel = 'month'
        
    else:
        if x_ticks is not None:
            ax.xaxis.set_major_locator(mtk.LinearLocator(x_ticks))
            ax.xaxis.set_major_formatter(mtk.FormatStrFormatter("%d"))
        if mode == 'hours':
            T_axis = process.timeline*24
            xlabel = 'Time, h.'
        elif mode == 'days':
            T_axis = process.timeline
            xlabel = 'Time, days'
        else:
            raise Exception("invalid mode!")
    
    if y_ticks is not None:
        ax.yaxis.set_major_locator(mtk.LinearLocator(y_ticks))
        ax.yaxis.set_major_formatter(mtk.FormatStrFormatter("%.2f"))
    
    ax.imshow(img.swapaxes(0, 1), aspect='auto', origin='lower',
              extent=[T_axis[0], T_axis[-1], Z_mesh[0], Z_mesh[-1]])
    ax.set_xlabel('t, hours', size=25)
    ax.set_ylabel('Z, m', size=25)
    ax.tick_params(labelsize=20)
    ax.axhline(lw=3, ls='--', color=color_waterline)

    cax_ice = fig.add_axes([ax.get_position().x0, ax.get_position().y0 - 0.15,
                            ax.get_position().width, 0.05])
    fig.colorbar(mcm.ScalarMappable(norm_ice, cmap_ice), cax=cax_ice,
                 orientation='horizontal', ticks=bounds_ice)
    cax_ice.set_xlabel(r'Ice temperature, $^{\circ}$C', size=25)
    cax_ice.tick_params(axis='x', labelsize=20)

    if process.snow_presence_history.any():
        cax_snow = fig.add_axes([ax.get_position().x0, ax.get_position().y0 - 0.3,
                                 ax.get_position().width, 0.05])
        fig.colorbar(mcm.ScalarMappable(norm_snow, cmap_snow), cax=cax_snow,
                     orientation='horizontal', ticks=bounds_snow)
        cax_snow.set_xlabel(r'Snow temperature, $^{\circ}$C', size=25)
        cax_snow.tick_params(axis='x', labelsize=20)
    
    if savepath is not None:
        fig.savefig(savepath, bbox_inches='tight', dpi=dpi)
        
        
def timeseries_err(process_sim, process_data,
                   rho_snow=snow_density, rho_water=water_density,
                   figsize=(30, 10), y_points=100,
                   mode='hours', year=1997, x_ticks=None, y_ticks=None,
                   label_data='data', label_sim='simulation',
                   legend_loc='best',
                   tmin_err=None, tmax_err=None, step_err=None,
                   levels_fill=None, levels_border = [-2, -1, -0.5, 0.5, 1, 2],
                   sigma_x=0, sigma_y=None,
                   cmap='seismic', rgb_background = [208, 245, 226],
                   savepath=None, dpi=200):

    Z_i_data, Z_s_data = get_Z(process_data, rho_snow, rho_water)
    Z_i_sim, Z_s_sim = get_Z(process_sim, rho_snow, rho_water)
    z_min = min(np.concatenate((Z_i_sim[:, 0], Z_i_data[:, 0])))
    z_max = max(np.concatenate((Z_s_sim[:, -1], Z_s_data[:, -1])))

    Z_mesh = np.linspace(z_min, z_max, y_points + 1)
    temp_mesh_data = np.array([np.interp(x=Z_mesh,
                                         xp=np.concatenate((Z_i_line, (Z_s_line[1:] if has_snow else []))),
                                         fp=np.concatenate(([T_oi], T_i_line, [T_is],
                                                           (T_s_line if has_snow else []),
                                                           ([T_sa] if has_snow else [])
                                                           )),
                                         left=np.nan, right=np.nan
                                         ) \
                          for Z_i_line, Z_s_line, T_oi, T_i_line, T_is, T_s_line, T_sa, has_snow \
                          in zip(Z_i_data, Z_s_data,
                                 process_data.oi_temp_history, process_data.ice_temp_history,
                                 process_data.is_temp_history, process_data.snow_temp_history,
                                 process_data.sa_temp_history, process_data.snow_presence_history
                                )
                               ])
    temp_mesh_sim = np.array([np.interp(x=Z_mesh,
                                        xp=np.concatenate((Z_i_line, (Z_s_line[1:] if has_snow else []))),
                                        fp=np.concatenate(([T_oi], T_i_line, [T_is],
                                                           (T_s_line if has_snow else []),
                                                           ([T_sa] if has_snow else [])
                                                           )),
                                        left=np.nan, right=np.nan
                                        ) \
                          for Z_i_line, Z_s_line, T_oi, T_i_line, T_is, T_s_line, T_sa, has_snow \
                          in zip(Z_i_sim, Z_s_sim,
                                 process_sim.oi_temp_history, process_sim.ice_temp_history,
                                 process_sim.is_temp_history, process_sim.snow_temp_history,
                                 process_sim.sa_temp_history, process_sim.snow_presence_history
                                )
                              ])

    data = (temp_mesh_sim - temp_mesh_data).T
        
    if levels_fill is None:
        if tmin_err is None:
            tmin_err = np.nanmin(data)
        if tmax_err is None:
            tmax_err = np.nanmax(data)
        levels_fill = get_bounds(tmin_err, tmax_err, step_err)
        
    max_abs_err = max(abs(tmin_err), tmax_err)
        
    if sigma_y is None:
        sigma_y = data.shape[1]/500.0

    fig = plt.figure(figsize=figsize)
    ax = plt.axes()
    
    if mode == 'month':
        daystart = (dt.datetime(year, 1, 1) - dt.datetime(1970, 1, 1)).days
        T_axis = daystart + process_data.timeline
        ax.xaxis.set_major_locator(mdt.MonthLocator())
        ax.xaxis.set_major_formatter(mdt.DateFormatter('%b'))
        xlabel = 'month'
        
    else:
        if x_ticks is not None:
            ax.xaxis.set_major_locator(mtk.LinearLocator(x_ticks))
            ax.xaxis.set_major_formatter(mtk.FormatStrFormatter("%d"))
        if mode == 'hours':
            T_axis = process_data.timeline*24
            xlabel = 'Time, h.'
        elif mode == 'days':
            T_axis = process_data.timeline
            xlabel = 'Time, days'
        else:
            raise Exception("invalid mode!")
    
    if y_ticks is not None:
        ax.yaxis.set_major_locator(mtk.LinearLocator(y_ticks))
        ax.yaxis.set_major_formatter(mtk.FormatStrFormatter("%.2f"))
        
    ax.set_facecolor(np.array(rgb_background)/256.0)
    contourf = ax.contourf(gauss_filter_with_nans(data, (sigma_x, sigma_y)),
                           levels=levels_fill, cmap=cmap, norm=clr.Normalize(-max_abs_err, max_abs_err),
                           extend='both', extent=[T_axis[0], T_axis[-1], Z_mesh[0], Z_mesh[-1]]
                          )
    contour = ax.contour(gauss_filter_with_nans(data, (sigma_x, sigma_y)),
                         levels=levels_border, cmap=cmap, norm=clr.Normalize(-max_abs_err/3, max_abs_err/3, clip=True),
                         linewidths=2.5, linestyles='--',
                         extent=[T_axis[0], T_axis[-1], Z_mesh[0], Z_mesh[-1]]
                        )
    ax.clabel(contour, fmt='%1.1f', fontsize=25)

    # ax.axhline(lw=3, ls='--', color='c')
    ax.plot(T_axis, Z_i_data[:, 0], lw=1.5, color='black')
    ax.plot(T_axis, Z_i_data[:, -1], lw=1.5, color='black')
    ax.plot(T_axis, Z_s_data[:, -1], lw=1.5, color='black', label=label_data)
    ax.plot(T_axis, Z_i_sim[:, 0], lw=3, ls=':', color='black')
    ax.plot(T_axis, Z_i_sim[:, -1], lw=3, ls=':', color='black')
    ax.plot(T_axis, Z_s_sim[:, -1], lw=3, ls=':', color='black', label=label_sim)

    ax.set_xlabel(xlabel, size=25)
    ax.set_ylabel('Z, m', size=25)
    ax.tick_params(labelsize=20)

    cax = fig.add_axes([ax.get_position().x0, ax.get_position().y0 - 0.15,
                        ax.get_position().width, 0.05])
    colorbar = fig.colorbar(contourf, cax=cax, orientation='horizontal', ticks=contourf.levels, drawedges=True)
    colorbar.outline.set_linewidth(3)
    colorbar.dividers.set_linewidth(3)
    colorbar.dividers.set_dashes((-0.5, (2.5, 6.5)))
    cax.tick_params(axis='x', labelsize=20)
    cax.set_xlabel(r'Temperature error, $^{\circ}$C', size=25)
    ax.legend(loc=legend_loc, prop={'size':20})
    
    if savepath is not None:
        fig.savefig(savepath, bbox_inches='tight', dpi=dpi)
        
        
def plot_characteristics(process_list, labels, savepath=None):
    fig, (ax_ice_th, ax_snow_th, ax_Tsu) = plt.subplots(nrows=3, figsize=(15, 35))
    
    for process, label in zip(process_list, labels):
        ax_ice_th.plot(process.timeline, process.ice_dz_history.sum(axis=1), label=label, lw=2.5)
        ax_snow_th.plot(process.timeline, process.snow_dz_history.sum(axis=1), label=label, lw=2.5)
        ax_Tsu.plot(process.timeline,
                    np.where(process.snow_presence_history, process.sa_temp_history, process.is_temp_history),
                    label=label, lw=2.5)
    
    ax_ice_th.set_title('Ice thickness', size=25)
    ax_ice_th.set_ylabel('Level, m.', size=20)
    ax_ice_th.set_xlabel('Time', size=20)
    ax_ice_th.tick_params(labelsize=15)
    ax_ice_th.legend(prop={'size': 20})
    ax_ice_th.grid()
    
    ax_snow_th.set_title('Snow thickness', size=25)
    ax_snow_th.set_ylabel('Level, m.', size=20)
    ax_snow_th.set_xlabel('Time', size=20)
    ax_snow_th.tick_params(labelsize=15)
    ax_snow_th.legend(prop={'size': 20})
    ax_snow_th.grid()
    
    ax_Tsu.set_title('Surface temperature', size=25)
    ax_Tsu.set_ylabel(r'Temperature, $^o C.$', size=20)
    ax_Tsu.set_xlabel('Time', size=20)
    ax_Tsu.tick_params(labelsize=15)
    ax_Tsu.legend(prop={'size': 20})
    ax_Tsu.grid()
    
    if savepath is not None:
        fig.savefig(savepath, bbox_inches='tight')
    else:
        plt.show()
        
        
def plot_errors(process_data, sim_process_list, labels, savepath=None):
    assert len(sim_process_list) == len(labels) - 1, \
    "Length of labels should be more than length of simulation processes list by one!"
    
    fig, (ax_ice_err, ax_snow_err, ax_Tsu_err) = plt.subplots(nrows=3, figsize=(15, 35))
    
    for process, label in zip(sim_process_list, labels):
        ax_ice_err.plot(process_data.timeline, process_data.ice_dz_history.sum(axis=1) \
                                             - process.ice_dz_history.sum(axis=1), 
                        label=label, lw=2.5)
        ax_snow_err.plot(process_data.timeline, process_data.snow_dz_history.sum(axis=1) \
                                              - process.snow_dz_history.sum(axis=1), 
                        label=label, lw=2.5)
        ax_Tsu_err.plot(process.timeline,
                    np.where(process_data.snow_presence_history,
                             process_data.sa_temp_history, process_data.is_temp_history) \
                  - np.where(process.snow_presence_history,
                             process.sa_temp_history, process.is_temp_history),
                    label=label, lw=2.5)
    
    ax_ice_err.set_title('Ice thickness error', size=25)
    ax_ice_err.set_ylabel('Level, m.', size=20)
    ax_ice_err.set_xlabel('Time', size=20)
    ax_ice_err.tick_params(labelsize=15)
    ax_ice_err.legend(prop={'size': 20})
    ax_ice_err.grid()
    
    ax_snow_err.set_title('Snow thickness error', size=25)
    ax_snow_err.set_ylabel('Level, m.', size=20)
    ax_snow_err.set_xlabel('Time', size=20)
    ax_snow_err.tick_params(labelsize=15)
    ax_snow_err.legend(prop={'size': 20})
    ax_snow_err.grid()
    
    ax_Tsu_err.set_title('Surface temperature error', size=25)
    ax_Tsu_err.set_ylabel(r'Temperature, $^o C.$', size=20)
    ax_Tsu_err.set_xlabel('Time', size=20)
    ax_Tsu_err.tick_params(labelsize=15)
    ax_Tsu_err.legend(prop={'size': 20})
    ax_Tsu_err.grid()
    
    if savepath is not None:
        fig.savefig(savepath, bbox_inches='tight')
    else:
        plt.show()

