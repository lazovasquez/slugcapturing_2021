#!/usr/bin/env python3

import json
import argparse
import numpy as np
import h5py
from mayavi import mlab
import matplotlib.pyplot as plt
import subprocess
import os

# extra
from json_extract import GetValue2
# warnings.filterwarnings('ignore')

if __name__ == '__main__':

    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='Erosion and deposition condition visualization')
    parser.add_argument('out_folder',
                        action='store',
                        metavar='OUT_FOLDER',
                        type=str,
                        help='Directory containing output files.')
    # parser.add_argument('--time_steps',
    #                     action='store',
    #                     metavar='TIME_STEPS',
    #                     type=int,
    #                     required=True,
    #                     help='Time steps [-].')
    # parser.add_argument('--time_step',
    #                     action='store',
    #                     metavar='TIME_STEP',
    #                     type=float,
    #                     required=True,
    #                     help='Time step [s].')
    # parser.add_argument('--simulation_time',
    #                     action='store',
    #                     metavar='SIMULATION_TIME',
    #                     type=float,
    #                     required=True,
    #                     help='Simulation time [s].')
    parser.add_argument('--axis',
                        action='store',
                        choices=['x', 'y', 'z'],
                        metavar='FLOW_AXIS',
                        type=str,
                        required=True,
                        help='Flow axis.')
    parser.add_argument('--show',
                        action='store_true',
                        default=False,
                        help='Show plot instead of saving to file.')
    parser.add_argument('--kappa_er',
                        action='store',
                        metavar='kappa_er',
                        type=float,
                        help='Erodibility coefficient [s/m].',
                        default=0.017)
    parser.add_argument('--rho_s',
                        action='store',
                        metavar='rho_s',
                        type=float,
                        help='Solid density of calcite[kg/m3].',
                        default=2710)
    parser.add_argument('--liquid_viscosity',
                        action='store',
                        metavar='liquid_viscosity',
                        type=float,
                        help='Liquid viscosity [Pa s].',
                        default=1.002e-3)
    parser.add_argument('--reaction_time',
                        action='store',
                        metavar='reaction_time',
                        type=float,
                        help='Reaction time [s].',
                        default=2) #1.7e-6
    parser.add_argument('--vmin',
                        action='store',
                        metavar='vmin_scale',
                        type=float,
                        help='Minimum scale value * 1e-6[m].',
                        default=None)
    parser.add_argument('--vmax',
                        action='store',
                        metavar='vmax_scale',
                        type=float,
                        help='Maximum scale value * 1e-6[m].',
                        default=None)
    arg = parser.parse_args()

    # Load centrelines input file
    with open(arg.out_folder + '/centerlines.json', mode='r') as file1:
        data = json.load(file1, encoding='utf-8')

    # Import parameters from config.json
    with open(arg.out_folder+'/config.json', mode='r') as file2:
        config_parameters = json.load(file2, encoding='utf-8')

    # Import absolute pressure
    # absolute_pressure = GetValue2(config_parameters).get_values('absolute_pressure')
    voxel_size = GetValue2(config_parameters).get_values('voxel_size')

    # Extract node geometry arrays from JSON
    nodes = sorted(data['graph']['nodes'], key=lambda node: int(node['id']))
    x = np.array([node['metadata']['node_coordinates']['x'] for node in nodes])
    y = np.array([node['metadata']['node_coordinates']['y'] for node in nodes])
    z = np.array([node['metadata']['node_coordinates']['z'] for node in nodes])

    # Extract link geometry arrays from JSON
    edges = sorted(data['graph']['edges'], key=lambda edge: int(edge['id']))
    source = np.array([int(edge['source']) for edge in edges])
    target = np.array([int(edge['target']) for edge in edges])
    link_squared_radius = np.array([edge['metadata']['link_squared_radius'] for edge in edges])
    link_radius = voxel_size*np.sqrt(link_squared_radius)

    link_length = voxel_size*np.array([edge['metadata']['link_length'] for edge in edges])

    with h5py.File(arg.out_folder + '/static_results.h5', "r") as file:
        # Load flow rate input file
        Q = np.squeeze(np.array(file['flow_rate_' + arg.axis], dtype=np.double))
        # q = np.reshape(Q, (3, int(Q.size/3)))
        # q = Q * 1e12  # convert to nL/s

        # Load flow speed input file
        V = np.squeeze(np.array(file['flow_speed_' + arg.axis], dtype=np.double))
        # v = V * 1e3  # convert to mm/s

        P = np.squeeze(np.array(file['pressures_' + arg.axis], dtype=np.double))

    # Calculate node origins (X, Y, Z) and node-to-node vectors (U[0], U[1], U[2])
    X = x[source]
    Y = y[source]
    Z = z[source]
    U = np.array([x[target] - x[source],
                  y[target] - y[source],
                  z[target] - z[source]])

    links_num = GetValue2(data).get_values('number_of_links')

    # Compute pressure gradient for each link
    pressure_gradients = np.linspace (0, links_num, num=links_num)
    for i in range(links_num):
        pressure_gradients[i] = (P[source[i]]-P[target[i]])+np.finfo(float).eps

    # Shear stresses
    tau_w = pressure_gradients*link_radius/(2*link_length)
    # tau_w = tau_w[0]
    # print(tau_w)
    tau_er = 2*arg.rho_s*arg.kappa_er*link_length/pressure_gradients
    erosion_rate = -arg.kappa_er*(tau_w-tau_er)

    # Create vectors
    cond_erosion = np.linspace (0, links_num, num=links_num)
    erosion_rate = np.linspace (0, links_num, num=links_num)
    t_er = np.linspace (0, links_num, num=links_num)
    modified_link_radius = np.linspace (0, links_num, num=links_num)

    for i in range(links_num):
        if tau_w[i] > tau_er[i]:
            cond_erosion[i] = 1
            erosion_rate[i] = -arg.kappa_er*(tau_w[i]-tau_er[i])
            t_er[i] = 2*arg.rho_s*link_length[i]/(arg.kappa_er*pressure_gradients[i])
            modified_link_radius[i] = link_radius[i]*np.exp(arg.reaction_time/t_er[i])
        else:
            cond_erosion[i] = 0
            erosion_rate[i] = 0
            modified_link_radius[i] = link_radius[i]

    link_radius_variation = modified_link_radius - link_radius
    link_radius_rate = modified_link_radius / link_radius

    # Print
    print('FLOWSIMULATOR::CMM SAYS:')

#     # POTENTIAL EROSION
#     # POINTS Creating 3D visualization for potential erosion along the flow axis
#     mlab.figure(size=(800, 700), bgcolor=(0.1, 0.1, 0.1))
#     mlab.points3d(X, Y, Z, cond_erosion, scale_mode='none', scale_factor=5, colormap='viridis', vmin=arg.vmin, vmax=arg.vmax)
#     mlab.outline()
#     mlab.orientation_axes()
#     mlab.scalarbar(title='Potential erosion [-]', orientation='vertical')
#     if (arg.show):
#         mlab.show()
#     else:
#         mlab.savefig(arg.out_folder + '/erosion_deposition_condition_plot_' + arg.axis + '.png')

#     # Plot CONDITION EROSION distribution
#     plt.figure(dpi=125)
#     plt.hist(cond_erosion, bins=5)
#     plt.xlabel(r'Potential erosion (0=none, erosion=1)', fontsize=16)
#     plt.ylabel(r'Number of capillary voxels', fontsize=16)
#     if (arg.show):
#         plt.show()
#     else:
#         plt.savefig(arg.out_folder + '/' + arg.filename + '_dist_links.png', dpi=125, bbox_inches='tight')

#     # # Plot CONDITION PHYSICAL distribution
#     # plt.figure(dpi=125)
#     # plt.hist(cond_erosion, bins=5)
#     # plt.xlabel(r'Potential physical process (deposition=1, 0=none, erosion=-1)', fontsize=16)
#     # plt.ylabel(r'Number of capillary voxels', fontsize=16)
#     # if (arg.show):
#     #     plt.show()
#     # else:
#     #     plt.savefig(arg.out_folder + '/' + arg.filename + '_dist_links.png', dpi=125, bbox_inches='tight')

#     # Number of links
#     count_erosion = np.count_nonzero(cond_erosion == 1)
#     print(count_erosion, 'links erosion of', links_num,'links')

#    # POINTS Creating 3D visualization for erosion rate along the flow axis
#     mlab.figure(size=(800, 700), bgcolor=(0.1, 0.1, 0.1))
#     mlab.points3d(X, Y, Z, erosion_rate, scale_mode='none', scale_factor=5, colormap='viridis', vmin=arg.vmin, vmax=arg.vmax)
#     mlab.outline()
#     mlab.orientation_axes()
#     mlab.scalarbar(title=r'm$_{er}$ [kg m-2]', orientation='vertical')
#     if (arg.show):
#         mlab.show()
#     else:
#         mlab.savefig(arg.out_folder + '/erosion_rate_plot_' + arg.axis + '.png')

#     # Maximum erosion rate
#     loc_max_erosion = np.where(erosion_rate == np.amin(erosion_rate))
#     print('Maximum erosion rate = ', min(erosion_rate), "link number: ", loc_max_erosion, "link diameter = ", 2.0 * modified_link_radius[loc_max_erosion]*voxel_size / 1.0e-6, '[1e-6 m]') # locate link, and find data

#     # LINK SQUARED RADIUS
#    # POINTS Creating 3D visualization for link squared radius along the flow axis
#     mlab.figure(size=(800, 700), bgcolor=(0.1, 0.1, 0.1))
#     mlab.points3d(X, Y, Z, 2.0 * np.sqrt(link_squared_radius) * voxel_size / 1.0e-6, colormap='viridis', vmin=arg.vmin, vmax=arg.vmax)
#     mlab.outline()
#     mlab.orientation_axes()
#     mlab.scalarbar(title='D [1e-6 m]', orientation='vertical')
#     if (arg.show):
#         mlab.show()
#     else:
#         mlab.savefig(arg.out_folder + '/' + arg.filename + '_plot.png')

#     # Plot pore size distribution
#     plt.figure(dpi=125)
#     plt.hist(2.0 * np.sqrt(link_squared_radius) * voxel_size / 1.0e-6, bins=20)
#     plt.xlabel(r'Capillary voxel diameter [$\mu\mathrm{m}$]', fontsize=16)
#     plt.ylabel(r'Number of capillary voxels', fontsize=16)
#     if (arg.show):
#         plt.show()
#     else:
#         plt.savefig(arg.out_folder + '/' + arg.filename + '_dist_links.png', dpi=125, bbox_inches='tight')

#     # MODIFIED LINK SQUARED RADIUS
#    # POINTS Creating 3D visualization for modified link squared radius along the flow axis
#     mlab.figure(size=(800, 700), bgcolor=(0.1, 0.1, 0.1))
#     mlab.points3d(X, Y, Z, 2.0 * modified_link_radius / 1.0e-6, colormap='viridis', vmin=arg.vmin, vmax=arg.vmax)
#     mlab.outline()
#     mlab.orientation_axes()
#     mlab.scalarbar(title=r'D$_2$ [1e-6 m]', orientation='vertical')
#     if (arg.show):
#         mlab.show()
#     else:
#         mlab.savefig(arg.out_folder + '/' + arg.filename + '_plot.png')

#     # Plot pore size distribution
#     plt.figure(dpi=125)
#     plt.hist(2.0 * modified_link_radius * voxel_size / 1.0e-6, bins=20)
#     plt.xlabel(r'Capillary voxel modified diameter [$\mu\mathrm{m}$]', fontsize=16)
#     plt.ylabel(r'Number of capillary voxels', fontsize=16)
#     if (arg.show):
#         plt.show()
#     else:
#         plt.savefig(arg.out_folder + '/' + arg.filename + '_dist_links.png', dpi=125, bbox_inches='tight')

# #####
#    # POINTS Creating 3D visualization for link squared radius variation along the flow axis
#     mlab.figure(size=(800, 700), bgcolor=(0.1, 0.1, 0.1))
#     mlab.points3d(X, Y, Z, 2.0 * link_radius_variation / 1.0e-6, colormap='viridis', vmin=0, vmax=max(link_radius_variation / 1.0e-6))
#     mlab.outline()
#     mlab.orientation_axes()
#     mlab.scalarbar(title=r'Variation D$_1$-D$_0$ [1e-6 m]', orientation='vertical')
#     if (arg.show):
#         mlab.show()
#     else:
#         mlab.savefig(arg.out_folder + '/' + arg.filename + '_plot.png')

#     # Plot pore size distribution
#     plt.figure(dpi=125)
#     plt.hist(2.0 * link_radius_variation * voxel_size / 1.0e-6, bins=20)
#     plt.xlabel(r'Link radius variation [$\mu\mathrm{m}$]', fontsize=16)
#     plt.ylabel(r'Number of capillary voxels', fontsize=16)
#     if (arg.show):
#         plt.show()
#     else:
#         plt.savefig(arg.out_folder + '/' + arg.filename + '_dist_links.png', dpi=125, bbox_inches='tight')

# #####
#    # POINTS Creating 3D visualization for link radius variation rate along the flow axis
#     mlab.figure(size=(800, 700), bgcolor=(0.1, 0.1, 0.1))
#     mlab.points3d(X, Y, Z, link_radius_rate, colormap='viridis', vmin=0, vmax=max(link_radius_rate))
#     mlab.outline()
#     mlab.orientation_axes()
#     mlab.scalarbar(title=r'Variation D$_1$/D$_0$', orientation='vertical')
#     if (arg.show):
#         mlab.show()
#     else:
#         mlab.savefig(arg.out_folder + '/' + arg.filename + '_plot.png')

#     # Plot pore size distribution
#     plt.figure(dpi=125)
#     plt.hist(2.0 * np.sqrt(link_radius_rate) * voxel_size / 1.0e-6, bins=30)
#     plt.xlabel(r'Variation D$_1$/D$_0$', fontsize=16)
#     plt.ylabel(r'Number of capillary voxels', fontsize=16)
#     if (arg.show):
#         plt.show()
#     else:
#         plt.savefig(arg.out_folder + '/' + arg.filename + '_dist_links.png', dpi=125, bbox_inches='tight')

    # Iterative process

    # Create directories
    os.mkdir('results_cmm')
    time_steps_test = 3
    for i in range(time_steps_test):
        os.mkdir('results_cmm/t_'+str(i))

    # # Update link squared radius
    # updated_link_squared_radius = (modified_link_radius*voxel_size)**2

    # # Update dict for link squared radius
    # for i in range(link_squared_radius.size):
    #     centerlines['graph']['edges'][i]['metadata']['link_squared_radius']=link_squared_radius[i]

    # # Replace link_square_radius in centerlines.json file
    # with open(output_folder+'/'+output_centerlines, mode='w') as fp:
    #     json.dump(centerlines, fp, sort_keys=False, indent=4)

    # Run simulator
    subprocess.run(['ls'])
