# -*- coding: utf-8 -*-
"""
Extended by Máté Szedlák

Used sources:
https://stackoverflow.com/questions/29188612/arrows-in-matplotlib-using-mplot3d
https://gist.github.com/jpwspicer/ea6d20e4d8c54e9daabbc1daabbdc027
"""
from copy import deepcopy
try:
    from matplotlib import pyplot
    from matplotlib.patches import FancyArrowPatch
    from mpl_toolkits.mplot3d import proj3d
except ImportError:
    print("Graphical libraries could not be loaded. GUI can not be used.")


class Arrow3D(FancyArrowPatch):
    """
    Vector drawer module based on the matplotlib library using external sources, like:
        https://stackoverflow.com/questions/29188612/arrows-in-matplotlib-using-mplot3d
        https://gist.github.com/jpwspicer/ea6d20e4d8c54e9daabbc1daabbdc027
    """
    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0, 0), (0, 0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        _xs3d, _ys3d, _zs3d = self._verts3d
        _xs, _ys, _zs = proj3d.proj_transform(_xs3d, _ys3d, _zs3d, renderer.M)
        self.set_positions((_xs[0], _ys[0]), (_xs[1], _ys[1]))
        FancyArrowPatch.draw(self, renderer)

    def plot_structure(truss, original=1, result=1, supports=1,
                      forces=1, reactions=1, scale_displacements=1, scale_forces=1, z_correction=1, show_values=1, save_plot=1, log=0):
        """
        General plotting method for trusses

        scale_displacements: Scale drawing of displacements
        scale_forces:   Scale force sign
        z_correction:    Scale z-axis
        """
        plotname = truss.title
        plot_width = 10.0              # Plot width in inches
        xframe = 0.0                     # Frame width at X direction | 0: auto
        yframe = 0.0                     # Frame width at Y direction | 0: auto
        zframe = 0.0                     # Frame width at Z direction | 0: auto
        scale_sup = scale_forces * 0.3        # Scale support sign  # All the others are input parameters

        # Stress coloring settings [R G B] - Examples considering pressure:
        # 0: [1, 0, 0]              Plain red
        # 1: [x, 0, 0]              Red to Black
        # 2: [1, 1-x, 0-x]          Red to White
        # 3: [1, (1-x)/2, (1-x)/2]  Red to MildRed - Distinct pressure and tension
        # 4: [x, 1-x, 0]            Red to Green
        _coloring = 3  # € [0, 1, 2, 3, 4]

        fig = pyplot.figure()
        _ax = fig.add_subplot(111, projection='3d')

        if truss.dof == 2:
            _ax.view_init(elev=90., azim=-90.)
            _ax.w_zaxis.line.set_lw(0.)
            _ax.set_zticklabels([])

        xmin = min(list(truss.nodal_coord[x][0] for x in range(truss.number_of_nodes())))
        xmax = max(list(truss.nodal_coord[x][0] for x in range(truss.number_of_nodes())))
        ymin = min(list(truss.nodal_coord[x][1] for x in range(truss.number_of_nodes())))
        ymax = max(list(truss.nodal_coord[x][1] for x in range(truss.number_of_nodes())))
        zmin = min(list(truss.nodal_coord[x][2] for x in range(truss.number_of_nodes())))
        zmax = max(list(truss.nodal_coord[x][2] for x in range(truss.number_of_nodes())))

        deltax = xmax - xmin
        deltay = ymax - ymin

        if (xframe == 0):
            xframe = max(deltax * 0.05, 2)
        if (yframe == 0):
            yframe = max(deltay * 1.5, 2)
        if (zframe == 0):
            zframe = max(deltay * 1.5, 2)

        if truss.dof == 3:
            plot_height = plot_width * ((deltay + yframe*2)/(deltax + xframe*2)) * 0.3
        else:
            plot_height = plot_width * 0.5 * 0.5
        fig.set_size_inches(plot_width, plot_height)

        _ax.set_xlim3d(xmin - xframe, xmax + xframe)
        _ax.set_ylim3d(ymin - yframe, ymax + yframe)
        _ax.set_zlim3d(zmin - zframe, zmax + zframe)

        if original == result:
            _coloring = 0

        # Giving plot names
        if original and (not result) and supports and (not reactions):
            plotname += ' - 01 Initial structure'
            if forces:
                plotname += ' with forces'
        elif original and result:
            plotname += ' - 02 Deformation'
            if not reactions:
                plotname += ' with reactions'
        elif (not original) and result:
            plotname += ' - 03 Stresses'
            if not reactions:
                plotname += ' with reactions'
        else:
            plotname += ' - Unnamed'

        print(plotname + ": ")
        if result:
            dipslay_displacement = deepcopy(truss.nodal_coord_def)
            if scale_displacements != 1.0:
                if log:
                    print('Displacements are scaled with factor: ' + str(scale_displacements))
                for i in range(truss.number_of_nodes()):
                    for j in range(3):
                        dipslay_displacement[i][j] = (truss.nodal_coord_def[i][j] -
                                                      truss.nodal_coord[i][j]) * scale_displacements + truss.nodal_coord[i][j]

        for i in range(truss.number_of_elements()):
            # Plot original structure
            if original:
                _ax.plot([truss.nodal_coord[truss.nodal_connections[i][1]][0], truss.nodal_coord[truss.nodal_connections[i][0]][0]],
                         [truss.nodal_coord[truss.nodal_connections[i][1]][1], truss.nodal_coord[truss.nodal_connections[i][0]][1]],
                         zs=[truss.nodal_coord[truss.nodal_connections[i][1]][2], truss.nodal_coord[truss.nodal_connections[i][0]][2]], color='b')
            # Plot deformed structure
            if result:
                if True: #truss._post_processed:        TODO: post_processed state should be checked in a different way
                    if truss.stress_color[i] > 0:
                        if _coloring == 1:
                            rgb_col = [0, 0, abs(truss.stress_color[i])]
                        elif _coloring == 2:
                            rgb_col = [1-abs(truss.stress_color[i]),
                                       1-abs(truss.stress_color[i]), 1]
                        elif _coloring == 3:
                            rgb_col = [(1-abs(truss.stress_color[i]))/2,
                                       (1-abs(truss.stress_color[i]))/2, 1]
                        elif _coloring == 4:
                            rgb_col = [0, 1-abs(truss.stress_color[i]),
                                       abs(truss.stress_color[i])]
                        else:
                            rgb_col = [1, 0, 0]
                    else:
                        if _coloring == 1:
                            rgb_col = [abs(truss.stress_color[i]), 0, 0]
                        elif _coloring == 2:
                            rgb_col = [1, 1-abs(truss.stress_color[i]),
                                       1-abs(truss.stress_color[i])]
                        elif _coloring == 3:
                            rgb_col = [1, (1-abs(truss.stress_color[i]))/2,
                                       (1-abs(truss.stress_color[i]))/2]
                        elif _coloring == 4:
                            rgb_col = [abs(truss.stress_color[i]),
                                       1-abs(truss.stress_color[i]), 0]
                        else:
                            rgb_col = [1, 0, 0]
                else:
                    print('Stresses are not calculated')
                    rgb_col = [1, 0, 0]
                _ax.plot([dipslay_displacement[truss.nodal_connections[i][1]][0], dipslay_displacement[truss.nodal_connections[i][0]][0]],
                         [dipslay_displacement[truss.nodal_connections[i][1]][1], dipslay_displacement[truss.nodal_connections[i][0]][1]],
                         zs=[dipslay_displacement[truss.nodal_connections[i][1]][2], dipslay_displacement[truss.nodal_connections[i][0]][2]],
                         color=rgb_col)

        if forces:
            for force in truss.known_f_not_zero:
                if truss.force[force] < 0:
                    value = -1.0
                else:
                    value = 1.0
                if force % 3 == 0:
                    f_dir = [value*scale_forces, 0, 0]
                elif force % 3 == 1:
                    f_dir = [0, value*scale_forces, 0]
                else:
                    f_dir = [0, 0, value*scale_forces*z_correction]
                if original:
                    f_arrow = Arrow3D([truss.nodal_coord[force//3][0], truss.nodal_coord[force//3][0] + f_dir[0]],
                                      [truss.nodal_coord[force//3][1], truss.nodal_coord[force//3][1] + f_dir[1]],
                                      [truss.nodal_coord[force//3][2], truss.nodal_coord[force//3][2] + f_dir[2]],
                                      mutation_scale=20, lw=1, arrowstyle="-|>", color="k")
                else:
                    f_arrow = Arrow3D([truss.nodal_coord_def[force // 3][0], truss.nodal_coord_def[force // 3][0] + f_dir[0]],
                                      [truss.nodal_coord_def[force // 3][1], truss.nodal_coord_def[force // 3][1] + f_dir[1]],
                                      [truss.nodal_coord_def[force // 3][2], truss.nodal_coord_def[force // 3][2] + f_dir[2]],
                                      mutation_scale=20, lw=1, arrowstyle="-|>", color="k")
                _ax.add_artist(f_arrow)

        if reactions:
            e_previous = -100
            for reaction in truss.known_displacement_a:
                value = 0.0             # Maybe this is useless <XXX>
                if truss.force[reaction] < 0:
                    value = -1.0
                elif truss.force[reaction] > 0:
                    value = 1.0
                if reaction % 3 == 0:
                    f_dir = [value*scale_forces, 0, 0]
                elif reaction % 3 == 1:
                    f_dir = [0, value*scale_forces, 0]
                else:
                    f_dir = [0, 0, value*scale_forces*z_correction]
                if abs(truss.force[reaction]) > 0:
                    f_arrow = Arrow3D([truss.nodal_coord[reaction//3][0], truss.nodal_coord[reaction//3][0] + f_dir[0]],
                                      [truss.nodal_coord[reaction//3][1], truss.nodal_coord[reaction//3][1] + f_dir[1]],
                                      [truss.nodal_coord[reaction//3][2], truss.nodal_coord[reaction//3][2] + f_dir[2]],
                                      mutation_scale=20, lw=1, arrowstyle="-|>", color="darkolivegreen")
                    _ax.add_artist(f_arrow)
                    if show_values:
                        _ax.set_xticklabels([])
                        _ax.set_yticklabels([])
                        _ax.set_zticklabels([])
                        if not reaction//3 == e_previous//3:
                            if truss.dof == 3:
                                _ax.text(truss.nodal_coord[reaction//3][0],
                                         truss.nodal_coord[reaction//3][1],
                                         truss.nodal_coord[reaction//3][2],
                                         "{:10.2f}".format(truss.force[(reaction//3) * 3 + 0])+'\n' +
                                         "{:10.2f}".format(truss.force[(reaction//3) * 3 + 1])+'\n' +
                                         "{:10.2f}".format(truss.force[(reaction//3) * 3 + 2]),
                                         fontsize=12, horizontalalignment='right')
                            elif truss.dof == 2:
                                _ax.text(truss.nodal_coord[reaction//3][0],
                                         truss.nodal_coord[reaction//3][1],
                                         truss.nodal_coord[reaction//3][2],
                                         "{:10.2f}".format(truss.force[(reaction//3) * 3 + 0])+'\n' +
                                         "{:10.2f}".format(truss.force[(reaction//3) * 3 + 1]),
                                         fontsize=12, horizontalalignment='right')
                e_previous = reaction

        if supports:
            for support_index in truss.known_displacement_a:
                if support_index % 3 == 0:
                    f_dir = [-1.0 * scale_sup, 0, 0]
                    col = 'g'
                elif support_index % 3 == 1:
                    f_dir = [0, -1.0 * scale_sup, 0]
                    col = 'y'
                else:
                    f_dir = [0, 0, -1.0 * scale_sup * z_correction]
                    col = 'brown'

                #if support_index % 3 == 2 and truss.dof == 3:
                _ax.plot([truss.nodal_coord[support_index // 3][0], truss.nodal_coord[support_index // 3][0]+f_dir[0]],
                         [truss.nodal_coord[support_index // 3][1], truss.nodal_coord[support_index // 3][1]+f_dir[1]],
                         color=col, linewidth=4.0)
                #else:
                #    _ax.plot(
                #        [truss.nodal_coord[support_index // 3][0], truss.nodal_coord[support_index // 3][0] + f_dir[0]],
                #        [truss.nodal_coord[support_index // 3][1], truss.nodal_coord[support_index // 3][1] + f_dir[1]],
                #        color=col, linewidth=4.0)

        path = None

        # pyplot.show()
        if save_plot:
            path = './Results/%s.png' % plotname
            fig.savefig(path)
            print("'%s.png' is saved." % plotname)
            print('------------------------------------')

        return path

