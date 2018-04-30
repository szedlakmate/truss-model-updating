# -*- coding: utf-8 -*-
"""
Created on Wed Dec 27 10:02:39 2017
"""
from copy import deepcopy
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d


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

    def plotstructure(struct, show_orig, show_result, show_supports,
            show_forces, show_reactions, scale_displapcements, scale_forces, z_correction, show_values, save_plot, log=0):
        """
        General plotting method for structures

        scale_displapcements: Scale drwaing of displacements
        scale_forces:   Scale force sign
        z_correction:    Scale z-axis
        """
        plotname = struct.name
        plot_width = 10.0              # Plot width in inches
        xframe = 0                     # Frame width at X direction
        yframe = 0                     # Frame width at Y direction
        zframe = 0                     # Frame width at Z direction
        scale_sup = scale_forces*0.3        # Scale support sign  # All the others are input parameters

        # Stress coloring settings [R G B] - Examples considering pressure:
        # 0: [1, 0, 0]              Plain red
        # 1: [x, 0, 0]              Red to Black
        # 2: [1, 1-x, 0-x]          Red to White
        # 3: [1, (1-x)/2, (1-x)/2]  Red to MildRed - Distincts pressure and tension
        # 4: [x, 1-x, 0]            Red to Green
        _coloring = 3  # € [0, 1, 2, 3, 4]

        fig = plt.figure()
        _ax = fig.add_subplot(111, projection='3d')

        if struct.DOF == 2:
            _ax.view_init(elev=90., azim=-90.)
            _ax.w_zaxis.line.set_lw(0.)
            _ax.set_zticklabels([])

        xmin = min(list(struct.nodal_coord[x][0] for x in range(struct.node_num)))
        xmax = max(list(struct.nodal_coord[x][0] for x in range(struct.node_num)))
        ymin = min(list(struct.nodal_coord[x][1] for x in range(struct.node_num)))
        ymax = max(list(struct.nodal_coord[x][1] for x in range(struct.node_num)))
        zmin = min(list(struct.nodal_coord[x][2] for x in range(struct.node_num)))
        zmax = max(list(struct.nodal_coord[x][2] for x in range(struct.node_num)))

        deltax = xmax - xmin
        deltay = ymax - ymin

        xframe = max(deltax * 0.05, 2)
        yframe = max(deltay * 1.5, 2)

        if struct.DOF == 3:
            plot_height = plot_width * ((deltay + yframe*2)/(deltax + xframe*2)) * 0.3
        else:
            plot_height = plot_width * 0.5 * 0.5
        fig.set_size_inches(plot_width, plot_height)

        _ax.set_xlim3d(xmin - xframe, xmax + xframe)
        _ax.set_ylim3d(ymin - yframe, ymax + yframe)
        _ax.set_zlim3d(zmin - zframe, zmax + zframe)

        if show_orig == show_result:
            _coloring = 0

        # Giving plot names
        if show_orig == 1 and show_result == 0 and show_supports == 1 and show_reactions == 0:
            plotname += ' - Initial structure'
            if show_forces:
                plotname += ' with forces'
        elif show_orig == 1 and show_result == 1:
            plotname += ' - Deformation'
            if show_reactions == 0:
                plotname += ' with reactions'
        elif show_orig == 0 and show_result == 1:
            plotname += ' - Stresses'
            if show_reactions == 0:
                plotname += ' with reactions'
        else:
            plotname += ' - Unnamed'

        print(plotname + ": ")
        if show_result:
            dipslaydisplacement = deepcopy(struct.nodal_coord_def)
            if scale_displapcements != 1.0:
                if log:
                    print('Displacements are scaled with factor: ' + str(scale_displapcements))
                for i in range(struct.node_num):
                    for j in range(3):
                        dipslaydisplacement[i][j] = (struct.nodal_coord_def[i][j] -
                        struct.nodal_coord[i][j]) * scale_displapcements + struct.nodal_coord[i][j]

        for i in range(struct.element_num):
            # Plot undeformed structure
            if show_orig:
                _ax.plot([struct.nodal_coord[struct.node[i][1]][0], struct.nodal_coord[struct.node[i][0]][0]],
                    [struct.nodal_coord[struct.node[i][1]][1], struct.nodal_coord[struct.node[i][0]][1]],
                    zs=[struct.nodal_coord[struct.node[i][1]][2], struct.nodal_coord[struct.node[i][0]][2]], color='b')
            # Plot deformed structure
            if show_result:
                if struct._post_processed:
                    if struct.stress_color[i] > 0:
                        if _coloring == 1:
                            rgb_col = [0, 0, abs(struct.stress_color[i])]
                        elif _coloring == 2:
                            rgb_col = [1-abs(struct.stress_color[i]),
                                       1-abs(struct.stress_color[i]), 1]
                        elif _coloring == 3:
                            rgb_col = [(1-abs(struct.stress_color[i]))/2,
                                       (1-abs(struct.stress_color[i]))/2, 1]
                        elif _coloring == 4:
                            rgb_col = [0, 1-abs(struct.stress_color[i]),
                                       abs(struct.stress_color[i])]
                        else:
                            rgb_col = [1, 0, 0]
                    else:
                        if _coloring == 1:
                            rgb_col = [abs(struct.stress_color[i]), 0, 0]
                        elif _coloring == 2:
                            rgb_col = [1, 1-abs(struct.stress_color[i]),
                                       1-abs(struct.stress_color[i])]
                        elif _coloring == 3:
                            rgb_col = [1, (1-abs(struct.stress_color[i]))/2,
                                       (1-abs(struct.stress_color[i]))/2]
                        elif _coloring == 4:
                            rgb_col = [abs(struct.stress_color[i]),
                                       1-abs(struct.stress_color[i]), 0]
                        else:
                            rgb_col = [1, 0, 0]
                else:
                    print('Stresses are not calculated')
                    rgb_col = [1, 0, 0]
                _ax.plot([dipslaydisplacement[struct.node[i][1]][0], dipslaydisplacement[struct.node[i][0]][0]],
                        [dipslaydisplacement[struct.node[i][1]][1], dipslaydisplacement[struct.node[i][0]][1]],
                        zs=[dipslaydisplacement[struct.node[i][1]][2], dipslaydisplacement[struct.node[i][0]][2]],
                        color=rgb_col)

        if show_forces:
            for i in struct.known_f_not_zero:
                if struct.force[i] < 0:
                    value = -1.0
                else:
                    value = 1.0
                if i % 3 == 0:
                    f_dir = [value*scale_forces, 0., 0.]
                elif i % 3 == 1:
                    f_dir = [0., value*scale_forces, 0.]
                else:
                    f_dir = [0., 0., value*scale_forces*z_correction]
                f_arrow = Arrow3D([struct.nodal_coord[i//3][0], struct.nodal_coord[i//3][0] + f_dir[0]],
                                  [struct.nodal_coord[i//3][1], struct.nodal_coord[i//3][1] + f_dir[1]],
                                  [struct.nodal_coord[i//3][2], struct.nodal_coord[i//3][2] + f_dir[2]],
                                  mutation_scale=20, lw=1, arrowstyle="-|>", color="k")
                _ax.add_artist(f_arrow)

        if show_reactions:
            e_previous = -100
            for i in struct.known_dis_a:
                value = 0.0             # Maybe this is useless <XXX>
                if struct.force[i] < 0:
                    value = -1.0
                elif struct.force[i] > 0:
                    value = 1.0
                if i % 3 == 0:
                    f_dir = [value*scale_forces, 0., 0.]
                elif i % 3 == 1:
                    f_dir = [0., value*scale_forces, 0.]
                else:
                    f_dir = [0., 0., value*scale_forces*z_correction]
                if abs(struct.force[i]) > 0:
                    f_arrow = Arrow3D([struct.nodal_coord[i//3][0], struct.nodal_coord[i//3][0] + f_dir[0]],
                                      [struct.nodal_coord[i//3][1], struct.nodal_coord[i//3][1] + f_dir[1]],
                                      [struct.nodal_coord[i//3][2], struct.nodal_coord[i//3][2] + f_dir[2]],
                                      mutation_scale=20, lw=1, arrowstyle="-|>", color="darkolivegreen")
                    _ax.add_artist(f_arrow)
                    if show_values:
                        _ax.set_xticklabels([])
                        _ax.set_yticklabels([])
                        _ax.set_zticklabels([])
                        if not i//3 == e_previous//3:
                            if struct.DOF == 3:
                                _ax.text(struct.nodal_coord[i//3][0],
                                    struct.nodal_coord[i//3][1],
                                    struct.nodal_coord[i//3][2],
                                    "{:10.2f}".format(struct.force[(i//3)*3+0])+'\n' +
                                    "{:10.2f}".format(struct.force[(i//3)*3+1])+'\n' +
                                    "{:10.2f}".format(struct.force[(i//3)*3+2]),
                                    fontsize=12, horizontalalignment='right')
                            elif struct.DOF == 2:
                                _ax.text(struct.nodal_coord[i//3][0],
                                    struct.nodal_coord[i//3][1],
                                    struct.nodal_coord[i//3][2],
                                    "{:10.2f}".format(struct.force[(i//3)*3+0])+'\n' +
                                    "{:10.2f}".format(struct.force[(i//3)*3+1]),
                                    fontsize=12, horizontalalignment='right')
                e_previous = i

        if show_supports:
            for i in struct.known_dis_a:
                if i % 3 == 0:
                    f_dir = [-1.0 * scale_sup, 0., 0.]
                    col = 'g'
                elif i % 3 == 1:
                    f_dir = [0., -1.0 * scale_sup, 0.]
                    col = 'y'
                else:
                    f_dir = [0., 0., -1.0 * scale_sup * z_correction]
                    col = 'brown'
                if i % 3 != 2 or struct.DOF == 3:
                    _ax.plot([struct.nodal_coord[i//3][0], struct.nodal_coord[i//3][0]+f_dir[0]],
                        [struct.nodal_coord[i//3][1], struct.nodal_coord[i//3][1]+f_dir[1]],
                        zs=[struct.nodal_coord[i//3][2], struct.nodal_coord[i//3][2]+f_dir[2]],
                        color=col, linewidth=4.0)
        plt.show()
        if save_plot:
            fig.savefig("./Structures/" + plotname + '.png')
            print("'" + plotname + ".png' is saved.")
            print('------------------------------------')
        return
