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

    def plotstructure(struct, showorig, showresult, showsupports,
            showforces, showreactions, scaledisp, scale_f, z_corr, showvalues, saveplot, log=0):
        """
        General plotting method for structures

        scaledisp: Scale drwaing of displacements
        scale_f:   Scale force sign
        z_corr:    Scale z-axis
        """
        plotname = struct.name
        plot_width = 18.0              # Plot width in inches
        xframe = 0                     # Frame width at X direction
        yframe = 0                     # Frame width at Y direction
        zframe = 0                     # Frame width at Z direction
        scale_sup = scale_f*0.3        # Scale support sign  # All the others are input parameters

        # Stress coloring settings [R G B] - Examples considering pressure:
        # 0: [1, 0, 0]              Plain red
        # 1: [x, 0, 0]              Red to Black
        # 2: [1, 1-x, 0-x]          Red to White
        # 3: [1, (1-x)/2, (1-x)/2]  Red to MildRed - Distincts pressure and tension
        # 4: [x, 1-x, 0]            Red to Green
        _coloring = 3  # € [0, 1, 2, 3, 4]

        fig = plt.figure()
        _ax = fig.add_subplot(111, projection='3d')

        if struct.dof == 2:
            _ax.view_init(elev=90., azim=-90.)
            _ax.w_zaxis.line.set_lw(0.)
            _ax.set_zticklabels([])

        xmin = min(list(struct.nodalcoord[x][0] for x in range(struct.nodenum)))
        xmax = max(list(struct.nodalcoord[x][0] for x in range(struct.nodenum)))
        ymin = min(list(struct.nodalcoord[x][1] for x in range(struct.nodenum)))
        ymax = max(list(struct.nodalcoord[x][1] for x in range(struct.nodenum)))
        zmin = min(list(struct.nodalcoord[x][2] for x in range(struct.nodenum)))
        zmax = max(list(struct.nodalcoord[x][2] for x in range(struct.nodenum)))

        deltax = xmax - xmin
        deltay = ymax - ymin

        xframe = max(deltax * 0.05, 2)
        yframe = max(deltay * 1.5, 2)

        if struct.dof == 3:
            plot_height = plot_width * ((deltay + yframe*2)/(deltax + xframe*2)) * 0.3
        else:
            plot_height = plot_width * 0.5
        fig.set_size_inches(plot_width, plot_height)

        _ax.set_xlim3d(xmin - xframe, xmax + xframe)
        _ax.set_ylim3d(ymin - yframe, ymax + yframe)
        _ax.set_zlim3d(zmin - zframe, zmax + zframe)

        if showorig == showresult:
            _coloring = 0

        # Giving plot names
        if showorig == 1 and showresult == 0 and showsupports == 1 and showreactions == 0:
            plotname += ' - Initial structure'
            if showforces:
                plotname += ' with forces'
        elif showorig == 1 and showresult == 1:
            plotname += ' - Deformation'
            if showreactions == 0:
                plotname += ' with reactions'
        elif showorig == 0 and showresult == 1:
            plotname += ' - Stresses'
            if showreactions == 0:
                plotname += ' with reactions'
        else:
            plotname += ' - Unnamed'

        print(plotname + ": ")
        if showresult:
            dipslaydisplacement = deepcopy(struct.nodalcoord_def)
            if scaledisp != 1.0:
                if log:
                    print('Displacements are scaled with factor: ' + str(scaledisp))
                for i in range(struct.nodenum):
                    for j in range(3):
                        dipslaydisplacement[i][j] = (struct.nodalcoord_def[i][j] -
                        struct.nodalcoord[i][j]) * scaledisp + struct.nodalcoord[i][j]

        for i in range(struct.elenum):
            # Plot undeformed structure
            if showorig:
                _ax.plot([struct.nodalcoord[struct.node[i][1]][0], struct.nodalcoord[struct.node[i][0]][0]],
                    [struct.nodalcoord[struct.node[i][1]][1], struct.nodalcoord[struct.node[i][0]][1]],
                    zs=[struct.nodalcoord[struct.node[i][1]][2], struct.nodalcoord[struct.node[i][0]][2]], color='b')
            # Plot deformed structure
            if showresult:
                if struct.postprocessed():
                    if struct.stresscolor[i] > 0:
                        if _coloring == 1:
                            rgb_col = [0, 0, abs(struct.stresscolor[i])]
                        elif _coloring == 2:
                            rgb_col = [1-abs(struct.stresscolor[i]),
                                       1-abs(struct.stresscolor[i]), 1]
                        elif _coloring == 3:
                            rgb_col = [(1-abs(struct.stresscolor[i]))/2,
                                       (1-abs(struct.stresscolor[i]))/2, 1]
                        elif _coloring == 4:
                            rgb_col = [0, 1-abs(struct.stresscolor[i]),
                                       abs(struct.stresscolor[i])]
                        else:
                            rgb_col = [1, 0, 0]
                    else:
                        if _coloring == 1:
                            rgb_col = [abs(struct.stresscolor[i]), 0, 0]
                        elif _coloring == 2:
                            rgb_col = [1, 1-abs(struct.stresscolor[i]),
                                       1-abs(struct.stresscolor[i])]
                        elif _coloring == 3:
                            rgb_col = [1, (1-abs(struct.stresscolor[i]))/2,
                                       (1-abs(struct.stresscolor[i]))/2]
                        elif _coloring == 4:
                            rgb_col = [abs(struct.stresscolor[i]),
                                       1-abs(struct.stresscolor[i]), 0]
                        else:
                            rgb_col = [1, 0, 0]
                else:
                    print('Stresses are not calculated')
                    rgb_col = [1, 0, 0]
                _ax.plot([dipslaydisplacement[struct.node[i][1]][0], dipslaydisplacement[struct.node[i][0]][0]],
                        [dipslaydisplacement[struct.node[i][1]][1], dipslaydisplacement[struct.node[i][0]][1]],
                        zs=[dipslaydisplacement[struct.node[i][1]][2], dipslaydisplacement[struct.node[i][0]][2]],
                        color=rgb_col)

        if showforces:
            for i in struct.known_f_notzero:
                if struct.force[i] < 0:
                    value = -1.0
                else:
                    value = 1.0
                if i % 3 == 0:
                    f_dir = [value*scale_f, 0., 0.]
                elif i % 3 == 1:
                    f_dir = [0., value*scale_f, 0.]
                else:
                    f_dir = [0., 0., value*scale_f*z_corr]
                f_arrow = Arrow3D([struct.nodalcoord[i//3][0], struct.nodalcoord[i//3][0] + f_dir[0]],
                                  [struct.nodalcoord[i//3][1], struct.nodalcoord[i//3][1] + f_dir[1]],
                                  [struct.nodalcoord[i//3][2], struct.nodalcoord[i//3][2] + f_dir[2]],
                                  mutation_scale=20, lw=1, arrowstyle="-|>", color="k")
                _ax.add_artist(f_arrow)

        if showreactions:
            e_previous = -100
            for i in struct.known_dis_a:
                value = 0.0             # Maybe this is useless <XXX>
                if struct.force[i] < 0:
                    value = -1.0
                elif struct.force[i] > 0:
                    value = 1.0
                if i % 3 == 0:
                    f_dir = [value*scale_f, 0., 0.]
                elif i % 3 == 1:
                    f_dir = [0., value*scale_f, 0.]
                else:
                    f_dir = [0., 0., value*scale_f*z_corr]
                if abs(struct.force[i]) > 0:
                    f_arrow = Arrow3D([struct.nodalcoord[i//3][0], struct.nodalcoord[i//3][0] + f_dir[0]],
                                      [struct.nodalcoord[i//3][1], struct.nodalcoord[i//3][1] + f_dir[1]],
                                      [struct.nodalcoord[i//3][2], struct.nodalcoord[i//3][2] + f_dir[2]],
                                      mutation_scale=20, lw=1, arrowstyle="-|>", color="darkolivegreen")
                    _ax.add_artist(f_arrow)
                    if showvalues:
                        _ax.set_xticklabels([])
                        _ax.set_yticklabels([])
                        _ax.set_zticklabels([])
                        if not i//3 == e_previous//3:
                            if struct.dof == 3:
                                _ax.text(struct.nodalcoord[i//3][0],
                                    struct.nodalcoord[i//3][1],
                                    struct.nodalcoord[i//3][2],
                                    "{:10.2f}".format(struct.force[(i//3)*3+0])+'\n' +
                                    "{:10.2f}".format(struct.force[(i//3)*3+1])+'\n' +
                                    "{:10.2f}".format(struct.force[(i//3)*3+2]),
                                    fontsize=12, horizontalalignment='right')
                            elif struct.dof == 2:
                                _ax.text(struct.nodalcoord[i//3][0],
                                    struct.nodalcoord[i//3][1],
                                    struct.nodalcoord[i//3][2],
                                    "{:10.2f}".format(struct.force[(i//3)*3+0])+'\n' +
                                    "{:10.2f}".format(struct.force[(i//3)*3+1]),
                                    fontsize=12, horizontalalignment='right')
                e_previous = i

        if showsupports:
            for i in struct.known_dis_a:
                if i % 3 == 0:
                    f_dir = [-1.0 * scale_sup, 0., 0.]
                    col = 'g'
                elif i % 3 == 1:
                    f_dir = [0., -1.0 * scale_sup, 0.]
                    col = 'y'
                else:
                    f_dir = [0., 0., -1.0 * scale_sup * z_corr]
                    col = 'brown'
                if i % 3 != 2 or struct.dof == 3:
                    _ax.plot([struct.nodalcoord[i//3][0], struct.nodalcoord[i//3][0]+f_dir[0]],
                        [struct.nodalcoord[i//3][1], struct.nodalcoord[i//3][1]+f_dir[1]],
                        zs=[struct.nodalcoord[i//3][2], struct.nodalcoord[i//3][2]+f_dir[2]],
                        color=col, linewidth=4.0)
        plt.show()
        if saveplot:
            fig.savefig("./Structures/" + plotname + '.png')
            print("'" + plotname + ".png' is saved.")
            print('------------------------------------')
        return
