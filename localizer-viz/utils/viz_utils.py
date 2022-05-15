import mayavi
from mayavi import mlab
import os
import numpy as np
import matplotlib as mpl
from matplotlib import cm
from sklearn.metrics.pairwise import rbf_kernel


class MTLabel:
    def __init__(self, hemi, vertices, name):
        self.hemi = hemi
        self.vertices = vertices
        self.name = name


def use_phong_shading(surf, **kwargs):
    """
    Change the shading method to Phong to get a shiny surface
    Parameters.

    ----------
    :param kwargs:
    :return: None
    """
    ambiant = kwargs.pop('ambient', 0.4225)
    specular = kwargs.pop('specular', 0.333)
    specular_power = kwargs.pop('specular_power', 66)
    diffuse = kwargs.pop('diffuse', 0.6995)
    interpolation = kwargs.pop('interpolation', 'phong')
    surf.actor.property.ambient = ambiant
    surf.actor.property.specular = specular
    surf.actor.property.specular_power = specular_power
    surf.actor.property.diffuse = diffuse
    surf.actor.property.interpolation = interpolation


def wrap_freeview_command(patient, overlay_fname=None, overlay_max_value=None, kind='elecs'):
    assert kind in ['elecs', 'bold', 'MT-Localizer']
    if kind == 'bold' or kind == 'MT-Localizer':
        cmd = 'freeview -f %s:curvature_method=binary:overlay=%s:overlay_color=heat:overlay_threshold=2,8' % (
            os.path.abspath(patient.inflated_path),
            os.path.abspath(patient.fmri_mb_sig_path))
    elif kind == 'elecs':
        cmd = 'freeview -f %s:curvature_method=binary:overlay=%s:overlay_color=colorwheel:overlay_threshold=0.1,%d' % (
            os.path.abspath(patient.inflated_path),
            os.path.abspath(overlay_fname),
            overlay_max_value)
    else:
        raise KeyError
    return cmd

# # MT-Localizer overlay
# # thresholding
# bold_upth_ind, bold_upth_value, bold_th = bpt.fmri_thresholding(patient.bold_overlay, th=2.)
# bold_upth_coor = patient.inflated['vert'][bold_upth_ind]
#
# # maya plot brain
# mesh, mlab = bpt.ctmr_gauss_plot(patient.inflated['faces'],
#                                  patient.inflated['vert'],
#                                  elecs=bold_upth_coor,
#                                  weights=bold_upth_value,
#                                  gsp=5.,
#                                  opacity=0.8,
#                                  vmin=-8,
#                                  vmax=8)
# bpt.el_add(elocs_inflated, color=plt.get_cmap('tab10')(0)[:-1], fontsize=2, label_offset=-1.5)
#
# mlab.show()

# ctmr_brain_plot.py


def fmri_thresholding(fmri_overlay, ratio=1e-3, th=None):
    fmri_overlay = fmri_overlay.squeeze()
    ind = np.argsort(fmri_overlay)
    # largest 1e-3
    if th is None:
        th = fmri_overlay[int(len(ind) * ratio)]
    ind_th = np.nonzero(fmri_overlay > th)[0]
    fmri_overlay_th = fmri_overlay[ind_th]
    return ind_th, fmri_overlay_th, th


def ctmr_gauss_plot(tri, vert, color=(0.8, 0.8, 0.8), elecs=None, weights=None,
                    opacity=1.0, representation='surface', line_width=1.0, gsp=10,
                    cmap=mpl.cm.get_cmap('RdBu_r'), show_colorbar=True, new_fig=True, vmin=None, vmax=None,
                    ambient=0.4225, specular=0.333, specular_power=66, diffuse=0.6995, interpolation='phong'):
    ''' This function plots the 3D brain surface mesh

    Parameters
    ----------
        color : tuple
            (n,n,n) tuple of floats between 0.0 and 1.0, background color of brain
        elecs : array-like
            [nchans x 3] matrix of electrode coordinate values in 3D
        weights : array-like
            [nchans x 1] - if [elecs] is also given, this will color the brain vertices
            according to these weights
        msize : float
            size of the electrode.  default = 2
        opacity : float (0.0 - 1.0)
            opacity of the brain surface (value from 0.0 - 1.0)
        cmap : str or mpl.colors.LinearSegmentedColormap
            colormap to use when plotting gaussian weights with [elecs]
            and [weights]
        representation : {'surface', 'wireframe'}
            surface representation
        line_width : float
            width of lines for triangular mesh
        gsp : float
            gaussian smoothing parameter, larger makes electrode activity
            more spread out across the surface if specified

    Returns
    -------
    mesh : mayavi mesh (actor)
    mlab : mayavi mlab scene
    '''
    # if color is another iterable, make it a tuple.
    color = tuple(color)

    brain_color = []
    # c = np.zeros(vert.shape[0],)

    if elecs is not None:
        brain_color = np.zeros(vert.shape[0], )
        for i in range(elecs.shape[0]):
            gauss_wt = weights[i] * rbf_kernel(vert, elecs[[i]], gamma=1 / gsp).squeeze()
            brain_color += gauss_wt

        # scale the colors so that it matches the weights that were passed in
        brain_color *= (np.abs(weights).max() / np.abs(brain_color).max())
        if vmin is None and vmax is None:
            vmin, vmax = -np.abs(brain_color).max(), np.abs(brain_color).max()

    # plot cortex and begin display
    if new_fig:
        mlab.figure(fgcolor=(0, 0, 0), bgcolor=(1, 1, 1), size=(1200, 900))

    if elecs is not None:
        kwargs = {}
        if type(cmap) == str:
            kwargs.update(colormap=cmap)

        mesh = mlab.triangular_mesh(vert[:, 0], vert[:, 1], vert[:, 2], tri,
                                    representation=representation, opacity=opacity,
                                    line_width=line_width, scalars=brain_color,
                                    vmin=vmin, vmax=vmax, **kwargs)

        if type(cmap) == mpl.colors.LinearSegmentedColormap:
            mesh.module_manager.scalar_lut_manager.lut.table = (cmap(np.linspace(0, 1, 255)) * 255).astype('int')
    else:
        mesh = mlab.triangular_mesh(vert[:, 0], vert[:, 1], vert[:, 2], tri,
                                    color=color, representation=representation,
                                    opacity=opacity, line_width=line_width)

    # cell_data = mesh.mlab_source.dataset.cell_data
    # cell_data.scalars = brain_color
    # cell_data.scalars.name = 'Cell data'
    # cell_data.update()

    # mesh2 = mlab.pipeline.set_active_attribute(mesh, cell_scalars = 'Cell data')
    # mlab.pipeline.surface(mesh)
    if weights is not None and show_colorbar:
        mlab.colorbar()

    # change OpenGL mesh properties for phong point light shading
    mesh.actor.property.ambient = ambient
    mesh.actor.property.specular = specular
    mesh.actor.property.specular_power = specular_power
    mesh.actor.property.diffuse = diffuse
    mesh.actor.property.interpolation = interpolation
    # mesh.scene.light_manager.light_mode = 'vtk'
    if opacity < 1.0:
        mesh.scene.renderer.set(use_depth_peeling=True)  # , maximum_number_of_peels=100, occlusion_ratio=0.0

    # Make the mesh look smoother
    for child in mlab.get_engine().scenes[0].children:
        poly_data_normals = child.children[0]
        poly_data_normals.filter.feature_angle = 80.0  # Feature angle says which angles are considered hard corners

    return mesh, mlab


def el_add(numbers, elecs, color=(1., 0., 0.), msize=2, label_offset=-1.0, ambient=0.3261, specular=1,
           specular_power=16, diffuse=0.6995, interpolation='phong', fontsize=5, **kwargs):
    '''This function adds the electrode matrix [elecs] (nchans x 3) to
    the scene.

    Parameters
    ----------
        elecs : numpyarray {number: (x, y, z)}
            [nchans x 3] matrix of electrode coordinate values in 3D
        color : tuple (triplet) or numpy array
            Electrode color is either a triplet (r, g, b),
            or a numpy array with the same shape as [elecs] to plot one color per electrode
        msize : float
            size of the electrode.  default = 2
        label_offset : float
            how much to move the number labels out by (so not blocked by electrodes)
        **kwargs :
            any other keyword arguments that can be passed to points3d
    '''

    # Get the current keyword arguments
    cur_kwargs = dict(color=color, scale_factor=msize, resolution=25)

    # Allow the user to override the default keyword arguments using kwargs
    cur_kwargs.update(kwargs)

    # plot the electrodes as spheres
    # If we have one color for each electrode, color them separately
    if type(color) is np.ndarray:
        if color.shape[0] == elecs.shape[0]:
            # for e in np.arange(elecs.shape[0]):
            #     points = mlab.points3d(elecs[e,0], elecs[e,1], elecs[e,2], scale_factor = msize,
            #                        color = tuple( color[e,:] ) , resolution=25)
            unique_colors = np.array(list(set([tuple(row) for row in color])))
            for individual_color in unique_colors:
                indices = np.where((color == individual_color).all(axis=1))[0]
                cur_kwargs.update(color=tuple(individual_color))
                points = mlab.points3d(elecs[indices, 0], elecs[indices, 1], elecs[indices, 2],
                                       **cur_kwargs)
        else:
            print('Warning: color array does not match size of electrode matrix')

    # Otherwise, use the same color for all electrodes
    else:
        points = mlab.points3d(elecs[:, 0], elecs[:, 1], elecs[:, 2], **cur_kwargs)

    # Set display properties
    points.actor.property.ambient = ambient
    points.actor.property.specular = specular
    points.actor.property.specular_power = specular_power
    points.actor.property.diffuse = diffuse
    points.actor.property.interpolation = interpolation
    # points.scene.light_manager.light_mode = 'vtk'

    if numbers is not None:
        for ni, n in enumerate(numbers):
            mayavi.mlab.text3d(elecs[ni, 0] + label_offset, elecs[ni, 1], elecs[ni, 2], str(n),
                               orient_to_camera=True, scale=fontsize)  # line_width=5.0, scale=1.5)

    return points, mlab
