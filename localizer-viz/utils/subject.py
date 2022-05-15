import os
import nibabel as nib
import nibabel.freesurfer.io as fio
import numpy as np

from .io_freesurfer import read_electrodes_location
from .depth2surf import ElocsProjector


class Subject:
    def __init__(self, subj, f_elocs=None, fsfast_root=None, hemi='lh', load=True, **kwargs):
        assert hemi in ['lh', 'rh']
        self.hemi = hemi
        self.fsfast_root = fsfast_root
        self.subj = subj
        subject_path = os.path.join(os.getenv('FREESURFER_HOME'), 'subjects/%s' % subj)
        self.subject_path = subject_path
        self.pial_path = os.path.join(subject_path, 'surf/%s.pial' % hemi)
        self.inflated_path = os.path.join(subject_path, 'surf/%s.inflated' % hemi)
        self.curv_path = os.path.join(subject_path, 'surf/%s.curv' % hemi)

        # for MT only
        self.fsfast_root = fsfast_root
        self.fmri_mb_sig_path = os.path.join(fsfast_root,
                                             'mb/bold/mt.%s/mt/sig.nii.gz' % hemi) if fsfast_root is not None else None
        self.fmri_dot_sig_path = os.path.join(fsfast_root,
                                              'dot/bold/mt.%s/mt/sig.nii.gz' % hemi) if fsfast_root is not None else None
        if f_elocs is not None:
            elocs_name = f_elocs.split('/')
            self.elocs_path = '/'.join(elocs_name[:-1])
            self.elocs_name = elocs_name[-1]
        else:
            self.elocs_path = None
            self.elocs_name = None

        self._mb = None
        self._dot = None
        self._pial = None
        self._inflated = None
        self._curv = None
        self._elocs = None
        self._eloc_info = None
        self._projector = None

        if load:
            # read sig overlay
            self._mb = nib.load(self.fmri_mb_sig_path)

            # read pial surf
            vert_pial, faces_pial = fio.read_geometry(self.pial_path)
            self._pial = {'vert': vert_pial, 'faces': faces_pial}
            vert_inflated, faces_inflated = fio.read_geometry(self.inflated_path)
            self._inflated = {'vert': vert_inflated, 'faces': faces_inflated}

            # read curv
            self._curv = fio.read_morph_data(self.curv_path)

            # read point set
            self._elocs, self._eloc_info = read_electrodes_location(f_elocs)

            # build projector
            self._projector = ElocsProjector(surf_pial=self._pial['vert'],
                                             surf_inflated=self.inflated['vert'],
                                             proj_method='nearest',
                                             neighbor_method='radius', **kwargs)

        self.ind_elecs = None
        self.ind_surf_vertex = None
        self.coor_elecs_pial = None
        self.coor_elecs_inflated = None

    @property
    def pial(self):
        if self._pial is None:
            vert_pial, faces_pial = fio.read_geometry(self.pial_path)
            self._pial = {'vert': vert_pial, 'faces': faces_pial}
        return self._pial

    @property
    def inflated(self):
        if self._inflated is None:
            vert_inflated, faces_inflated = fio.read_geometry(self.inflated_path)
            self._inflated = {'vert': vert_inflated, 'faces': faces_inflated}
        return self._inflated

    @property
    def curv(self):
        if self._curv is None:
            self._curv = fio.read_morph_data(self.curv_path)
        return self._curv

    @property
    def mb_overlay(self):
        if self._mb is None:
            self._mb = nib.load(self.fmri_mb_sig_path)
        return self._mb.get_fdata()

    @property
    def dot_overlay(self):
        if self._dot is None:
            self._dot = nib.load(self.fmri_dot_sig_path)
        return self._dot.get_fdata()

    def load_label(self, label):
        label_name = label
        label_fname = ".".join([self.hemi, label_name, 'label'])

        filepath = os.path.join(self.subject_path, 'label', label_fname)

        if not os.path.exists(filepath):
            raise ValueError('Label file %s does not exist'
                             % filepath)
        ids = nib.freesurfer.read_label(filepath)
        label = np.zeros(self.inflated['vert'].shape[0])
        label[ids] = 1
        return ids, label

    def load_annot(self, annot, target_structure=None):
        """

        :param annot: str, name of the annotation file.
        :param target_structure: str, label name of what we need.
        :return:
        """
        annot_name = annot
        annot_fname = ".".join([self.hemi, annot_name, 'annot'])

        filepath = os.path.join(self.subject_path, 'label', annot_fname)

        if not os.path.exists(filepath):
            raise ValueError('Annot file %s does not exist'
                             % filepath)
        # ids = nib.freesurfer.read_label(filepath)
        labels, ctab, names = nib.freesurfer.read_annot(filepath)
        names = [n_.decode("utf-8") for n_ in names]
        if target_structure is not None:
            if self.hemi == 'lh':
                target_structure = target_structure + '_L'
            else:
                target_structure = target_structure + '_R'
            idx = names.index(target_structure)
            target_mask = np.zeros_like(labels)
            target_ids, = np.nonzero(labels == idx)
            target_mask[target_ids] = 1
            return target_ids, target_mask
        else:
            return labels, ctab, names

    def get_elecs_indices(self):
        return self.ind_elecs.copy(), self.ind_surf_vertex.copy()

    def eleoc_projector(self):
        ind_elecs, ind_surf_vertex = self._projector.elocs_projection(self._elocs)
        self.ind_elecs = ind_elecs
        self.ind_surf_vertex = ind_surf_vertex
        self.coor_elecs_pial = self._pial['vert'][self.ind_surf_vertex]
        self.coor_elecs_inflated = self.inflated['vert'][self.ind_surf_vertex]

    def save_eloc(self, save_dir=None):
        if save_dir is None:
            save_dir = self.elocs_path
        if self.ind_elecs is not None and self.ind_surf_vertex is not None:
            np.savez(os.path.join(save_dir, 'elocs.npz'),
                     ind_elecs=self.ind_elecs,
                     ind_surf_vertex=self.ind_surf_vertex,
                     coor_elecs_pial=self.coor_elecs_pial,
                     coor_elecs_inflated=self.coor_elecs_inflated)

    def load_eloc(self, save_dir=None):
        if save_dir is None:
            save_dir = self.elocs_path
        try:
            data = np.load(os.path.join(save_dir, 'elocs.npz'))
        except FileNotFoundError:
            print('Projected electrode vertex indices not found.')
            return -1
        self.ind_elecs = data['ind_elecs']
        self.ind_surf_vertex = data['ind_surf_vertex']
        self.coor_elecs_pial = data['coor_elecs_pial']
        self.coor_elecs_inflated = data['coor_elecs_inflated']
