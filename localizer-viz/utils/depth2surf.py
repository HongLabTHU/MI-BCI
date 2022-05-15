import numpy as np
from sklearn.neighbors import NearestNeighbors


class ElocsProjector:
    # c.r.s unit, i.e.voxel size of aseg volume
    # in a range of 3.5 mm
    default_radius = 3.5

    def __init__(self, surf_pial, surf_inflated, **kwargs):
        """

        :param surf_pial:
        :param kwargs:
            radius: epsilon-neighbor threshold
            k: k-neighbor threshold
            surf_inflated: inflated surface vertices
        """
        self.surf_pial = surf_pial
        self.surf_inflated = surf_inflated
        self.radius = kwargs.pop('radius', ElocsProjector.default_radius)
        k = kwargs.pop('k', 1)
        self.neigh_pial = NearestNeighbors(n_neighbors=k, radius=self.radius)
        # fit
        self.neigh_pial.fit(surf_pial)

    def elocs_projection(self, elocs):
        """

        :param elocs:
        :return:
        """
        nbrs_dist, nbrs_ind = self.neigh_pial.kneighbors(elocs, 1, return_distance=True)
        nbrs_dist = nbrs_dist.squeeze()
        nbrs_ind = nbrs_ind.squeeze()
        elecs_in_radius_indices_surf = nbrs_ind[nbrs_dist < self.radius]
        # 1-based
        elecs_in_radius_indices_elec = np.nonzero(nbrs_dist < self.radius)[0] + 1
        return elecs_in_radius_indices_elec, elecs_in_radius_indices_surf


def coordinate_check(elocs, elocs_info, **kwargs):
    if elocs_info['useRealRAS']:
        # elocs is in CT space, need transformation
        # TODO: load CT vox2ras_tkr and vox2ras matrix
        # TODO: load CT registration mat
        # TODO: transform
        reg_mat_path = kwargs.pop('reg_ct')
        ct_vol_path = kwargs.pop('vol_ct')
        raise NotImplementedError

    return elocs
