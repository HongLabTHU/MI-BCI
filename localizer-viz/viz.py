import argparse
import sys

import numpy as np
from mayavi import mlab
from surfer import Brain

from utils.subject import Subject
from utils.viz_utils import MTLabel


def parse_args():
    parser = argparse.ArgumentParser(
        description='Electrodes projection (depth to surface)'
    )
    parser.add_argument(
        '--subject',
        '-s',
        dest='subject',
        help='Freesurfer subject name',
        required=True,
        default=None,
        type=str
    )
    parser.add_argument(
        '--fsfast',
        dest='fsfast_root',
        help='Freesurfer subject fsfast root',
        required=True,
        default=None,
        type=str
    )
    parser.add_argument(
        '--hemi',
        dest='hemi',
        help='Which hemisphere the electrodes are located at? lh or rh',
        default=None,
        type=str
    )
    parser.add_argument(
        '--electrode-locations',
        '-eloc',
        dest='elocs',
        help='Path to electrode data file',
        required=True,
        type=str
    )
    if len(sys.argv) < 5:
        parser.print_help()
        sys.exit(1)

    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    patient = Subject(args.subject, args.elocs, fsfast_root=args.fsfast_root, hemi=args.hemi, load=False)
    patient.load_eloc()
    elecs_ind, elecs_vert_ind = patient.get_elecs_indices()

    # figure size
    size = (1200, 900)
    fig = mlab.figure(size=size)
    brain = Brain(args.subject, args.hemi, 'inflated',
                  title=args.subject,
                  figure=fig,
                  foreground='white',
                  background='black',
                  alpha=1,
                  show_toolbar=True,
                  # cortex=(0.6, 0.6, 0.6)
                  )

    # mt label
    label_mt_name = 'V5/MT+'
    mt_ids, mt_mask = patient.load_annot('BN_Atlas', target_structure=label_mt_name)
    brain.add_label(MTLabel(patient.hemi, mt_ids, label_mt_name),
                    borders=True,
                    alpha=0.7,
                    color=np.array((247, 221, 96)) / 255.)

    # add overlay
    brain.add_overlay(patient.mb_overlay.squeeze(), min=1.5, max=8, sign='pos')

    # add electrodes
    # color_non = '#efefef'
    # color_mt = '#96d7fa'
    # color_best = '#46a2f9'
    brain.add_foci(elecs_vert_ind, coords_as_verts=True, color='#96d7fa', scale_factor=0.3,
                   resolution=100)

    mlab.show()
