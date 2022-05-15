import argparse
import sys
import os

from utils.subject import Subject


def parse_args():
    parser = argparse.ArgumentParser(
        description='Electrodes projection (depth to surface)'
    )
    parser.add_argument(
        '--subject',
        '-s',
        dest='subject',
        help='Freesurfer subject name',
        default=None,
        type=str
    )
    parser.add_argument(
        '--fsfast',
        dest='fsfast_root',
        help='Freesurfer subject fsfast root',
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
        help='Path to electrode data file.',
        type=str
    )
    if len(sys.argv) < 4:
        parser.print_help()
        sys.exit(1)

    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    patient = Subject(args.subject, args.elocs, fsfast_root=args.fsfast_root, hemi=args.hemi)

    patient.eleoc_projector()
    # save eloc
    patient.save_eloc()

