import struct
import os
from copy import deepcopy

import numpy as np

TRIANGLE_FILE_MAGIC_NUMBER = 16777214
QUAD_FILE_MAGIC_NUMBER = 16777215
NEW_QUAD_FILE_MAGIC_NUMBER = 16777213

MRI_UCHAR = 0
MRI_INT = 1
MRI_LONG = 2
MRI_FLOAT = 3
MRI_SHORT = 4
MRI_BITMAP = 5
MRI_TENSOR = 6


def _fread3(fid):
    _b = fid.read(3)
    num = struct.unpack('>I', b'\0' + _b)[0]
    return num


def read_surf(fname):
    # read binary, big-endian
    fid = open(fname, 'rb')
    if not fid:
        raise Exception('could not open surface file %s' % fname)

    magic = _fread3(fid)

    if magic == QUAD_FILE_MAGIC_NUMBER or magic == NEW_QUAD_FILE_MAGIC_NUMBER:
        vnum = _fread3(fid)
        fnum = _fread3(fid)
        vertex_coords = np.fromfile(fid, dtype='>h', count=3 * vnum).astype(np.float64) / 100.

        faces = np.empty((fnum, 4), dtype=np.int32)
        for i in range(fnum):
            for n in range(4):
                faces[i, n] = _fread3(fid)
    elif magic == TRIANGLE_FILE_MAGIC_NUMBER:
        fid.readline()
        fid.readline()
        vnum = struct.unpack('>i', fid.read(4))[0]
        fnum = struct.unpack('>i', fid.read(4))[0]
        vertex_coords = np.fromfile(fid, dtype='>f', count=3 * vnum).astype(np.float64)
        faces = np.fromfile(fid, dtype='>i', count=3 * fnum).astype(np.int32)
        faces = np.reshape(faces, (fnum, 3))
    else:
        fid.close()
        raise Exception('ERROR: magic number %d unknown' % magic)

    vertex_coords = np.reshape(vertex_coords, (vnum, 3))
    fid.close()

    return vertex_coords, faces, magic


def read_annotation(filename, verbosity=False):
    def read_int(fp):
        b_ = fp.read(4)
        if len(b_) != 4:
            raise IOError
        num = struct.unpack('>i', b_)[0]
        return num

    # big-endian
    fp = open(filename, 'rb')

    if not fp:
        raise Exception('Annotation file cannot be opened')

    A = read_int(fp)

    tmp = np.fromfile(fp, dtype='>i', count=2 * A)
    vertices = tmp[::2]
    label = tmp[1::2]

    colortable = {}
    try:
        bool_ = read_int(fp)
    except IOError:
        if verbosity:
            print('No Colortable found.')
        fp.close()
        return vertices, label, colortable

    if bool_:
        numEntries = read_int(fp)

        if numEntries > 0:
            if verbosity:
                print('Reading from Original Version')
            colortable['numEntries'] = numEntries
            len_ = read_int(fp)
            colortable['orig_tab'] = fp.read(len_)[:-1].decode('utf-8')
            colortable['struct_names'] = dict()
            colortable['table'] = np.zeros((numEntries, 5), dtype=np.int32)
            for i in range(numEntries):
                len_ = read_int(fp)
                colortable['struct_names'][i] = fp.read(len_)[:-1].decode('utf-8')
                colortable['table'][i, :4] = np.array(struct.unpack('>4i', fp.read(4 * 4)))
                colortable['table'][i, 4] = colortable['table'][i, 0] + (colortable['table'][i, 1] << 8) + (
                        colortable['table'][i, 2] << 16) + (colortable['table'][i, 3] << 24)
            if verbosity:
                print('colortable with %d entries read (originally' % colortable['numEntries'],
                      colortable['orig_tab'], ')')
        else:
            version = -numEntries
            if verbosity:
                if version != 2:
                    print('Error! Does not handle version %d' % version)
                else:
                    print('Reading from version %d' % version)
            numEntries = read_int(fp)
            colortable['numEntries'] = numEntries
            len_ = read_int(fp)
            colortable['orig_tab'] = fp.read(len_)[:-1].decode('utf-8')

            colortable['struct_names'] = dict()
            colortable['table'] = np.zeros((numEntries, 5), dtype=np.int32)
            numEntriesToRead = read_int(fp)
            for i in range(numEntriesToRead):
                structure = read_int(fp)
                if structure < 0 and verbosity:
                    error_message = 'Error! Read entry, index %d' % structure
                    raise ValueError(error_message)
                if structure in colortable['struct_names'] and verbosity:
                    error_message = 'Error! Duplicate Structure %d' % structure
                    raise ValueError(error_message)
                len_ = read_int(fp)
                colortable['struct_names'][structure] = fp.read(len_)[:-1].decode('utf-8')
                colortable['table'][structure, :4] = np.array(struct.unpack('>4i', fp.read(4 * 4)))
                colortable['table'][structure, 4] = colortable['table'][structure, 0] + (
                        colortable['table'][structure, 1] << 8) + (colortable['table'][structure, 2] << 16) + (
                                                            colortable['table'][structure, 3] << 24)
            if verbosity:
                print('colortable with %d entries read (originally %s)' % (
                    colortable['numEntries'], colortable['orig_tab']))
    else:
        raise ValueError('Error! Should not be expecting bool = 0')

    fp.close()

    #  This makes it so that each empty entry at least has a string, even
    # if it is an empty string.This can happen with average subjects.
    for i in range(numEntries):
        if i not in colortable['struct_names']:
            colortable['struct_names'][i] = ''

    return vertices, label, colortable


def read_electrodes_location(filename):
    electrode_location = []
    electrode_info = {}
    with open(filename, 'r') as f:
        # data at first then info
        begin_info = False
        for line in f:
            # remove '\n'
            line = line.rstrip()
            if line == 'info':
                begin_info = True
                continue
            if not begin_info:
                eloc = [float(i) for i in line.split()]
                electrode_location.append(eloc)
            else:
                info_ = line.split()
                electrode_info[info_[0]] = int(info_[1])
    return np.array(electrode_location), electrode_info


def save_mgh(vol, fname, M, **mr_parms):
    """

    :param vol:
    :param fname:
    :param M: 4x4 vox2ras transform such that y(i1, i2, i3), xyz = M*[i1 i2 i3 1] where indicies are 0-based
    :param mr_parms: [tr flipangle te ti]
    :return:
    """
    # mri scan parameters
    tr = mr_parms.pop('tr', 0)
    flipangle = mr_parms.pop('flipangle', 0)
    te = mr_parms.pop('te', 0)
    ti = mr_parms.pop('ti', 0)

    # for 1 dim vol, pad 1 in the first dimension
    if vol.ndim == 1:
        vol = vol[None]

    # big endian
    fid = open(fname, 'wb')

    fid.write(struct.pack('>i', 1))  # magic

    ndim = []
    for i in range(4):
        if i < vol.ndim:
            ndim.append(vol.shape[i])
            fid.write(struct.pack('>i', vol.shape[i]))
        else:
            ndim.append(1)
            fid.write(struct.pack('>i', 1))

    if vol.ndim == 5:
        fid.write(struct.pack('>i', MRI_TENSOR))
    else:
        fid.write(struct.pack('>i', MRI_FLOAT))

    fid.write(struct.pack('>i', 1))  # dof (not used)

    UNUSED_SPACE_SIZE = 256
    USED_SPACE_SIZE = 3 * 4 + 4 * 3 * 4

    MdcD = M[:3, :3]
    delta = np.sqrt(np.sum(MdcD ** 2, axis=0))

    Mdc = MdcD / delta
    Pcrs_c = [ndim[0] / 2, ndim[1] / 2, ndim[2] / 2, 1]
    Pxyz_c = np.dot(M, Pcrs_c)[:3]

    fid.write(struct.pack('>h', 1))
    fid.write(delta.astype('>f').tostring('F'))
    fid.write(Mdc.astype('>f').tostring('F'))
    fid.write(Pxyz_c.astype('>f').tostring('F'))

    unused_space_size = UNUSED_SPACE_SIZE - 2
    unused_space_size -= USED_SPACE_SIZE
    fid.write(b'\0' * unused_space_size)

    fid.write(vol.astype('>f').tostring('F'))

    fid.write(struct.pack('>f', tr))
    fid.write(struct.pack('>f', flipangle))
    fid.write(struct.pack('>f', te))
    fid.write(struct.pack('>f', ti))

    fid.close()

    if fname.split('.')[-1].lower() in ['mgz', 'gz']:
        # compress file
        cmd = 'gzip -f %s; mv %s.gz %s' % (fname, fname, fname)
        retcode = os.system(cmd)
        if retcode:
            print("Child was terminated by signal", np.abs(retcode))


def save_electrodes_location(elocs, info, filename):
    """

    :param elocs: dict key: index, value: crs coordinates
    :param info: dict
    :param filename:
    :return:
    """
    # make new electrode info
    info = deepcopy(info)
    if info['numpoints'] != len(elocs):
        info['numpoints'] = len(elocs)

    with open(filename, 'w') as f:
        for key in elocs:
            eloc_ = elocs[key]
            f.write('%f %f %f\n' % (eloc_[0], eloc_[1], eloc_[2]))
        # append info
        f.write('info\n')
        f.write('%s %d\n' % ('numpoints', info['numpoints']))
        f.write('%s %d\n' % ('useRealRAS', info['useRealRAS']))


def read_transform_matrix(filename):
    with open(filename, 'r') as f:
        content_ = f.readlines()
        mat_content_ = content_[4:8]
        mat = []
        for l in mat_content_:
            mat.extend(l.rstrip().split())
        mat = ' '.join(mat)
        mat = np.fromstring(mat, dtype=np.float64, sep=' ')
        mat = mat.reshape((-1, 4))
    return mat
