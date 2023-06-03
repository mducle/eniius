import numpy as np
import nexusformat.nexus
import warnings
import json

VERSION = '0.1'

# Monkey patch the nexus write function to use fixed-width ASCII NX_class labels
# because IBEX [ISISICP] only supports this format (as it uses the old napi.c)
old_wg = nexusformat.nexus.tree.NXFile._writegroup
old_wd = nexusformat.nexus.tree.NXFile._writedata
def new_wg(self, group):
    links = old_wg(self, group)
    if group.nxpath != '' and group.nxpath != '/':
        path = self.nxpath + '/' + group.nxname
        if group.nxclass and group.nxclass != 'NXgroup':
            self[path].attrs['NX_class'] = np.array(group.nxclass, dtype='S')
    return links
def new_wd(self, data):
    if data._target is None and data._uncopied_data is None and data._memfile is None \
            and (not data.nxfile or data.nxfile.filename == self.filename) \
            and data.dtype is not None and data.is_string():
        data._shape = (1,)
        data._dtype = f'S{len(data._value)}'
        data._value = np.array(data._value, dtype='S')
    return old_wd(self, data)
nexusformat.nexus.tree.NXFile._writegroup = new_wg
nexusformat.nexus.tree.NXFile._writedata = new_wd

from nexusformat.nexus import *


def conv_types(obj):
    typ = type(obj)
    dtyp = np.dtype(typ)
    if dtyp.name == 'object':
        if isinstance(obj, np.ndarray):
            vl = obj.tolist()
            el = vl[0]
            for ii in range(1, len(obj.shape)):
                el = el[0]
            tp = np.dtype(type(el)).name
        elif obj is None:
            (tp, vl) = ('string', 'None')
        elif isinstance(obj, NXattr):
            val = obj.nxdata
            if isinstance(val, np.ndarray):
                val = val.tolist()
            if obj.dtype == 'object':
                (tp, vl) = (np.dtype(type(obj.nxdata)).name, val)
            else:
                (tp, vl) = (obj.dtype, val)
        else:
            raise RuntimeError(f'unrecognised type {typ} / {dtyp}')
    else:
        (tp, vl) = (dtyp.name, obj)
    if tp == 'str':
        tp = 'string'
    elif tp == 'float64':
        tp = 'float'
    elif tp == 'object':
        raise RuntimeError('Internal logical error')
    return tp, vl


class Writer:
    """
    Writes out files in various formats from a NeXus structure with instrument information
    """

    def __init__(self, nxobj):
        self.rootname = 'root'
        self.data = None
        self.sample = None
        self.nxobj = None
        if nxobj.nxclass == 'NXinstrument':
            self.inst = nxobj
            nxentry = NXentry()
            nxentry['instrument'] = self.inst
            self.sample = NXsample()
        else:
            if nxobj.nxclass == 'NXroot':
                entries = dir(nxobj)
                if len(entries) > 1:
                    warnings.warn('More than one entry in nxobject, using only first entry')
                self.rootname = entries[0]
                self.nxobj = nxobj
                nxentry = nxobj[self.rootname]
            elif nxobj.nxclass == 'NXentry':
                nxentry = nxobj
            else:
                raise RuntimeError('Input must be an NXroot, NXentry or NXinstrument')
            if 'instrument' not in dir(nxentry):
                raise RuntimeError('nxobj does not have an instrument entry')
            self.inst = nxentry['instrument']
            self.data = nxentry['data'] if 'data' in dir(nxentry) else None
            self.sample = nxentry['sample'] if 'sample' in dir(nxentry) else None
        if self.nxobj is None:
            self.nxobj = NXroot()
            self.nxobj[self.rootname] = nxentry


    def to_json(self, invar):
        # Converts a NeXus object to a JSON-compatible dictionary
        outfile, nxobj = (invar, self.nxobj) if isinstance(invar, str) else (None, invar)
        children = []
        for k, obj in nxobj.items():
            if hasattr(obj, 'nxclass'):
                if obj.nxclass == 'NXfield':
                    typ, val = conv_types(obj.nxdata)
                    entry = {'module':'dataset',
                             'config':{'name':k, 'values':val, 'type':typ}}
                    attrs = []
                else:
                    entry = {'name':k, 'type':'group'}
                    attrs = [{'name':'NX_class', 'dtype':'string', 'values':obj.nxclass}]
                    if len(obj._entries) > 0:
                        entry['children'] = self.to_json(obj)
            else:
                raise RuntimeError(f'unrecognised object key {k}')
            for n, v in obj.attrs.items():
                typ, val = conv_types(v)
                attrs.append({'name':n, 'dtype':typ, 'values':val})
            if len(attrs) > 0:
                entry['attributes'] = attrs
            children.append(entry)
        if outfile is not None:
            if not outfile.endswith('.json'):
                outfile += '.json'
            with open(outfile, 'w') as f:
                f.write(json.dumps({'children':children}, indent=4))
        else:
            return children


    def to_nxspe(self, outfile, ei=25, det_file=None):
        if not outfile.endswith('.nxspe'):
            outfile += '.nxspe'
        with nxopen(outfile, 'w') as root:
            root['w1'] = NXentry()
            root['w1/definition'] = NXfield('NXSPE', version='1.3')
            root['w1/NXSPE_info'] = NXcollection(fixed_energy=NXfield(ei, units='meV'),
                                                           ki_over_kf_scaling=True,
                                                           psi=NXfield(np.nan, units='degrees'))
            root['w1/instrument'] = self.inst
            if self.sample is not None:
                root['w1/sample'] = self.sample
            if self.data is None:
                n = 5
                dmat = np.random.rand(n, n)
                emat = np.random.rand(n, n) / 10.
                th = np.linspace(-25, 120, n)
                phi = np.zeros(n)
                en = np.linspace(ei/20, ei*0.9, n)
                wd = np.ones(n) / 10.
                dd = np.ones(n) * 4.
                dat = NXdata(data=NXfield(dmat, axes='polar:energy', signal=1), error=emat,
                             errors=emat, energy=NXfield(en, units='meV'),
                             azimuthal=th, azimuthal_width=wd, polar=phi, polar_width=wd, distance=dd)
                root['w1/data'] = dat
            else:
                root['w1/data'] = self.data
            if det_file:
                for ky, val in self._parse_det(det_file).items():
                    root[f'w1/{ky}'] = val


    def to_icp(self, outfile, det_file=None):
        if not outfile.endswith('.nxs'):
            outfile += '.nxs'
        with nxopen(outfile, 'w') as root:
            root['mantid_workspace_1'] = NXentry()
            root['mantid_workspace_1/program_name'] = NXfield('eniius', version=VERSION)
            root['mantid_workspace_1/instrument'] = self.inst
            if det_file:
                for ky, val in self._parse_det(det_file).items():
                    root[f'mantid_workspace_1/{ky}'] = val


    def _parse_det(self, det_file):
        # Assumes ISIS format; cols=(det#, delta, L2, code, theta, phi, W_xyz, a_xyz, det_123)
        with open(det_file, 'r') as f:
            titles = [next(f) for x in range(3)][2].split()
            if titles[0].startswith('det') and titles[1].startswith('no'):
                titles = titles[1:]
            titles = ','.join(titles[6:])
        detdat = np.loadtxt(det_file, skiprows=3)
        def _getsubset(detdat, fnm, idx):
            fd = {f'number_of_{fnm}': NXfield(np.array(len(idx), dtype='uint64'))}
            fd['detector_number'] = NXfield(detdat[idx,0].astype(np.int32))
            fd['detector_offset'] = NXfield(detdat[idx,1])
            fd['distance'] = NXfield(detdat[idx,2], units='metre')
            fd['polar_angle'] = NXfield(detdat[idx,4], units='degree')
            fd['azimuthal_angle'] = NXfield(detdat[idx,5], units='degree')
            fd['user_table_titles'] = NXfield(titles)
            fd['code'] = NXfield(detdat[idx,3])
            for j in range(6, detdat.shape[1]):
                fd[f'user_table_{j-5}'] = NXfield(detdat[idx,j])
            assert len(titles.split(',')) == detdat.shape[1]-6, "Number of titles not commensurate with table width"
            return fd
        return {'instrument/physical_monitors': NXdetector(**_getsubset(detdat, 'monitors', np.where(detdat[:,3] == 1)[0])),
                'instrument/physical_detectors': NXdetector(**_getsubset(detdat, 'detectors', np.where(detdat[:,3] != 1)[0]))}
