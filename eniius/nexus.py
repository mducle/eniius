from nexusformat.nexus import *
import nexusformat.nexus as nexus
import mcstasscript
import numpy as np
import warnings
import json
import sys
import os

from .mcstas import NX2COMP_MAP, AffineRotate, NXoff

comps_path = os.path.abspath(os.path.join(os.path.dirname(__file__), 'mcstas-comps'))

def get_nx_component(nxobj, nxtype=None, nxname=None):
    # Recursively looks through a NeXus object for a specific component type or name
    for k, v in nxobj.items():
        if nxtype is not None and isinstance(v, nxtype):
            return v
        elif nxname is not None and v.nxname == nxname:
            return v
        elif hasattr(v, 'keys') and len(v.keys()) > 0:
            v = get_nx_component(v, nxtype, nxname)
            if v is not None:
                return v


class NXinst2McStas():
    # Class to convert a NeXus component to a McStas one

    def __init__(self, instname, nx_inst):
        self.instname = instname
        self.mc_inst = mcstasscript.McStas_instr(self.instname, package_path=comps_path)
        self.nx_inst = nx_inst
        self.comps = []

        for label, comp in self.nx_inst.items():
            if not hasattr(comp, 'entries'):
                continue
            if 'mcstas' in comp.entries:
                comp_pars, comp_name, comp_ord = self._nx2mc_previous(label, json.loads(comp.mcstas.nxvalue))
            else:
                comp_ord = None
                try:
                    comp_pars, comp_name, comp_pos = self._nx2mc_general(label, comp)
                except RuntimeError as err:
                    warnings.warn(str(err))
                    continue
            nxtransform = get_nx_component(comp, nxtype=nexus.NXtransformations)
            if nxtransform:  # NXtransformations overwrite parameters defined by component
                comp_pos = self._get_pos_from_transform(nxtransform)
            self.comps.append([label, comp_pars, comp_name, comp_pos, comp_ord]) 
        for idc in self._get_order():
            label, comp_pars, comp_name, comp_pos = tuple(self.comps[idc][:4])
            for idx, cp in enumerate(comp_pars):
                mc_comp = self.mc_inst.add_component(label if idx==0 else f'{label}{idx}', comp_name)
                mc_comp.set_parameters(**cp)
                for posdat in comp_pos:
                    # posdat is of the form [method, [values]], method is e.g. set_AT, set_ROTATED
                    getattr(mc_comp, posdat[0])(*posdat[1])

    def _get_order(self):
        order = []
        for comp in self.comps:
            # Ignore cases where mcstas_order is negative (eniius internal components)
            if comp[4] is not None and comp[4] >= 0:
                order.append(comp[4])
            elif comp[4] is None:
                order.append(self._get_distance_of_comp(comp[2]))
        return np.argsort(order)

    def _get_distance_of_comp(self, comp_name):
        comp_pos = [cc[3] for cc in self.comps if cc[2] == comp_name][0]
        if len(comp_pos) == 0:
            return 0
        distvec = comp_pos[0][1][0]
        sgn = np.sign(distvec[-1] if isinstance(distvec, list) else distvec)
        dist = np.sqrt(np.sum(np.square(distvec))) * sgn
        rel_name = comp_pos[0][1][1] if (len(comp_pos[0][1]) > 1) else None
        while rel_name is not None:
            comp_pos = [cc[3] for cc in self.comps if cc[2] == rel_name][0]
            distvec = comp_pos[0][1][0]
            sgn = np.sign(distvec[-1] if isinstance(distvec, list) else distvec)
            dist = dist + (np.sqrt(np.sum(np.square(distvec))) * sgn)
            rel_name = comp_pos[0][1][1] if (len(comp_pos[0][1]) > 1) else None
        return dist

    def _get_affinelist_from_transform(self, nxtransform):
        dep_dict = {v.depends_on:k for k, v in nxtransform.entries.items()}
        # Construct chain on the "depends_on" keyword.
        init_node = [k for k,v in dep_dict.items() if k == '.' or k not in nxtransform]
        # Prefer absolute, other use a node which is outside this set of transformations
        init_node = '.' if '.' in init_node else init_node[0]
        order = [dep_dict[init_node]]
        while order[-1] in dep_dict:
            order.append(dep_dict[order[-1]])
        transform = np.eye(4)
        transform_list, transform_added, dep = ([], False, init_node)
        for name in order:
            transform = np.matmul(AffineRotate.from_nxfield(nxtransform[name]).transform, transform)
            try:
                tr = AffineRotate(transformation_matrix=transform, depends_on=dep)
                transform_added = False
            except AssertionError:
                transform_list.append(tr)
                transform_added = True
                dep = nxtransform[name].depends_on
        if not transform_added:
            transform_list.append(tr)
        return transform_list

    def _get_pos_from_transform(self, nxtransform):
        transform_list = self._get_affinelist_from_transform(nxtransform)
        if len(transform_list) > 1:
            raise NotImplementedError("TODO: implement multiple ARMs for complex transformations")
        tr = transform_list[0]
        relative = tr.depends_on if tr.depends_on != '.' else None
        comp_pos = [['set_AT', [list(tr.transform[:3, 3]), relative]]]
        if tr.is_rotation:
            comp_pos.append(['set_ROTATED', [list(AffineRotate.get_euler_angles(tr.transform[:3, :3])), relative]])
        return comp_pos

    def _nx2mc_previous(self, label, compdict):
        # Recreate a previously saved McStas component from the "mcstas" field
        comp_name = compdict.pop('mcstas_component')
        comp_order = compdict.pop('mcstas_order')
        return [compdict], comp_name, comp_order

    def _nx2mc_general(self, label, comp):
        # Tries to convert a general NX component to McStas (without McStas field)
        comp_pos = []
        try:
            comp_pars, comp_name, comp_pos = getattr(self, comp.nxclass)(**comp.entries)
        except AttributeError:
            if comp.nxclass in NX2COMP_MAP:
                compinfo = NX2COMP_MAP[comp.nxclass]
                comp_name = compinfo[0]
                comp_pars = [{v:comp.entries[k].nxdata for k, v in compinfo[1].items() if k in comp.entries}]
                if len(compinfo) > 2:
                    comp_pos = [[v, [comp.entries[k].nxdata]] for k, v in compinfo[2].items() if k in comp.entries]
            else:
                raise RuntimeError(f'NeXus component "{comp}" type "{comp.nxclass}" could not be converted to McStas')
        return comp_pars, comp_name, comp_pos

    @classmethod
    def NXdisk_chopper(cls, **kw):
        params = {'nu':kw['rotation_speed'].nxdata, 'radius':kw['radius'].nxdata}
        params['yheight'] = kw['slit_height'].nxdata if 'slit_height' in kw else params['radius'] / 10
        if all([p in kw.keys() for p in ['slit_edges', 'slits']]):
            warnings.warn('NXdisk_chopper defines both "slits" and "slit_edges". '
                          'Will use "slit_edges" and ignore "slit" and/or "slit_angle"')
        if 'slit_edges' in kw:
            sl = np.array(kw['slit_edges'].nxdata)
            slit_wd = np.diff(sl)
            slit_cen = (sl[1:] + sl[:-1]) / 2.
            pars = [{'theta_0':slit_wd[i], 'phase':slit_cen[i], 'nslit':1, **params} for i in range(len(slit_cen))]
        elif 'slit' in kw:
            slit_wd = kw['slit_angle'].nxdata if 'slit_angle' in kw else 5
            pars [{'theta_0':slit_wd, 'nslit':kw['slit'].nxdata, **params}]
        return pars, 'DiskChopper', []

    def show_components(self):
        self.mc_inst.show_components()

