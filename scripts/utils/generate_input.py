#!/usr/bin/env python3
import ruamel.yaml
ruamel.yaml.representer.RoundTripRepresenter.ignore_aliases = lambda x, y: True
yaml = ruamel.yaml.YAML(typ='rt')  # Round trip loading and dumping
yaml.default_flow_style = None
yaml.explicit_start = True
yaml.explicit_end = True
yaml.preserve_quotes = True
yaml.version = '1.1'

###############################################################################
if __name__ == "__main__":
    # Inputs
    ctools_input = {}
    ctools_input['Profiling Output'] = {}
    ctools_input['Profiling Output']['filename_prefix'] = 'cldera_stats'
    ctools_input['Profiling Output']['Flush Frequency'] = 100
    ctools_input['Profiling Output']['Enable Output'] = 'true'
    ctools_input['Profiling Output']['Save Geometry Fields'] = 'true'
    #mask_loc = '/projects/cldera/cldera-tools/input/masks/' # Sandia HPC
    #mask_loc = '/nscratch/jwatkin/e3sm/cldera/model_input/cldera-tools/masks/' # Sandia HPC
    #mask_loc = '' # Sandia HPC (copy)
    #mask_loc = '/sems-data-store/ACME/cldera/cldera-tools/input/masks/' # Mappy
    mask_loc = '/global/cfs/cdirs/m4014/cldera-tools/input/masks/' # Perlmutter
    global_mask = mask_loc + 'global.ne30.nc'
    zonal_mask = mask_loc + 'zonal.regions.ne30.nc'
    ar6_mask = mask_loc + 'ar6.regions.ne30.nc'
    masks2d = {'_glb':global_mask, '_znl':zonal_mask, '_ar6':ar6_mask}
    masks3d = {'_glb':global_mask, '_znl':zonal_mask, '_ar6':ar6_mask}
    #masks2d = {'_ar6':ar6_mask}
    #masks3d = {'_ar6':ar6_mask}
    level_bounds = {'':[0,71], '_ltropo':[54,71], '_strato':[15,21], '_53':[19,19], '_998':[71,71]}
    # vars3d = ['T','u','v','SO201','SO202','SO203','H2SO401','H2SO402','H2SO403','Mass_so401','Mass_so402','Mass_so403',
    #           'E90j','QRL','QRS','ST80_25j']
    #vars3d = ['T','u','v','SO201','SO202','SO203','H2SO401','H2SO402','H2SO403','Mass_so401','Mass_so402','Mass_so403',
    #          'QRL','QRS']
    vars3d = ['T','u','v','SO2','H2SO4','Mass_so4','QRL','QRS']
    #vars2d_aod = ['AEROD_v','ABSORB','AODVIS','AODALL','AODABS','AODSO401','AODSO402','AODSO403']
    vars2d_aod = ['AEROD_v','ABSORB','AODVIS','AODALL','AODABS','AODSO4']
    #vars2d = ['BURDENSO401','BURDENSO402','BURDENSO403','TS','TREFHT','QREFHT','QFLX','SHFLX','LHFLX','FLDS','FLDSC','FLNS','FLNSC','FLNT','FLNTC','FLUT',
    vars2d = ['BURDENSO4','TS','TREFHT','QREFHT','QFLX','SHFLX','LHFLX','FLDS','FLDSC','FLNS','FLNSC','FLNT','FLNTC','FLUT',
              'FLUTC','FSDS','FSDSC','FSDS_d2','FSNS','FSNSC','FSNT','FSNTC','FSNTOA','FSNTOAC','FSUTOA','FSUTOAC',
              'LWCF','SOLIN','SOLL','SOLLD','SOLLD_d2','SOLL_d2','SOLS','SOLSD','SOLSD_d2','SOLS_d2','SWCF']
    ctools_input['Fields To Track'] = vars3d + vars2d_aod + vars2d

    # 3D Inputs
    for var in vars3d:
        ctools_input[var] = {}
        ctools_input[var]['Compute Stats'] = []
        for lvl_reg, lvl_bnd in level_bounds.items():
            for mask_reg, mask_filename in masks3d.items():
                name = var + mask_reg + lvl_reg
                ctools_input[var]['Compute Stats'].append(name)
                ctools_input[var][name] = {}
                ctools_input[var][name]['type'] = 'pipe'
                ctools_input[var][name]['inner'] = {}
                ctools_input[var][name]['inner']['type'] = 'vertical_contraction'
                ctools_input[var][name]['inner']['level_bounds'] = lvl_bnd
                ctools_input[var][name]['inner']['weight_field'] = 'pdel'
                ctools_input[var][name]['outer'] = {}
                ctools_input[var][name]['outer']['type'] = 'masked_integral'
                ctools_input[var][name]['outer']['mask_file_name'] = mask_filename
                ctools_input[var][name]['outer']['mask_field'] = 'mask' + mask_reg
                ctools_input[var][name]['outer']['weight_field'] = 'area'

    # 2D AOD Inputs
    for var in vars2d_aod:
        ctools_input[var] = {}
        ctools_input[var]['Compute Stats'] = []
        for mask_reg, mask_filename in masks2d.items():
            name = var + mask_reg
            ctools_input[var]['Compute Stats'].append(name)
            ctools_input[var][name] = {}
            ctools_input[var][name]['type'] = 'bounded_masked_integral'
            ctools_input[var][name]['mask_file_name'] = mask_filename
            ctools_input[var][name]['mask_field'] = 'mask' + mask_reg
            ctools_input[var][name]['weight_field'] = 'area'
            ctools_input[var][name]['valid_bounds'] = [0.0, 1.0e10]

    # 2D Inputs
    for var in vars2d:
        ctools_input[var] = {}
        ctools_input[var]['Compute Stats'] = []
        for mask_reg, mask_filename in masks2d.items():
            name = var + mask_reg
            ctools_input[var]['Compute Stats'].append(name)
            ctools_input[var][name] = {}
            ctools_input[var][name]['type'] = 'masked_integral'
            ctools_input[var][name]['mask_file_name'] = mask_filename
            ctools_input[var][name]['mask_field'] = 'mask' + mask_reg
            ctools_input[var][name]['weight_field'] = 'area'

    # Write yaml file
    with open('cldera_profiling_config.yaml', 'w') as file:
        yaml.dump(ctools_input, file)
