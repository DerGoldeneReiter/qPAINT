
props_dir='/fs/pool/pool-schwille-paint/Data/p04.lb-FCS/19-06-05_N=12/id63_5nM_p35uW_2/19-06-05_FS'
props_name='id63_5nM_p35uW_2_MMStack_Pos0.ome_locs_render_picked_props_ig1.hdf5'

align_dir='/fs/pool/pool-schwille-paint/Data/p04.lb-FCS/19-06-05_N=12/id63_5nM_p1700uW_1/19-06-05_FS'
align_name='id63_5nM_p1700uW_1_MMStack_Pos0.ome_locs_render_align_picked.hdf5'

#################################################### Load packages
import os #platform independent paths
import importlib
# Load user defined functions
import var_io as io
import picasso.io

importlib.reload(io)

#%%

############################################################# Read locs, apply props & save locs
###  Define paths
props_path=os.path.join(props_dir,props_name)
align_path=os.path.join(align_dir,align_name)

### File read in
props,props_info=io.read_locs(props_path)
align,align_info=io.read_locs(align_path)

### Limit align to filter groups in props
align.set_index('group',inplace=True)
props_groups=props.group.values
align_filter=align.loc[props_groups,:]

### Reset index
align_filter.reset_index(inplace=True)

### Convert to DataFrame to numpy.recarray
align_filter_rec=align_filter.to_records(index=False)

### Save filtered picked file
savepath=align_path.replace('.hdf5','_filter.hdf5')
picasso.io.save_locs(savepath,align_filter_rec,align_info)