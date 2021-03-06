############################################################# Load packages
import os
import matplotlib.pyplot as plt
import picasso.io
import picasso.render
import autopick

plt.style.use('~/qPAINT/styles/paper.mplstyle')

############################################################## Define data
render_dir='/fs/pool/pool-schwille-paint/Data/p04.lb-FCS/19-06-18_N=48/id133_2-5nM_p35uW_1/19-07-01_FS'
render_name='id133_2-5nM_p35uW_1_MMStack_Pos0.ome_locs_render.hdf5'

############################################################## Load data
#### Paths
render_path=os.path.join(render_dir,render_name)

#### Load rendered locs
locs_render,info_render=picasso.io.load_locs(render_path)

#### Render image
oversampling=5
min_blur_width=0.02
print('... rendering with oversampling of %i'%(oversampling))
n_locs,image_render=picasso.render.render(
        locs_render,
        info_render,
        oversampling=oversampling,
        viewport=None,
        blur_method="gaussian",
        min_blur_width=min_blur_width
        )

#%%
############################################################# Set parameters
render_box = 9
render_mng = 300

pick_diameter=2

############################################################## Automated pick detection
picklocs_render=autopick._autopick(image_render,render_mng,render_box,fit=True)

#%%
############################################################## View picks
#### Set view in real pixels (no oversampling)
view_x=0
view_y=-10
view_width=40

f=plt.figure(num=12,figsize=[7,7])
f.subplots_adjust(bottom=0.,top=1.,left=0.,right=1.)
f.clear()
ax=f.add_subplot(111)
ax.imshow(image_render,cmap='jet',vmin=0,vmax=15,interpolation='nearest')
ax.scatter(picklocs_render.x,picklocs_render.y,s=200,marker='+',alpha=1,color='w',linewidths=1)
ax.set_xlim(view_x*oversampling,(view_x+view_width)*oversampling)
ax.set_ylim(view_y*oversampling,(view_y+view_width)*oversampling)
ax.grid(False)

#%%
############################################################### Save picks
save_path=render_path.replace('.hdf5','_autopicks.yaml')
picklocs=picklocs_render.copy()
picklocs.x=(picklocs_render.x)/oversampling
picklocs.y=(picklocs_render.y)/oversampling
autopick._locs2picks(picklocs,pick_diameter,save_path)

############################################################### Save yaml
info_autopick={'box_size':render_box,'mng':render_mng,'pick_diameter':pick_diameter}
import var_io
info_list=[info_render,info_autopick]
var_io.create_yaml(info_list,save_path.replace('.yaml','_log.yaml'))

