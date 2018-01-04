#%% Filter groups
filter_mean_x=[100.,400.]
filter_mean_y=[100.,400.]
filter_mean_frame=[29000,33000]
filter_std_frame=[14000,20000]

# Spatially filter (center)
istrue_filter_mean_x=(groupprops['mean_x']>filter_mean_x[0])&(groupprops['mean_x']<filter_mean_x[1])
istrue_filter_mean_y=(groupprops['mean_y']>filter_mean_y[0])&(groupprops['mean_y']<filter_mean_y[1])
istrue_filter_spatial=istrue_filter_mean_x&istrue_filter_mean_y
istrue_filter_spatial_inverse=~np.array(istrue_filter_spatial)

# Filter according to signal stability
istrue_filter_mean_frame=(groupprops['mean_frame']>filter_mean_frame[0])&(groupprops['mean_frame']<filter_mean_frame[1])
istrue_filter_std_frame=(groupprops['std_frame']>filter_std_frame[0])&(groupprops['std_frame']<filter_std_frame[1])
istrue_filter_stable=istrue_filter_mean_frame&istrue_filter_std_frame

istrue_all=istrue_filter_spatial&istrue_filter_stable
#istrue_all=istrue_filter_spatial_inverse&istrue_filter_stable

groupprops_filter=groupprops[:][istrue_filter_stable]
