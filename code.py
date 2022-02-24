try:
    from matplotlib import pyplot as plt # visualization
    from termcolor import colored # print colored text
    import cartopy.crs as ccrs # projected visualizations
    import matplotlib.colors as colors # colors for visualizations
    import xarray as xr # process netCDF s-5p data
    import numpy as np # data manupulation
    import cartopy # improved visualizations
    import matplotlib.gridspec as gridspec # create subplot
    from glob import iglob # data access in file manager
    from os.path import join # same
    from functools import reduce # string manipulation
    import pandas as pd # data manipulation
    import harp # preprocess L2 to L3 Sentinel5p data
    from harp._harppy import NoDataError, CLibraryError # import specific error that may happen when processing L2 to L3
    import itertools # dict manipulation
    import cartopy.feature as cf
    import cartopy.io.shapereader as shpreader
    from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
    import matplotlib.patches as mpatches
    from shapely.ops import cascaded_union
    import shapely.vectorized
    from os.path import exists
    import imageio # create gif
    import ipywidgets as widgets
    import h5py
except ModuleNotFoundError:
    print ('Module import error')
else:
    print(colored('\nAll libraries properly loaded. Ready to start!!!','green'),'\n')



    # look for NO2 products in target folder
product_path = "E:\Mini Project\Data"
input_files_OFFL = sorted(list(iglob(join(product_path, '**','*OFFL*NO2*.nc'), recursive=True)))
print(colored('NO2 OFFL products detected:','blue'), len(input_files_OFFL))

#select level 2 NO2 Pproduct to be explored
s5p_file = input_files_OFFL[0]
print(colored('product selected for analysis:','blue'),s5p_file)



# open - Global attributes
with xr.open_dataset(s5p_file) as s5p_img_GA:
    print (colored('Global attributes of product:\n' , 'blue'), s5p_img_GA)
    
# open product-GROUP METADATA/GRANULE_DESCRIPTION
with xr.open_dataset(s5p_file, group='METADATA/GRANULE_DESCRIPTION') as s5p_img_MT:
    print(colored('\nMETADATA/GRANULE_DESCRIPTION Group:\n' , 'blue'), s5p_img_MT)
    
# open product-GROUP PRODUCT
with xr.open_dataset(s5p_file, group = 'PRODUCT') as s5p_img_PRD:
    print(colored('\n PRODUCT Group:\n' , 'blue'), s5p_img_PRD)



    # open product - GROUP PRODUCT
NO2 = s5p_img_PRD['nitrogendioxide_tropospheric_column']
print(colored('Dimension: names for each axis(e.g.,("x", "y", "z"):','blue'), NO2.dims)
print(colored('\ncoordinates: dick-like container of arrays that label each point: ','blue'), NO2.coords)
print(colored('\nAttributes: dict to hold arbitary metadata (attributes):\n','blue'), NO2.attrs)
print(colored('\nValues: a numpy.ndarray holding the array values:\n','blue'), NO2.values)



#convert values to molecules/cm2
NO2 = NO2 * NO2. multiplication_factor_to_convert_to_molecules_percm2

# First visualization.Simple plot
plt.figure(figsize=(14,6))
ax= plt.axes()
NO2[0].plot.pcolormesh(ax=ax, x='longitude', y='latitude', add_colorbar=True, cmap='magma_r', vmin=0);



# Second visualization using Cartopy and Orthographic projection
plt.figure(figsize=(17,12))
ax = plt.axes(projection=ccrs.Orthographic(11, 42))
NO2[0].plot.pcolormesh(ax=ax, x='longitude', y='latitude', add_colorbar=True, cmap='magma_r',transform=ccrs.PlateCarree(), vmin=0)
ax.add_feature(cartopy.feature.RIVERS)
ax.set_title('S-5p L2 NO$_2$ (2020-04-01) | qa_value > 0')
ax.coastlines('10m')
#ax.set_global()
ax.stock_img()
ax.gridlines()
                       
# Save figure to file
plt.savefig("E:\Mini Project\Processing\one", bbox_inches='tight', dpi=300)



NO2_filter = NO2.where(s5p_img_PRD['qa_value'] > 0.75, drop=True)

# Plot original data vs quality flag
plt.figure(figsize=(30,10))

# Plot qa_value image
ax1 = plt.subplot(121, projection=ccrs.Orthographic(11, 42))
s5p_img_PRD['qa_value'][0].plot.pcolormesh(ax=ax1, x='longitude', y='latitude', add_colorbar=True, cmap='Spectral', transform=ccrs.PlateCarree())
ax1.set_title('S-5p L2 NO$_2$ (2021-12-01) | Quality Flag Layer')
ax1.add_feature(cartopy.feature.LAND, edgecolor='black')
ax1.add_feature(cartopy.feature.OCEAN)
ax1.add_feature(cartopy.feature.RIVERS)
ax1.coastlines('10m')
ax1.gridlines()

# Plot masked NO2 data
ax2 = plt.subplot(122, projection=ccrs.Orthographic(11, 42))
NO2_filter[0].plot.pcolormesh(ax=ax2, x='longitude', y='latitude', add_colorbar=True, cmap='Spectral', transform=ccrs.PlateCarree(), vmin=0)
ax2.set_title('Filtered S-5p L2 NO$_2$ (2020-04-01) | Quality Flag Layer > 0.75')
ax2.add_feature(cartopy.feature.RIVERS)
ax2.coastlines('10m')
ax2.stock_img()
ax2.gridlines()

# Save figure to file
plt.savefig("E:\Mini Project\Processing\QualityFilter.png", bbox_inches='tight', dpi=300)




# Define AQI (coordinates UR and LL corners lat/lon order) and subset product
ur = (13, 77.5)
ll = (8, 74.80)

NO2_subset = NO2_filter.where((NO2_filter.longitude < ur[1]) & (NO2_filter.longitude > ll[1]) & (NO2_filter.latitude > ll[0]) & (NO2_filter.latitude < ur[0]), drop=True)

# Plot subset
plt.figure(figsize=(12,12))
ax = plt.axes(projection=ccrs.Orthographic(11, 42))

stt_prv = cf.NaturalEarthFeature(category='cultural',
    name='admin_1_states_provinces_lines',
    scale='10m',facecolor='none')

NO2_subset[0].plot.pcolormesh(ax=ax, x='longitude', y='latitude', add_colorbar=True, cmap='magma_r',transform=ccrs.PlateCarree(), vmin=0)
ax.add_feature(cartopy.feature.LAND, edgecolor='black')
ax.add_feature(stt_prv, linewidth=0.2, edgecolor='black')
ax.set_title('NO$_2$ Concentration over Kerala (2021-12-01) | qa_value > 0.75')

# Add POI to map
ax.text(10.10, 76.25, 'Ernakulam', transform=ccrs.Geodetic());

# Save figure to file
plt.savefig("E:\Mini Project\Processing\Ernakulam.png", bbox_inches='tight', dpi=300)






# Export path
export_path = "E:\Mini Project\Processing"

for i in input_files_OFFL:
    harp_L2_L3 = harp.import_product(i, operations=" \
                             tropospheric_NO2_column_number_density_validity>75; \
                             derive(tropospheric_NO2_column_number_density [Pmolec/cm2]); \
                             derive(datetime_stop {time}); \
                             latitude > 13 [degree_north] ; latitude < 16.6 [degree_north] ; longitude > 77.5 [degree_east] ; longitude < 83.6 [degree_east]; \
                             bin_spatial(360, 13, 0.01, 610, 77.5, 0.01); \
                             derive(latitude {latitude}); derive(longitude {longitude});\
                             keep(NO2_column_number_density, tropospheric_NO2_column_number_density, stratospheric_NO2_column_number_density, NO2_slant_column_number_density, tropopause_pressure, absorbing_aerosol_index, cloud_fraction, sensor_altitude, sensor_azimuth_angle, sensor_zenith_angle, solar_azimuth_angle, solar_zenith_angle); \
                             ")
    export_folder = "{export_path}/{name}".format(export_path=export_path, name=i.split('/')[-1].replace('L2', 'L3'))
    harp.export_product(harp_L2_L3, export_folder, file_format='netcdf')
    
print(colored('All L2 products converted to L3','green'))




# Save two original attributes from the the L2 product (time_coverage_start and time_coverage_end)
export_path = "E:\Mini Project\Data"
attributes = {
         i.split('/')[-1]: {
             'time_coverage_start': xr.open_dataset(i).attrs['time_coverage_start'],
             'time_coverage_end': xr.open_dataset(i).attrs['time_coverage_end'],
         } for i in input_files_OFFL}
# Print only first items in dictionary as example
dict(itertools.islice(attributes.items(), 1))




# Create list of L3 products
product_path = "E:\Mini Project\Processing\L2 to L3"
filename_L3 = sorted(list(iglob(join(export_path, '**', '*L3__NO2*.nc'), recursive=True)))

# Create a time coordinate with np datatype determine64. Important to allow time indexing later
def preprocess(ds):
    ds['time'] = pd.to_datetime(np.array([attributes[ds.attrs['source_product']]['time_coverage_start']])).values
    return ds

L3_APRL_20_21 = xr.open_mfdataset(filename_L3, combine='nested', concat_dim='time', preprocess=preprocess, chunks={'time': 100})
L3_APRL_20_21



product_path = "E:\Mini Project\Processing\L2 to L3"
L3_APRL_20 = L3_APRL_20_21.sel(time=slice('2020-04-01', '2020-04-07'))
L3_APRL_20 = L3_APRL_20.resample(time='ID').mean(dim='time', skipna=None)

L3_APRL_21 = L3_APRL_20_21.sel(time=slice('2021-04-01', '2021-04-07'))
L3_APRL_21 = L3_APRL_21.resample(time='ID').mean(dim='time', skipna=None)
                                           
print(colored('Dataset for April 2020:\n', 'blue'),L3_APRL_20, colored('\n\nDataset for April 2021:\n', 'blue'),L3_APRL_21)




# Mean value for April 2020
L3_APRL_20_mean = L3_APRL_20.mean(dim='time')
NO2_APRL_20_mean = L3_APRL_20_mean['tropospheric_NO2_column_number_density']

# Mean value for April 2021
L3_APRL_21_mean = L3_APRL_21.mean(dim='time')
NO2_APRL_21_mean = L3_APRL_21_mean['tropospheric_NO2_column_number_density']

print(colored('Dataset for April 2020:\n','blue'),NO2_APRL_20_mean, colored('\n\nDataset for April 2021:\n','blue'),NO2_APRL_21_mean)




# Select a month to display
data = NO2_APRL_21_mean # change by 'NO2_APRL_20_mean' accordingly
year = 2021 # change by '2020' accordingly



# Image desing
fig = plt.figure(figsize=(18, 6))

# Main map
ax = fig.add_subplot(1,1,1, projection=ccrs.PlateCarree())
ax.set_extent([6, 15, 43, 47.5])

states_provinces = cf.NaturalEarthFeature(
    category='cultural',
    name='admin_0_countries',
    scale='10m',
    facecolor='#DEDEDE')

im = data.plot.pcolormesh(ax=ax, transform=ccrs.PlateCarree(), cmap='magma_r', vmin=0, vmax=15, x='longitude', y='latitude', zorder=3)
im.colorbar.remove()

# Add text
ax.text(0, 1.07, 'Average NO2 concentrations', fontsize = 17, transform=ax.transAxes)
ax.text(0, 1.02, 'Kerala, December {}'.format, fontsize = 13, transform=ax.transAxes)
ax.text(0.61, -0.21, "Data: ESA Sentinel-5p / TROPOMI\nCredits: Contains Copernicus data (2021) processed bu RUS Copernicus", fontsize=12, color='gray', multialignment='right', transform=ax.transAxes)

# Add countries
ax.add_feature(states_provinces, edgecolor='black')

ax.coastlines("10m", zorder=3);
ax.add_feature(cartopy.feature.BORDERS.with_scale('10m'),zorder=3)

# Set colorbar properties
cbar_ax = fig.add_axes([0.25, -0.015, 0.25, 0.01])
cbar = plt.colorbar(im, cax=cbar_ax, orientation='horizontal', ticks=[0,1,5,10,15])
cbar.set_label(r"$10^{15}$ molecules / cm$_2$)", labelpad=-50, fontsize=14)
cbar.outline.set_visible(False)
cbar.ax.set_yticklabels(['0','1','5','10','15'])

# Set plot frame color
gl = ax.gridlines(draw_labels=True, linewidth=1, color='gray', alpha=0.3, linestyle=':')
gl.xlabels_top = False
gl.ylabels_right = False
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER

# Kerala map
ax = fig.add_subplot(2,5,10, projection=ccrs.PlateCarree()) #nrwo ncol index
ax.set_extent([5,19,36,49])
ax.add_feature(states_provinces, linewidth=0.5, edgecolor='black')
ax.add_patch(mpatches.Rectangle(xy=[13, 77.5], width=5.9, height=2.9, facecolor='blue', alpha=0.3, transform=ccrs.PlateCarree(), zorder=3))

# Save figure to file
plt.savefig("E:\Mini Project\processing.png".format(year), bbox_inches='tight', dpi=600);




# Create a gif to show the difference
filenames = []
images = []
for filename in filenames:
    images.append(imageio.imread(filename))
imageio.mimsave("E:\Mini Project\Processing\gif", images, fps=1)



