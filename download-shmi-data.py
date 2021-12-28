#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 22 17:46:26 2021

@author: cecilia
"""
# imports
## bulk mkdir for flow sites: mkdir site_{1..20044}
## find and remove empty dirs in folder: find . -type d -empty -delete

import requests
import pandas as pd
import geopandas as gpd

from shapely.geometry import Point


#%% Download all flow data and save to csv

url = 'https://opendata-download-hydroobs.smhi.se/api/version/latest/parameter/1.json'

resp = requests.get(url=url)
data = resp.json()
sites = [ str(x['id']) for x in data['station']]

sites_df = pd.json_normalize(data['station'])
sites_df = sites_df[['id', 'name', 'latitude', 'longitude']]
sites_df['id'] = sites_df['id'].astype(str)

flow_df = pd.DataFrame(columns=['Datum (svensk sommartid)', 'Vattenföring', 'Kvalitet'])

for site_id in sites:
    url='https://opendata-download-hydroobs.smhi.se/api/version/latest/parameter/1/station/' + site_id + '/period/corrected-archive/data.csv'
    print(url)
    df=pd.read_csv(url, 
                   names=['Datum (svensk sommartid)', 'Vattenföring', 'Kvalitet'], 
                   sep=';',
                   skiprows=7,
                   index_col=False)
    df['site_id'] = site_id
    flow_df = flow_df.append(df)
    print(len(flow_df))

flow_df = flow_df.merge(sites_df, how='left', left_on='site_id', right_on='id')
flow_df.to_csv('flow.csv', index=False)

#%% Only keep sites with continuous data for 10+ years

flow_df = pd.read_csv('./flow.csv').iloc[:, 1:]

flow_df.columns = ['date', 'flow', 'quality', 'site_id', 'id',
       'name', 'latitude', 'longitude']

flow_df = flow_df.drop('id', axis=1)
flow_df['date'] = pd.to_datetime(flow_df['date'])

## Are there other quality codes that I can keep?
flow_df = flow_df[flow_df['quality'] == 'G']

freq_df = flow_df.groupby('site_id')['date'].agg([('year_count', lambda x: x.dt.year.max() - x.dt.year.min()), 
                ('obs_per_year', lambda x: x.count() /(x.dt.year.max() - x.dt.year.min()))])
freq_df = freq_df.reset_index()
freq_df['site_id'] = freq_df['site_id'].astype(str)

sites_df.columns = ['site_id', 'name', 'latitude', 'longitude']
flow_sites_df = pd.merge(sites_df, freq_df, on='site_id')
final_flow_sites_df = flow_sites_df[(flow_sites_df['year_count'] >= 10) & (flow_sites_df['obs_per_year'] >= 200)]

final_flow_sites_df.to_csv('flow_sites.csv', index=False)

#%% Create shapefile of suitable flow sites

final_flow_sites_df['geometry'] = final_flow_sites_df.apply(lambda row: Point(row.longitude, row.latitude),axis=1)
gdf = gpd.GeoDataFrame(final_flow_sites_df, geometry = final_flow_sites_df.geometry)
gdf.to_file(driver = 'ESRI Shapefile', crs="EPSG:4326", filename= "./flow_sites/flow_sites.shp")

#%% Download snow depth data and save to csv

url = 'https://opendata-download-metobs.smhi.se/api/version/latest/parameter/8.json'

resp = requests.get(url=url)
data = resp.json()
sites = [ str(x['id']) for x in data['station']]

sites_df = pd.json_normalize(data['station'])
sites_df = sites_df[['id', 'name', 'latitude', 'longitude']]
sites_df['id'] = sites_df['id'].astype(str)

snow_df = pd.DataFrame(columns=['date', 'time', 'snow_depth', 'quality']) 

for site_id in sites:
    url='https://opendata-download-metobs.smhi.se/api/version/latest/parameter/8/station/' + site_id + '/period/corrected-archive/data.csv'
    print(url)
    df=pd.read_csv(url, 
                   names=['date', 'time', 'snow_depth', 'quality'], 
                   sep=';',
                   skiprows=11,
                   index_col=False)
    df['site_id'] = site_id
    snow_df = snow_df.append(df)
    print(len(snow_df))

snow_df = snow_df.merge(sites_df, how='left', left_on='site_id', right_on='id')
snow_df.to_csv('snow.csv', index=False)

#%% Only keep sites with data for 10+ years

snow_df = pd.read_csv('./snow.csv', low_memory=False)

snow_df = snow_df.drop('id', axis=1)

snow_df = snow_df[pd.notna(snow_df['snow_depth'])]
snow_df = snow_df[snow_df['date'] != 'Datum']

snow_df['date'] = pd.to_datetime(snow_df['date'])

## Are there other quality codes that I can keep?
snow_df = snow_df[snow_df['quality'] == 'G']

freq_df = snow_df.groupby('site_id')['date'].agg([('year_count', lambda x: x.dt.year.max() - x.dt.year.min()), 
                ('obs_per_year', lambda x: x.count() /(x.dt.year.max() - x.dt.year.min()))])
freq_df = freq_df.reset_index()
freq_df['site_id'] = freq_df['site_id'].astype(str)

sites_df.columns = ['site_id', 'name', 'latitude', 'longitude']
snow_sites_df = pd.merge(sites_df, freq_df, on='site_id')
final_snow_sites_df = snow_sites_df[(snow_sites_df['year_count'] >= 10) & (snow_sites_df['obs_per_year'] >= 100)]

final_snow_sites_df.to_csv('snow_sites.csv', index=False)

#%% Create snow sites shapefile

final_snow_sites_df['geometry'] = final_snow_sites_df.apply(lambda row: Point(row.longitude, row.latitude),axis=1)
gdf = gpd.GeoDataFrame(final_snow_sites_df, geometry = final_snow_sites_df.geometry)
gdf.to_file(driver = 'ESRI Shapefile', crs="EPSG:4326", filename= "./snow_sites/snow_sites.shp")

#%% Calculating drainage area from waterbody polygons, raster is too computationally intensive

flow_sites = gpd.read_file("./flow_sites/flow_snapped_joined.gpkg")
catchments = gpd.read_file("./catchments/zhyd_fixed.gpkg")

basin_list = []

for ind, row in flow_sites.iterrows():
    
    target_id = row.ZHYD
    basin_id = row.Bas0_ID
    
    #Filter shp by basin id
    basin = catchments[catchments.Bas0_ID == basin_id]
                
    # Select all upstream polygons
    us_wbs = basin[basin.Is_Upstrea == 1]

    valid_poly_ids = set()
    valid_poly_ids.add(target_id)
    potential_poly_ids = []
    
    # Go downstream until end of basin or target polygon found, add ids to set of potential polygons
    for i, row2 in us_wbs.iterrows():
        upstream_poly_ids = set()
        upstream_poly_ids.add(row2.ZHYD)
        current = row2
        while True:
            current_list = [x for i,x in basin[basin.ZHYD == current.NextDown].iterrows()]
            if len(current_list) == 0:
                break
            if len(current_list) > 1:
                raise ValueError('too many ds polygons')
            upstream_poly_ids.add(current.ZHYD)
            if current.Exutoire == 1 or current.NextDown == '' or current.ZHYD == target_id:
                break
            current = current_list[0]
        potential_poly_ids.append(upstream_poly_ids)
    
    # Increase set of valid polygons if any of the existing valid ids is in the potential set
    old_len = 0
    current_len = len(potential_poly_ids)
    while True:
        for i in potential_poly_ids:
            if bool(i & valid_poly_ids):
                valid_poly_ids.update(i)
                potential_poly_ids.remove(i)
        old_len = current_len
        current_len = len(potential_poly_ids)
        if current_len == old_len:
            break
    us_catchment = basin[basin.ZHYD.isin(valid_poly_ids)].dissolve()
    us_catchment['site_id'] = row.site_id
    basin_list.append(us_catchment)
    print(ind)

upstream_gdf = pd.concat([x for x in basin_list]).pipe(gpd.GeoDataFrame)
upstream_gdf.plot()
upstream_gdf.to_file(driver = 'ESRI Shapefile', crs="EPSG:4326", filename= "./catchments/us_catchments.shp")

#%% Download all basin shapefiles and save to shp

url = 'https://opendata-download-hydrography.smhi.se/api/version/2.0/spec/inspire-hyp.json'

resp = requests.get(url=url)
data = resp.json()
basins = [ str(x['key']) for x in data['drainageBasin']]

gdf_list = []

for basin in basins:
    url='https://opendata-download-hydrography.smhi.se/api/version/2.0/spec/inspire-hyp/basin/' + basin + '/smhi.inspire-hyp_rb_' + basin + '-se.v2.4258.gml.zip'
    print(url)
    try:
        gdf=gpd.read_file(url)
        gdf_list.append(gdf)
    except:
        print(basin, 'not found')

basins_gdf = pd.concat([x for x in gdf_list]).pipe(gpd.GeoDataFrame)

basins_gdf.plot()
basins_gdf.to_file(driver = 'ESRI Shapefile', crs="EPSG:4258", filename= "./basins/basins.shp")

#%% Find snow depth sites in each site catchment, and have one row per flow/snow site pair (in QGIS)
## Then find distance between snow and flow site for each polygon and export to csv

us_catchments = gpd.read_file("./catchments/us_catchments_snow.gpkg")
flow_sites = gpd.read_file("./flow_sites/flow_snapped_joined.gpkg")
snow_sites = gpd.read_file("./snow_sites/snow_sites.shp")

us_catchments = us_catchments[us_catchments.snow_site_id.values != None]

## Calculate distance between flow and snow site
for i, row in us_catchments.iterrows():
    flow_site = flow_sites[flow_sites.site_id == row.site_id].geometry.iloc[0]
    snow_site = snow_sites[snow_sites.site_id == row.snow_site_id].geometry.iloc[0]   
    points_df = gpd.GeoDataFrame({'geometry': [flow_site, snow_site]}, crs='EPSG:4326')
    points_df = points_df.to_crs('EPSG:3006')
    points_df2 = points_df.shift() 
    us_catchments.loc[i, 'distance'] = points_df.distance(points_df2).iloc[1]

us_catchments.to_csv('site_pairs.csv', index=False)