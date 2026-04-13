'''
GRASS GIS based hydrological analysis
=====================================
Author: Shubham Kumar
Date: April 2026
-------------------------------------
Usage: 
    python hydrological_analysis.py --shp path_to_shapefile \
        --output path_to_output_directory \
        --threshold flow_accumulation_threshold

Description:
    Given a shapefile of a hydrological unit (e.g., a watershed), this script performs the following steps:
    1. Sets up a GRASS GIS environment.
    2. Imports the shapefile and a Digital Elevation Model (DEM) into GRASS GIS.
    3. Fills sinks in the DEM to ensure proper flow direction.
    4. Calculates flow direction and flow accumulation.
    5. Extracts stream networks based on a specified flow accumulation threshold.
    6. Exports the resulting stream network as a GeoJSON file for visualization and further analysis.

'''

import geopandas as gpd
import matplotlib.pyplot as plt
import elevation
from pathlib import Path
import argparse
import shutil
import os
import sys
import subprocess
import time

import rasterio

MIN_WATERSHED_SIZE = 5555 # Translats to roughly 500 hectares at 30m resolution (500 * 10000 / (30*30) = 5555 cells)
# 0. Argument parser

def parse_args():
    parser = argparse.ArgumentParser(
        description="GRASS GIS Hydrological Analysis Pipeline"
    )
    parser.add_argument(
        "--shp", required=True,
        help="Path to input shapefile"
    )
    parser.add_argument(
        "--output", required=True,
        help="Directory for output files"
    )
    parser.add_argument(
        "--grassdb", default=str(Path.home() / "grassdata"),
        help="Directory for GRASS GIS database (default: ~/grassdata)"
    )
    parser.add_argument(
        "--threshold", type=int, default=500,
        help="Flow-accumulation threshold for stream extraction (default: 500 cells)"
    )
    parser.add_argument(
        "--min_watershed_size", type=int, default=1000,
        help="Minimum watershed size in cells (default: 1000)"
    )

    return parser.parse_args()

# 1. GRASS Setup
def setup_grass_session(grassdb: str, epsg: int, project_name: str = "hydro_project"):
    grassdb_path = Path(grassdb).resolve()
    location = project_name
    mapset = "PERMANENT"

    try:
        from grass_session import Session
        session = Session()
        
        # Create location FIRST if it doesn't exist (grass_session can't auto-create)
        loc_path = grassdb_path / location
        if not loc_path.exists():
            grassdb_path.mkdir(parents=True, exist_ok=True)
            grass_bin = shutil.which("grass")
            subprocess.run([
                grass_bin, "-c", f"EPSG:{epsg}", "-e", str(loc_path)
            ], check=True)
            print(f"[INFO] GRASS location created: {loc_path}")
        
        # Now open the existing location
        session.open(
            gisdb=str(grassdb_path),
            location=location,
            mapset=mapset,
        )
        print("grass_session initialised successfully.")
        return session, None
    except ImportError:
        pass

    # Fallback: manual gsetup (your existing code is fine here)
    grass_bin = shutil.which("grass") 
    if grass_bin is None:
        sys.exit("ERROR: GRASS GIS not found. Install it and ensure 'grass' is on PATH.")

    result = subprocess.run([grass_bin, "--config", "python_path"], capture_output=True, text=True)
    grass_python = result.stdout.strip()
    if grass_python and grass_python not in sys.path:
        sys.path.insert(0, grass_python)

    loc_path = grassdb_path / location
    if not loc_path.exists():
        grassdb_path.mkdir(parents=True, exist_ok=True)
        subprocess.run([grass_bin, "-c", f"EPSG:{epsg}", "-e", str(loc_path)], check=True)
        print(f"GRASS location created: {loc_path}")

    import grass.script.setup as gsetup
    gsetup.init(str(grassdb_path), location, mapset)
    print("GRASS environment initialised via grass.script.setup")
    return None, grass_bin

def get_utm_epsg_for_bbox(bbox):
    west, south, east, north = bbox
    
    from pyproj.aoi import AreaOfInterest
    from pyproj.database import query_utm_crs_info
    
    utm_info = query_utm_crs_info(
        datum_name="WGS 84",
        area_of_interest=AreaOfInterest(
            west_lon_degree=west,
            south_lat_degree=south,
            east_lon_degree=east,
            north_lat_degree=north
        )
    )

    if utm_info:
        return utm_info[0].code
    else:
        raise ValueError("No suitable UTM CRS found for the given AOI.")

def fill_sinks(dem_raster):
    print("Filling sinks in DEM...")
    import grass.script as gs
    filled_dem = "dem_filled"
    flow_dir = "flow_dir"
    gs.run_command("r.fill.dir", input=dem_raster, output=filled_dem, direction=flow_dir, overwrite=True)
    return filled_dem, flow_dir

def natural_depressions(dem_filled, dem_original):
    import grass.script as gs
    natural_depressions = "natural_depressions"
    gs.run_command("r.mapcalc", expr=f"{natural_depressions} = {dem_filled} - {dem_original}", overwrite=True)
    return natural_depressions

def calculate_flow_accumulation(dem_filled, threshold):
    import grass.script as gs
    gs.run_command("r.watershed", 
                elevation=dem_filled, 
                accumulation="flow_acc", 
                threshold=MIN_WATERSHED_SIZE,
                basin="micro_watersheds",
                #stream="streams_raw",
                flags="as",          # -a: positive accumulation; -s: single-flow (D8) 
                overwrite=True)
    return "flow_acc", "micro_watersheds"

def export_outputs(output_dir, rasters_to_export: dict, vectors_to_export: dict):
    import grass.script as gs
    for name, raster in rasters_to_export.items():
        output_path = Path(output_dir) / f"{name}.tif"
        gs.run_command("r.out.gdal",
                   input=raster,
                   output=str(output_path),
                   format="GTiff",
                   createopt="COMPRESS=LZW",
                   overwrite=True)
        print(f"Exported raster: {output_path}")

    for name, vector in vectors_to_export.items():
        output_path = Path(output_dir) / f"{name}.geojson"
        gs.run_command("v.out.ogr",
                   input=vector,
                   output=str(output_path),
                   format="GeoJSON",
                   overwrite=True)
        print(f"Exported vector: {output_path}")
    return


# Main function
def main():
    # 0. Arguments parsing
    args = parse_args()

    os.makedirs(args.output, exist_ok=True)

    # Reading shapefile and determining UTM zone
    watershed_gdf = gpd.read_file(args.shp)
    watershed_gdf = watershed_gdf.to_crs(epsg=4326)  

    watershed_gdf.plot(color='white', edgecolor='gray', figsize=(15,12))
    plt.title("Input Watershed Boundary")
    plt.savefig(Path(args.output) / "watershed_boundary.png", dpi=300, bbox_inches='tight')
    plt.show()
    

    minx, miny, maxx, maxy = watershed_gdf.dissolve().total_bounds
    buffer = 0.01 # Nearly 1 km buffer at regions near equator(in degrees)
    bbox = (minx -buffer, miny -buffer, maxx +buffer, maxy +buffer)
    epsg = get_utm_epsg_for_bbox(bbox)

    location_of_dem = Path(args.output) / "dem_raw.tif"
    location_of_dem = location_of_dem.resolve()

    elevation.clip(bounds=bbox, output=str(location_of_dem), product='SRTM1')
    print()
    print(f'DEM raster downloaded and saved to: {location_of_dem}')
    elevation.clean()

    fig, ax = plt.subplots(figsize=(10,8))

    import rasterio
    from rasterio.plot import show
    rasterin = rasterio.open(location_of_dem)

    show(rasterin, ax=ax, cmap='terrain', title='SRTM 30m DEM')
    watershed_gdf.boundary.plot(color="black", ax=ax, linewidth=1.5)
    plt.colorbar(ax.images[0], ax=ax, label='Elevation(m)')
    plt.savefig(Path(args.output) / "dem_with_watershed.png", dpi=300, bbox_inches='tight')
    plt.show()
    # 1. Setup GRASS session
    
    name_of_proj = Path(args.output).resolve().name
    session, grass_bin = setup_grass_session(args.grassdb, epsg, name_of_proj)

    import grass.script as gs

    subprocess.run([
    "gdalwarp",
    "-t_srs", f"EPSG:{epsg}",
    str(location_of_dem),
    str(Path(args.output) / f"dem_{epsg}.tif")
    ])

    gs.run_command("r.in.gdal",
               input=str(Path(args.output) / f"dem_{epsg}.tif"),
               output="dem_utm",
               overwrite=True)

    gs.run_command("g.region", raster="dem_utm", flags="p")
    
    watershed_utm_path = Path(args.output) / "watershed_utm.shp"
    watershed_gdf.to_crs(epsg=epsg).to_file(watershed_utm_path)
    gs.run_command("v.in.ogr",
                   input=str(watershed_utm_path),
                   output="watershed",
                   overwrite=True)
    gs.run_command("r.mask",
                   vector="watershed")
    
    print("DEM imported into GRASS and region set to DEM extent.")

    fig, ax = plt.subplots(figsize=(10, 8))

    dem_path = Path(args.output) / f"dem_{epsg}.tif"
    with rasterio.open(dem_path) as src:
        show(src, ax=ax, cmap='terrain')

    plt.title("DEM in GRASS GIS (UTM)")
    plt.axis('off')
    plt.show()

    
    dem_filled, flow_dir = fill_sinks("dem_utm")
    
    depressions = natural_depressions(dem_filled, "dem_utm")
    
    
    flow_accumulation, micro_watersheds = calculate_flow_accumulation(dem_filled, args.threshold)

    gs.run_command("r.stream.extract",
               elevation=dem_filled,
               accumulation=flow_accumulation,
               stream_raster="streams_rast",
               stream_vector="streams_vect",
               threshold=args.threshold, 
               overwrite=True)

    gs.run_command("r.to.vect",
               input="micro_watersheds",
               output="watersheds_vect",
               type="area",
               overwrite=True)
    
    gs.run_command("r.mask", flags="r") 
    
    rasters_to_export = {
            "dem_filled":           dem_filled,
            "flow_direction":       flow_dir,
            "flow_accumulation":    flow_accumulation,
            "natural_depressions":  depressions
        }
    vectors_to_export = {
            "streams":              "streams_vect",
            "micro_watersheds":     "watersheds_vect"
        }
    export_outputs(args.output, rasters_to_export, vectors_to_export)

if __name__ == "__main__":
    main()
