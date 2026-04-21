# Hydrological Analysis
Generic pipeline that runs standard hydrology algorithms for a given area of interest

## Setup
First clone/download this Github repository. 

```
git clone https://github.com/mastershubham/hydrology-pipeline
```
Move to the appropriate directory.
```
cd hydrology-pipeline
```

We are using Docker for working in a containerized environment. Ensure that docker is installed on your device and has network connection. 

```
docker build -t <name_of_the_image> .
```
Put up the data somewhere preferably under hydrology-pipeline or some folder under it.

Now run the command as:
```
python hydrological_analysis.py --shp <path_to_shapefile> --output <output_directory> --grassdb <location_of_grass_database>
```
