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
docker run --rm \
    -v ./data:/app/data \
    -v ~/grassdata:/root/grassdata \
    -w /app \
    <name_of_the_image> \
    python3 hydrological_analysis.py \
        --shp ./data/masalia_tehsil_boundary.shp \
        --output masalia \
        --grassdb /root/grassdata
```
