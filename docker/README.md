# gpusimilarity Docker Container

The preferred method of using Docker is to just to download the publicly
released pre-built docker images.  In addition, the enclosed Dockerfile can be
built to produce an image capable of running GPUSimilarity.

## Requirements
* nvidia-docker (regular docker will not expose GPUs)
* Newest versions of nvidia drivers
* Recommended for production:  systemd / nginx

## How to use pre-built container
Both releases and master are always built and uploaded publicly already.  You
can access them via:
For master:
```
docker pull klorton/gpusimilarity:latest
```

For releases:

```
docker pull klorton/gpusimilarity:v1.0
```

The underlying requirements will change very little, so after the initial
build/pull new docker layers will only require a few megabytes for
gpusimilarity changes.

## How to build (recommended only for development)
From the docker directory:
```
docker build -t gpusim:internal .
```

## How to use to build a fingerprint file
```
docker run -v /PATH/TO/SMIGZ_DIR:/data -it gpusim python3 \
/gpusimilarity/bld/python/gpusim_createdb.py /data/INPUT.smi.gz \
/data/OUTPUT.fsim
```

This will run single threaded.  If you're building huge databases you should
enter the container interactively with the data directory mounted, and follow
the regular instructions from the parent ReadMe to run multithreaded.

## How to start a gpusimilarity server interactively
```
/bin/nvidia-docker run --net=host -v /path/to/fsim/files:/mnt/fsim -it \
klorton/gpusimilarity:latest python3 /gpusimilarity/bld/python/gpusim_server.py \
/mnt/fsim/1.fsim --port 8080 --http_interface
```

Once the server says "Ready for Searches", you should be able to test it by
either using a graphical web browser or links to access "http://localhost:8080"

## Example of a systemd configuration file
See [the provided service file](gpusimilarity.service) to be placed in `/etc/systemd/system/`

## Additional steps you'll want to take for production
* The `--http_interface` option is only for testing **never** use it in
  production as it hasn't been hardened at all
* You should use a systemd script to start/stop the service, including restarts
* **Never** expose the http from this service directly, use nginx or equivalent with certificates to expose via ssl/https.
