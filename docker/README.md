# gpusimilarity Docker Container

The preferred method of using Docker is to just to download the publicly
released pre-built docker images.  In addition, the enclosed Dockerfile can be
built to produce an image capable of running GPUSimilarity.

##Requirements
* nvidia-docker (regular docker will not expose GPUs)
* Newest versions of nvidia drivers
* Recommended for production:  systemd / nginx

## How to use pre-built container
Both releases and master are always built and uploaded publicly already.  You
can access them via:
For master:
`docker pull klorton/gpusimilarity:latest`
For releases:
`docker pull klorton/gpusimilarity:v1.0`

The underlying requirements will change very little, so after the initial
build/pull new docker layers will only require a few megabytes for
gpusimilarity changes.

## How to build (recommended only for development)
From the docker directory:
```
docker build -t gpusim:internal .
```

## How to use to build a fingerprint file
TODO

## How to start a gpusimilarity server interactively
```
/bin/nvidia-docker run --net=host -v /path/to/fsim/files:/mnt/fsim -it \
klorton/gpusimilarity:latest bash -c 'python3 /gpusimilarity/bld/python/gpusim_server.py \
/path/to/fsim/files/1.fsim --port 8080 --http_interface'
```

## Example of a systemd configuration file
```
[Unit]
Description=Service to return GPUSimilarity results
After=docker.service

[Service]
Type=simple
Restart=always
RestartSec=10
User=gpusim
ExecStart=/bin/nvidia-docker run --net=host -v /path/to/fsim:/mnt/fsim klorton/gpusimilarity:latest python3 /gpusimilarity/bld/python/gpusim_server.py --port 8080 /mnt/fsim/your_file.fsim

[Install]
WantedBy=multi-user.target
```

## Additional steps you'll want to take for production
* The `--http_interface` option is only for testing **never** use it in production as it hasn't been hardened at all
* You should use a systemd script to start/stop the service, including restarts
* **Never** expose the http from this service directly, use nginx or equivalent with certificates to expose via ssl/https.
