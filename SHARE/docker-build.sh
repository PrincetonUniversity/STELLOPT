#!/bin/sh

# Run from STELLOPT root directory to build docker container
ln -sf SHARE/Dockerfile Dockerfile
docker build --tag=stellopt .
