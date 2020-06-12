#!/bin/bash

# Adapted from: https://denibertovic.com/posts/handling-permissions-with-docker-volumes/
# Add local user and group
# Either use the LOCAL_USER_ID and LOCAL_GROUP_ID if passed in at runtime or
# fallback

USER_ID=${LOCAL_USER_ID:-0}
GROUP_ID=${LOCAL_GROUP_ID:-0}
OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}
export OMP_NUM_THREADS
dname="/host/${LOCAL_DIR}"

if [[ ! -d "$dname" ]]; then
	echo "ERROR: directory $dname doesn't exist"
	echo "Make sure the variable LOCAL_DIR is set."
	exit 1
fi

if ! cd "$dname"; then
	echo "ERROR: could not cd to $dname."
	echo "Make sure the variable LOCAL_DIR is set."
	exit 2
fi

echo "Starting bgw-docker container with UID / GID = $USER_ID / $GROUP_ID at $dname"
groupadd -g $GROUP_ID -o bgw
useradd --shell /bin/bash -u $USER_ID -g $GROUP_ID -o -c "" -m bgw
echo 'ulimit -s unlimited' >> $HOME/.bashrc
ulimit -s unlimited

mpiexec xstelloptv2 -h