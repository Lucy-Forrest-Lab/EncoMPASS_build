#!/bin/bash

# Take source and destination, normalize any host:path forms so the host is always <computer name>, then run rsync over SSH with a standard set of flags and permissive rwX perms for everyone.
# You’d want to replace <computer name> with your real remote host

PATH_FROM=${1}
PATH_TO=${2}

IFS=":" read -a ARR_FROM <<< ${PATH_FROM}
if [[ "${#ARR_FROM[@]}" == "2" ]]
then
    PATH_FROM="<computer name>:${ARR_FROM[1]}"
fi

IFS=":" read -a ARR_TO <<< ${PATH_TO}
if [[ "${#ARR_TO[@]}" == "2" ]]
then
    PATH_TO="<computer name>:${ARR_TO[1]}"
fi

rsync -auzP --no-perms --no-g --chmod=ugo=rwX --rsh=ssh ${PATH_FROM} ${PATH_TO}
