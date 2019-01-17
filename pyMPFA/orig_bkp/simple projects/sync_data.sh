#!/bin/bash
declare -A outputmap
declare -A inputmap

outputroot='/media/anton/Elements/backup/'
inputroot='/home/anton/'
excludefile=$inputroot'data/simple\ projects/exclude_files.txt'
outputmap=([0]=$outputroot'HDD.Backup/MCBDataCloud/'
           [1]=$outputroot'HDD.Doc/'
           [2]=$outputroot'SSD.Data/'
           [3]=$outputroot'HDD.Home/')
inputmap=([0]=$inputroot'backup/MCBDataCloud/'
          [1]=$inputroot'Doc/'
          [2]=$inputroot'data/'
          [3]=$inputroot
          )
indexlist=${!inputmap[*]}
for indexitem in $indexlist; do
    options='-aruh --progress'
    inputdir=${inputmap[$indexitem]}
    outputdir=${outputmap[$indexitem]}
    if [[ $inputdir = $inputroot ]]; then
        options=$options '--exclude-from='$excludefile
    fi
    printf "Sync folder" "$inputdir" "to" "$outputdir"
    rsync $options $inputdir $outputdir
done
