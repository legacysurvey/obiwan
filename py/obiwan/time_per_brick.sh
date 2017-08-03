#!/bin/bash

# Does string occur in file?
# 0 for true, 1 for false
inFile() {
	local fn=$1
	local text=$2
	if grep -e "$text" --quiet $fn; then
		return 0
	else
		return 1
	fi
}

# Wall time for given legacypipe tractor for a log file
getWAll() {
	local log=$1
	local stage=$2
	# If text doesn't exist then legacypipe error occured
	if inFile $log "Resources for stage $stage"; then
		local a=`grep -e "Resources for stage $stage" -A 13 $log|grep "Grand total Wall:"|awk '{print $4,$5}'`
	else
		local a=NotExist
	fi
	# return
	echo $a
}

## Main ##
export logs=$1

for stage in tims mask_junk srcs fitblobs coadds wise_forced writecat; do
	export outfn=wall_$stage.txt
	for log in `cat $logs`;do 
		wall=$(getWAll $log $stage)
		brick=`basename $log|sed s#.log##g`
		echo $brick wall $wall >> $outfn
	done
done


