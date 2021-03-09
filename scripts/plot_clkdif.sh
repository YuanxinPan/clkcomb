#!/bin/bash

[ $# -ne 2 ] && echo "usage: plot_clkdif.sh dif_file log_file" && exit 1

# Input
dif_file=$1
log_file=$2

WORK_DIR=`pwd`
TEMP_DIR=`mktemp -d`

ACs=(` awk '/satclk/ {print $2}' $dif_file | sort -u`)
PRNs=(`awk '/satclk/ {print $3}' $dif_file | sort -u`)

echo ${ACs[*]} > $TEMP_DIR/ac_list && sed -i 's/ /\n/g' $TEMP_DIR/ac_list

for prn in ${PRNs[*]}
do
    echo $prn...

    for ac in ${ACs[*]}
    do
        grep "$ac $prn"        "$dif_file" > $TEMP_DIR/dif_${ac}_$prn
        grep "del .. $ac $prn" "$log_file" > $TEMP_DIR/del_${ac}_$prn
    done

    # plot
    plot_clkdif.py $TEMP_DIR $prn
done

rm -rf $TEMP_DIR
