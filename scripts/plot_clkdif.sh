#!/bin/bash

# clkcomb - Clock and phase bias products Combination
# Copyright (C) 2021 Yuanxin Pan
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice,
#    this list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
# 3. Neither the name of the copyright holder nor the names of its contributors
#    may be used to endorse or promote products derived from this software
#    without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.

[ $# -ne 2 ] && echo "usage: plot_clkdif.sh dif_file log_file" && exit 1

# Check python installation
python -V > /dev/null 2>&1 || { echo "python3 (with matplotlib) needs to be installed"; exit 1; }

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
        grep "$ac $prn"         "$dif_file" > $TEMP_DIR/dif_${ac}_$prn
        grep "del .. $ac  $prn" "$log_file" > $TEMP_DIR/del_${ac}_$prn
    done

    # plot
    plot_clkdif.py $TEMP_DIR $prn || break
done

rm -rf $TEMP_DIR
