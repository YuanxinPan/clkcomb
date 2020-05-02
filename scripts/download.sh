#!/bin/bash

######################################################################
##                        Message Colors                            ##
######################################################################
NC='\033[0m'
RED='\033[0;31m'
CYAN='\033[0;36m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'

MSGERR="${RED}error:$NC"
MSGWAR="${YELLOW}warning:$NC"
MSGINF="${BLUE}::$NC"
MSGSTA="${BLUE}===>$NC"

######################################################################
##                     Funciton definations                         ##
######################################################################
main()
{
    [ $# -ne 3 ] && Print_Help && return 1

    ## Set date
    year=$1 doy=`echo $2 | awk '{printf("%03d\n",$1)}'`
    [ ${#year} -eq 2 ] && year=20$1

    local gweek=(`gweek $doy $year`)
    week=${gweek[0]}; dow=${gweek[1]}
    echo "$year $doy | $week $dow"

    [ $3 -eq 0 ] && DownloadGPSFloat  $week $dow
    [ $3 -eq 1 ] && DownloadGNSSFloat $week $dow
    [ $3 -eq 2 ] && DownloadGPSFixed  $week $dow
}

Print_Help() { # purpose: print usage for download.sh
               # usage  : Print_Help
    echo " -----------------------------------------------------------------------"
    echo "  Purpose  :    download orbit/clock/bias for combination"
    echo "  Usage    :    download.sh year doy opt(0-2)"
    echo "                   -- 0: float GPS  orbit/clock"
    echo "                   -- 1: float GNSS orbit/clock"
    echo "                   -- 2: fixed GPS  orbit/clock/bias"
    echo "  Example  :    download.sh 2020 1 0"
    echo " -----------------------------------------------------------------------"
}

DownloadGPSFloat() { # purpose: download float GPS clock/orbit for combination
                     # usage  : DownloadGPSFloat
	local ac sp3 clk url
    for ac in cod esa gfz jpl emr grg igs
    do
        # echo $ac
        sp3="${ac}${week}${dow}.sp3"
        clk="${ac}${week}${dow}.clk"
        [ "$ac" = 'cod' ] && sp3=${sp3/sp3/eph}
        [ "$ac" = 'igs' ] && clk=${clk}_30s

        url="ftp://cddis.gsfc.nasa.gov/pub/gps/products/${week}/"
        WgetDownload ${url}/${sp3}.Z && uncompress ${sp3}.Z
        WgetDownload ${url}/${clk}.Z && uncompress ${clk}.Z

        [ "$ac" = 'cod' -a -f "$sp3" ] && mv ${sp3} ${sp3/eph/sp3}
        [ "$ac" = 'igs' -a -f "$clk" ] && mv ${clk} ${clk%_30s}
    done
}

DownloadGNSSFloat() { # purpose: download float multi-GNSS clock/orbit for combination
                      # usage  : DownloadGNSSFloat
    local ac sp3 clk url
    for ac in cod gfz grg wum
    do
        # echo $ac
        sp3="${ac^^}0MGXFIN_${year}${doy}0000_01D_05M_ORB.SP3.gz"
        clk="${ac^^}0MGXFIN_${year}${doy}0000_01D_30S_CLK.CLK.gz"
        [ "$ac" = 'gfz' ] && sp3=${sp3/FIN/RAP} && clk=${clk/FIN/RAP}
        [ $ac = wum -o $ac = grg ] && sp3=${sp3/05M/15M}

        url="ftp://igs.gnsswhu.cn/pub/gnss/products/mgex/${week}/"
        WgetDownload ${url}/${sp3} && gunzip -f ${sp3} && mv ${sp3%.gz} ${ac}${week}${dow}.sp3
        WgetDownload ${url}/${clk} && gunzip -f ${clk} && mv ${clk%.gz} ${ac}${week}${dow}.clk
    done
}

DownloadGPSFixed() { # purpose: download GPS IRC bias/clock/orbit for combination
                     # usage  : DownloadGPSFixed
    local url sp3 clk bia
    # grg
    url="ftp://cddis.gsfc.nasa.gov/pub/gps/products/${week}/"
    sp3="grg${week}${dow}.sp3"
    clk="grg${week}${dow}.clk"
    WgetDownload ${url}/${sp3}.Z && uncompress ${sp3}.Z
    WgetDownload ${url}/${clk}.Z && uncompress ${clk}.Z

    # COD
    url="ftp://ftp.aiub.unibe.ch/CODE/${year}/"
    sp3="COD${week}${dow}.EPH"
    clk="COD${week}${dow}.CLK"
    bia="COD${week}${dow}.BIA"
    WgetDownload ${url}/${sp3}.Z && uncompress ${sp3}.Z && cp $sp3 whu${week}${dow}.sp3 && mv $sp3 cod${week}${dow}.sp3
    WgetDownload ${url}/${clk}.Z && uncompress ${clk}.Z && mv $clk cod${week}${dow}.clk
    WgetDownload ${url}/${bia}.Z && uncompress ${bia}.Z && mv $bia cod${week}${dow}.bia

    # whu
    url="ftp://igs.gnsswhu.cn/pub/whu/phasebias/${year}/"
    clk="WHU5IGSFIN_${year}${doy}0000_01D_30S_CLK.CLK"
    bia="WHU0IGSFIN_${year}${doy}0000_01D_01D_ABS.BIA"
    WgetDownload ${url}/clock/${clk}.Z && uncompress ${clk}.Z && mv $clk whu${week}${dow}.clk
    WgetDownload ${url}/bias/${bia}.Z  && uncompress ${bia}.Z && mv $bia whu${week}${dow}.bia
}

WgetDownload() { # purpose: download a file with wget
                 # usage  : WgetDownload url
    local url="$1"
    local args="-nv -nc -c -t 3 --connect-timeout=10 --read-timeout=60"
    cmd="wget ${args} ${url}"
    $cmd
    [ -e $(basename "${url}") ] && return 0 || return 1
}

######################################################################
##                               Entry                              ##
######################################################################
main "$@"

