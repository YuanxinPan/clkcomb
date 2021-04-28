#!/bin/bash

#####################################################
##                                                 ##
##  Purpose: Download precise products for pppx    ##
##                                                 ##
##  Author : Yuanxin Pan    yxpan@whu.edu.cn       ##
##                                                 ##
##  Version: 1.0                                   ##
##                                                 ##
##  Date   : Nov-06, 2020                          ##
##                                                 ##
##    Copyright (C) 2020 by Yuanxin Pan            ##
##                                                 ##
#####################################################

# Message prefix
MSGERR="\033[0;31merror:\033[0m"
MSGWAR="\033[1;33mwarning:\033[0m"
MSGINF="\033[1;34m::\033[0m"
MSGSTA="\033[1;34m===>\033[0m"

# HOST
HOST="https://cddis.nasa.gov/archive"
HOST="ftps://gdc.cddis.eosdis.nasa.gov"
WUM_HOST="ftp://igs.ign.fr/pub/igs"
WUM_HOST="ftp://igs.gnsswhu.cn/pub/gnss"

# wget
WGET_ARGS="--auth-no-challenge"                                     # https
WGET_ARGS="--ftp-user anonymous --ftp-password yxpan@whu.edu.cn"    # ftps

######################################################################
##                     Funciton definations                         ##
######################################################################
main()
{
    [ $# -ne 3 ] && Help && return 1

    ## Set date
    year=$1 doy=`echo $2 | awk '{printf("%03d\n",$1)}'`
    [ ${#year} -eq 2 ] && year=20$1

    local gweek=(`gweek $doy $year`)
    week=${gweek[0]}; dow=${gweek[1]}
    echo "$year $doy | $week $dow"

    [ $3 -eq 0 ] && DownloadGPSFloat  $week $dow
    [ $3 -eq 1 ] && DownloadGNSSFloat $week $dow
    [ $3 -eq 2 ] && DownloadGPSFixed  $week $dow
    [ $3 -eq 3 ] && DownloadGNSSFixed $week $dow
}

Help() { # purpose: print usage for download.sh
               # usage  : Help
    echo " ------------------------------------------------------------"
    echo "  Purpose  :    download orbit/clock/bias for combination"
    echo "  Usage    :    download.sh year doy opt(0-2)"
    echo "                   -- 0: float GPS  orbit/clock"
    echo "                   -- 1: float GNSS orbit/clock"
    echo "                   -- 2: fixed GPS  orbit/clock/bias"
    echo "                   -- 3: fixed GNSS orbit/clock/bias"
    echo "  Example  :    download.sh 2020 1 0"
    echo " ------------------------------------------------------------"
}

DownloadGPSFloat() { # purpose: download float GPS clock/orbit for combination
                     # usage  : DownloadGPSFloat
	local ac sp3 clk url
    for ac in cod emr esa gfz grg jpl igs
    do
        # echo $ac
        sp3="${ac}${week}${dow}.sp3"
        clk="${ac}${week}${dow}.clk"
        [ "$ac" = 'cod' ] && sp3=${sp3/sp3/eph}
        [ "$ac" = 'igs' ] && clk=${clk}_30s

        url="$HOST/gps/products/${week}"
        WgetDownload $url/${sp3}.Z "$WGET_ARGS" && uncompress ${sp3}.Z
        WgetDownload $url/${clk}.Z "$WGET_ARGS" && uncompress ${clk}.Z

        [ "$ac" = 'cod' -a -f "$sp3" ] && mv ${sp3} ${sp3/eph/sp3}
        [ "$ac" = 'igs' -a -f "$clk" ] && mv ${clk} ${clk%_30s}
    done
}

DownloadGNSSFloat() { # purpose: download float multi-GNSS clock/orbit for combination
                      # usage  : DownloadGNSSFloat
    local ac sp3 clk url
    for ac in cod esa gfz grg wum
    do
        # echo $ac
        sp3="${ac^^}0MGXFIN_${year}${doy}0000_01D_05M_ORB.SP3.gz"
        clk="${ac^^}0MGXFIN_${year}${doy}0000_01D_30S_CLK.CLK.gz"
        [ "$ac" = gfz ] && sp3=${sp3/FIN/RAP} && clk=${clk/FIN/RAP}
        [ $ac = esa ] && sp3=${sp3/MGX/MGN} && clk=${clk/MGX/MGN}
        [ $ac = wum -o $ac = grg ] && sp3=${sp3/05M/15M}

        url="$HOST/gps/products/mgex/${week}"
        [ $ac = 'esa' ] && url="http://navigation-office.esa.int/products/gnss-products/${week}"
        [ $ac = 'wum' ] && url="$WUM_HOST/products/mgex/${week}"

        WgetDownload $url/${sp3} "$WGET_ARGS" && gunzip -f ${sp3} && mv ${sp3%.gz} ${ac}${week}${dow}.sp3
        WgetDownload $url/${clk} "$WGET_ARGS" && gunzip -f ${clk} && mv ${clk%.gz} ${ac}${week}${dow}.clk
    done
}

DownloadGPSFixed() { # purpose: download GPS IRC bias/clock/orbit for combination
                     # usage  : DownloadGPSFixed
    local url sp3 clk bia
    # grg
    url="$HOST/gps/products/${week}"
    sp3="grg${week}${dow}.sp3"
    clk="grg${week}${dow}.clk"
    WgetDownload $url/${sp3}.Z "$WGET_ARGS" && uncompress ${sp3}.Z
    WgetDownload $url/${clk}.Z "$WGET_ARGS" && uncompress ${clk}.Z

    # COD
    url="ftp://ftp.aiub.unibe.ch/CODE/${year}"
    sp3="COD${week}${dow}.EPH"
    clk="COD${week}${dow}.CLK"
    bia="COD${week}${dow}.BIA"
    WgetDownload ${url}/${sp3}.Z && uncompress ${sp3}.Z && cp $sp3 whu${week}${dow}.sp3 && mv $sp3 cod${week}${dow}.sp3
    WgetDownload ${url}/${clk}.Z && uncompress ${clk}.Z && mv $clk cod${week}${dow}.clk
    WgetDownload ${url}/${bia}.Z && uncompress ${bia}.Z && mv $bia cod${week}${dow}.bia

    # whu
    url="ftp://igs.gnsswhu.cn/pub/whu/phasebias/${year}"
    clk="WHU5IGSFIN_${year}${doy}0000_01D_30S_CLK.CLK"
    bia="WHU0IGSFIN_${year}${doy}0000_01D_01D_ABS.BIA"
    WgetDownload ${url}/clock/${clk}.Z && uncompress ${clk}.Z && mv $clk whu${week}${dow}.clk
    WgetDownload ${url}/bias/${bia}.Z  && uncompress ${bia}.Z && mv $bia whu${week}${dow}.bia
}

DownloadGNSSFixed() { # purpose: download GNSS IRC bias/clock/orbit for combination
                      # usage  : DownloadGNSSFixed
    local url sp3 clk bia
    # grg
    url="$HOST/gps/products/mgex/${week}"
    sp3="GRG0MGXFIN_${year}${doy}0000_01D_15M_ORB.SP3.gz"
    clk="GRG0MGXFIN_${year}${doy}0000_01D_30S_CLK.CLK.gz"
    WgetDownload $url/${sp3} "$WGET_ARGS" && gunzip -f ${sp3} && mv ${sp3%.gz} grg${week}${dow}.sp3
    WgetDownload $url/${clk} "$WGET_ARGS" && gunzip -f ${clk} && mv ${clk%.gz} grg${week}${dow}.clk

    # COD
    url="ftp://ftp.aiub.unibe.ch/CODE_MGEX/CODE/${year}"
    sp3="COM${week}${dow}.EPH.Z"
    clk="COM${week}${dow}.CLK.Z"
    bia="COM${week}${dow}.BIA.Z"
    WgetDownload ${url}/${sp3} && uncompress ${sp3} && mv ${sp3%.Z} cod${week}${dow}.sp3
    WgetDownload ${url}/${clk} && uncompress ${clk} && mv ${clk%.Z} cod${week}${dow}.clk
    WgetDownload ${url}/${bia} && uncompress ${bia} && mv ${bia%.Z} cod${week}${dow}.bia

    # whu
    url="ftp://igs.gnsswhu.cn/pub/whu/phasebias/${year}"
    clk="WHU5MGXFIN_${year}${doy}0000_01D_30S_CLK.CLK.Z"
    bia="WHU0MGXFIN_${year}${doy}0000_01D_01D_ABS.BIA.Z"
    WgetDownload ${url}/clock/${clk} && uncompress ${clk} && mv ${clk%.Z} whu${week}${dow}.clk
    WgetDownload ${url}/bias/${bia}  && uncompress ${bia} && mv ${bia%.Z} whu${week}${dow}.bia

    url="$WUM_HOST/products/mgex/${week}"
    sp3="WUM0MGXFIN_${year}${doy}0000_01D_15M_ORB.SP3.gz"
    WgetDownload ${url}/${sp3} && gunzip ${sp3} && mv ${sp3%.gz} whu${week}${dow}.sp3
}

WgetDownload() { # purpose: download a file with wget
                 # usage  : WgetDownload url
    local url="$1"
    # local args="$2 -nv -N -t 3 --connect-timeout=10 --read-timeout=60"
    local args="$2 -nv -nc -c -t 3 --connect-timeout=10 --read-timeout=60"
    cmd="wget ${args} ${url}"
    $cmd
    [ -e $(basename "${url}") ] && return 0 || return 1
}

######################################################################
##                               Entry                              ##
######################################################################
main "$@"

