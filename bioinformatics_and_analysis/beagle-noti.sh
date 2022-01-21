#!/bin/sh

ifile=island_thrush_angsd.beagle.gz
ofile=island_thrush_angsd.beagle.noti.gz

zcat $ifile \
 | awk '$2$3!="13" && $2$3!="31" && $2$3!="20" && $2$3!="02"' \
 | gzip -c - \
 > $ofile
