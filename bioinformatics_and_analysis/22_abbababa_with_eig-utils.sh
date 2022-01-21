#!/bin/sh

source ./env.sh

pfx=island_thrush_all_filter

~/eig-utils/abba-baba \
	-p fourpop_shorter_names.txt \
	${pfx}.ind ${pfx}.geno ${pfx}.snp \
	> ${pfx}.abba-baba.txt
