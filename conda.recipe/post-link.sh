#!/usr/bin/env bash

cd ${SP_DIR}/${PKG_NAME}
wget ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.52-86.0/interproscan-5.52-86.0-64-bit.tar.gz
wget ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.52-86.0/interproscan-5.52-86.0-64-bit.tar.gz.md5
md5sum -c interproscan-5.52-86.0-64-bit.tar.gz.md5
tar --exclude-from=${RECIPE_DIR}/exclude.txt -pxvzf interproscan-5.52-86.0-*-bit.tar.gz
rm interproscan-5.52-86.0-*-bit.tar.gz
ln -s ${PREFIX}/bin ${SP_DIR}/${PKG_NAME}/interproscan-5.52-86.0/interproscan.sh

$PYTHON ${SP_DIR}/${PKG_NAME}/interproscan-5.52-86.0/initial_setup.py 





#emerald_download_data 
#ln -s $SP_DIR/emeraldbgc/interproscan-5.52-86.0/interproscan.sh $PREFIX/bin/interproscan.sh
#ln -s ~/interproscan.sh $PREFIX/bin/interproscan.sh
