#!/usr/bin/env bash

env > ~/env.test.txt
wget ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.52-86.0/interproscan-5.52-86.0-64-bit.tar.gz
wget ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.52-86.0/interproscan-5.52-86.0-64-bit.tar.gz.md5
md5sum -c interproscan-5.52-86.0-64-bit.tar.gz.md5
PY=$CONDA_PREFIX/bin/python
PKG_FILE=$($PY -c "import ${PKG_NAME};print(${PKG_NAME}.__file__)")
PKG_DIR=$(dirname $PKG_FILE)
echo $PKG_DIR >> ~/env.test.txt
tar --exclude-from=${PKG_DIR}/exclude.txt -pxvzf interproscan-5.52-86.0-*-bit.tar.gz -C ${PKG_DIR}
rm interproscan-5.52-86.0-*-bit.tar.gz
ln -s ${PKG_DIR}/interproscan-5.52-86.0/interproscan.sh ${CONDA_PREFIX}/bin

cd ${PKG_DIR}/interproscan-5.52-86.0/
$PY initial_setup.py 





#emerald_download_data 
#ln -s $SP_DIR/emeraldbgc/interproscan-5.52-86.0/interproscan.sh $PREFIX/bin/interproscan.sh
#ln -s ~/interproscan.sh $PREFIX/bin/interproscan.sh
