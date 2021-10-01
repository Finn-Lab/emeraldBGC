#!/usr/bin/env bash

set -e

rsync -avh --progress rsync://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.52-86.0/interproscan-5.52-86.0-64-bit.tar.gz .

rsync -avh --progress rsync://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.52-86.0/interproscan-5.52-86.0-64-bit.tar.gz.md5 .
