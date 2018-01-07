#! /bin/bash

mkdir tmp
cd tmp

wget ftp://ftp1.nci.nih.gov/pub/CGAP/mitelman.tar.gz

tar zxvf mitelman.tar.gz

mv molclingene.dat ../

cd ../

rm -rf tmp


