#!/bin/bash

modulename=new-module

echo 'Extracting...'
mkdir -p $modulename
tar -xzf template.tar.gz -C $modulename
#cp ~/Downloads/JAGS-WIENER-MODULE-1.1/* $modulename -r

echo 'Editing...'
cp src-ct.cc                     $modulename/src/ct.cc
cp src-distributions-HEADER.h    $modulename/src/distributions/DCT.h
cp src-distributions-SOURCE.cc   $modulename/src/distributions/DCT.cc
cp src-distributions-Makefile.am $modulename/src/distributions/Makefile.am
cp src-Makefile.am               $modulename/src/Makefile.am
cp Makefile.am                   $modulename/Makefile.am
cp configure.ac                  $modulename/configure.ac

cd $modulename

echo 'Compiling...'
autoreconf -fvi && ./configure --prefix=/usr

make && sudo make install
