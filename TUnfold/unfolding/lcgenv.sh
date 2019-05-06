# setup new env
# Using gcc with newer binutils
# Using globaly provided gcc from /cvmfs/sft.cern.ch/lcg/contrib
source /cvmfs/sft.cern.ch/lcg/contrib/gcc/6.2.0binutils/x86_64-slc6/setup.sh;
export GCC__HOME="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/gcc/6.2.0/x86_64-slc6";
# Using gcc with newer binutils
# Using globaly provided gcc from /cvmfs/sft.cern.ch/lcg/contrib
# found package Python-2cba0
# found package sqlite-9c81e
source /cvmfs/sft.cern.ch/lcg/contrib/gcc/6.2.0binutils/x86_64-slc6/setup.sh;
export LD_LIBRARY_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/Boost/1.66.0/x86_64-slc6-gcc62-opt/lib:$LD_LIBRARY_PATH";
export PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/Python/2.7.15/x86_64-slc6-gcc62-opt/bin:$PATH";
export LD_LIBRARY_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/Python/2.7.15/x86_64-slc6-gcc62-opt/lib:$LD_LIBRARY_PATH";
export PKG_CONFIG_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/Python/2.7.15/x86_64-slc6-gcc62-opt/lib/pkgconfig:$PKG_CONFIG_PATH";
export PYTHONPATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/Python/2.7.15/x86_64-slc6-gcc62-opt/lib/python2.7/site-packages:$PYTHONPATH";
export PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/sqlite/3210000/x86_64-slc6-gcc62-opt/bin:$PATH";
export LD_LIBRARY_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/sqlite/3210000/x86_64-slc6-gcc62-opt/lib:$LD_LIBRARY_PATH";
export PKG_CONFIG_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/sqlite/3210000/x86_64-slc6-gcc62-opt/lib/pkgconfig:$PKG_CONFIG_PATH";
export SQLITE__HOME="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/sqlite/3210000/x86_64-slc6-gcc62-opt";
export PYTHON__HOME="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/Python/2.7.15/x86_64-slc6-gcc62-opt";
cd "/cvmfs/sft.cern.ch/lcg/releases/LCG_94/Python/2.7.15/x86_64-slc6-gcc62-opt"
export PYTHONHOME="${PYTHON__HOME}"
cd - 1>/dev/null # from /cvmfs/sft.cern.ch/lcg/releases/LCG_94/Python/2.7.15/x86_64-slc6-gcc62-opt
export BOOST__HOME="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/Boost/1.66.0/x86_64-slc6-gcc62-opt";
cd "/cvmfs/sft.cern.ch/lcg/releases/LCG_94/Boost/1.66.0/x86_64-slc6-gcc62-opt"
export CPLUS_INCLUDE_PATH=${BOOST__HOME}/include:$CPLUS_INCLUDE_PATH
export C_INCLUDE_PATH=${BOOST__HOME}/include:$C_INCLUDE_PATH
cd - 1>/dev/null # from /cvmfs/sft.cern.ch/lcg/releases/LCG_94/Boost/1.66.0/x86_64-slc6-gcc62-opt
# Using gcc with newer binutils
# Using globaly provided gcc from /cvmfs/sft.cern.ch/lcg/contrib
# found package libxml2-830a9
# found package graphviz-dc59a
# found package pango-1771b
# found package harfbuzz-c8713
# found package glib-53b12
# found package gettext-f094d
# found package pcre-9b589
# found package pkg_config-588e0
# found package libffi-26487
# found package freetype-4f115
# found package zlib-da225
# found package cairo-1af41
# found package fontconfig-63041
# found package expat-fe6d4
# found package pkg_config-588e0
# found package gperf-699d7
# found package freetype-4f115
# found package png-1a97c
# found package freetype-4f115
# found package pixman-bdd51
# found package pkg_config-588e0
# found package fontconfig-63041
# found package cairo-1af41
# found package freetype-4f115
# found package glib-53b12
# found package libtool-9ad34
# found package fontconfig-63041
# found package expat-fe6d4
# found package blas-bff5d
# found package numpy-d30c2
# found package setuptools-c7ed7
# found package Python-2cba0
# found package sqlite-9c81e
# found package Python-2cba0
# found package Davix-56a1c
# found package Boost-4fdfa
# found package Python-2cba0
# found package libxml2-830a9
# found package fftw-a8420
# found package tbb-d3621
# found package CASTOR-f6b00
# found package mysql-d1585
# found package gfal-6fc75
# found package xrootd-c201f
# found package Python-2cba0
# found package R-48a34
# found package cairo-1af41
# found package zeromq-183c3
# found package zlib-da225
# found package srm_ifce-be254
# found package Python-2cba0
# found package dcap-cdd28
# found package Vc-7fbe0
# found package oracle-e33b7
# found package GSL-32fc5
source /cvmfs/sft.cern.ch/lcg/contrib/gcc/6.2.0binutils/x86_64-slc6/setup.sh;
export PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/ROOT/6.14.04/x86_64-slc6-gcc62-opt/bin:$PATH";
export LD_LIBRARY_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/ROOT/6.14.04/x86_64-slc6-gcc62-opt/lib:$LD_LIBRARY_PATH";
export LD_LIBRARY_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/ROOT/6.14.04/x86_64-slc6-gcc62-opt/lib/JupyROOT:$LD_LIBRARY_PATH";
export LD_LIBRARY_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/ROOT/6.14.04/x86_64-slc6-gcc62-opt/lib/JsMVA:$LD_LIBRARY_PATH";
export PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/libxml2/2.9.7/x86_64-slc6-gcc62-opt/bin:$PATH";
export LD_LIBRARY_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/libxml2/2.9.7/x86_64-slc6-gcc62-opt/lib:$LD_LIBRARY_PATH";
export PKG_CONFIG_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/libxml2/2.9.7/x86_64-slc6-gcc62-opt/lib/pkgconfig:$PKG_CONFIG_PATH";
export LD_LIBRARY_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/libxml2/2.9.7/x86_64-slc6-gcc62-opt/lib/cmake:$LD_LIBRARY_PATH";
export LIBXML2__HOME="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/libxml2/2.9.7/x86_64-slc6-gcc62-opt";
export PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/graphviz/2.28.0/x86_64-slc6-gcc62-opt/bin:$PATH";
export LD_LIBRARY_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/graphviz/2.28.0/x86_64-slc6-gcc62-opt/lib:$LD_LIBRARY_PATH";
export PKG_CONFIG_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/graphviz/2.28.0/x86_64-slc6-gcc62-opt/lib/pkgconfig:$PKG_CONFIG_PATH";
export LD_LIBRARY_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/graphviz/2.28.0/x86_64-slc6-gcc62-opt/lib/graphviz:$LD_LIBRARY_PATH";
export PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/pango/1.40.13/x86_64-slc6-gcc62-opt/bin:$PATH";
export LD_LIBRARY_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/pango/1.40.13/x86_64-slc6-gcc62-opt/lib:$LD_LIBRARY_PATH";
export PKG_CONFIG_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/pango/1.40.13/x86_64-slc6-gcc62-opt/lib/pkgconfig:$PKG_CONFIG_PATH";
export PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/harfbuzz/1.6.3/x86_64-slc6-gcc62-opt/bin:$PATH";
export LD_LIBRARY_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/harfbuzz/1.6.3/x86_64-slc6-gcc62-opt/lib:$LD_LIBRARY_PATH";
export PKG_CONFIG_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/harfbuzz/1.6.3/x86_64-slc6-gcc62-opt/lib/pkgconfig:$PKG_CONFIG_PATH";
export PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/glib/2.52.2/x86_64-slc6-gcc62-opt/bin:$PATH";
export LD_LIBRARY_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/glib/2.52.2/x86_64-slc6-gcc62-opt/lib:$LD_LIBRARY_PATH";
export PKG_CONFIG_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/glib/2.52.2/x86_64-slc6-gcc62-opt/lib/pkgconfig:$PKG_CONFIG_PATH";
export LD_LIBRARY_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/glib/2.52.2/x86_64-slc6-gcc62-opt/lib/gio:$LD_LIBRARY_PATH";
export LD_LIBRARY_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/glib/2.52.2/x86_64-slc6-gcc62-opt/lib/glib-2.0:$LD_LIBRARY_PATH";
export PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/gettext/0.19.8.1/x86_64-slc6-gcc62-opt/bin:$PATH";
export LD_LIBRARY_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/gettext/0.19.8.1/x86_64-slc6-gcc62-opt/lib:$LD_LIBRARY_PATH";
export LD_LIBRARY_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/gettext/0.19.8.1/x86_64-slc6-gcc62-opt/lib/gettext:$LD_LIBRARY_PATH";
export GETTEXT__HOME="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/gettext/0.19.8.1/x86_64-slc6-gcc62-opt";
export PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/pcre/8.38/x86_64-slc6-gcc62-opt/bin:$PATH";
export LD_LIBRARY_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/pcre/8.38/x86_64-slc6-gcc62-opt/lib:$LD_LIBRARY_PATH";
export PKG_CONFIG_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/pcre/8.38/x86_64-slc6-gcc62-opt/lib/pkgconfig:$PKG_CONFIG_PATH";
export PCRE__HOME="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/pcre/8.38/x86_64-slc6-gcc62-opt";
export PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/pkg_config/0.28/x86_64-slc6-gcc62-opt/bin:$PATH";
export PKG_CONFIG__HOME="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/pkg_config/0.28/x86_64-slc6-gcc62-opt";
export LD_LIBRARY_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/libffi/3.2.1/x86_64-slc6-gcc62-opt/lib:$LD_LIBRARY_PATH";
export LD_LIBRARY_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/libffi/3.2.1/x86_64-slc6-gcc62-opt/lib/libffi-3.2.1:$LD_LIBRARY_PATH";
export PKG_CONFIG_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/libffi/3.2.1/x86_64-slc6-gcc62-opt/lib/pkgconfig:$PKG_CONFIG_PATH";
export LD_LIBRARY_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/libffi/3.2.1/x86_64-slc6-gcc62-opt/lib64:$LD_LIBRARY_PATH";
export LIBFFI__HOME="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/libffi/3.2.1/x86_64-slc6-gcc62-opt";
export GLIB__HOME="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/glib/2.52.2/x86_64-slc6-gcc62-opt";
export PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/freetype/2.6.3/x86_64-slc6-gcc62-opt/bin:$PATH";
export LD_LIBRARY_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/freetype/2.6.3/x86_64-slc6-gcc62-opt/lib:$LD_LIBRARY_PATH";
export PKG_CONFIG_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/freetype/2.6.3/x86_64-slc6-gcc62-opt/lib/pkgconfig:$PKG_CONFIG_PATH";
export LD_LIBRARY_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/zlib/1.2.11/x86_64-slc6-gcc62-opt/lib:$LD_LIBRARY_PATH";
export PKG_CONFIG_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/zlib/1.2.11/x86_64-slc6-gcc62-opt/lib/pkgconfig:$PKG_CONFIG_PATH";
export ZLIB__HOME="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/zlib/1.2.11/x86_64-slc6-gcc62-opt";
export FREETYPE__HOME="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/freetype/2.6.3/x86_64-slc6-gcc62-opt";
export PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/cairo/1.15.8/x86_64-slc6-gcc62-opt/bin:$PATH";
export LD_LIBRARY_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/cairo/1.15.8/x86_64-slc6-gcc62-opt/lib:$LD_LIBRARY_PATH";
export PKG_CONFIG_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/cairo/1.15.8/x86_64-slc6-gcc62-opt/lib/pkgconfig:$PKG_CONFIG_PATH";
export LD_LIBRARY_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/cairo/1.15.8/x86_64-slc6-gcc62-opt/lib/cairo:$LD_LIBRARY_PATH";
export PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/fontconfig/2.12.6/x86_64-slc6-gcc62-opt/bin:$PATH";
export LD_LIBRARY_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/fontconfig/2.12.6/x86_64-slc6-gcc62-opt/lib:$LD_LIBRARY_PATH";
export PKG_CONFIG_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/fontconfig/2.12.6/x86_64-slc6-gcc62-opt/lib/pkgconfig:$PKG_CONFIG_PATH";
export PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/expat/2.2.5/x86_64-slc6-gcc62-opt/bin:$PATH";
export LD_LIBRARY_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/expat/2.2.5/x86_64-slc6-gcc62-opt/lib:$LD_LIBRARY_PATH";
export PKG_CONFIG_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/expat/2.2.5/x86_64-slc6-gcc62-opt/lib/pkgconfig:$PKG_CONFIG_PATH";
export EXPAT__HOME="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/expat/2.2.5/x86_64-slc6-gcc62-opt";
export PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/gperf/3.1/x86_64-slc6-gcc62-opt/bin:$PATH";
export GPERF__HOME="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/gperf/3.1/x86_64-slc6-gcc62-opt";
export FONTCONFIG__HOME="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/fontconfig/2.12.6/x86_64-slc6-gcc62-opt";
export PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/png/1.6.17/x86_64-slc6-gcc62-opt/bin:$PATH";
export LD_LIBRARY_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/png/1.6.17/x86_64-slc6-gcc62-opt/lib:$LD_LIBRARY_PATH";
export PKG_CONFIG_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/png/1.6.17/x86_64-slc6-gcc62-opt/lib/pkgconfig:$PKG_CONFIG_PATH";
export PNG__HOME="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/png/1.6.17/x86_64-slc6-gcc62-opt";
export LD_LIBRARY_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/pixman/0.34.0/x86_64-slc6-gcc62-opt/lib:$LD_LIBRARY_PATH";
export PKG_CONFIG_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/pixman/0.34.0/x86_64-slc6-gcc62-opt/lib/pkgconfig:$PKG_CONFIG_PATH";
export PIXMAN__HOME="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/pixman/0.34.0/x86_64-slc6-gcc62-opt";
export CAIRO__HOME="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/cairo/1.15.8/x86_64-slc6-gcc62-opt";
export HARFBUZZ__HOME="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/harfbuzz/1.6.3/x86_64-slc6-gcc62-opt";
export PANGO__HOME="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/pango/1.40.13/x86_64-slc6-gcc62-opt";
export PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/libtool/2.4.2/x86_64-slc6-gcc62-opt/bin:$PATH";
export LD_LIBRARY_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/libtool/2.4.2/x86_64-slc6-gcc62-opt/lib:$LD_LIBRARY_PATH";
export LIBTOOL__HOME="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/libtool/2.4.2/x86_64-slc6-gcc62-opt";
export GRAPHVIZ__HOME="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/graphviz/2.28.0/x86_64-slc6-gcc62-opt";
export LD_LIBRARY_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/blas/0.2.20.openblas/x86_64-slc6-gcc62-opt/lib64:$LD_LIBRARY_PATH";
export PKG_CONFIG_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/blas/0.2.20.openblas/x86_64-slc6-gcc62-opt/lib64/pkgconfig:$PKG_CONFIG_PATH";
export BLAS__HOME="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/blas/0.2.20.openblas/x86_64-slc6-gcc62-opt";
export PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/numpy/1.14.2/x86_64-slc6-gcc62-opt/bin:$PATH";
export LD_LIBRARY_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/numpy/1.14.2/x86_64-slc6-gcc62-opt/lib:$LD_LIBRARY_PATH";
export PYTHONPATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/numpy/1.14.2/x86_64-slc6-gcc62-opt/lib/python2.7/site-packages:$PYTHONPATH";
export PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/setuptools/36.0.1/x86_64-slc6-gcc62-opt/bin:$PATH";
export LD_LIBRARY_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/setuptools/36.0.1/x86_64-slc6-gcc62-opt/lib:$LD_LIBRARY_PATH";
export PYTHONPATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/setuptools/36.0.1/x86_64-slc6-gcc62-opt/lib/python2.7/site-packages:$PYTHONPATH";
export PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/Python/2.7.15/x86_64-slc6-gcc62-opt/bin:$PATH";
export LD_LIBRARY_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/Python/2.7.15/x86_64-slc6-gcc62-opt/lib:$LD_LIBRARY_PATH";
export PKG_CONFIG_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/Python/2.7.15/x86_64-slc6-gcc62-opt/lib/pkgconfig:$PKG_CONFIG_PATH";
export PYTHONPATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/Python/2.7.15/x86_64-slc6-gcc62-opt/lib/python2.7/site-packages:$PYTHONPATH";
export PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/sqlite/3210000/x86_64-slc6-gcc62-opt/bin:$PATH";
export LD_LIBRARY_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/sqlite/3210000/x86_64-slc6-gcc62-opt/lib:$LD_LIBRARY_PATH";
export PKG_CONFIG_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/sqlite/3210000/x86_64-slc6-gcc62-opt/lib/pkgconfig:$PKG_CONFIG_PATH";
export SQLITE__HOME="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/sqlite/3210000/x86_64-slc6-gcc62-opt";
export PYTHON__HOME="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/Python/2.7.15/x86_64-slc6-gcc62-opt";
cd "/cvmfs/sft.cern.ch/lcg/releases/LCG_94/Python/2.7.15/x86_64-slc6-gcc62-opt"
export PYTHONHOME="${PYTHON__HOME}"
cd - 1>/dev/null # from /cvmfs/sft.cern.ch/lcg/releases/LCG_94/Python/2.7.15/x86_64-slc6-gcc62-opt
export SETUPTOOLS__HOME="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/setuptools/36.0.1/x86_64-slc6-gcc62-opt";
export NUMPY__HOME="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/numpy/1.14.2/x86_64-slc6-gcc62-opt";
export PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/Davix/0.6.7/x86_64-slc6-gcc62-opt/bin:$PATH";
export LD_LIBRARY_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/Davix/0.6.7/x86_64-slc6-gcc62-opt/lib64:$LD_LIBRARY_PATH";
export PKG_CONFIG_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/Davix/0.6.7/x86_64-slc6-gcc62-opt/lib64/pkgconfig:$PKG_CONFIG_PATH";
export LD_LIBRARY_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/Boost/1.66.0/x86_64-slc6-gcc62-opt/lib:$LD_LIBRARY_PATH";
export BOOST__HOME="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/Boost/1.66.0/x86_64-slc6-gcc62-opt";
cd "/cvmfs/sft.cern.ch/lcg/releases/LCG_94/Boost/1.66.0/x86_64-slc6-gcc62-opt"
export CPLUS_INCLUDE_PATH=${BOOST__HOME}/include:$CPLUS_INCLUDE_PATH
export C_INCLUDE_PATH=${BOOST__HOME}/include:$C_INCLUDE_PATH
cd - 1>/dev/null # from /cvmfs/sft.cern.ch/lcg/releases/LCG_94/Boost/1.66.0/x86_64-slc6-gcc62-opt
export DAVIX__HOME="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/Davix/0.6.7/x86_64-slc6-gcc62-opt";
export PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/fftw3/3.3.4/x86_64-slc6-gcc62-opt/bin:$PATH";
export LD_LIBRARY_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/fftw3/3.3.4/x86_64-slc6-gcc62-opt/lib:$LD_LIBRARY_PATH";
export PKG_CONFIG_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/fftw3/3.3.4/x86_64-slc6-gcc62-opt/lib/pkgconfig:$PKG_CONFIG_PATH";
export FFTW__HOME="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/fftw3/3.3.4/x86_64-slc6-gcc62-opt";
export LD_LIBRARY_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/tbb/2018_U1/x86_64-slc6-gcc62-opt/lib:$LD_LIBRARY_PATH";
export TBB__HOME="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/tbb/2018_U1/x86_64-slc6-gcc62-opt";
export PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/castor/2.1.13-6/x86_64-slc6-gcc62-opt/bin:$PATH";
export LD_LIBRARY_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/castor/2.1.13-6/x86_64-slc6-gcc62-opt/lib64:$LD_LIBRARY_PATH";
export CASTOR__HOME="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/castor/2.1.13-6/x86_64-slc6-gcc62-opt";
export PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/mysql/5.7.20/x86_64-slc6-gcc62-opt/bin:$PATH";
export LD_LIBRARY_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/mysql/5.7.20/x86_64-slc6-gcc62-opt/lib:$LD_LIBRARY_PATH";
export PKG_CONFIG_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/mysql/5.7.20/x86_64-slc6-gcc62-opt/lib/pkgconfig:$PKG_CONFIG_PATH";
export LD_LIBRARY_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/mysql/5.7.20/x86_64-slc6-gcc62-opt/lib/plugin:$LD_LIBRARY_PATH";
export MYSQL__HOME="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/mysql/5.7.20/x86_64-slc6-gcc62-opt";
export PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/Grid/gfal/1.13.0-0/x86_64-slc6-gcc62-opt/bin:$PATH";
export LD_LIBRARY_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/Grid/gfal/1.13.0-0/x86_64-slc6-gcc62-opt/lib64:$LD_LIBRARY_PATH";
export PYTHONPATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/Grid/gfal/1.13.0-0/x86_64-slc6-gcc62-opt/lib64/python2.6/site-packages:$PYTHONPATH";
export GFAL__HOME="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/Grid/gfal/1.13.0-0/x86_64-slc6-gcc62-opt";
export PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/xrootd/4.8.4/x86_64-slc6-gcc62-opt/bin:$PATH";
export LD_LIBRARY_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/xrootd/4.8.4/x86_64-slc6-gcc62-opt/lib:$LD_LIBRARY_PATH";
export PYTHONPATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/xrootd/4.8.4/x86_64-slc6-gcc62-opt/lib/python2.7/site-packages:$PYTHONPATH";
export LD_LIBRARY_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/xrootd/4.8.4/x86_64-slc6-gcc62-opt/lib64:$LD_LIBRARY_PATH";
export XROOTD__HOME="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/xrootd/4.8.4/x86_64-slc6-gcc62-opt";
export PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/R/3.2.5/x86_64-slc6-gcc62-opt/bin:$PATH";
export LD_LIBRARY_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/R/3.2.5/x86_64-slc6-gcc62-opt/lib64:$LD_LIBRARY_PATH";
export PKG_CONFIG_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/R/3.2.5/x86_64-slc6-gcc62-opt/lib64/pkgconfig:$PKG_CONFIG_PATH";
export LD_LIBRARY_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/R/3.2.5/x86_64-slc6-gcc62-opt/lib64/R:$LD_LIBRARY_PATH";
export PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/zeromq/4.2.5/x86_64-slc6-gcc62-opt/bin:$PATH";
export LD_LIBRARY_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/zeromq/4.2.5/x86_64-slc6-gcc62-opt/lib:$LD_LIBRARY_PATH";
export PKG_CONFIG_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/zeromq/4.2.5/x86_64-slc6-gcc62-opt/lib/pkgconfig:$PKG_CONFIG_PATH";
export ZEROMQ__HOME="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/zeromq/4.2.5/x86_64-slc6-gcc62-opt";
export R__HOME="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/R/3.2.5/x86_64-slc6-gcc62-opt";
cd "/cvmfs/sft.cern.ch/lcg/releases/LCG_94/R/3.2.5/x86_64-slc6-gcc62-opt"
export R_HOME=${R__HOME}/lib64/R
cd - 1>/dev/null # from /cvmfs/sft.cern.ch/lcg/releases/LCG_94/R/3.2.5/x86_64-slc6-gcc62-opt
export PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/Grid/srm-ifce/1.13.0-0/x86_64-slc6-gcc62-opt/bin:$PATH";
export LD_LIBRARY_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/Grid/srm-ifce/1.13.0-0/x86_64-slc6-gcc62-opt/lib64:$LD_LIBRARY_PATH";
export PKG_CONFIG_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/Grid/srm-ifce/1.13.0-0/x86_64-slc6-gcc62-opt/lib64/pkgconfig:$PKG_CONFIG_PATH";
export SRM_IFCE__HOME="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/Grid/srm-ifce/1.13.0-0/x86_64-slc6-gcc62-opt";
export PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/Grid/dcap/2.47.7-1/x86_64-slc6-gcc62-opt/bin:$PATH";
export LD_LIBRARY_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/Grid/dcap/2.47.7-1/x86_64-slc6-gcc62-opt/lib64:$LD_LIBRARY_PATH";
export LD_LIBRARY_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/Grid/dcap/2.47.7-1/x86_64-slc6-gcc62-opt/lib64/dcap:$LD_LIBRARY_PATH";
export DCAP__HOME="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/Grid/dcap/2.47.7-1/x86_64-slc6-gcc62-opt";
export LD_LIBRARY_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/Vc/1.3.2/x86_64-slc6-gcc62-opt/lib:$LD_LIBRARY_PATH";
export LD_LIBRARY_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/Vc/1.3.2/x86_64-slc6-gcc62-opt/lib/cmake:$LD_LIBRARY_PATH";
export VC__HOME="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/Vc/1.3.2/x86_64-slc6-gcc62-opt";
export PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/oracle/11.2.0.3.0/x86_64-slc6-gcc62-opt/bin:$PATH";
export LD_LIBRARY_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/oracle/11.2.0.3.0/x86_64-slc6-gcc62-opt/lib:$LD_LIBRARY_PATH";
export ORACLE__HOME="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/oracle/11.2.0.3.0/x86_64-slc6-gcc62-opt";
export PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/GSL/2.5/x86_64-slc6-gcc62-opt/bin:$PATH";
export LD_LIBRARY_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/GSL/2.5/x86_64-slc6-gcc62-opt/lib:$LD_LIBRARY_PATH";
export PKG_CONFIG_PATH="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/GSL/2.5/x86_64-slc6-gcc62-opt/lib/pkgconfig:$PKG_CONFIG_PATH";
export GSL__HOME="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/GSL/2.5/x86_64-slc6-gcc62-opt";
export ROOT__HOME="/cvmfs/sft.cern.ch/lcg/releases/LCG_94/ROOT/6.14.04/x86_64-slc6-gcc62-opt";
cd "/cvmfs/sft.cern.ch/lcg/releases/LCG_94/ROOT/6.14.04/x86_64-slc6-gcc62-opt"
test -s $ROOT__HOME/bin/thisroot.sh && source $ROOT__HOME/bin/thisroot.sh
cd - 1>/dev/null # from /cvmfs/sft.cern.ch/lcg/releases/LCG_94/ROOT/6.14.04/x86_64-slc6-gcc62-opt
