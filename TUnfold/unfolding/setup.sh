#!/bin/zsh

export LCGENV_PATH=/cvmfs/sft.cern.ch/lcg/releases
#cat $LCGENV_PATH/HEAD/Readme.md
export LCG_VERSION=LCG_94
export PLATFORM=x86_64-slc6-gcc62-opt
#export PLATFORM=x86_64-slc6-gcc49-opt

echo "# setup new env" > lcgenv.sh

$LCGENV_PATH/lcgenv/latest/lcgenv -p ${LCG_VERSION} ${PLATFORM} gcc >> lcgenv.sh                           
$LCGENV_PATH/lcgenv/latest/lcgenv -p ${LCG_VERSION} ${PLATFORM} Boost >> lcgenv.sh
#$LCGENV_PATH/lcgenv/latest/lcgenv -p ${LCG_VERSION} ${PLATFORM} lhapdf 6.2.0 >> lcgenv.sh
#$LCGENV_PATH/lcgenv/latest/lcgenv -p ${LCG_VERSION} ${PLATFORM} lhapdf 6.1.6 >> lcgenv.sh
$LCGENV_PATH/lcgenv/latest/lcgenv -p ${LCG_VERSION} ${PLATFORM} ROOT >> lcgenv.sh
#$LCGENV_PATH/lcgenv/latest/lcgenv -p ${LCG_VERSION} ${PLATFORM} fastjet >> lcgenv.sh

#$LCGENV_PATH/lcgenv/latest/lcgenv -p ${LCG_VERSION} ${PLATFORM} automake >> lcgenv.sh
#$LCGENV_PATH/lcgenv/latest/lcgenv -p ${LCG_VERSION} ${PLATFORM} autoconf >> lcgenv.sh
#$LCGENV_PATH/lcgenv/latest/lcgenv -p ${LCG_VERSION} ${PLATFORM} libtool >> lcgenv.sh
#$LCGENV_PATH/lcgenv/latest/lcgenv -p ${LCG_VERSION} ${PLATFORM} m4 >> lcgenv.sh
#$LCGENV_PATH/lcgenv/latest/lcgenv -p ${LCG_VERSION} ${PLATFORM} lhapdfsets >> lcgenv.sh
#$LCGENV_PATH/lcgenv/latest/lcgenv -p ${LCG_VERSION} ${PLATFORM} doxygen >> lcgenv.sh
#$LCGENV_PATH/lcgenv/latest/lcgenv -p ${LCG_VERSION} ${PLATFORM} swig >> lcgenv.sh
#$LCGENV_PATH/lcgenv/latest/lcgenv -p ${LCG_VERSION} ${PLATFORM} zlib >> lcgenv.sh
#$LCGENV_PATH/lcgenv/latest/lcgenv -p ${LCG_VERSION} ${PLATFORM} yoda 1.6.7 >> lcgenv.sh
#$LCGENV_PATH/lcgenv/latest/lcgenv -p ${LCG_VERSION} ${PLATFORM} Cairo >> lcgenv.sh

source lcgenv.sh

#asetup Athena,21.0.61,here

export CC=$(which gcc)
export CXX=$(which g++)
#export LD=$CC
