#!/bin/bash

cd $HOME/tracking/code/Ctrax
stdeb_run_setup --default-distribution karmic-ads
version=`cat version.txt`
cd deb_dist/Ctrax-$version
echo "**********************************************************"
dpkg-buildpackage -rfakeroot -uc -us
cd ..
echo "**********************************************************"
dput private Ctrax_$version-1*changes
