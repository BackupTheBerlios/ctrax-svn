#!/bin/bash

sudo echo "" >> /dev/null && \

sudo rm -r /Library/Frameworks/Python.framework/Versions/2.6/lib/python2.6/site-packages/Ctrax_mac* && \
sudo rm -r build/lib.macosx-10.3-fat-2.6/Ctrax_mac && \

python setup.py build && \
sudo python setup.py install && \
Ctrax-script-mac.py
