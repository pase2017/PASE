#!/bin/bash

cmake -Duse_hypre=1 -Dhypre_root=~/Software/hypre-2.11.2/src/hypre \
      ../

#echo ""
#make VERBOSE=1
