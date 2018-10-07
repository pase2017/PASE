#!/bin/bash

cmake -Duse_hypre=1 -Dhypre_root=~/Software/hypre/hypre-2.12.0/src/hypre \
      ../

#echo ""
#make VERBOSE=1
