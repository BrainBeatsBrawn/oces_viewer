#!/bin/bash

# Command to view velox-head.gltf:

pushd build && ninja && popd  &&  ./build/oces_viewer -f data3/velox-head.gltf -r 0.0002 -c -0.00054,-0.00009,-0.00002 -x1,0,0 -o0.2 -s -y -v

# Argument description:
#
# -f path to the glTF
#
# These options are all used to project the ommatidial values onto a 2D plane for convenient visualization:
#
# -r Radius of the sphere used to project the ommatidia onto a 2D plane
# -c Location of sphere centre for the projection (this is mirrored for the second eye)
# -x An axis for the additional rotation of the ommatidia before they are projected using the default equirectangular projection
# -o The amount in radians for the additional rotation about the axis
# -s (Optional) show the projection sphere
# -y (Optional) show the projection rays
