# oces_viewer

This repository provides a (draft) library module `oces::reader` and a couple of programs for reading Open Compound Eye Standard (OCES) files.
Open Compound Eye Standard was created by Blayze Millward during his PhD work with Mike Mangan at The University of Sheffield.
This reader does not yet provide a full implementation of the standard.

* **oces2cray** Opens an OCES file and outputs the ommatidia data in libcompound-ray's CSV 'eye file format' to stdout
* **oces_viewer** Opens an OCES file and visualizes the eye and head.

![Screenshot of Velox head and eyes](https://github.com/sebsjames/oces_viewer/blob/main/data/velox-head.png?raw=true)


## Getting the code

Clone this repo from github with submodules.

```bash
git clone git@github.com:sebsjames/oces_viewer
cd oces_viewer
git submodule init
git submodule update
```

The submodules are: [sebsjames/mathplot](https://github.com/sebsjames/mathplot), [sebsjames/maths](https://github.com/sebsjames/maths), [args](https://github.com/Taywee/args) and [tinygltf](https://github.com/sebsjames/tinygltf).


## Build tools

To compile oces_viewer, you need:

* clang-20 or higher.
* cmake version 3.28.5 or higher. Either `apt install cmake` on Ubuntu 25+ or download and build cmake from the cmake.org download page (it's an easy, reliable compile).
* ninja (`apt install ninja-build`)


## Build

Do a cmake build

```bash
cd oces_viewer
mkdir build
pushd build
CC=clang-20 CXX=clang++-20 cmake .. -GNinja
ninja
popd
```

## Run the program

To see command line options:
```bash
./build/oces_viewer -h
```

To view the velox-head data:
```bash
./build/oces_viewer -f data/velox-head.gltf -r 0.0002 -c -0.00056,0.00005,-0.00005 -x1,0,0 -o0.2 -s -y
```
