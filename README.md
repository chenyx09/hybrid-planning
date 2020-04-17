Introduction
============

This project is about robot exploration using Lidar point cloud and incrementally construct an abstraction for high level control.

Prerequisites
============
IRIS from MIT Tedrake's group
<https://github.com/rdeits/iris-distro>

Installation is a bit tricky, if you install with python 2, you need to add
YOUR_IRIS_DIR/build/install/lib/python3.X/dist-packages to your python path where 3.X is your python3 version.

IRIS also need other prerequisites such as Mosek, check the manual on git for instructions.

OSQP
Check <https://osqp.org/docs/interfaces/python.html>

Pycddlib
Check <https://pypi.org/project/pycddlib/>

TuLip
Check <https://github.com/tulip-control/tulip-control>
