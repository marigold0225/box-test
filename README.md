## kinetic approach in a box

## how to install 

```
git clone xxx

./setup.ps1 (this is a script for windows powershell)
```
## Attention! 
The project is written in C++ and built using Ninja. Make sure you have the necessary dependencies installed.
You can modify the CMake configuration file CMakeLists.txt, to choose a compilation method that suits you.
Of course, you can also compile on a Linux system.

## how to use 
project structure:

```
│   CMakeLists.txt
│   README.md
│   setup.ps1
│
├───build
├───include
│       Box.h
│       Particle.h
│       Simulation.h
│
├───out
│       Box_test.exe
│       box_test.png
│       initial_state.csv
│       input.txt
│       particles.csv
│       plot.py
│
└───src
        Box.cpp
        main.cpp
        Particle.cpp
        Simulation.cpp
```
The input file 'input.txt' should be placed alongside the executable 'Box_test.exe'.

