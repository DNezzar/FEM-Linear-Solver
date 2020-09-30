# FEM-Linear-Solver

**Files:**

- *fem.py* : It is the Python input file where you have to specified the input parameters such as coordinates, element connectivity, ..., and where you will have the ouput (displacment ,strain, stress and graph)

- *libmod.dll* : It is the library that comes with the *module.pyd* file

- *module.pyd* : File imported at the beginning of the *fem.py* file. It is obtained from the conversion and compilation of Fortran file to Python using F2PY module (Numpy)

- *solver.f90* : Fortran file with the main program 

**Manual:**

- Download only the *fem.py*, *libmod.dll* and *module.pyd*. 

- Modify the input parameters in the *fem.py* file as you want

- In command console (in files directory), run the *fem.py* like: ```$ python fem.py```

**Requierement:**

- Python3 (any version) but it has to be 64 bit

- Numpy and Matplotlib have to be installed

**Limitation:**

- ***Cannot impose diplacments***
- ***Cannot impose distributed forces***
