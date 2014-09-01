IMPORTANT: ALPHA STATUS

symodesys
=========
Systems of ordinary differential equations as easy as Python.

# About
symodesys is an open-source Python framework for manipulation and
analysis of systems of ordinary differential equations.  It is tightly
built upon Sympy and can be seen as a set of convenience classes and
templates for code generation which turns Sympy into a poverful ODE
system analysis tool.

The code is generated using templates, and hence no restriction
exists on what langauges are used (C/Fortran/C++ used for GSL/ODEPACK/odeint).

Systems of ODE's are appearing in all of the natural sciences. Sympy
together with numpy + scipy + matplotlib etc makes python an excellent
language to build a versatile common framework for modelling of ODE systems.

Central to symodesys is NOT to reinvent the wheel. The less code put
into symodesys the better. Instead it aims to make sympy (as the most
promising free-software CAS in the views of the author) an efficient
tool for dealing with systems of ODE's.

Symodesys offers (for IVP's):
* General Interactive plotting through Enthought's Chaco library
* Code generation support for: GSL, Sundials and ODEPACK

Other software built using symodesys:
* [TODO] symchemkin - Flexible framework for analysis of systems of
  chemical reactions

# TODO
* More test cases where analytic solution is known.
* [TODO] Add options of dynamic step size GSL integrator
* [TODO] Extend ODEPACK support LSODE dense and banded.
* [TODO] Add code generation support for: RODAS (look at assimulo), ODEINT, RKC
* [TODO] Improve code generation to write loops when possible (let compiler unroll if optimal, decalare global const paramters), param structs etc.
* [TODO] Look into generating Jacobian as a sympy.SparseMatrix
* [TODO] Visualization of sensitivity of system by accepting uncertain input parameters (normal distr for starters..)
* [TODO] Optimization of paramters (fiting) to match e.g. numerical data.
* [TODO] Add routines for scaling (and automatic rescaling) in e.g. the IVP class. Note: this is a simple special case of variable transformation and may hence be handled by that code, only convenience routines are needed to be added.
* [TODO] Export of trajectories. (easy)
* [TODO] Add python code writer. (low priority)


# Prerequisities
(Tested against)
Python 2.7
Sympy 0.7.3
Cython 0.19.1

# Optional packages
Gnu sciencific library (GSL v. 1.15)
Sundials 
ODEPACK
[TODO] RODAS
[TODO] rkc
[TODO] vode (low priority)

# Installation
Installing and building symodesys can be a little challenging due to
external dependencies. For example, when debugging the setup script I
invoke it like:

``env DEBUG_DISTUTILS=1 COMPILER_VENDOR=gnu C_INCLUDE_PATH=~/.local/include:$C_INCLUDE_PATH python -m pudb setup.py build``

But if you have Sundials 2.5.0 on your systems standard include path
and you don't have Intel's fortran compiler installed (and my setup
script is clever enough) even pip should do the right thing.

To install Chaco etc. I needed to do:

```
sudo pip install --upgrade enable --allow-external PIL --allow-unverified PIL
sudo pip install chaco traitsui
apt-get install python-qt4-gl
```

And prepend ``python`` command as e.g. ``QT_API=pyqt python ode.py``.


# Status
symodesys is currently in its prototyping stage and hence big changes
to API and frequent refactoring is expected. If you are fine with this
you may of course fork and I'll happily accept pull requests.

# Roadmap
1. Refactor and refine API until convergence is reached (current
expected time frame for this: months)
2. Write at least one non-trivial python package based on symodesys
3. Write tests and expand doc strings and comments in code
4. Write proper documentation and doctests (in Sympy style)

# Known issues
## ODEPACK
opkda1.f fails using ifort 13.0.1 with the following error
``
opkda1.f(9498): error #6633: The type of the actual argument differs from the type of the dummy argument.   [RWORK]
     1   RWORK(LACOR), IA, JA, IC, JC, RWORK(LWM), RWORK(LWM), IPFLAG,
---------------------------------------------------^
compilation aborted for opkda1.f (code 1)
``
I have not looked into whether this is a real bug detected by the Intel compiler.

# Similar projects
* symneqsys
## Common dependencies
* symvarsub
* pycompilation

## Author
Written by Björn Dahlgren. Copyright 2012-2014.

## License
Open Source. Released under the very permissive simplified
(2-clause) BSD license. See LICENCE.txt for further details.

## Install SUNDIALS 2.5.0
Note that you need sundials 2.5.0 (too new for Ubuntu 12.04 LTS)
Download from http://computation.llnl.gov/casc/sundials/main.html

e.g.:

``
   mkdir /tmp/sundials_build/
   cd /tmp/sundials_build/
   cmake-gui ~/Downloads/sundials-2.5.0/
``

There are many ways to configure Sundials depending on your
system. I use Intel C compiler (icc) and ~/.local as prefix,
and ``CMAKE_C_FLAGS=-fPIC``. (static linking)

For running a symodesys simulation you need to make sure the include
and library dirs are searched, e.g.:

``env C_INCLUDE_PATH=~/.local/include:$C_INCLUDE_PATH LIBRARY_PATH=~/.local/lib python codegen.py -i sundials``

if you build a shared sundials library at a custom prefix path you must
ensure ``LD_LIBRARY_PATH`` is also set accordingly.
