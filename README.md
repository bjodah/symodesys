IMPORTANT: NOT READY FOR PUBLIC ADOPTION: PRE ALPHA STATUS

symodesys
=========
Systems of ordinary differential equations as easy as Python.

# About
symodesys is an open-source Python framework for manipulation and
analysis of systems of ordinary differential equations.  It is tightly
built upon Sympy and can be seen as a set of convenience classes and
templates for code generation which turns Sympy into a poverful ODE
system analysis tool.

Systems of ODE's are appearing in all of the natural sciences. Sympy
together with numpy + scipy + matplotlib etc makes python an excellent
language to build a well tested and versatile common framework for
models of ODE systems.

Central to symodesys is NOT to reinvent the wheel. The less code put
into symodesys the better. Instead it aims to make sympy (as the most
promising free-software CAS in the views of the author) an efficient
tool for dealing with systems of ODE's.

Symodesys offers (for IVP's):
* [TODO] Add options of dynamic step size GSL integrator
* [TODO] Add code generation support for: Sundials, ODEPACK and RODAS (look at assimulo)
* [TODO] Improve code generation to write loops when possible (let compiler unroll if optimal, decalare global const paramters)
* [TODO] Look into generating Jacobian as a sympy.SparseMatrix
* [TODO] General Interactive plotting through Enthought's Chaco library
* [TODO] Optimization of paramters (fiting) to match e.g. numerical data.
* [TODO] Add routines for scaling (and automatic rescaling) in e.g. the IVP class.
* [TODO] Export of trajectories

Other software built using symodesys:
* [TODO] symchemkin - Flexible framework for analysis of systems of
  chemical reactions

# Prerequisities
Python 2.x (TODO: determine x, probably x=7)
Sympy 0.y (TODO: determine y, probably y=6)
Cython 0.z (TODO: determine z, probably z=18)

# Optional packages
Sundials u (TODO: determine u, probably u=)

# Installation
Clone git-repo. Add to $PYTHONPATH

# Status
symodesys is currently in its prototyping stage and hence big changes
to API and frequent refactoring is expected. If you are fine with this
you may of course fork and I'll happily accept pull requests.

# Roadmap
1. Refactor and refine API until convergence is reached (current
expected time frame for this: months)
2. Write at least one non-trivial python package based on symodesys
3. Write tests and expande doc strings and comments in code
4. Write proper documentation and doctests (in Sympy style)
5. Announce on Sympy maillist for comments.
6. If received well, and as the project matures, it might be included
in Sympy in some way.

# Philosophy
* Small codebase to aid future maintainence
* Write programs that do one thing and do it well (Unix) - let's avoid feature creep
* Code reuse - Leverage other peoples code as much as possible (external dependencies
  is _not_ a problem thanks to PyPI) and provide means good class hierarchy for 3rd party extensions.
* KISS (Keep it simple stupid) - Minimum amount of magic code in the project
* Try to be in line with best practices (code reviews are welcome).

## Author
Written by Björn Dahlgren. Copyright 2012-2013.

## License
Open Soucrce. Released under the very permissive simplified
(2-clause) BSD license. See LICENCE.tx for further details.
