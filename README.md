symodesys
=========

symodesys is an open-source Python framework for manipulation and analysis of systems of ordinary differential equations.
It is tightly built upon sympy and can be seen as a set of convenience functions and templates for code generation which turns
Sympy into a poverful ODE system analysis tool.

Systems of ODE's are appearing in all of the natural sciences. Sympy together with numpy + scipy + matplotlib etc makes python an excellent
language to build a well tested and versatile common framework for models of ODE systems.

Central to symodesys is NOT to reinvent the wheel. The less code put into symodesys the better. Instead it aims to make sympy (as the most
promising free-software CAS in the views of the author) an efficient tool for dealing with systems of ODE's.

Symodesys offers (for IVP's):
* [TODO] Interfacing and _code generation_ for state of the art numerical packages (GSL, Sundials, ODEPACK)
* [TODO] Interactive plotting through Enthought's Chaco library
* Static plot generation through matplotlib
* [TODO] Optimization of paramters (fiting) to match e.g. numerical data.
* [TODO] Export of trajectories

Other software built using symodesys:
* [TODO] symchemkin - Flexible framework for analysis of systems of chemical reactions

## About
Written by Björn Dahlgren. Copyright 2012.

## License
Open Soucrce. Released under the very permissive simplified (2-clause) BSD license. See LICENCE.tx for further details.
