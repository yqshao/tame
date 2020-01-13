===========================================
Molecular Dynamics Post-Processing Programm
===========================================

Introduction
============

This is a **work in progress** set of tools to post process MD
trajectories, mainly aiming for calculating transport properties for
electrolyte systems. **The code currently works only for LAMMPS output
files.**

Code structure
==============

FrameData and FrameArray
------------------------

``FrameData`` is the core data type in mdppp. It represents a set of
named variables which can be loaded by a generator. Children of a
``FrameData`` are ``FrameArray`` s.

``FrameArray`` is a subclass of numpy ndarray. It accepts usual
operations of numpy arrays (simple arithmetic operations, numpy
ufuncs, etc.; **NOTE** ``np.expand_dims`` is **not** usable yet, but
one can use ``array[np.newaxis,:]`` instead), with the difference that
the arrays are automatically updated when the parent ``FrameData``
updates.

``FrameArray`` can be extened to take different operations during
frame update, for instance, a time correlation function which requires
some intermediate variable to be computed and cached whenever a new
frame is read. This can be done by implementing a customzied
``update()`` method.

A set of common operations: e.g. time average ``tavg``, time
autocorrelation function ``tacf``, mean square displacement ``tmsd``
are implemented in mdppp, they can be called as functions and returns
a customized ``FrameArray``.

Modules
-------

- ``io``: functions for loading ``FrameData`` from trajectories or log files.
- ``analysis``: functions or ``FrameOp`` for ``FrameArray`` of given
  tasks, e.g. computing RDF, MSD or autocorrelations.
- ``recipes``: functions for specialized tasks (calculating N-E/G-K
  conductivity, viscosity, etc.)

Usage
=====

Installation
------------

.. code-block:: bash
		
   pip install git+https://github.com/yqshao/mdppp.git

Using a recipe
--------------

Recipes are organzied in subcommands of the mdppp command, for
example, the following command computes the Green-Kubo conductivity
taking element 1 and 2 as cation and anion species:

.. code-block:: bash

   mdppp cond jacf prod.dump --c-type 1 --a-type 2

Writing a recipe
----------------

**Example is yet to be written**		
