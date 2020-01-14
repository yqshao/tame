Introduction
============

mdppp is a library written for processing molecular dynamisc
trajectories. It attemps to provide a easy way to analyze long
trajectories by introducing the ``FrameArray``: a numpy-array-like
data structure which is promptly updated when the trajectory is
continuously loaded.

Analysis protcols are expressed as recipes: documented scripts that
perform specific tasks. The recipes shipped with mdppp can be easily
invoked with the ``mdppp`` shell command, for instance: the following
command computes the conductivity of a electrolyte system using the
Green-Kubo relation:

``mdppp cond jacf prod.dump --c-type 1 --a-type 2``

Recipes can be written as individual Python scripts (welcome to
contribute!). Writting recipes in mdppp is intended to be easy and
intuitive. See the following chapter for a tutorials of writting such
analyses or recipes.
