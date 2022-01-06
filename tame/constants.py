# Length
Ang = 1e-10
cm = 1e-2

# Time
ns = 1e-9
ps = 1e-12
fs = 1e-15

# Pressure
atm = 101325

# Constants
NA = 6.0221409
kB = 1.38064852e-23
e = 1.60217662e-19

# Collections of pre_defines unit sets
class unit_collection():
    pass

units = unit_collection
units.length   = Ang
units.time     = fs
units.velocity = Ang/fs
