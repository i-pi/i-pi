""" Harmonic potential """

def harm_driver(cell, pos):
    """ Silly harmonic potential, with unit frequency in a.u."""
    pot = (pos ** 2).sum() * 0.5
    force = -pos  # makes a zero force with same shape as pos
    vir = cell * 0.0  # makes a zero virial with same shape as cell
    extras = "nada"
    return pot, force, vir, extras
