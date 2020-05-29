# ------------------------------------------------------------------
# q-TIP4P/f dipole moment example (real part)
# ------------------------------------------------------------------

import numpy as np

name = "qTIP4P-cmumu-re"

method = "0thorder_re"

Ashape = (3,)

Bshape = (3,)

qm = -1.1128
qh = -0.50 * qm
gam = 0.73612


def Afunc0(qsum, psum):
    nndof, npl = qsum.shape
    natoms = nndof / 3
    nmolec = natoms / 3
    ans = np.zeros((npl, 3))
    x = np.zeros((nmolec, 3, 3))
    for i in range(npl):
        x[:] = qsum[:, i].reshape(nmolec, 3, 3)
        for j in range(nmolec):
            ans[i] += (
                (gam * x[j, 0, :] + 0.5 * (1.0 - gam) * (x[j, 1, :] + x[j, 2, :])) * qm
                + x[j, 1, :] * qh
                + x[j, 2, :] * qh
            )
    return ans


def Bfunc(qsum, psum):
    return Afunc0(qsum, psum)


# ------------------------------------------------------------------
# q-TIP4P/f dipole moment example (imaginary part)
# ------------------------------------------------------------------

# import numpy as np
#
# name = "qTIP4P-cmumu-im"
#
# method = "1storder_im"
#
# Ashape = (3,)
#
# Bshape = (3,)
#
# Change as appropriate #
# nmolec = 128
# natoms = nmolec*3
# nndof = natoms*3
# npl = 64
# ------------------------
#
# qm = -1.1128
# qh = -0.50*qm
# gam = 0.73612
# grad = np.zeros((3, nndof))
# x component, oxygen
# grad[0].reshape((nmolec,3,3))[:,0,0] = gam*qm
# y component, oxygen
# grad[1].reshape((nmolec,3,3))[:,0,1] = gam*qm
# z component, oxygen
# grad[2].reshape((nmolec,3,3))[:,0,2] = gam*qm
# x component, H1
# grad[0].reshape((nmolec,3,3))[:,1,0] = 0.5*(1.0-gam)*qm + qh
# y component, H1
# grad[1].reshape((nmolec,3,3))[:,1,1] = 0.5*(1.0-gam)*qm + qh
# z component, H1
# grad[2].reshape((nmolec,3,3))[:,1,2] = 0.5*(1.0-gam)*qm + qh
# x component, H2
# grad[0].reshape((nmolec,3,3))[:,2,0] = 0.5*(1.0-gam)*qm + qh
# y component, H2
# grad[1].reshape((nmolec,3,3))[:,2,1] = 0.5*(1.0-gam)*qm + qh
# z component, H2
# grad[2].reshape((nmolec,3,3))[:,2,2] = 0.5*(1.0-gam)*qm + qh
# grad = np.array([grad]*npl)
#
# def Afunc1(qsum):
#    return grad
#
# def Afunc0(qsum, psum):
#    nndof, npl = qsum.shape
#    natoms = nndof / 3
#    nmolec = natoms / 3
#    ans = np.zeros((npl, 3))
#    x = np.zeros((nmolec, 3, 3))
#    for i in xrange(npl):
#        x[:] = qsum[:,i].reshape(nmolec, 3, 3)
#        for j in xrange(nmolec):
#            ans[i] += (gam*x[j, 0, :] + 0.5*(1.0-gam)*(x[j, 1, :] + x[j, 2, :]))*qm + \
#                 x[j, 1, :]*qh + x[j, 2, :]*qh
#    return ans
#
# def Bfunc(qsum, psum):
#    qm = -1.1128
#    qh = -0.50*qm
#    gam = 0.73612
#    nndof, npl = qsum.shape
#    natoms = nndof / 3
#    nmolec = natoms / 3
#    ans = np.zeros((npl, 3))
#    x = np.zeros((nmolec, 3, 3))
#    for i in xrange(npl):
#        x[:] = qsum[:,i].reshape(nmolec, 3, 3)
#        for j in xrange(nmolec):
#            ans[i] += (gam*x[j, 0, :] + 0.5*(1.0-gam)*(x[j, 1, :] + x[j, 2, :]))*qm + \
#                 x[j, 1, :]*qh + x[j, 2, :]*qh
#    return ans
