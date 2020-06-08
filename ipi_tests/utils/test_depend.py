#!/usr/bin/env python3
from ipi.utils.depend import *
import numpy as np
import threading
import time


class A(dobject):
    def __init__(self):

        dself = dd(self)

        dself.scalar = depend_value(name="a_scalar", value=1)
        dself.vector = depend_array(name="a_vector", value=np.zeros(10))

        dself.dscalar = depend_value(
            name="d_scalar", func=self.get_scalar, dependencies=[dself.scalar]
        )
        dself.dvector = depend_value(
            name="d_vector",
            func=self.get_vector,
            dependencies=[dself.dscalar, self.vector],
        )

    def get_scalar(self):
        return self.scalar * 2

    def get_vector(self):
        return self.vector * self.dscalar


class B(dobject):
    def __init__(self):

        dself = dd(self)

        dself.scalar = depend_value(name="b_scalar", value=1)
        dself.vector = depend_array(name="b_vector", value=np.zeros(10))

    def bind(self, A):

        self.A = A
        dself = dd(self)
        dself.dscalar = depend_value(
            name="db_scalar",
            func=self.get_scalar,
            dependencies=[dself.scalar, dd(A).dscalar],
        )
        dself.dvector = depend_value(
            name="db_vector",
            func=self.get_vector,
            dependencies=[dself.dscalar, self.vector],
        )

    def get_scalar(self):
        return self.scalar - self.A.scalar

    def get_vector(self):
        return self.vector * self.dscalar


myA = A()

mas = 5
mav = 0.3
mbs = 2
mbv = 3

myA.scalar = mas
myA.vector[:] = mav

myB = B()
myB.scalar = mbs
myB.vector[:] = mbv
myB.bind(myA)


def threadA(Aobj):
    for i in range(10000):
        Aobj.scalar = np.sqrt(i)
        time.sleep(0.0001)


def threadB(Aobj, Bobj):
    for i in range(10000):
        Aobj.scalar = i
        time.sleep(0.0001)
        Bobj.scalar = 4 + i
        Bobj.vector = np.sqrt(i)
        print(Bobj.dvector[0])


# myB.dvector should always contain B.vector*B.scalar- A.scalar
print(myB.dvector - mbv * (mbs - mas))

# now try with threads
ta = threading.Thread(target=threadA, name="TA", kwargs={"Aobj": myA})
ta.daemon = True
tb = threading.Thread(target=threadB, name="TB", kwargs={"Bobj": myB, "Aobj": myA})
tb.daemon = True
ta.start()
tb.start()
ta.join()
tb.join()

print(myB.dvector - myB.vector * (myB.scalar - myA.scalar))
