import numpy as np
import ctypes as ct

def ctype_2D_double_pointer(arr):
    assert isinstance(arr,np.ndarray)
    ctarr = arr.shape[0]*[None]
    for i in range(arr.shape[0]):
        ctarr[i] = np.ctypeslib.as_ctypes(arr[i,:])
    arrp = (ct.POINTER(ct.c_double)*arr.shape[0])(*ctarr)
    return arrp
