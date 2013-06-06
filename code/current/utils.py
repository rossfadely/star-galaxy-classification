import numpy as np
import ctypes as ct

def ctype_2D_double_pointer(arr):
    assert isinstance(arr,np.ndarray)
    ctarr = arr.shape[0]*[None]
    for i in range(arr.shape[0]):
        ctarr[i] = np.ctypeslib.as_ctypes(arr[i,:])
    arrp = (ct.POINTER(ct.c_double)*arr.shape[0])(*ctarr)
    return arrp

def get_train_and_test(Ndata,truth,label_types,fraction,
                       rel_fractions=None):

    # shuffle
    ind = np.random.permutation(Ndata)

    # get training indicies
    Ntrain = np.int(np.round(Ndata * fraction))
    Ntest  = Ndata - Ntrain
    train_ind = ind[:Ntrain]
    train_labels = truth[train_ind]

    # organize training data into classes
    class_train_ind = {}
    for l in label_types:
        class_train_ind[l] = np.sort(train_ind[
                np.where(train_labels==l)[0]])
    
    # relative fractions?
    if rel_fractions!=None:
        assert isinstance(rel_fractions,np.ndarray)
        N = np.array([class_train_ind[l].shape[0] for l in label_types])
        idx = np.where(N==N.min())[0]
        rel_fractions /= rel_fractions[idx]
        N = np.round(rel_fractions * N.min())
        train_ind = np.array([],dtype=np.int)
        for i,l in enumerate(label_types):
            class_train_ind[l] = class_train_ind[l][:N[i]]
            train_ind = np.append(train_ind,class_train_ind[l])

    # sort and get test indicies
    test_ind = np.sort(ind)
    train_ind = np.sort(train_ind)
    test_ind = np.delete(test_ind,train_ind)
    test_labels = truth[test_ind]

    # organize test data into classes
    class_test_ind = {}
    for l in label_types:
        class_test_ind[l] = np.sort(test_ind[
                np.where(test_labels==l)[0]])

    return class_train_ind,class_test_ind,train_ind,test_ind

def quick_ML_assess(star_chi2,gal_chi2,truth):

    correct = 0.
    cstar = 0.
    cgalaxy = 0.
    for i in range(truth.shape[0]):
        if ((star_chi2[i].min()<gal_chi2[i].min()) & (truth[i]==1)):
            correct += 1.
            cstar += 1.
        if ((star_chi2[i].min()>gal_chi2[i].min()) & (truth[i]==0)):
            correct += 1.
            cgalaxy += 1.
    ind = np.where(truth!=-99)
    print 'ML total'
    print correct, truth[ind].shape[0],correct/truth[ind].shape[0]
    ind = np.where(truth==1)[0]
    print 'ML star'
    print cstar,truth[ind].shape[0],cstar/truth[ind].shape[0]
    ind = np.where(truth==0)[0]
    print 'ML gal'
    print cgalaxy,truth[ind].shape[0],cgalaxy/truth[ind].shape[0],'\n'

def quick_HB_assess(star_prob,gal_prob,truth):
    correct = 0.
    cstar = 0.
    cgalaxy = 0.
    for i in range(truth.shape[0]):
        if ((star_prob[i]>gal_prob[i]) & (truth[i]==1)):
            correct += 1.
            cstar+=1.
        if ((star_prob[i]<gal_prob[i]) & (truth[i]==0)):
            correct += 1.
            cgalaxy += 1.
    ind = np.where(truth!=-99)
    print 'HB total'
    print correct, truth[ind].shape[0],correct/truth[ind].shape[0]
    ind = np.where(truth==1)[0]
    print 'HB star'
    print cstar,truth[ind].shape[0],cstar/truth[ind].shape[0]
    ind = np.where(truth==0)[0]
    print 'HB gal'
    print cgalaxy,truth[ind].shape[0],cgalaxy/truth[ind].shape[0],'\n'
