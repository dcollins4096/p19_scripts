"""I hate having to open a file, unpickle from it, close the file.
This just unpickles from the filename given"""


import cPickle

def load(filename, *args, **kwargs):
    inputfile = open(filename,'r')
    output = cPickle.load(inputfile,*args,**kwargs)
    inputfile.close()
    return output

def dump(object,filename, *args, **kwargs):
    file = open(filename,'w')
    output = cPickle.dump(object,file,*args,**kwargs)
    file.close()

def bload(filename, *args, **kwargs):
    file = open(filename,'rb')
    output = cPickle.load(file,*args,**kwargs)
    file.close()
    return output

def bdump(object,filename, *args, **kwargs):
    print "bdump"
    file = open(filename,'wb')
    output = cPickle.dump(object,file,protocol=-1,*args,**kwargs)
    file.close()
