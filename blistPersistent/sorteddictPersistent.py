from blist import sorteddict
from durus.persistent import Persistent

class sorteddict(sorteddict, Persistent):
    def __del__(self, key):
        super(sorteddict,self).__del__(key)
        super(Persistent,self)._p_note_change()
    
    def __setitem__(self, key, value):
        super(sorteddict,self).__setitem__(key, value)
        super(Persistent,self)._p_note_change()
    
    def clear(self):
        super(sorteddict,self).clear()
        super(Persistent,self)._p_note_change()
    
    def update(self, *args, **kw):
        super(sorteddict,self).update(self, *args, **kw)
        super(Persistent,self)._p_note_change()