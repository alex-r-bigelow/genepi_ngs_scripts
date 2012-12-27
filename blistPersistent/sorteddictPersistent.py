from blist import sorteddict
from durus.persistent import Persistent

class sorteddict(sorteddict, Persistent):
    def __del__(self, key):
        super(sorteddict,self).__del__(key)
        super(Persistent,self)._p_note_change(self)
    
    def __setitem__(self, key, value):
        super(sorteddict,self).__setitem__(self, key, value)
        super(Persistent,self)._p_note_change(self)
    
    def clear(self):
        super(sorteddict,self).clear(self)
        super(Persistent,self)._p_note_change(self)
    
    def update(self, other):
        super(sorteddict,self).update(self, other)
        super(Persistent,self)._p_note_change(self)