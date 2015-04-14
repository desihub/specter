"""
Implement a dictionary that caches the last N entries and throws the rest away
"""

class CacheDict(dict):
    """
    A dictionary that keeps only the last n items added
    """
    def __init__(self, n, d=None):
        """
        Create CacheDict with n keys cached.

        If d is a dict, use that to initialize the key/value pairs
        """
        self._keys = [None,]*n
        self._current = -1
        self._n = n

        if type(d) == dict:
            for key, value in d.items():
                self[key] = value

    #- Needed by pickle to reconstruct the object.
    #- Pickle tries to reconstruct the dictionary via __setitem__ before
    #- filling in _keys, _current, _n.  This allows it to create a new
    #- object first to get those right, before calling __setitem__
    def __reduce__(self):
        return type(self), (self._n, dict(self))
        
    def __setitem__(self, key, value):
        """Sets the key/value, possibly dropping an earlier cached key/value"""
        if key in self:
            return
       
        i = self._current = (self._current + 1) % self._n
        if self._keys[i] is not None:
            del self[self._keys[i]]
            
        self._keys[i] = key
        dict.__setitem__(self, key, value)
        
