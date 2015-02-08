"""
Implement a dictionary that caches the last N entries and throws the rest away
"""

class CacheDict(dict):
    def __init__(self, n):
        self._keys = [None,]*n
        self._current = -1
        self._n = n
        
    def __setitem__(self, key, value):
        if key in self:
            return
            
        i = self._current = (self._current + 1) % self._n
        if self._keys[i] is not None:
            del self[self._keys[i]]
            
        self._keys[i] = key
        dict.__setitem__(self, key, value)
        
        