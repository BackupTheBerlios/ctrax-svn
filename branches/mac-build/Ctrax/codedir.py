import os
codedir = os.path.dirname(__file__)
(head,tail) = os.path.split(codedir)
if tail == 'library.zip':
    codedir = head
