"""The one exception type shared across this subtree.

It lives apart from bundle.py so that `python -m visualization serve` can report errors
without importing numpy, biopython and jsonschema. Serving a directory needs none of
them, and a fresh clone has no virtualenv to supply them.
"""


class BundleError(ValueError):
    pass
