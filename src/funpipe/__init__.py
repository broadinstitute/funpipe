# src
# └── package
#     ├── __init__.py
#     ├── moduleA.py
#     ├── moduleB.py
#     └── subpackage
#         ├── __init__.py
#         └── moduleC.py

# funpipe
# ├── .readthedocs.yml
# ├── CHANGELOG.md
# ├── CONDUCT.md
# ├── CONTRIBUTING.md
# ├── docs
# │   └── ...
# ├── LICENSE
# ├── README.md
# ├── pyproject.toml
# ├── src
# |   └── funpipe(package)
# |       ├── __init__.py
# |       ├── moduleA.py
# |       ├── moduleB.py
# |       └── variant(subpackage)
# |           ├── __init__.py
# |           └── moduleC.py
# └── tests
#     └── ...


# Absolute && Relative imports

# ===================================

# Import from moduleA in moduleB:

# from package.moduleA import XXX

# from .moduleA import XXX

# ====================================

# Import from moduleA in moduleC:

# from package.moduleA import XXX

# from ..moduleA import XXX

# ====================================

# Import from moduleC in moduleA:

# from package.subpackage.moduleC import XXX

# from .subpackage.moduleC import XXX