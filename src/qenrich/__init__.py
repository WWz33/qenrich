from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("qenrich")
except PackageNotFoundError:  # source checkout without an installed distribution
    __version__ = "0.1.0"
