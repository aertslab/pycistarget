from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("pycistarget")
except PackageNotFoundError:
    __version__ = "unknown"
