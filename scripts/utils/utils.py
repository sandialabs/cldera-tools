"""
Utilities
"""

###############################################################################
def expect(condition, error_msg, exc_type=SystemExit, error_prefix="ERROR:"):
###############################################################################
    """
    Similar to assert except doesn't generate an ugly stacktrace. Useful for
    checking user error, not programming error.

    >>> expect(True, "error1")
    >>> expect(False, "error2")
    Traceback (most recent call last):
        ...
    SystemExit: ERROR: error2
    """
    if not condition:
        msg = error_prefix + " " + error_msg
        raise exc_type(msg)
