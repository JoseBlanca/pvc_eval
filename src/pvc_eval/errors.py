class PvcEvalInputError(ValueError):
    """A fault of the input that stops the run (doc/spec.md §3).

    Its message is what the user reads: the file, the line or the position,
    and how to put it right. The command line prints it without a traceback.
    """
