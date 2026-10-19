"""What the user reads when something in the input goes wrong.

Public: PvcEvalInputError, the error that stops the run, and format_count, which writes a
count with its noun for the messages of errors and warnings.
"""


class PvcEvalInputError(ValueError):
    """A fault of the input that stops the run (doc/spec.md §3).

    Its message is what the user reads: the file, the line or the position,
    and how to put it right. The command line prints it without a traceback.
    """


def format_count(number: int, noun: str) -> str:
    """A number and its noun, in the singular for 1: "1 sample", "3 samples"."""
    return f"{number} {noun}" if number == 1 else f"{number} {noun}s"
