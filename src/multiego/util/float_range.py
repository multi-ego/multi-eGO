class FloatRange(object):
    """
    Represents a range of floating-point values.

    Attributes:
    - start (float): The starting value of the range.
    - end (float): The ending value of the range.
    """

    def __init__(self, start, end):
        """
        Initializes a FloatRange object with a start and end value.

        Args:
        - start (float): The starting value of the range.
        - end (float): The ending value of the range.
        """
        self.start = start
        self.end = end

    def __eq__(self, other):
        """
        Checks if a value falls within the defined range.

        Args:
        - other (float): The value to check against the range.

        Returns:
        - bool: True if the value falls within the range, False otherwise.
        """
        return self.start <= other <= self.end
