from bhr.pipe import Pipe


class Borehole:

    def __init__(self, inputs: dict) -> None:

        self.length = inputs["depth"]
        self.pipe = Pipe({**inputs["pipe"], **inputs["fluid"], "length": self.length * 2})
        self.filling = None
