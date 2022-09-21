class Ring():
    def __init__(self, x, y, z, pos):
        self.x = x      # x position of the ring
        self.y = y      # y position of the ring
        self.z = z      # z position of the ring
        self.pos = pos  # orientation of the ring - "xy" or "xz" or "zy"

    def __str__(self):
        return f"x = {self.x} y = {self.y} z = {self.z} orientation: {self.pos}"
