#   Object class to simplify working with rings orientation and coordinates

class Ring():
    def __init__(self, x, y, z, pos, radius):
        self.x = x      # x position of the ring
        self.y = y      # y position of the ring
        self.z = z      # z position of the ring
        self.pos = pos  # orientation of the ring - "xy" or "xz" or "zy"
        self.r = radius # Radius of the ring

#   Important to make parameters visible in console

    def __str__(self):
        return f"x = {self.x} y = {self.y} z = {self.z} orientation: {self.pos} Radius: {self.r}"
