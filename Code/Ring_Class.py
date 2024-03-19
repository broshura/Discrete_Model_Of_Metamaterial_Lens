# Object class for rings
class Ring():
    def __init__(self, x, y, z, pos, radius, width):
        self.x = x      # x position of the ring
        self.y = y      # y position of the ring
        self.z = z      # z position of the ring
        self.pos = pos  # orientation of the ring - "x" or "y" or "z"
        self.r = radius # Radius of the ring
        self.w = width  # Width of the strip

#   Important to make parameters visible in console

    def __repr__(self):
        return f"x = {self.x} y = {self.y} z = {self.z} orientation: {self.pos} Radius: {self.r}"

    def __str__(self):
        return f"x = {self.x} y = {self.y} z = {self.z} orientation: {self.pos} Radius: {self.r}"
