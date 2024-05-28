# Object class for rings
class Ring():
    def __init__(self, x, y, z, pos, radius, width, L, C, R):
        self.x = x      # x position of the ring
        self.y = y      # y position of the ring
        self.z = z      # z position of the ring

        self.pos = pos  # orientation of the ring - "x" or "y" or "z"

        self.r = radius # Radius of the ring
        self.w = width  # Width of the strip

        self.L = L      # Self-inductance
        self.C = C      # Capacitance
        self.R = R      # Resistance

    def M(self, w):
        return self.R/1j/w - self.L + 1/(w ** 2 * self.C)

    def Z(self, w):
        return self.R - 1j * w * self.L + 1j/(w * self.C)

    def sigma(self, w):
        return 1/(self.R - 1j * w * self.L + 1j/(w * self.C))

#   Important to make parameters visible in console

    def __repr__(self):
        return f"x = {self.x} y = {self.y} z = {self.z} orientation: {self.pos} Radius: {self.r}"

    def __str__(self):
        return f"x = {self.x} y = {self.y} z = {self.z} orientation: {self.pos} Radius: {self.r}"
