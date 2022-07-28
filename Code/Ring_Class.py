class Ring():
    def __init__(self, x, y, z, pos):
        self.x = x      # x position of the ring
        self.y = y      # y position of the ring
        self.z = z      # z position of the ring
        self.pos = pos  # orientation of the ring - "xy" or "xz" or "zy"

    def __repr__(self):
        return "x = " + str(self.x) + " y = " + str(self.y) + " z = " + str(self.z) + " orientation: " + self.pos
