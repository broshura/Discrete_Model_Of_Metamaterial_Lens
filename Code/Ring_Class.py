class Ring():
    def __init__(self, x, y, z, pos):
        self.x = x      # x position of the ring
        self.y = y      # y position of the ring
        self.z = z      # z position of the ring
        self.pos = pos  # orientation of the ring - "xy" or "xz" or "zy"

    def __repr__(self):
        return "x = " + str(round(self.x, 1)) + " y = " + str(round(self.y, 1)) + " z = " + str(round(self.y, 1)) + " orientation: " + self.pos
