# Electrode class: defines the attributes of an Electrode object
class Electrode:

    # initialization method
    def __init__(self, r_e, x, y, index):
        self.r_e = r_e # electrode radius
        self.x = x # x position on the grid
        self.y = y # y position on the grid
        self.index = index
        # number of wire crossings?

    # printing method
    def __str__(self):
        ret = "Electrode details:" + "\n\tIndex: " + str(self.index) + "\n\tRadius: " \
                + str(self.r_e) + "\n\tPosition: (" + str(self.x) + ", " + \
                str(self.y) + ")"
        return ret

    # returns array with coordinates of center
    def getCenter(self):
        return [self.x, self.y]

    # returns radius
    def getRadius(self):
        return self.r_e

    # returns index number
    def getIndex(self):
        return self.index
