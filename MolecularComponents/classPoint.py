import string
from math import sqrt

class Point:
    def __init__(self, x, y, z):
        self.update(x, y, z)
    def update(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z
    def pdb_print(self):
        print "%7.3f %7.3f %7.3f" % (self.x, self.y, self.z)
    def get_coordinate_list(self):
        new_list = []
        new_list.append(self.x)
        new_list.append(self.y)
        new_list.append(self.z)
        return new_list
    def get_coordinates(self):
        return self.x,self.y,self.z
    def get_distance_from(self, other_point):
        return self.dist(other_point)
    def dist(self, other_point):
        return sqrt(pow(self.x-other_point.x,2) + pow(self.y-other_point.y,2) + pow(self.z-other_point.z,2))
    def point_line_distance(self, p1, p2):
        """
        calculate p4, the point on the line between p1 and p2 to which p3 is perpendicular
        then return the distance between p3 and p4.
           1_____3_____2         is relevant, but, 
        3  1___________2         and
           1___________2    3    are not.
        
        returns -1 if p3 is not perpendicular to p1-p2 between them. 
        
        however, if 3 is the alpha carbon of 1 or 2, shielding is surely not as bad when
        the alpha carbon is just to the inside of the line. so, a simple weighting function
        attempts to fix the problem by multiplying the actual distance by two. Therefore
        less shielding is reported when considering the alpha carbons of the amino acids
        being investigated.
        """
        dot = 0
        dot = dot + (( p2.x - p1.x )*( self.x - p1.x))	# dot product of a2a1 and a2a3, x component
        dot = dot + (( p2.y - p1.y )*( self.y - p1.y))	# y component
        dot = dot + (( p2.z - p1.z )*( self.z - p1.z))	# z component
        dst = p1.dist(p2)
        if dst == 0:
            # dst is only zero when alpha and beta atoms are being compared from the same residue
            # This should only occur for GLY or alpha carbon models. In these cases, the alpha carbon
            # is not shielding the beta carbon, so return a high number.
            return 10.0
        else:
            dot = dot / p1.dist(p2)
        if ((dot > dst) or (dot < 0)):
            return -1
        x = ((dot / dst) * (p2.x - p1.x)) + p1.x
        y = ((dot / dst) * (p2.y - p1.y)) + p1.y
        z = ((dot / dst) * (p2.z - p1.z)) + p1.z
        p3 = Point(x,y,z)
        dst = p2.dist(p3)
        return dst;
