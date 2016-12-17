from math import *
from scipy import *

"""
   Returns the angle made by the triangle P-Q-R in degrees
"""
def angle (P,Q,R):
  P_T=P-Q
  Q_T=r_[0,0,0] # Origin
  R_T=R-Q
  # Fix for a bug that if the two vectors are the same
  # then temp should equal 1.0, but it comes out slightly
  # larger than 1.0 because of a loss in precision.
  if (P_T[0] == R_T[0] and P_T[1] == R_T[1] and P_T[2] == R_T[2]):
    temp = 1.0
  else:
    temp = dot(P_T,R_T)/(mag(P_T)*mag(R_T))

  return acos (temp)*180/pi

"""
   Returns the distance between two vectors with an arbitrary
   number of components.
"""
def distance (v1,v2):
  if (len(v1) != len(v2)):
    print "ERROR: Vectors are not the same length"
    return None
  sumSqr=0
  for i in range(len(v1)):
    sumSqr+=(v1[i]-v2[i])*(v1[i]-v2[i])
  return sqrt(sumSqr)

"""
   Return the magnitude of a vector with an arbitrary number
   of components.
"""
def mag (v):
  sumSqr=0.0
  for comp in v:
    sumSqr+=comp*comp
  return sqrt(sumSqr)

"""
   Compute the cross product of two 3D vectors.
   Returns a scipy vector.
"""
def cross_product (a,b):
  return r_[a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]]

"""
   The rotation matrix is calculated after all of the vectors are translated
   to point P (i.e. P is the new origin). Therefore, you will have to translate
   the new position back to the original coordinate system. I've provided a function
   called new2old, which takes a point in the new coordinate system and returns it
   as a vector in the original coordinate system. The rotation matrix that is returned
   will transfer a point in the new coordinate system back to the original coordinate
   system (with respect to P).
"""
def findRotationMatrix (P,Q,R):
  # Make scipy vectors out of P,Q,R
  Pv=r_[P.x,P.y,P.z]
  Qv=r_[Q.x,Q.y,Q.z]
  Rv=r_[R.x,R.y,R.z]    
  """ 1.) Find a vector that is normal to the plane (P,Q,R) by
          finding two vectors (a,b) and finding a point orthogonal
          to both. Vectors a and b lie in the plane translated to
          origin. """
  a=Qv-Pv
  b=Rv-Pv
  n=cross_product(a,b)
  n=n/mag(n) # Normalize
  """ 2.) Create the new axes. These will be used to transform the
          original axis into the plane """
  y_new=a/mag(a)
  z_new=n # Now all vectors in the new axis with a 0 z component will lie
          # in the plane
  x_new=cross_product(y_new,z_new)/mag(cross_product(y_new,z_new))
  """ 3.) Finding euler's angle requires a little bit of geometry. Please refer to
          M.E. Rose, Elementary Theory of Angular Momemtum, Wiley: New York 1957 for
          an explanations of the angles. Unfortunately there are a lot of different
          conventions, and I randomly picked the one mentioned above. """
  # Beta is measured between the new z axis and the old z axis. All angles are measured
  # counter-clockwise, so if the the y value of the new z axis is < 0 then we have to substract
  # the angle from 360 degrees.
  #line=cross_product([0,0,1],z_new)
  if (z_new[1] < 0):
    # line is a vector that specifies the intersection of the two x,y planes
    # alpha and gamma around found using this line
    line=cross_product(z_new,[0,0,1])
    beta=2*pi-acos(dot(z_new,[0,0,1])/(mag(z_new)*mag([0,0,1])))
  else:
    # line is a vector that specifies the intersection of the two x,y planes
    # alpha and gamma around found using this line
    line=cross_product([0,0,1],z_new)
    beta=acos(dot(z_new,[0,0,1])/(mag(z_new)*mag([0,0,1])))
  # Alpha is the angle between 'line' and the original y axis.
  if (line[0] > 0):
    alpha=2*pi-acos(dot(line,[0,1,0])/(mag(line)*mag([0,1,0])))
  else:
    alpha=acos(dot(line,[0,1,0])/(mag(line)*mag([0,1,0])))
  # Gamma is the angle between 'line' and the new y axis
  if ((beta < pi and y_new[2] < 0) or (beta > pi and y_new[2] > 0)):
    gamma=2*pi-acos(dot(line,y_new)/(mag(line)*mag(y_new)))
  else:
    gamma=acos(dot(line,y_new)/(mag(line)*mag(y_new)))

  """ 4.) Now all we need to do is find the rotation matrix. This can be
          found by multiplying each of the three rotation matrixes (one for each angle). """
  rot11=cos(alpha)*cos(beta)*cos(gamma)-sin(alpha)*sin(gamma)
  rot12=-cos(alpha)*cos(beta)*sin(gamma)-sin(alpha)*cos(gamma)
  rot13=cos(alpha)*sin(beta)
  rot21=sin(alpha)*cos(beta)*cos(gamma)+cos(alpha)*sin(gamma)
  rot22=-sin(alpha)*cos(beta)*sin(gamma)+cos(alpha)*cos(gamma)
  rot23=sin(alpha)*sin(beta)
  rot31=-sin(beta)*cos(gamma)
  rot32=sin(beta)*sin(gamma)
  rot33=cos(beta) 
  rot=r_[[[rot11,rot12,rot13]],
         [[rot21,rot22,rot23]],
         [[rot31,rot32,rot33]]]
  return rot

"""
"""
def new2old (target_new,rot,P):
  Pv=r_[P.x,P.y,P.z]
  """ Multiply the rotation matrix by the target point
      in the new coordinate system.
      ** Note to go the other direction you need to invert or
      transpose the rotation matrix. """
  target=dot(rot,target_new) # Matrix multiplication
  return target+Pv

"""
   Returns the position of a point (target) in the same plane as
   P, Q, and R with the angle Q-P-Target => angle. The algorithm
   find the plane described by the points P, Q, and R, and translates
   and rotates the coordinate system to line up with the plane. This makes
   it easy to place the target point in the plane.
"""
def findPlanarPosition (distance,angle,P,Q,R,check=False):
  rot=findRotationMatrix(P,Q,R)
  target_new=[distance*sin(radians(angle)),distance*cos(radians(angle)),0]
  target=new2old (target_new,rot,P)

  # I've included some geometric checks to make sure everything is working
  # correctly. If these tests pass then the algorithm should be working
  # correctly.
  if (check == True):
    Pv=r_[P.x,P.y,P.z]
    Qv=r_[Q.x,Q.y,Q.z]
    Rv=r_[R.x,R.y,R.z]    
    a=Qv-Pv
    b=Rv-Pv
    n=cross_product(a,b)
    n=n/mag(n)
    y_new=a/mag(a)
    z_new=n
    target_p=target-Pv
    angle_gen=acos(dot(target_p,y_new)/(mag(target_p)*mag(y_new)))*180/pi
    angle_ref=angle
    if (angle > 180): angle_ref=360-angle
    if (angle_gen-angle_ref > 1.0):
      print "ERROR: angle generated does not match specified angle: (gen)%f != (ref spec)%f"%(angle_gen,angle_ref)
    if (dot(target_p,z_new) > 0.0001): print "ERROR: target is not in the plane"

  return target

"""
"""
def findCirclePosition (distance,angle,circle_angle,P,Q,R,check=False):
  rot=findRotationMatrix(P,Q,R)
  # Find the position in the new coordinate system
  r=abs(distance*sin(radians(angle)))
  y=distance*cos(radians(angle))
  z=r*sin(radians(circle_angle))
  x=r*cos(radians(circle_angle))
  target_new=[x,y,z]
  target=new2old (target_new,rot,P)

  # I've included some geometric checks to make sure everything is working
  # correctly. If these tests pass then the algorithm should be working
  # correctly.
  if (check == True):
    Pv=r_[P.x,P.y,P.z]
    Qv=r_[Q.x,Q.y,Q.z]
    Rv=r_[R.x,R.y,R.z]    
    a=Qv-Pv
    b=Rv-Pv
    n=cross_product(a,b)
    n=n/mag(n)
    y_new=a/mag(a)
    z_new=n
    target_p=target-Pv
    angle_gen=acos(dot(target_p,y_new)/(mag(target_p)*mag(y_new)))*180/pi
    angle_ref=angle
    if (angle > 180): angle_ref=360-angle
    if (angle_gen-angle_ref > 1.0):
      print "ERROR: angle generated does not match specified angle: (gen)%f != (ref spec)%f"%(angle_gen,angle_ref)

  return target
