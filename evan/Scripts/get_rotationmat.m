function Rmat = get_rotationmat(point1,point2)
%Rmat = get_rotationmat(point1,point2)
%
%Calculate rotation matrix Rmat to go from point1 to point2, such that
%point1 * Rmat = point2

radius = sqrt(point1(1)^2 + point1(2)^2 + point1(3)^2);

point1 = point1 / radius;
point2 = point2 / radius;

rotangle = acos(dot(point2,point1));
rotaxis = cross(point2,point1);

Rmat = rotationmat3D(rotangle,rotaxis);
    