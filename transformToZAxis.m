function T = transformToZAxis(point, direction)
    % Normalize the direction vector
    direction = direction / norm(direction);
    
    % Create the z-axis unit vector
    zAxis = [0; 0; 1];
    
    % Calculate the rotation axis using cross product
    rotationAxis = cross(zAxis, direction);
    
    % Calculate the angle between the z-axis and the direction vector
    angle = acos(dot(zAxis, direction));
    
    % Create the skew-symmetric matrix for the rotation axis
    K = [0 -rotationAxis(3) rotationAxis(2);
         rotationAxis(3) 0 -rotationAxis(1);
         -rotationAxis(2) rotationAxis(1) 0];
    
    % Calculate the rotation matrix using Rodrigues' rotation formula
    R = eye(3) + sin(angle) * K + (1 - cos(angle)) * (K * K);
    
    % Create the translation matrix
    T_translate = eye(4);
    T_translate(1:3, 4) = point;
    
    % Create the rotation matrix in homogeneous coordinates
    T_rotate = eye(4);
    T_rotate(1:3, 1:3) = R;
   
    % Combine the translation and rotation matrices
    T = T_translate * T_rotate;
end