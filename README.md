# quat_m
Matlab package for Quaternion arithmetic  
Written by Dan Oates (WPI Class of 2020)

### Description
This package contains the Quat class for Quaternions. It implements the
following arithetic operations:

- Euclydian norm [norm(q)]
- Conjugation [conj(q)]
- Inversion [inv(q)]
- Normalization [unit(q)]
- Addition [q1 + q2]
- Subtraction [q1 - q2]
- Unary plus [+q]
- Unary minus [-q]
- Left division [q1 \ q2]
- Right division [q1 / q2]

It also has functions for rotation matrices, vector rotations, and axis-angle
decomposition for unit quaternions.

### Cloning and Submodules
Clone this repo as '+quat' and add the containing dir to the Matlab path.