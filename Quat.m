classdef Quat
    %QUAT Class for quaternion vectors
    %   
    %   Author: Dan Oates (WPI Class of 2020)
    
    properties
        w;  % w - component
        x;  % x - component
        y;  % y - component
        z;  % z - component
    end
    
    methods
        function q = Quat(varargin)
            %QUAT Construct quaternion vector
            %   q = QUAT(w, x, y, z) Element-based form
            %   q = QUAT([ax; ay; az], t) Axis-angle form
            %       [ax; ay; az] = Axis of rotation
            %       t = Rotation angle [rad]
            %   q = QUAT([w; x; y; z]) Vector form
            %   q = QUAT() Union quaternion [1, 0, 0, 0]
            if nargin == 4
                % Element-based form
                q.w = varargin{1};
                q.x = varargin{2};
                q.y = varargin{3};
                q.z = varargin{4};
            elseif nargin == 2
                % Axis-angle form
                ax = varargin{1};
                ax = ax / norm(ax);
                t_2 = 0.5 * varargin{2};
                ct_2 = cos(t_2);
                st_2 = sin(t_2);
                q.w = ct_2;
                q.x = ax(1) * st_2;
                q.y = ax(2) * st_2;
                q.z = ax(3) * st_2;
            elseif nargin == 1
                % Vector form
                v = varargin{1};
                q.w = v(1);
                q.x = v(2);
                q.y = v(3);
                q.z = v(4);
            elseif nargin  == 0
                % Union quaternion
                q.w = 1;
                q.x = 0;
                q.y = 0;
                q.z = 0;
            else
                error('Invalid arguments.')
            end
        end
        
        function n = norm(q)
            %n = NORM(q) Euclydian magnitude
            n = norm(vector(q));
        end
        
        function q = conj(q)
            %q = CONJ(q) Conjugation [w, -x, -y, -z]
            import('quat.Quat');
            q = Quat(q.w, -q.x, -q.y, -q.z);
        end
        
        function q = inv(q)
            %q = INV(q) Inverse [conj(q) / norm(q)^2]
            import('quat.Quat');
            s = 1 / norm(q)^2;
            q = Quat(q.w*s, -q.x*s, -q.y*s, -q.z*s);
        end
        
        function q = unit(q)
            %q = UNIT(q) Unit quaternion of q
            import('quat.Quat');
            s = 1 / norm(q);
            q = Quat(q.w*s, q.x*s, q.y*s, q.z*s);
        end
        
        function q = pos_w(q)
            %q = POS_W(q) Quaternion with positive w-component
            if q.w < 0
                q = -q;
            end
        end
        
        function [ax, t] = axis(q)
            %[ax, t] = AXIS(q) Axis-angle form of quaternion
            %   ax = Unit vector axis of rotation
            %   t = Angle of rotation [rad]
            q.w = min(max(-1, q.w), +1);
            t_2 = acos(q.w);
            t = 2 * t_2;
            ax_s = 1 / sin(t_2);
            ax = zeros(3, 1);
            ax(1) = q.x * ax_s;
            ax(2) = q.y * ax_s;
            ax(3) = q.z * ax_s;
        end
        
        function v = rotate(q, v)
            %v = ROTATE(q, v) Rotates vector v by unit quaternion q
            v = mat_rot(q) * v;
        end
        
        function v = vector(q)
            %v = VECTOR(obj) Vector form [w; x; y; z]
            v = [q.w; q.x; q.y; q.z];
        end
        
        function M = mat_int(q)
            %M = MAT_INT(q) Intrinsic matrix [q'*q = M*q']
            M = [...
                [+q.w, -q.x, -q.y, -q.z]; ...
                [+q.x, +q.w, +q.z, -q.y]; ...
                [+q.y, -q.z, +q.w, +q.x]; ...
                [+q.z, +q.y, -q.x, +q.w]];
        end
        
        function M = mat_ext(q)
            %M = MAT_EXT(q) Extrinsic matrix [q*q' = M*q']
            M = [...
                [+q.w, -q.x, -q.y, -q.z]; ...
                [+q.x, +q.w, -q.z, +q.y]; ...
                [+q.y, +q.z, +q.w, -q.x]; ...
                [+q.z, -q.y, +q.x, +q.w]];
        end
        
        function R = mat_rot(q)
            %R = MAT_ROT(q) Rotation matrix of unit quaternion [q*v*q^=1 = R*v]
            xx = 2 * q.x * q.x;
            yy = 2 * q.y * q.y;
            zz = 2 * q.z * q.z;
            wx = 2 * q.w * q.x;
            wy = 2 * q.w * q.y;
            wz = 2 * q.w * q.z;
            xy = 2 * q.x * q.y;
            xz = 2 * q.x * q.z;
            yz = 2 * q.y * q.z;
            R = [...
                [1 - yy - zz, xy - wz, xz + wy]; ...
                [xy + wz, 1 - xx - zz, yz - wx]; ...
                [xz - wy, yz + wx, 1 - xx - yy]];
        end
        
        function q = mtimes(q1, q2)
            %q = MTIMES(q1, q2) Quaternion multiplication
            import('quat.Quat');
            w_ = q1.w*q2.w - q1.x*q2.x - q1.y*q2.y - q1.z*q2.z;
            x_ = q1.w*q2.x + q1.x*q2.w + q1.y*q2.z - q1.z*q2.y;
            y_ = q1.w*q2.y - q1.x*q2.z + q1.y*q2.w + q1.z*q2.x;
            z_ = q1.w*q2.z + q1.x*q2.y - q1.y*q2.x + q1.z*q2.w;
            q = Quat(w_, x_, y_, z_);
        end
        
        function q = mrdivide(q1, q2)
            %q = MRDIVIDE(q1, q2) Returns q such that q * q2 = q1
            q = q1 * inv(q2);
        end
        
        function q = mldivide(q1, q2)
            %q = MLDIVIDE(q1, q2) Returns q such that q1 * q = q2
            q = inv(q1) * q2;
        end
        
        function q = plus(q1, q2)
            %q = PLUS(q1, q2) Quaternion addition
            import('quat.Quat');
            w_ = q1.w + q2.w;
            x_ = q1.x + q2.x;
            y_ = q1.y + q2.y;
            z_ = q1.z + q2.z;
            q = Quat(w_, x_, y_, z_);
        end
        
        function q = minus(q1, q2)
            %q = MINUS(q1, q2) Quaternion subtraction
            import('quat.Quat');
            w_ = q1.w - q2.w;
            x_ = q1.x - q2.x;
            y_ = q1.y - q2.y;
            z_ = q1.z - q2.z;
            q = Quat(w_, x_, y_, z_);
        end
        
        function q = uminus(q)
            %q = UMINUS(q) Unary minus [-w, -x, -y, -z]
            import('quat.Quat');
            q = Quat(-q.w, -q.x, -q.y, -q.z);
        end
        
        function q = uplus(q)
            %q = UPLUS(q) Unary plus [+w, +x, +y, +z]
        end
    end
end