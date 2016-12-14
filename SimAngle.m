function U=SimAngle(I,R,angle,typ)
    % from the code of Nico Vervliet and Lieven De Lathauwer
    % (fixedAngleVect in tensorlab 3.0)
    v = zeros(I, R);
    for k = 1:R
        v(k, k) = sqrt(1 - sum(v(1:k-1,k).^2));
        v(k, k+1:end) = (cos(angle) - sum(v(1:k-1,k).^2))./ v(k, k);
    end
    v(R, R) = sqrt(1 - sum(v(1:R-1,R).^2));
    rot = orth(typ(I));
    v = rot*v;
    U=v;
    clear v k rot;
    
    