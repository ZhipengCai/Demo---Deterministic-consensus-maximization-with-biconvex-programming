function theta = enforceSO3Constraint(theta0)
    R0 = vec2rotMat(theta0);
    R = closestRotMat(R0);
    theta = rotMat2Vec(R);
end