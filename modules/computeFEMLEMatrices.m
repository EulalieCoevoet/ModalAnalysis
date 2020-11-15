function [K,M] = computeFEMLEMatrices(V,T,young,nu,density)
    [volume,g1,g2,g3,g4] = computeVolumeAndGradient(V,T);
    K = computeStiffnessMatrix(V,T,young,nu,volume,g1,g2,g3,g4);
    M = computeMassMatrix(V,T,density,volume);
end
