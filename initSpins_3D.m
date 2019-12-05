function spin = initSpins(numSpins_xDim, numSpins_yDim, numSpins_zDim, p)
%INITSPINS Initialize a configuration of spins.
%   spin = INITSPINS(numSpinsPerDim, p) returns a configuration of spins
%   with |numSpinsPerDim| spins along each dimension and a proportion |p|
%   of them pointing upwards. |spin| is a matrix of +/- 1's.

%   Copyright 2017 The MathWorks, Inc.

spin = sign(p - rand(numSpins_xDim, numSpins_yDim, numSpins_zDim));
