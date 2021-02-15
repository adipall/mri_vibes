% function boolean = safeNumberCoefficients( n, nmax, ncoeffs )
function boolean = safeNumberCoefficients( n, nmax, ncoeffs )
if ( nmax >= max(n) ) && (nmax <= ncoeffs )
  boolean = true;
else
  boolean = false;
end
