%function value=series(angles,coefficient,indices,trigonometry)
function value=series(angles,coefficient,indices,trigonometry)
switch lower(trigonometry)
    case {'cosine'}
        matrix = cos(angles*indices');
    case {'sine'}
        matrix = sin(angles*indices');
    otherwise
        disp('Unknown method.');
        value = [];
        return;
end
value = matrix*coefficient;
