% function testFailed = checkNodes(numberSubdomains)
%if ~testFailed, then gid==node_num_map 
function testFailed = checkNodes(numberSubdomains)
wrongValue = ones(numberSubdomains,1,'int32');
for subdomain =1:numberSubdomains,
     fromExodus = getSubdomainNodes(subdomain);
     fromSierra = getSubdomainNodes2(subdomain);
     rowBased = size(fromExodus,2) + size(fromSierra,2) - 2;
     if ~rowBased,
         wrongSize = size(fromExodus,1) - size(fromSierra,1);
         if ~wrongSize,
             wrongValue(subdomain) = max(abs(fromExodus - fromSierra ));
        end
    end
end
testFailed = max( abs( wrongValue ));
