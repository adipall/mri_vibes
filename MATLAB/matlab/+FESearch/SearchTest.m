classdef SearchTest < matlab.unittest.TestCase

    methods (Test)
        function testSearch(testCase)
            actualSolution = twoBlocks();
            expectedSolution = 1;
            testCase.verifyEqual(actualSolution,expectedSolution);
        end
    end
end
%%
function state = twoBlocks()
    [nodexyz, idx, disp, ms] = FESearch.Example.setTwoBlockMesh();
    obj=FESearch.Search;
    inp=obj.ParticleSearch(nodexyz(idx,:),ms,nodexyz(idx,:),max(abs(disp(idx,:)))); % disp could be dt*v
    % For each surfaceElement, inp{surfaceElement,:} is the nodes closest to ms(surfaceElement,:)

    numElement = size(ms,1);
    for element=1:numElement,
        neighbor = setdiff(  inp{element,:},  ms(element,:));
        assert( size( inp{element,:},1 ) - size( neighbor,1 ) == 4 );
    end
    state=1;
end
