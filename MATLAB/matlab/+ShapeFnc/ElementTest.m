classdef ElementTest < matlab.unittest.TestCase

    methods (Test)
        function testQuad(testCase)
            actualSolution = twoBlocks();
            expectedSolution = 1;
            testCase.verifyEqual(actualSolution,expectedSolution);
        end
    end
end
%%
function state = twoBlocks()
    [coord, idx, disp, connectivity] = FESearch.Example.setTwoBlockMesh();
    obj=ShapeFnc.Qua4( connectivity, coord );
    state=1;
end
