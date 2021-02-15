classdef MeshTest < matlab.unittest.TestCase

    methods (Test)
        function testMesh(testCase)
            actualSolution = twoBlocks();
            expectedSolution = 1;
            testCase.verifyEqual(actualSolution,expectedSolution);
        end
        function testTriangle(testCase)
            actualSolution = oneBlock();
            expectedSolution = 1;
            testCase.verifyEqual(actualSolution,expectedSolution);
        end
    end
end
%%
function state = twoBlocks()
    [coord, idx, disp, connectivity] = FESearch.Example.setTwoBlockMesh();
    elementType = 'Qua4';
    fexo = FEMesh.Exodus( elementType, connectivity(1:12), coord(idx,:) );
    fexo.Blocks(2).ID=2;
    fexo.Blocks(2).Connectivity = connectivity(13:24);
    fexo.Blocks(2).ElementType = elementType;
    state=1;
end
function state = oneBlock()
    % cd code/Salinas/tools/matlab
    obj = FEMesh.icosahedron();
    elementType = 'Tri3';
    fexo = FEMesh.Exodus( elementType, obj.connectivity, obj.coordinates);
    state=1;
end

