function conclusion=testMasterSlaveOrder(mpcOrder,mpcVector)
conclusion= 'No, the slave node is not the first node in each constraint equation';
numMpc = size(mpcOrder,1);
slaveLocation=cumsum([1;mpcOrder(2:numMpc)]);
residual=norm( mpcVector(slaveLocation,3) - ones(numMpc,1) , 1);
if residual == 0,
    conclusion= 'Yes, the slave node is the first node in each constraint equation';
end
