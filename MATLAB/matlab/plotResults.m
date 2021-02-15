solverNames = ["Esmond", "NPT", "NPTHTS", "Pardiso"];
numRuns = 5;
for i=1:length(solverNames)
  a = strcat('load timings_report_', solverNames(i));
  eval(a);
  a = ['threadCount' num2str(i) ' = threadCount;'];
  eval(a);
  for j=1:numRuns
    a = ['ddTime = run' num2str(j) '_dat;'];
    eval(a);
    b = ['rsltTime = run' num2str(j) '_rslt;'];
    eval(b);
    if (j == 1)
      a = ['ddTimeBestSol' num2str(i) ' = ddTime;'];
      eval(a);
      b = ['rsltTimeBestSol' num2str(i) ' = rsltTime;'];
      eval(b)
    else
      a = ['ddTimeBestSol' num2str(i) ' = getMin(ddTimeBestSol' num2str(i) ', ddTime);'];
      eval(a);
      b = ['rsltTimeBestSol' num2str(i) ' = getMin(rsltTimeBestSol' num2str(i) ', rsltTime);'];
      eval(b);
    end
  end
end
colors = ['b' 'g' 'k' 'm'];
markers = ['o' 's' 'd' 'x'];
q = char(39);
plotCount = 1;
for i=1:length(MPICount)
  for k=1:2
    if (k == 1)
      dataName = 'ddTimeBestSol';
      titleName = 'dd\_solver.dat solve times for ';
    else
      dataName = 'rsltTimeBestSol';
      titleName = 'rslt file Elapsed Times for ';
    end
    a = ['figure(' num2str(plotCount) ')'];
    eval(a)
    a = 'plot(';
    b = 'h_legend = legend(';
    for j=1:length(solverNames)
      a = [a 'threadCount' num2str(j) ', ' dataName num2str(j) ', ' q colors(j) markers(j) '-' q ', '];
      b = strcat(b, q, solverNames(j), q);
      if (j == length(solverNames))
        a = [a q 'LineWidth' q ', 2, ' q 'MarkerSize' q ', 10)'];
        b = strcat(b, ');');
      else
        b = strcat(b, ', ');
      end
    end
    eval(a)
    eval(b)
    xlabel('number of threads', 'FontSize', 20);
    ylabel ('time (sec)', 'FontSize', 20);
    set(gca, 'fontsize', 20)
    a = [titleName num2str(MPICount(i)) ' MPI ranks'];
    title(a, 'FontSize', 14)
    plotCount = plotCount + 1;
  end  
end

function B = getMin(currentMin, A)
numRows = size(A,1);
numCols = size(A,2);
B = currentMin;
for i=1:numRows
  for j=1:numCols
    if (A(i,j) < currentMin(i,j))
      B(i,j) = A(i,j);
    end
  end
end
end
