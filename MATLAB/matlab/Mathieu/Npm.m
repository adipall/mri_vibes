% function normalizationFactor= Npm(category,coefficients)
function normalizationFactor= Npm(category,coefficients)
normalizationFactor= pi * sum( coefficients.* coefficients )';
if category == 1  %-even-even
    normalizationFactor = normalizationFactor + (coefficients(1,:).* coefficients(1,:))';
end