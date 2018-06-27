function Exp = Expectation(n, numAgents, Exp, D, U, state, storedNum)

for k = 1:numAgents
    Exp(n,k) = sum((D(:,k) - U(:,:,k) * state(:,k)).^2)/storedNum;
end
    
    
