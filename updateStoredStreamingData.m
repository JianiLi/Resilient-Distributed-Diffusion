function [D, U] = updateStoredStreamingData( d,u, D, U)

D(1:(end-1),:) = D(2:end,:);
D(end,:) = d;

U(1:(end-1),:,:) = U(2:end,:,:);
U(end,:,:) = u;