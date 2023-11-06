function [Y]=function_ds(X)
    WW1=X;
    DWW1= (max(sum(X),eps)).^(-0.5); 
    WW1=  diag(DWW1)*WW1*diag(DWW1);
    Y=SK_normalize(WW1);
end