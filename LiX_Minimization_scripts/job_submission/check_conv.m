function stop = check_conv(~,optimValues,~,Gradient_Tol_RMS,Gradient_Tol_Max)

if ( rms(optimValues.gradient) < Gradient_Tol_RMS ) && max(abs(optimValues.gradient)) < Gradient_Tol_Max
    stop = true;
else
    stop = false;
end

% disp(['Gradient: ' num2str(optimValues.gradient',' %.5f')])
% disp(['Lattice E: ' num2str(optimValues.fval,' %.5f')])
% disp(['Lattice params: ' num2str(x,' %.5f')])

end