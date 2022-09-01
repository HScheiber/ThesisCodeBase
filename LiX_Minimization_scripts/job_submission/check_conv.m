function stop = check_conv(~,optimValues,state,E_Unphys)

switch state
    case 'iter'
        if optimValues.fval < E_Unphys
            stop = true;
       else
            stop = false;
       end
    otherwise
        stop = false;
end

% disp(['Gradient: ' num2str(optimValues.gradient',' %.5f')])
% disp(['Lattice E: ' num2str(optimValues.fval,' %.5f')])
% disp(['Lattice params: ' num2str(x,' %.5f')])

end