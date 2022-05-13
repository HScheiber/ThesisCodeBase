function PlotState = Update_PlotState_Version(RefState,PlotState)

    N = length(PlotState);

    if isfield(PlotState,'ViewTrajectoryh')
        % Set SliceTrajectoryh to ViewTrajectoryh
        [PlotState(1:N).SliceTrajectoryh] = deal(PlotState(1:N).ViewTrajectoryh);
        
        % Remove ViewTrajectoryh field
        PlotState = rmfield(PlotState,'ViewTrajectoryh');
    end

    if ~isfield(PlotState,'compareQlButton1')       
        for idx = 1:N
            PlotState(idx).compareQlButton1 = RefState.compareQlButton1;
        end
    end

    if ~isfield(PlotState,'compareQlButton2')
        for idx = 1:N
            PlotState(idx).compareQlButton2 = RefState.compareQlButton2;
        end
    end

    if ~isfield(PlotState,'Modelfilter_b1')
        for idx = 1:N
            PlotState(idx).Modelfilter_b1.Value = false;
            PlotState(idx).Modelfilter_b2.Value = true;
        end
    end
end