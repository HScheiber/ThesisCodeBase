function output = unique_plot(app)

    Unique_plot = (app.data_type == 0);
    if isempty(Unique_plot)
        output = [];
        return
    end
    
    if ~isnan(app.energy_type_plot)
        Unique_plot( find(app.data_type == 1, 1, 'first') ) = true;
    end
    if ~isnan(app.rdf_type_plot)
        Unique_plot( find(app.data_type == 2, 1, 'first') ) = true;
    end
    if ~isnan(app.msd_type_plot)
        Unique_plot( find(app.data_type == 3, 1, 'first') ) = true;
    end
    if ~isnan(app.AvQl_type_plot)
        Unique_plot( find(app.data_type == 4, 1, 'first') ) = true;
    end
    if ~isnan(app.NN_type_plot)
        Unique_plot( find(app.data_type == 5, 1, 'first') ) = true;
    end
    output = Unique_plot;
end



