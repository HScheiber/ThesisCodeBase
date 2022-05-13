% Script for plotting energy components of BO models compared with others

%% Analysis Parameters
fs = 24;
Salts = {'LiF' 'LiCl' 'LiBr' 'LiI'};
Theory = 'TF';
Basenum = 'C';
Midnum = 'G';
reps = {'1' '2' '3' '4' '5'}; %  {'1' '2' '3' '4' '5'}{'a' 'b' 'c' 'd' 'e'}
Structure = 'Rocksalt';
Plot_diff = false; % If true, plot energy differences w.r.t. reference theory

% Which comparisons to plot
plot_Pot = true;
plot_CoulSR = true;
plot_CoulLR = true;
plot_LJSR = true;
plot_CoulSR_MM = false;
plot_CoulSR_MX = false;
plot_CoulSR_XX = false;
plot_LJSR_MM = true;
plot_LJSR_MX = true;
plot_LJSR_XX = true;

plot_switches = [plot_CoulSR plot_CoulLR plot_LJSR plot_Pot ...
    plot_CoulSR_MM plot_CoulSR_MX plot_CoulSR_XX ...
    plot_LJSR_MM plot_LJSR_MX plot_LJSR_XX];

Models = cell(length(reps),1);
for idx = 1:length(reps)
    Models{idx} = [Basenum Midnum reps{idx}];
end

filename = [Theory '_Models_' Models{1} '-' Models{end} '_Energy_Breakdown.png'];

if Plot_diff
    N_Models = length(Models);
    Colours = cbrewer('qual','Dark2',N_Models);
    Model_Labels = Models;
else
    N_Models = length(Models)+1;
    Colours = [cbrewer('qual','Dark2',N_Models-1); 0 0 0];
    Model_Labels = [Models; Theory];
end
N_Salts = length(Salts);
N_Plots = sum(plot_switches);


%% Other parameters
savefile = false; % switch to save the final plots to file

% Find model in results
ML_results_dir = 'C:\Users\Hayden\Documents\Patey_Lab\BO_Models';

ydat = zeros(N_Models,N_Plots,N_Salts);
ticklabels = cell(1,N_Plots);
for idx = 1:length(Salts)
    Salt = Salts{idx};
    [Metal,Halide] = Separate_Metal_Halide(Salt);
    
    % Load reference model
    sres = dir([ML_results_dir filesep Salt '_' Theory '_reference.mat']);
    % Load full opt data
    model_filename = fullfile(sres.folder,sres.name);
    data = load(model_filename);
    Minimization_Data = data.Minimization_Data;
    Structures = cell(1,length(Minimization_Data));
    for qq = 1:length(Minimization_Data)
        Structures{qq} = Minimization_Data{qq}.Structure;
    end
    Minimization_Data_ref = Minimization_Data{strcmpi(Structure,Structures)};
    
    for jdx = 1:N_Models
        
        % Find the final optimized model data
        if jdx == N_Models && ~Plot_diff
            Minimization_Data_Choice = Minimization_Data_ref;
        else
            Model = Models{jdx};
            sres = dir([ML_results_dir filesep '*' Salt '*' Theory '*' 'Model_' Model '_*fullopt.mat']);

            if length(sres) > 1
                warning(['Multiple results found for: ' Salt ', ' Theory ', Model ' Model '. Using first result.']);
                sres = sres(1);
            elseif isempty(sres)
                error(['No results found for: ' Salt ', ' Theory ', Model ' Model '.']);
            end

            % Load full opt data
            model_filename = fullfile(sres.folder,sres.name);
            data = load(model_filename);
            Minimization_Data = data.Minimization_Data;
            Structures = cell(1,length(Minimization_Data));
            for qq = 1:length(Minimization_Data)
                Structures{qq} = Minimization_Data{qq}.Structure;
            end
            Minimization_Data_Choice = Minimization_Data{strcmpi(Structure,Structures)};
        end
        
        pdx = 1;
        if plot_Pot
            if Plot_diff
                ydat(jdx,pdx,idx) = Minimization_Data_Choice.Pot - Minimization_Data_ref.Pot;
            else
                ydat(jdx,pdx,idx) = Minimization_Data_Choice.Pot;
            end
            ticklabels{pdx} = 'Total Potential';
            pdx = pdx+1;
        end
        
        if plot_CoulSR
            if Plot_diff
                ydat(jdx,pdx,idx) = Minimization_Data_Choice.CoulSR - Minimization_Data_ref.CoulSR;
            else
                ydat(jdx,pdx,idx) = Minimization_Data_Choice.CoulSR;
            end
            ticklabels{pdx} = 'Total Coulomb (S.R.)';
            pdx = pdx+1;
        end
        
        if plot_CoulLR
            if Plot_diff
                ydat(jdx,pdx,idx) = Minimization_Data_Choice.CoulLR - Minimization_Data_ref.CoulLR;
            else
                ydat(jdx,pdx,idx) = Minimization_Data_Choice.CoulLR;
            end
            ticklabels{pdx} = 'Total Coulomb (L.R.)';
            pdx = pdx+1;
        end

        if plot_LJSR
            if Plot_diff
                ydat(jdx,pdx,idx) = Minimization_Data_Choice.LJSR - Minimization_Data_ref.LJSR;
            else
                ydat(jdx,pdx,idx) = Minimization_Data_Choice.LJSR;
            end
            ticklabels{pdx} = 'Total vdW';
            pdx = pdx+1;
        end
        
        if plot_CoulSR_MM
            if Plot_diff
                ydat(jdx,pdx,idx) = Minimization_Data_Choice.CoulSR_MM - Minimization_Data_ref.CoulSR_MM;
            else
                ydat(jdx,pdx,idx) = Minimization_Data_Choice.CoulSR_MM;
            end
            ticklabels{pdx} = ['Coulomb: ' Metal '-' Metal ' (S.R.)'];
            pdx = pdx+1;
        end
        
        if plot_CoulSR_MX
            if Plot_diff
                ydat(jdx,pdx,idx) = Minimization_Data_Choice.CoulSR_MX - Minimization_Data_ref.CoulSR_MX;
            else
                ydat(jdx,pdx,idx) = Minimization_Data_Choice.CoulSR_MX;
            end
            ticklabels{pdx} = ['Coulomb: ' Metal '-' Halide ' (S.R.)'];
            pdx = pdx+1;
        end
        
        if plot_CoulSR_XX
            if Plot_diff
                ydat(jdx,pdx,idx) = Minimization_Data_Choice.CoulSR_XX - Minimization_Data_ref.CoulSR_XX;
            else
                ydat(jdx,pdx,idx) = Minimization_Data_Choice.CoulSR_XX;
            end
            ticklabels{pdx} = ['Coulomb: ' Halide '-' Halide ' (S.R.)'];
            pdx = pdx+1;
        end
        
        if plot_LJSR_MM
            if Plot_diff
                ydat(jdx,pdx,idx) = Minimization_Data_Choice.LJSR_MM - Minimization_Data_ref.LJSR_MM;
            else
                ydat(jdx,pdx,idx) = Minimization_Data_Choice.LJSR_MM;
            end
            ticklabels{pdx} = ['vdW: ' Metal '-' Metal];
            pdx = pdx+1;
        end
        
        if plot_LJSR_MX
            if Plot_diff
                ydat(jdx,pdx,idx) = Minimization_Data_Choice.LJSR_MX - Minimization_Data_ref.LJSR_MX;
            else
                ydat(jdx,pdx,idx) = Minimization_Data_Choice.LJSR_MX;
            end
            ticklabels{pdx} = ['vdW: ' Metal '-' Halide];
            pdx = pdx+1;
        end
        
        if plot_LJSR_XX
            if Plot_diff
                ydat(jdx,pdx,idx) = Minimization_Data_Choice.LJSR_XX - Minimization_Data_ref.LJSR_XX;
            else
                ydat(jdx,pdx,idx) = Minimization_Data_Choice.LJSR_XX;
            end
            ticklabels{pdx} = ['vdW: ' Halide '-' Halide];
            pdx = pdx+1;
        end
    end
end

% Create figure and axis
figh = figure('WindowState','maximized','NumberTitle','off',...
    'Name','','Visible','On');
t = tiledlayout(figh,N_Salts,1,'TileSpacing','tight');        

for idx = 1:N_Salts
    yy = ydat(:,:,idx); % Select the plots for each Salt
    axobj = nexttile(t);
    hold(axobj,'on')
    
    p = bar(axobj,1:N_Plots,yy,'FaceColor','flat','Visible','on','BarWidth',1);
    
    % Apply colours
    for jdx = 1:N_Models
        p(jdx).CData = Colours(jdx,:);
        p(jdx).LineWidth = 2;
    end
    
    % Add plot title and y axis
    set(axobj,'box','on','TickLabelInterpreter','latex');
    set(axobj,'XMinorTick','Off','YMinorTick','on','FontSize',fs);
    grid(axobj,'minor');
    ylim(axobj,'padded')
    set(axobj,'XLim',[0.5 N_Plots+0.5])
    xticks(axobj,1:N_Plots)
    axobj.XAxis.TickLength = [0,0];
    xticklabels(axobj,[])
    title(axobj, Salts{idx},'Fontsize',fs,'Interpreter','latex')
end
if Plot_diff
    ylabel(t, ['$E_{\mathrm{Model}}$ - $E_{\mathrm{' Theory '}}$ [kJ mol$^{-1}$]'],'Fontsize',fs,'Interpreter','latex');
else
    ylabel(t, 'Energy [kJ mol$^{-1}$]','Fontsize',fs,'Interpreter','latex');
end
xticklabels(axobj,ticklabels)

legh = legend(axobj,p,Model_Labels,'Position',[0.9347 0.4005 0.0622 0.2442],'Orientation','Vertical',...
    'Interpreter','latex','Box','off','fontsize',fs,'NumColumns', 1);

title(legh,['\bf{\underline{Models}}'],'Fontsize',fs,'Interpreter','latex')


if savefile
    set(figh,'renderer','opengl')
    exportgraphics(figh,filename,'Resolution',300)
end