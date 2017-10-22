function matlab_plots

    function main
        corbanFigureDefaults
        
        %plotFigure1
        plotFigure2
        saveAllFigures
    end

%% Global Variables
load('score_timecourses.mat');
load('score_timecourses_2.mat');

original_final_scores = original_score_timecourse(end,:);
relaxed_final_scores = relaxed_score_timecourse(end,:);

initial = -762.59235;

original_delta_scores = original_final_scores - initial;
relaxed_delta_scores = relaxed_final_scores - initial;

pick = 7;

main
%% Plots
    function plotFigure1
        setupFigure(1,"Changes In Score Compared");
        ax1 = subplot(1,2,1); hold on; box off;
        b = bar(original_delta_scores,...
            'LineStyle','None');
        xlabel('Decoy Number');
        ylabel('\DeltaScore (REU)');
        title('Beginning With Initial Complex');
        set(gca,'XTick',1:16);
        ax2 = subplot(1,2,2); hold on; box off;
        bar(relaxed_delta_scores,...
            'LineStyle','None');
        bar(pick,relaxed_delta_scores(pick),...
            'LineStyle','None');
        xlabel('Decoy Number');
        title('Beginning with Relaxed Complex');
        set(gca,'XTick',1:16);
        linkaxes([ax1, ax2],'xy')
    end


    function plotFigure2
        numsteps = length(relaxed_score_timecourse);
        setupFigure(2,"Minimization Timecourses");
        plot([0, numsteps],[0 0],'k-','LineWidth',1)
        p = plot(0:numsteps, [0, 0, 0; ...
            (relaxed_score_timecourse(:,[1, 7, 10]) - initial)],...
            '.-');
        xlim([0, numsteps])
        xlabel('Algorithm Step');
        set(gca,'XTick',1:7:numsteps);
        labels = strcat("Mutation ", num2str((1:6)'));
        labels = [labels; "Fast Relax"];
        set(gca,'XTickLabel',labels);
        ylabel('\DeltaScore (REU)');
        legend(p,'Decoy 1','Decoy 7','Decoy 10')
        grid on;
    end

end
