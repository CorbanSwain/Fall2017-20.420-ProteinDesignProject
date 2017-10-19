corbanFigureDefaults

initial = -762.59235;

delta_scores = scores - initial;

figure(1); clf; box off;
bar(delta_scores)
xlabel('Decoy Number')
ylabel('\DeltaScore (REU)')
set(gca,'XTick',1:16)
saveAllFigures();