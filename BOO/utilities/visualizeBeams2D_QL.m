function visualizeBeams2D_QL(selectedfpangles,beamlog)
load('4pi_angles.mat')

activeBeams = [];
for ii = 1:size(selectedfpangles,1)
    [minval,ind] = min(sum(abs(selectedfpangles(ii,:)-angles),2));
    activeBeams(ii) = ind;
end

markerSize = 50;
pos_leg = [0.3 0.95 0.4 0.03];
fontsize = 15;

figure();
set(gcf,'units','normalized','outerposition',[0.2 0.2 0.4 0.6])
beam_infeasible = find(beamlog==0);
beam_feasible = find(beamlog==1);
scatter(theta_Varian(beam_infeasible, 2), theta_Varian(beam_infeasible, 1),markerSize,  'kx', 'linewidth', 1.5)
hold on
scatter(theta_Varian(beam_feasible, 2), theta_Varian(beam_feasible, 1),markerSize, 'bo', 'linewidth', 1)
hold on
if(nnz(activeBeams)~=0)
scatter(theta_Varian(activeBeams,2),theta_Varian(activeBeams,1),markerSize,'ro','filled')
end

hold off
axis([85 275 -5 360])
x_tick_vec = 90:15:270;
y_tick_vec = 0:30:360;
set(gca, 'XDir','reverse')
set(gca, 'XTick', x_tick_vec, 'YTick', y_tick_vec, 'YTickLabel', rem(540-y_tick_vec, 360), 'XTickLabel', rem(540-flip(x_tick_vec), 360), 'fontsize', fontsize)
ylabel('Gantry (VarianIEC)', 'fontsize', fontsize)
xlabel('Couch (VarianIEC)', 'fontsize', fontsize)

hl= legend('infeasible', 'feasible', 'selected beam','orientation', 'horizontal');
set(hl, 'position', pos_leg, 'units', 'normalized', 'fontsize', 15, 'fontweight', 'bold')

ch = get(legend, 'Children');
textCh = ch(strcmp(get(ch, 'Type'), 'text'));
for iText = 1:numel(textCh) 
    set(textCh(iText), 'Position', get(textCh(iText), 'Position') + [-0.015 0 0])
end
legend boxoff



