clear; close all; clc;
Cost = @(X) rastrigin(X);

alphaVal = 0.3;
xmin = -5.2; xmax = 5.2;
ymin = -5.2; ymax = 5.2;

figure
hold on
limits = repmat([xmin xmax], 2, 1);
[X,Y] = meshgrid(linspace(limits(1,1),limits(1,2),100),...
                   linspace(limits(2,1),limits(2,2),100));
Z = rastrigin([X(:)'; Y(:)']);
Z = reshape(Z,size(X));
surf(X,Y,Z,'linestyle','none','FaceAlpha',1);
colormap("parula")
%title({'Rastrigin''s Function','-5.2 < x < 5.2 and -5.2 < y < 5.2','The minimum of the function is z = 0 at (x, y) = (0, 0).'})
xlabel('x')
ylabel('y')
zlabel('z')
grid on

az=120;
el=20;
view(az,el);
c = colorbar;

savefig('figs\rastrigin.fig');
print(gcf,'figs\rastrigin','-depsc','-r300')
print(gcf,'figs\rastrigin','-dpng','-r300')

%%
figure
set(gcf, 'WindowState', 'maximized');
% Since Z is a matrix, max(Z) is a vector, so we need max(max(Z)) to maximize Z.
MaxValue = max(max(Z));
DeltaV = 50;
V1 = [0:DeltaV:4*DeltaV];
V2 = [10*DeltaV:10*DeltaV:MaxValue/2];
V = [V1 V2];

ax = gca;
[M, c] = contour(X, Y, Z, 'ShowText', 'off', 'HandleVisibility', 'off');
contourcmap("parula",5,"Colorbar","on","Location","horizontal", ...
    "TitleString","Fitness values for Rastrigin''s function")
hold on
plot(0,0,'m.','markersize',24)
text(0.05,0.05,'x*','Color','black','FontSize',14,'FontWeight','bold')
grid off
%title({'Contour plot of Rastrigin''s function'})
xlabel('x')
ylabel('y')

ticks=[-3:1:3];
set(c,'linewidth',1);

%% Nelder Mead
% Coefficients
global calpha cbeta cgamma1 cgamma2 cdelta iter BestCost VarMin VarMax op
calpha = 0; cbeta = 0; cgamma1 = 0; cgamma2 = 0; cdelta = 0;


Simplex = [-2.5 2.5; -2.2 0; 1.1 -1.2]; % Initial Simplex
X = Simplex(:,1); Y = Simplex(:,2); 
for i=1:size(Simplex,1)
    F(i) = rastrigin(Simplex(i,:));
end
[~,sortInd] = sort(F);

Simplex = Simplex(sortInd,:);
F = F(sortInd);
B = Simplex(1,:);
G = Simplex(2,:);
W = Simplex(3,:);
X = Simplex(:,1); Y = Simplex(:,2); 
fill(X,Y,'white','FaceAlpha',0.1,'LineWidth',1.25)
legend('Global minima', 'Initial Simplex','Location','northoutside','Orientation','horizontal');

savefig('figs\rastrigin_contour.fig');
print(gcf,'figs\rastrigin_contour','-depsc','-r300')
print(gcf,'figs\rastrigin_contour','-dpng','-r300')

fill(NaN,NaN,'blue','FaceAlpha',alphaVal,'LineWidth',1.25, 'HandleVisibility', 'on', 'DisplayName','Expansion');
fill(NaN,NaN,'magenta','FaceAlpha',alphaVal,'LineWidth',1.25, 'HandleVisibility', 'on', 'DisplayName','Reflection');
fill(NaN,NaN,'yellow','FaceAlpha',alphaVal,'LineWidth',1.25, 'HandleVisibility', 'on', 'DisplayName','Outside Contraction');
fill(NaN,NaN,'cyan','FaceAlpha',alphaVal,'LineWidth',1.25, 'HandleVisibility', 'on', 'DisplayName','Inside Contraction');
fill(NaN,NaN,'red','FaceAlpha',alphaVal,'LineWidth',1.25, 'HandleVisibility', 'on', 'DisplayName','Shrink');
lgd = legend;
lgd.NumColumns = 4;
disp(' k Best            Good            Worst          fbest     Operation         ');
disp('-----------------------------------------------------------------------------------')
fprintf('%2d (%5.3f, %5.3f) (%5.3f, %5.3f) (%5.3f, %5.3f) %6.3f %s\n', 0, B(1), B(2), ...
    G(1), G(2), W(1), W(2), F(1), 'Initialize Simplex');
% Second Simplex
maksIter = 100;
for iter=1:maksIter
    [Simplex,F,clrCode] = nelderMead(Cost,Simplex,F);
    X = Simplex(:,1); Y = Simplex(:,2);
    if (iter<=5) || (mod(iter,10)==0)
       fprintf('%2d (%5.3f, %5.3f) (%5.3f, %5.3f) (%5.3f, %5.3f) %6.3f %s\n', iter, Simplex(1,1), Simplex(1,2), ...
       Simplex(2,1), Simplex(2,2), Simplex(3,1), Simplex(3,2), F(1), op);
       fill(X,Y,clrCode,'FaceAlpha',alphaVal,'LineWidth',1.25,'HandleVisibility', 'off');
       savefig(['figs\nms\rastrigin_nms_' num2str(iter) '.fig']);
       print(gcf,['figs\nms\rastrigin_nms_' num2str(iter)],'-depsc','-r300')
       print(gcf,['figs\nms\rastrigin_nms_' num2str(iter)],'-dpng','-r300')
    end
end
disp('-----------------------------------------------------------------------------------')
plot(Simplex(1,1),Simplex(1,2),'bs','MarkerSize',12,'HandleVisibility', 'on','DisplayName','NMS Solution','MarkerFaceColor','green')
plot(0,0,'m.','markersize',24,'HandleVisibility', 'off');
hold off
savefig(['figs\nms\rastrigin_nms_' num2str(maksIter) '.fig']);
print(gcf,['figs\nms\rastrigin_nms_' num2str(maksIter)],'-depsc','-r300')
print(gcf,['figs\nms\rastrigin_nms_' num2str(maksIter)],'-dpng','-r300')