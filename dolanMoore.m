%Generates figures for optimization performance profiling
% 
% Syntax:
%   dolanMoore(T)
%   dolanMoore(T, algs, funId, type)
%   dolanMoore(T, algs, funId, type, varargin)
%   [scaledT, stairsT] = dolanMoore(T, algs, funId, type, varargin)
% 
% Description:
%
%   dolanMoore(T) calls the performance profiling algorithm presented by
%   Elizabeth D. Dolan and Jorge J. More.
%   It takes the performance measure as an input and generate the
%   associated figures.
%
%   [scaledT, stairsT] = oppa(T, algs, settings) also enables usage of
%   pre-defined algorithm names and also pre-defined plot properties.
%   Function returns the scaled values of the measures and their
%   corresponding values in the performance profile.
%   
% Input Arguments
%   T         - Performance measure matrix. Rows of T should represent
%               experiments and columns represent algorithms/software to be
%               compared.
%   algs      - Names of the algorithms for graphic legend
%   varargin  - Plot(stairs) function properties  
%
% Examples:
%   
% Output Arguments
%   
%   See also STAIRS.
% 
% Sertalp B. Cay, Pelin Cay 2014
% Modified by Korhan Gunel, 2022


function [scaledX, stairsT] = dolanMoore(T, algs, funId, type, varargin)
    % Control algs
    n = size(T,2); % # of algorithms
    legendX = {};
    if(exist('algs','var'))
        %for i=1:n
            legendX = algs;
        %end
    else
        for i=1:n
            legendX{i} = strcat('Algorithm', int2str(i));
        end
    end
    % Eliminate entries with zero
    T(T(:,:)==0)=inf;
    % Scale every row by dividing its minimum
    scaledX = bsxfun(@rdivide,T,min(T')');
    % Select a border
    upperborder = ceil(max(scaledX(scaledX(:,:)<Inf)));
    if upperborder <2
        scaledX = scaledX + 1;
        upperborder = ceil(max(scaledX(scaledX(:,:)<Inf)));
    end
    if upperborder > 10^10
        scaledX = 1+scaledX./max(max(scaledX));
        upperborder = ceil(max(scaledX(scaledX(:,:)<Inf)));
    end
    % Sort the columns
    stairsT = sort(scaledX);
    % Value array
    stairsT = [ones(1,n); stairsT; upperborder*ones(1,n)];
    % Generate the graphic
    hold all
    linestyles = {'--',':','-.','-'};
    markers = {'s','o','+','*','x','d','v','^','<','>','p','h','.','+','*','none'};
    colors = [12, 20, 25, 1:11 13:14];
    %set(gca, 'LineStyleOrder', {'-', ':', '--', '-.'}); % different line styles
    set(0,'DefaultAxesColorOrder','remove');
    % Handle the inf
    stairsT(stairsT(:,:)==inf)=upperborder;
    for i=1:n
        A = (cumsum(stairsT(:,i)<upperborder)-1)/size(T,1);
        B = stairsT(:,i);
        stairs(B, A,'Color',myColorCodes(colors(mod(i,n)+1)),'LineStyle',linestyles{mod(i,n)+1},'Marker',markers{mod(i,n)+1}, 'MarkerSize',10,'LineWidth', 2.5);
    end
    % Correct the limits
    xlim([1 upperborder]);
    ylim([0 1]);
    set(gca,'FontSize',14)
    set(gca,'xscale','log');
    set(gca,'xtick',2.^(1:log2(upperborder)));
    %legend(legendX,'Location','southoutside','Orientation','horizontal');
    legend(legendX,'Location','southeast','Orientation','horizontal');
%     if type == 1
%          title(['Log_2 Scaled Performance Profile of f_{' num2str(funId) '} based on best fitness values']);
%     else
%          title(['Log_2 Scaled Performance Profile of f_{' num2str(funId) '} based on running time']);
%     end
    
    ylabel('P((log_2 (r_{p,s}) \leq \tau: 1 \leq s \leq n_s )');
    xlabel('\tau');
    box on;
    set(gcf,'units','points','position',[10,10,800,600]);
    
    hold off
    
return



