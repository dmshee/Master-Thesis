function [hgroup]=gplot_mean(x,group,wl,varargin)
% function to plot mean and standard deviation of groups of spectra (x)
% color-coded by group.
% input:
% x: spectra, samples in columns, variables in rows, as in normal plot
% function
%
% group: group number used for color coding; one label per group number from 1:max(group); even if there is
% no samples for a certain group number. 'group' should be a natural
% number!
% varargin{1}: legend in a cell format: {'h1','h2','h3','h4'};
%
% output:
% h: cell with figure axis for each group 
%
%[h]=gplot_mean(x,group,{'h1','h2','h3','h4'})
%
%nbins=size(unique(y),1);colors=lines(nbins);
%gplot_mean(x,y,legend0,colors)
%
% Example to change the legend
% legend ([h(1),h(2),h(4)],{'h1','h2','h4'})% only for groups with
% data
%
%%[h]=gplot(x,group,{'h1','h2','h3','h4'})
[C, ia, ic] = unique(group);
if size(x,1)~=size(ic,1)
   error('"x" and "group" are not the same size; review function inputs')
end
nz=size(x,2);
nbins=max(ic);

if nargin>4&& ~isempty(varargin{2}) % select colors
colours=varargin{2};
if size(colours,2)&& ~3
    error ('input colours should be in RGB matlab format: 3 columns per color with 0 to 1 values')
end
else
colours=jet(nbins);
end

hgroup=plot([]); %inicializes figure handle for the legend
figure,

for i=1:nbins
h{i}=plot(wl,mean(x(ic==i,:))','color',colours(i,:));
std1=mean(x(ic==i,:))+std(x(ic==i,:));
std2=mean(x(ic==i,:))-std(x(ic==i,:));
hold on, plot(wl,std1','color',colours(i,:),'LineStyle','--')
hold on, plot(wl,std2','color',colours(i,:),'LineStyle','--')
hold on, fill([wl flip(wl)],[std1 flip(std2)],colours(i,:),'LineStyle','none'); alpha(.25)    
if ~isempty (h{i}) % to include the legend of non empty classes only
    hgroup=[hgroup,h{i}(1)];
    end
hold on,
end

if nargin>2 && ~isempty(varargin{1}) %plots the legend
lg=varargin{1};    
text_legend=lg(unique(ic));
legend(hgroup,text_legend);
end



end