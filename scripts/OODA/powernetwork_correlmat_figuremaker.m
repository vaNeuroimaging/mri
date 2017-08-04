function h = powernetwork_correlmat_figuremaker(matrix,colors_new,varargin)
%
% Name:figure_corrmat_powernetwork.m
% $Revision: 1.1 $
% $Date: 2015/08/04 15:05:28 $
% FROM Tim - works for his network structure

networks = {'Unlabelled';'SM';'SM (lat)';'CO';'Auditory';'DMN';'Memory';'Visual';'FP';'Salience';'Sub-cortex';'Ventral Attn';'Dorsal Attn';'Cerebellum'};
%load('/data/cn5/caterina/TaskConn_Methods/all_data/Power_module_colors.mat');
%colors_new(1,:) = 0;
h = figure('Color',[0.8275 0.8275 0.8275],'Position',[56 143 1295 807]); %[56 143 1095 807]
%h = figure('Color',[0.8275 0.8275 0.8275],'Position',[2061 162 1095 807]);
if nargin>2
    if ~isempty(varargin{1}) & ~isempty(varargin{2})
        climlow = varargin{1};
        climhigh = varargin{2};
        imagesc(matrix,[climlow climhigh]);
    else
        imagesc(matrix);
    end
else
    imagesc(matrix);
end

vline_new([29 59 64 78 91 149 154 185 210 228 241 250 261],'r',3);
hline_new([29 59 64 78 91 149 154 185 210 228 241 250 261],'r',3);
tickpos = [14 44 61 70 84 119 151 168 197 218 234 244 255 262];

netpos = [0.5 29 59 64 78 91 149 154 185 210 228 241 250 261 264.5];
ax = axis;

set(gca,'XTick',tickpos,'Xlim',[ax(1) ax(2)]);

% for i = 1:length(networks)
%     rx = rectangle('Position',[netpos(i) 264.5 netpos(i+1)-netpos(i) 50],'FaceColor',[colors_new(i,1) colors_new(i,2) colors_new(i,3)],'Clipping','off','LineWidth',2);
%     ry = rectangle('Position',[-50 netpos(i) 50.5 netpos(i+1)-netpos(i)],'FaceColor',[colors_new(i,1) colors_new(i,2) colors_new(i,3)],'Clipping','off','LineWidth',2);
% end

set(gca,'XTicklabel','');
set(gca,'YTicklabel','');    

%tx_outline = text(tickpos,ones(1,length(tickpos))*265,networks);
tx= text(tickpos,ones(1,length(tickpos))*265,networks);
%set(tx_outline,'HorizontalAlignment','right','VerticalAlignment','top','Rotation',45)
set(tx,'HorizontalAlignment','right','VerticalAlignment','top','Rotation',45);

for i = 1:length(tx)
     %set(tx_outline(i),'Color','k','FontSize',15,'FontWeight','bold'); 
     set(tx(i),'Color',[colors_new(i,1) colors_new(i,2) colors_new(i,3)],'FontName','Helvetica','FontSize',14,'FontWeight','bold');   
end
set(gca,'FontWeight','bold','FontSize',14);


%tickpos = [10 41 59 68 82 117 148 165 195 215 231 242 251 258.5];
tickpos = [9 39 56 67 80 116 145 164 194 214 229 241 250 256.5];
ty= text(-1*ones(1,length(tickpos)),tickpos,networks);
set(ty,'HorizontalAlignment','right','VerticalAlignment','top')

for i = 1:length(ty)
    set(ty(i),'Color',[colors_new(i,1) colors_new(i,2) colors_new(i,3)],'FontName','Helvetica','FontSize',16,'FontWeight','bold');   
end
colorbar;
set(gca,'FontWeight','bold','FontSize',16);

if nargin>4
    if ~isempty(varargin{4})
        title(titletext,'FontWeight','bold','FontSize',14);
    end
end

