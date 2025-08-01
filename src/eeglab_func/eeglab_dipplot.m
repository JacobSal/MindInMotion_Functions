
function [fig] = eeglab_dipplot(STUDY, ALLEEG, clusters, varargin)
%   Detailed explanation goes here
%   IN: 
%   OUT: 
%   IMPORTANT: 
%     %## HIRES EXAMPLE
%     
%     HIRES_TEMPLATE = 'M:\jsalminen\GitHub\par_EEGProcessing\src\_data\_resources\mni_icbm152_nlin_sym_09a\mni_icbm152_t1_tal_nlin_sym_09a.nii';
%     if ~ispc
%         HIRES_TEMPLATE = convertPath2UNIX(HIRES_TEMPLATE);
%     else
%         HIRES_TEMPLATE = convertPath2Drive(HIRES_TEMPLATE);
%     end
%     %- assign hires_template default
%     tmp = strsplit(HIRES_TEMPLATE,filesep);
%     fpath = strjoin(tmp(1:end-1),filesep);
%     fname = tmp{end};
%     ext = strsplit(fname,'.');
%     fname = ext{1};
%     ext = ext{end};
%     hires_mesh = [fpath filesep fname '_mesh.mat'];
%     hires_mri = [fpath filesep fname '_mri.mat'];
%     DIPPLOT_STRUCT.mri = hires_mri;
%     DIPPLOT_STRUCT.meshdata = hires_mesh;
%     DIPPLOT_STRUCT.mri = HIRES_TEMPLATE;
%     DIPPLOT_STRUCT.meshdata = [];
%     
% CAT CODE
%  _._     _,-'""`-._
% (,-.`._,'(       |\`-/|
%     `-.-' \ )-`( , o o)
%           `-    \`_`"'-
% Code Designer: Chang Liu, Jacob Salminen
% Code Date: 04/28/2023, MATLAB R2020b
% Copyright (C) Jacob Salminen, jsalminen@ufl.edu
% Copyright (C) Chang Liu, liu.chang1@ufl.edu
cat_logo();
%  See also  pop_clustedit(), dipplot()        
%
% Authors:  Edits (02/15/2024) Jacob Salminen
%% FUNCTION SETUP
tt = tic();
fprintf('Setting up function...\n');
eeglab_fpath = get_eeglab_fpath();
%- set default paths for boundary element head model
PATH_EEGLAB_BEM  = [eeglab_fpath filesep 'plugins' filesep 'dipfit' filesep 'standard_BEM' filesep];
MNI_MRI = [PATH_EEGLAB_BEM filesep 'standard_mri.mat'];
MNI_VOL = [PATH_EEGLAB_BEM filesep 'standard_vol.mat'];
% MNI_CHAN_1005 = [PATH_EEGLAB_BEM filesep 'elec' filesep 'standard_1005.elc'];
% COORD_TRANSFORM_MNI = [0 0 0 0 0 -1.5708 1 1 1];
MNI_MRI = load(MNI_MRI);
MNI_MRI = MNI_MRI.mri;
fig = [];

%## INPUT
%-
PLOT_TYPE = 'all_nogroup';
GROUP_MARKERS = {'o','^','hexagram','diamond','pentagram','hexagram'};
%-
% COLOR_PALETTE = {[1 1 1],...        % White
%             [1 1 0]...             % Yellow
%             [221,52,151]/255,...    % Pink
%             [1 0 0],...             % Red
%             [250 140 0]/255,...     % oragne
%             [210 173 255]/255,...   % purple
%             [0.5 0.5 0],...         % Olive
%             [0.5 0 0.5],...         % Purple
%             [0.5 0 0],...           % Maroon
%             [0 1 1],...             % Aqua
%             [0 1 0],...             % Lime
%             [0 0.5 0.5],...         % Teal
%             [0 0.5 0],...           % Green
%             [0 0 1],...             % Blue
%             [0 0 0.5],...           % Navy
%             [0.8 0.8 0.8]};          % Gray
% Choosing and sorting 13 colors for clusters: 1white, 2yellow, 3pink, 4Red, 5orange, 6purple, 
% 7olive, 8purple, 9maroon, 10aqua, 11lime, 12teal, 13green, 14blue, 15navy, 16gray
% 1red, 2lime, 3blue, 4yellow, 5green, 6aqua, 7
% COLOR_PALETTE = COLOR_PALETTE([4 11 14 2 13 10 5 6 15 16 1 7 9 3]);
%-- alternative
COLOR_PALETTE = num2cell(linspecer(16),2)';

%## DIPPLOT STRUCT
DEFAULT_DIPPLOT_STRUCT = struct('rvrange',[0,15],... % this is a value from 0 to 100 (e.g., rv = 0.15 is 15)
    'summary','off',...
    'mri',MNI_MRI,...
    'coordformat','MNI',...
    'transform',[],...
    'image','mri',...
    'verbose','on',...
    'plot','on',...
    'color',{{[0,0,1]}},...
    'view',[0.5,0.5,0.5],...
    'mesh','off',...
    'meshdata',MNI_VOL,...
    'axistight','off',...
    'gui','off',...
    'num','off',...
    'cornermri','off',...
    'drawedges','off',...
    'projimg','off',...
    'projlines','off',...
    'projwidth',1,...
    'projcol',{{[0,0,0]}},...
    'projalpha',1, ...
    'dipolesize',30,...
    'dipolelength',0,...
    'pointout','off',...
    'sphere',1,...
    'spheres','off',...
    'normlen','off',...
    'dipnames',{{}},...
    'holdon','on',...
    'camera','set',...
    'density','off');
%     'spheresize',1,... not an actual argument
DEFAULT_PLOT_STRUCT = struct( ...
    'figure_position',[100,100,720,720], ...
    'marker_face_alpha',0.6, ...
    'marker_linew_mult',0.04, ...
    'marker_sized_mult',4, ...
    'marker_type','o' ...
    );
%## Define Parser
p = inputParser;
%## REQUIRED
addRequired(p,'STUDY',@isstruct);
addRequired(p,'ALLEEG',@isstruct);
addRequired(p,'clusters',@isnumeric);
%## PARAMETER
addParameter(p,'PLOT_TYPE',PLOT_TYPE,@ischar);
addParameter(p,'COLOR_PALETTE',COLOR_PALETTE,@iscell);
addParameter(p,'DIPPLOT_STRUCT',DEFAULT_DIPPLOT_STRUCT,@(x) validate_struct(x,DEFAULT_DIPPLOT_STRUCT));
addParameter(p,'PLOT_STRUCT',DEFAULT_PLOT_STRUCT,@(x) validate_struct(x,DEFAULT_PLOT_STRUCT));
parse(p,STUDY,struct.empty,clusters,varargin{:});
%## DEFAULTS
PLOT_STRUCT = set_defaults_struct(p.Results.PLOT_STRUCT,DEFAULT_PLOT_STRUCT);
DIPPLOT_STRUCT = set_defaults_struct(p.Results.DIPPLOT_STRUCT,DEFAULT_DIPPLOT_STRUCT);

PLOT_TYPE = p.Results.PLOT_TYPE;
COLOR_PALETTE = p.Results.COLOR_PALETTE;
% PLOT_STRUCT = p.Results.PLOT_STRUCT;
% DIPPLOT_STRUCT = p.Results.DIPPLOT_STRUCT;
fprintf('done. %0.2f\n',toc(tt));
%% SET REST OF DEFAULTS (NOT FUNCTIONAL FOR THIS)
% PLOT_STRUCT = set_defaults_struct(PLOT_STRUCT,DEFAULT_PLOT_STRUCT);
% DIPPLOT_STRUCT = set_defaults_struct(DIPPLOT_STRUCT,DEFAULT_DIPPLOT_STRUCT);
%% (PARAM SETUP)
try
    tmp = strsplit(DIPPLOT_STRUCT.mri,filesep);
    fpath = strjoin(tmp(1:end-1),filesep);
    fname = tmp{end};
    ext = strsplit(fname,'.');
    fname = ext{1};
    ext = ext{end};
    MRI_FPATH = [fpath filesep fname '_dipplotmri.mat'];
    VOL_FPATH = [fpath filesep fname '_dipplotvol.mat'];
catch e
    fprintf('Tried assigning DIPPLOT_STRUCT.mri...\n');
    fprintf('\n%s\n',getReport(e));
end
tt = tic;
fprintf('Assigning MRI & Volume for plotting...\n');
try
    if isstruct(DIPPLOT_STRUCT.mri)
        fprintf('Not loading mri, using inputted struct...\n');
    else
        fprintf('Loading %s & and creating .mat file for mri...\n',DIPPLOT_STRUCT.mri);
        if strcmp(ext,'nii')
            fprintf('Creating mri.mat and meshdata.mat\n');
            ftmri_tmp = ft_read_mri([fpath filesep fname '.' ext]);
            ftmri_tmp.coordsys = 'acpc';
            
            %- reformat to eeglab style...
            mri = [];
            mri.dim = ftmri_tmp.dim;
            mri.xgrid = 1:ftmri_tmp.dim(1);
            mri.ygrid = 1:ftmri_tmp.dim(2);
            mri.zgrid = 1:ftmri_tmp.dim(3);
            mri.anatomy = uint8(ftmri_tmp.anatomy*1.5); %uint8(tmp.anatomy); %uint8(tmp.anatomy);
%             mri.transform = [1,0,0,-91;...
%                             0,1,0,-126;...
%                             0,0,1,-73;...
%                             0,0,0,1];
            mri.transform = ftmri_tmp.transform;
            mri.hdr = ftmri_tmp.hdr;
%             mri.unit = ftmri_tmp.unit;
%             mri.coordsys = 'mni152'; %tmp.coordsys;
            save(MRI_FPATH,'mri');
%             mri = ftmri_tmp;
        elseif strcmp(ext,'mat')
            tmp = load(DIPPLOT_STRUCT.mri);
            DIPPLOT_STRUCT.mri = tmp{1};
            mri = tmp;
        end
    end
    %- check meshdata
    if ~exist(DIPPLOT_STRUCT.meshdata,'file') || isempty(DIPPLOT_STRUCT.meshdata)
        %--
        cfg             = [];
        cfg.output      = {'scalp', 'skull', 'brain'};
        cfg.spmmethod   = 'new';
        cfg.spmversion  = 'spm12';
        segmentation    = ft_volumesegment(cfg, ftmri_tmp);
        %--
        cfg             = [];
        cfg.tissue      = {'scalp', 'skull', 'brain'};
        cfg.numvertices = [500,1000,1500];
        mesh            = ft_prepare_mesh(cfg, segmentation);
%                 save([fpath filesep fname '_mesh.mat'],'mesh');
        cfg             = [];
        cfg.method      = 'bemcp';
        ftvol_tmp       = ft_prepare_headmodel(cfg, mesh);
        %- convert to eeglab format....
        vol = [];
        %-
%         vol.ft_vol = ftvol_tmp;
%         vol.bnd = ftvol_tmp.bnd; %ftvol_tmp.bnd(3:-1:1);
%         vol.cond = ftvol_tmp.cond;
%         vol.mat = ftvol_tmp.mat;
%         vol.type = 'bemcp';
        %- 
        vol.ft_vol = ftvol_tmp;
        vol.bnd = ftvol_tmp.bnd(3:-1:1);
        vol.cond = ftvol_tmp.cond;
        tmp_mat = ftvol_tmp.mat;
        tmp_mat = padarray(tmp_mat,[2500,0],0,'post');
        tmp_mat = sparse(tmp_mat);
        vol.mat = tmp_mat;
        vol.type = 'bemcp';
        %-
%         vol = ftvol_tmp;
        %-
        save(VOL_FPATH,'vol');
        DIPPLOT_STRUCT.meshdata = VOL_FPATH;
        %## DEBUG
        %{
        %-
        figure;
        cfg              = [];
        cfg.funparameter = 'scalp'; 'brain','skull'
        cfg.anaparameter = 'anatomy';
        cfg.method       = 'ortho';
        ft_sourceplot(cfg, segmentation, mri);
        %}

        %## EEGLAB & MINE VALIDATION
        %-- volume
        %{
        vol_eeglab = load('M:\jsalminen\GitHub\MIND_IN_MOTION_PRJ\MindInMotion_YoungerOlderAdult_KinEEGCorrs\submods\eeglab\plugins\dipfit\standard_BEM\standard_vol.mat');
        vol_eeglab = vol_eeglab.vol;
        vol = load(DIPPLOT_STRUCT.meshdata);
        mesh = vol.vol.ft_vol;
        tmp = mesh.bnd(1);
        figure;
        ft_plot_mesh(tmp,'facecolor','none','edgecolor',[0,0,0]);
        hold on;            
        ft_plot_mesh(vol_eeglab.bnd(3),'facecolor','none','edgecolor',[1,0,0]);
        %(03/21/2025) JS, a little bit of none overlap, but essentially the
        same
        %}
        %-- mri
        %{
        mri_eeglab = load('M:\jsalminen\GitHub\MIND_IN_MOTION_PRJ\MindInMotion_YoungerOlderAdult_KinEEGCorrs\submods\eeglab\plugins\dipfit\standard_BEM\standard_mri.mat');
        mri_eeglab = mri_eeglab.mri;
        mri = DIPPLOT_STRUCT.mri;
        disp(mri.transform); disp(mri_eeglab.transform);
        % figure;
        % ft_plot_mesh(tmp,'facecolor','none','edgecolor',[0,0,0]);
        % hold on;            
        % ft_plot_mesh(vol_eeglab.bnd(3),'facecolor','none','edgecolor',[1,0,0]);
        %}
    else
        fprintf('Using mesh/vol file: %s\n',DIPPLOT_STRUCT.meshdata);
    end
catch e
    fprintf('Error. %s\n',getReport(e));
    fprintf('\nUsing deafult mri and volume...\n\n');
    DIPPLOT_STRUCT.mri = MNI_MRI;
    DIPPLOT_STRUCT.meshdata = MNI_VOL;
end
fprintf('done. %0.2f\n',toc(tt));

%% ===================================================================== %%
%## AXIS TIGHT PARAMETERS
%-- getting center slices
sz = size(DIPPLOT_STRUCT.mri.anatomy);
sz1 = floor(sz(1)/2);
sz2 = floor(sz(2)/2);
sz3 = floor(sz(3)/2);

%## INITIATE FIGURES
fig = figure('color','black');
set(fig,'Position',PLOT_STRUCT.figure_position)
set(gca,'box','off', ...
    'xtick',[], ...
    'ytick',[],...
    'ztick',[], ...
    'xcolor',[0,0,0], ...
    'ycolor',[0,0,0]);

if length(clusters) > length(COLOR_PALETTE)
    COLOR_PALETTE = linspecer(length(clusters)+1);
end

%## LOOP THROUGH DIPOLES
%--
hold on;
for cc = 1:length(clusters)
    fprintf('Displaying cluster %i...\n',clusters(cc));
    tt = tic;
    %--
    tmp_dipplot_struct = DIPPLOT_STRUCT;
    cl_i = clusters(cc);       
    sets_i = STUDY.cluster(cl_i).sets;
    comps_i = STUDY.cluster(cl_i).comps;
    dips = cell(length(sets_i),1);
    %--
    for i = 1:length(sets_i)
        tmp = [];
        ind = sets_i(i) == sets_i;
        tmp.grp = ALLEEG(sets_i(i)).group;
        tmp.posxyz = ALLEEG(sets_i(i)).dipfit.model(comps_i(ind)).posxyz;
        tmp.momxyz = ALLEEG(sets_i(i)).dipfit.model(comps_i(ind)).momxyz;
        tmp.rv = ALLEEG(sets_i(i)).dipfit.model(comps_i(ind)).rv;
        tmp.pot = ALLEEG(sets_i(i)).dipfit.model(comps_i(ind)).datapot;
%         tmp.component = comps_i(ind);
        dips{i} = tmp;
    end
    %- concatenate to struct & remove empty dips.
    dips = dips(~cellfun(@isempty,dips));
    dips = cat(2,dips{:});
    %- colors
    try
        tmp_dipplot_struct.color = COLOR_PALETTE(cl_i); %repmat(colors{cl_i},[length(dips),1]);
    catch e
        fprintf('Error. Color couldn''t be selected\n');
        fprintf('%s\n',getReport(e));
        tmp_dipplot_struct.color = [1,0,0];
    end
    %## DIPPLOT SETTINGS
    switch PLOT_TYPE
        case 'average_nogroup'
            %- average dip
            mean_dip = mean(cat(1,dips.posxyz),1);
            dist = cellfun(@(x) norm(x-mean_dip),{dips.posxyz});
            mean_dist = mean(dist);
            std_dist = std(dist);
            avgdip = struct('posxyz',mean(cat(1,dips.posxyz),1),...
                'momxyz',mean(cat(1,dips.momxyz),1),...
                'rv',mean(cat(1,dips.rv),1));
            %## DEBUG ORIGIN
            %{
            avgdip = struct('posxyz',[0,0,0],...
                'momxyz',[0,0,0],...
                'rv',0.1);
            %}
            % DIPPLOT_STRUCT.dipolesize =  5.8376*(mean_dist+std_dist*3)*3*1.5/10; 
            tmp_dipplot_struct.dipolesize =  5.8376*360*DIPPLOT_STRUCT.dipolesize;
            % (02/20/2024) JS, no good way of defining the dipole size, but
            % can provide approximations with the right scaling
            % 1150 (with default zoom) ~197mm? => 5.8376 dipsize/mm            
            % cfg = tmp_dipplot_struct;
            %## SCATTER 3
            % color_hold = tmp_dipplot_struct.color;
            scatter3_dipplot(avgdip,tmp_dipplot_struct);   
            %## 
            % h = findobj(fig.Children(end).Children,'type','Scatter');
            % for ii = 1:length(h)
            %     hold on;
            %     tmp_dipplot_struct.color = color_hold;
            %     h(ii).MarkerFaceColor = tmp_dipplot_struct.color{1};
            %     h(ii).Marker = PLOT_STRUCT.marker_type;
            %     h(ii).ColorVariable = [];
            %     h(ii).CData = [0,0,0];
            %     h(ii).LineWidth = h(ii).SizeData*PLOT_STRUCT.marker_linew_mult;
            %     h(ii).MarkerFaceAlpha = PLOT_STRUCT.marker_face_alpha;
            %     h(ii).SizeData = h(ii).SizeData*PLOT_STRUCT.marker_sized_mult;
            % end   
           
        case 'all_nogroup'
            %- plot all dipoles
            tmp_dipplot_struct.dipolesize = 5.8376*360*DIPPLOT_STRUCT.dipolesize;
            % dipplot(dips,DIPPLOT_STRUCT);
            % color_hold = tmp_dipplot_struct.color;
            scatter3_dipplot(dips,tmp_dipplot_struct);
            
            %##
            % h = findobj(fig.Children(end).Children,'type','Scatter');
            % for ii = 1:length(h)
            %     hold on;
            %     % tmp_dipplot_struct.color = color_hold;
            %     % h(ii).MarkerFaceColor = tmp_dipplot_struct.color{1};
            %     h(ii).Marker = PLOT_STRUCT.marker_type;
            %     h(ii).ColorVariable = [];
            %     h(ii).CData = [0,0,0];
            %     h(ii).LineWidth = h(ii).SizeData*PLOT_STRUCT.marker_linew_mult;
            %     h(ii).MarkerFaceAlpha = PLOT_STRUCT.marker_face_alpha;
            %     h(ii).SizeData = h(ii).SizeData*PLOT_STRUCT.marker_sized_mult;
            % end     
        case 'all_group'
            %- plot all dipoles
            tmp_dipplot_struct.dipolesize = 5.8376*(20)*3*1.5/10*DIPPLOT_STRUCT.dipolesize;
            grps = unique({dips.grp});
            if length(grps) > 1
                color_chng = (1/(length(grps)))+0.1;
            else
                color_chng = 1;
            end
            color_hold = tmp_dipplot_struct.color;
            % cnt = 1;
            % dipplot(dips,DIPPLOT_STRUCT);
            % cnt = length(fig.Children(end).Children);
            scatter3_dipplot(dips,tmp_dipplot_struct);            
            
            %##
            h = findobj(fig.Children(end).Children,'type','Scatter');
            for ii = 1:length(h)
                hold on;
                tmp_dipplot_struct.color = color_hold;
                gg = find(strcmp(dips(ii).grp,grps));
                kk = color_chng*(gg-1);
                if kk ~= 0
                    tmp_dipplot_struct.color = {COLOR_PALETTE{cl_i}*(kk)};
                end
                h(ii).Marker = GROUP_MARKERS{gg};
                h(ii).ColorVariable = [];
                h(ii).CData = [0,0,0];
                h(ii).MarkerFaceColor = tmp_dipplot_struct.color{1};
                h(ii).LineWidth = h(ii).SizeData*PLOT_STRUCT.marker_linew_mult;
                h(ii).MarkerFaceAlpha = PLOT_STRUCT.marker_face_alpha;
                h(ii).SizeData = h(ii).SizeData*PLOT_STRUCT.marker_sized_mult;
            end      
        otherwise
            fprintf('You didn''t choose a viable plotting option...\n')
            % %- plot all dipoles
            % DIPPLOT_STRUCT.dipolesize = 5.8376*(10)*3/10;
            % % dipplot(dips,DIPPLOT_STRUCT);
            % scatter3_dipplot(dips,DIPPLOT_STRUCT);
    end
    hold on;
end
%## CHANGE GENERAL PROPERTIES 
switch PLOT_TYPE
    case {'all_nogroup','average_nogroup'}
        h = findobj(fig.Children(end).Children,'type','Scatter');
        for ii = 1:length(h)
            hold on;
            % tmp_dipplot_struct.color = color_hold;
            % h(ii).MarkerFaceColor = tmp_dipplot_struct.color{1};
            h(ii).Marker = PLOT_STRUCT.marker_type;
            h(ii).ColorVariable = [];
            h(ii).CData = [0,0,0];
            h(ii).LineWidth = h(ii).SizeData*PLOT_STRUCT.marker_linew_mult;
            h(ii).MarkerFaceAlpha = PLOT_STRUCT.marker_face_alpha;
            h(ii).SizeData = h(ii).SizeData*PLOT_STRUCT.marker_sized_mult;
        end
end
%## OVERRIDE TEXTURES && TIGHT AXES
% for ci = 1:3
%     fig.Children(end).Children(ci).EdgeColor = 'texturemap';
%     fig.Children(end).Children(ci).FaceColor = 'flat';
% end
if strcmp(DIPPLOT_STRUCT.axistight,'off')
    for ci = 1:3
        switch ci
            case 1
                xlim = fig.Children(end).Children(ci).XData(1,2)-fig.Children(end).Children(ci).XData(1,1);
                ylim = fig.Children(end).Children(ci).YData(2,1)-fig.Children(end).Children(ci).YData(1,1);
                im = squeeze(DIPPLOT_STRUCT.mri.anatomy(1:xlim,1:ylim,sz3))';
            case 2
                xlim = fig.Children(end).Children(ci).XData(1,2)-fig.Children(end).Children(ci).XData(1,1);
                ylim = fig.Children(end).Children(ci).ZData(2,1)-fig.Children(end).Children(ci).ZData(1,1);
                im = squeeze(DIPPLOT_STRUCT.mri.anatomy(1:xlim,sz2,1:ylim))';                            
            case 3
                xlim = fig.Children(end).Children(ci).YData(1,2)-fig.Children(end).Children(ci).YData(1,1);
                ylim = fig.Children(end).Children(ci).ZData(2,1)-fig.Children(end).Children(ci).ZData(1,1);
                im = squeeze(DIPPLOT_STRUCT.mri.anatomy(sz1,1:xlim,1:ylim))';                            
        end
        fig.Children(end).Children(ci).CData = im;        
    end
    colormap(fig,'gray');
end


%## INCLUDE 3D BRAINMESH
if strcmp(DIPPLOT_STRUCT.mesh,'on')
    hold on;
    % tcmap = linspecer(40);
    FACE_COLOR = 'none'; %tcmap(20,:);
    EDGE_COLOR = [0,0,0];
    %--
    vol = load(DIPPLOT_STRUCT.meshdata);
    mesh = vol.vol.ft_vol;
    tmp = mesh.bnd(1);
    %--
    ft_plot_mesh(tmp, ...
        'facecolor',FACE_COLOR, ...
        'edgecolor',EDGE_COLOR, ...
        'facealpha',0.3, ...
        'edgealpha',0.1);
    drawnow;
end
hold off;
drawnow;
fprintf('done. %0.2f\n',toc(tt));
end

%%
function [eeglab_fpath] = get_eeglab_fpath()
    %## GET EEGLAB PATH
    if ~ispc
        tmp = strsplit(path(),':');
    else
        tmp = strsplit(path(),';');
    end
    b1 = regexp(tmp,'eeglab','end');
    b2 = tmp(~cellfun(@isempty,b1));
    eeglab_fpath = b2{1};
    tmp = strsplit(eeglab_fpath,filesep);
    eeglab_ind = find(strcmp(tmp,'eeglab'));
    if ~ispc
        eeglab_fpath = [filesep strjoin(tmp(1:eeglab_ind),filesep)];
    else
        eeglab_fpath = strjoin(tmp(1:eeglab_ind),filesep);
    end
    fprintf('EEGLAB path: %s\n',eeglab_fpath);
end
