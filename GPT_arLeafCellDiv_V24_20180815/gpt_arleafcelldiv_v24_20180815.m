function m = gpt_arleafcelldiv_v24_20180815( m )
%m = gpt_arleafcelldiv_v24_20180815( m )
%   Morphogen interaction function.
%   Written at 2018-08-15 17:05:00.
%   GFtbox revision 6021, 2018-07-18 15:25.
%   Model last saved to SVN as revision 4004, 2012-02-08 14:49:55.085950.

% The user may edit any part of this function lying between lines that
% begin "%%% USER CODE" and "%%% END OF USER CODE".  Those lines themselves
% delimiters themselves must not be moved, edited, deleted, or added.

    if isempty(m), return; end

    fprintf( 1, '%s found in %s\n', mfilename(), which(mfilename()) );

    setGlobals();
    realtime = m.globalDynamicProps.currenttime;
    dt = m.globalProps.timestep;

%%% USER CODE: INITIALISATION

% This is the interaction function for the GFtbox model associated with the
% paper "Spatiotemporal coordination of cell division and growth during
% organ morphogenesis", Figures 8 and 9.
%
% Due to a change in our terminology while the paper was in preparation,
% some of the morphogens and cell factors have different names here and in
% the paper.  Everything here with "mdgf" in its name relates to PMF in the
% paper, and likewide "cdiv" relates to COMP.

m = leaf_set_userdatastatic( m, 'collectsplitdata', true );

%% model init
m = leaf_setproperty( m, ...
    'bioApresplitproc', @bioApresplitproc );

PLOTTING_LINEAGE = false;
if false % PLOTTING_LINEAGE
    CELL_EDGE_COLOR = [0.95 0.95 0.95];  % Almost white.
else
    CELL_EDGE_COLOR = [0 0 0];  % Black.
end

if Steps(m) == 0
    
    % Model Visualisation and Setup Parameters
    m = leaf_setproperty( m, 'mingradient', 0, 'timestep', 1, 'userpolarisation',false, ...
        'bioAsplitcells', true, 'allowSplitBio',true, ...
        'bioApullin', 0.09 );
    m = leaf_plotoptions( m, 'FEthinlinesize', 1, 'layeroffset', 0.45, 'bioAlinecolor', CELL_EDGE_COLOR );
    m = leaf_fix_vertex( m, 'vertex', find(m.nodes(:,2)<=min(m.nodes(:,2)+0.001)), 'dfs', 'y' );
    
    % Model output times to plot images for the paper.
    %m.userdata.outputtimes = [86 89 100 124 140 147 153 162 175 182 205 216 240 264 313 338 362 386 412];
    m.userdata.outputtimes = [86 115 132 140 147 156 164 178 204 230 256 282 308 334 360 386 412]; %linspace(178,412,10)
    m.userdata.output = 0; % make high resolution images when time = m.userdata.outputtimes; 0 = no output.
    
    if m.globalProps.IFsetsoptions
        % Setting up options.
        % Most of the modelnames correspond to figures in the paper, as
        % annotated.  Some do not, and were for other investigations:
        %   EPIMODEL_NO_MID_COMP
        %   SPCH_ON_PLATES_NO_MID_COMP
        modelnames = { 'SUBEPIMODEL', ... % Figure 8 E,F - spch subepidermis model
                       'EPIMODEL', ... % Figure 8 G,H - spch epidermis model and Figure 9B and S17B - later stage mature spch epidermis model
                       'EPIMODEL_NO_POLARISER', ... % Figure S17C - Like EPIMODEL, but with no polariser, to observe the statistics on cell neighbour numbers.
                       'EPIMODEL_NO_MID_COMP', ... % (Not used.) Like EPIMODEL, but without the turning off of id_cdiv on the midline.
                       'SPCH_ON_PLATES_NO_MID_COMP', ... % Not used.) Like SPCH_ON_PLATES, but without the turning off of id_cdiv on the midline.
                       'EPIMODEL-FACTOR-DISTALGROWTH', ... % Figure 9 C - later stage spch feedback-free epidermis model
                       'EPIMODEL-CELLS-DISTALGROWTH', ... % Figure 9 J - later stage spch cell-feedback epidermis model
                       'SHIFT_LATE_FACTOR_EARLY', ... % Figure 9 H
                       'SHIFT_LATE_FACTOR_LATE', ... % Figure 9 I
                       'CHANGE_GROWTH_RATE_FACTOR_LOW', ... % Figure 9 F
                       'CHANGE_GROWTH_RATE_FACTOR_HIGH', ... % Figure 9 G
                       'CHANGE_CELL_DIVTHRESH_FACTOR_LOW', ... % Figure 9 E
                       'CHANGE_CELL_DIVTHRESH_FACTOR_HIGH', ... % Figure 9 D
                       'SHIFT_LATE_CELLS_EARLY', ... % Figure 9 O
                       'SHIFT_LATE_CELLS_LATE', ... % Figure 9 P
                       'CHANGE_GROWTH_RATE_CELLS_LOW', ... % Figure 9 M
                       'CHANGE_GROWTH_RATE_CELLS_HIGH', ... % Figure 9 N
                       'CHANGE_CELL_DIVTHRESH_CELLS_LOW', ... % Figure 9 L
                       'CHANGE_CELL_DIVTHRESH_CELLS_HIGH', ... % Figure 9 K
                       'SPCH_IN_CHAMBER', ... % Fig S14C - SPCH grown in the imaging chamber
                       'SPCH_ON_PLATES' ... % Fig 10B - SPCH grown on plates
                     };

        [m,options] = setUpModelOptions( m, ...
                'modelname',    modelnames,      'EPIMODEL', ...
                'modeldeco',    {'NONE', 'CELLSFORDIV'}, 'CELLSFORDIV', ...
                'cellIntercellSpaces', [],       true, ...
                'celldivnoise', [],              true, ...
                'bpgrad',       [],              0.195, ...
                'ppgrad',       [],              0.041, ...
                'plam',         [],              0.0235, ...
                'glate',        [0.0048 0.0096], 0.0048, ...
                'hlate',        [],              2.2, ...
                'bpol',         [],              0.1, ...
                'hmid',         [],              1.0, ...
                'plate',        [],              0.7, ...
                'BOUNDRYPOS',   [],              -0.005, ...
                'PgradTHRESH',  [],              0.628, ...
                'pmdgf',        [0.9 0],         0.9, ...
                'c_area_scaling', [1500 15000],  1500 ...
            );
    end
    
    %model csv output
    % Prepare the folder to extract growth data
    if exist('Extracted Growth', 'dir') == 0
        mkdir 'Extracted Growth';
    end
    
    % Prepare the csv file for leaf dimensions and average growth 
    logFileHeaders = {'realtime', 'LeafArea', 'LeafWidth', 'LeafLength', ...
        'late', 'plam', 'numberOfCells', 'meanCellArea', ...
        'numberOfCellInCDIV', 'meanCellAreaInCDIV', 'cellAreaInCDIV' ...
        'numberOfCellOutCDIV', 'meanCellAreaOutCDIV', 'cellAreaOutCDIV'};
    modelname = getModelOption( m, 'modelname' );
    createLogFile( logFileHeaders, [modelname,'_LeafData.csv'] );
end

OPTIONS = getModelOptions( m );
printModelOptions( m );
% m = leaf_setproperty( m, 'legendTemplate', ['%T: ' OPTIONS.modelname] );

GROWTH_START_TIME = 87;  % Everything before this is setup to initialise the simulation.

if realtime >= GROWTH_START_TIME
    switch OPTIONS.modelname
        case 'SPCH_IN_CHAMBER'
            GROWTH_RATIO = 1;
            PHYSIO_TIME_RATIO = 0.75;
        case { 'SPCH_ON_PLATES', 'SPCH_ON_PLATES_NO_MID_COMP' }
            GROWTH_RATIO = 0.6;
            PHYSIO_TIME_RATIO = 0.55; % growth rate wrt real time is therefore 0.6*0.55 = 0.33
        otherwise
            % Other versions of the model make no adjustment.
            GROWTH_RATIO = 1;
            PHYSIO_TIME_RATIO = 1;
    end
    physiotime = GROWTH_START_TIME + (realtime - GROWTH_START_TIME) * PHYSIO_TIME_RATIO;
    physiotimestep = dt * PHYSIO_TIME_RATIO;
else
    GROWTH_RATIO = 1;
    PHYSIO_TIME_RATIO = 1;
    physiotime = realtime;
    physiotimestep = dt;
end

BOUNDRYPOS = OPTIONS.BOUNDRYPOS;
PgradTHRESH = OPTIONS.PgradTHRESH;
celldivnoise = OPTIONS.celldivnoise;
cellIntercellSpaces = OPTIONS.cellIntercellSpaces;
modelname = OPTIONS.modelname;
modeldeco = OPTIONS.modeldeco;

% bpol = m.userdata.ranges.bpol.range(m.userdata.ranges.bpol.index)
bpol = OPTIONS.bpol;

m.plotdefaults.perelement = [];
m.plotdefaults.perelementaxes = [];
m.plotdefaults.perelementcomponents = [];


plotArabidopsisRun( m );
%%% END OF USER CODE: INITIALISATION

%%% SECTION 1: ACCESSING MORPHOGENS AND TIME.
%%% AUTOMATICALLY GENERATED CODE: DO NOT EDIT.

% Each call of getMgenLevels below returns four results:
% XXX_i is the index of the morphogen called XXX.
% XXX_p is the vector of all of its values.
% XXX_a is its mutation level.
% XXX_l is the "effective" level of the morphogen, i.e. XXX_p*XXX_a.
% In SECTION 3 of the automatically generated code, all of the XXX_p values
% will be copied back into the mesh.

    polariser_i = FindMorphogenRole( m, 'POLARISER' );
    P = m.morphogens(:,polariser_i);
    [kapar_i,kapar_p,kapar_a,kapar_l] = getMgenLevels( m, 'KAPAR' );  %#ok<ASGLU>
    [kaper_i,kaper_p,kaper_a,kaper_l] = getMgenLevels( m, 'KAPER' );  %#ok<ASGLU>
    [kbpar_i,kbpar_p,kbpar_a,kbpar_l] = getMgenLevels( m, 'KBPAR' );  %#ok<ASGLU>
    [kbper_i,kbper_p,kbper_a,kbper_l] = getMgenLevels( m, 'KBPER' );  %#ok<ASGLU>
    [knor_i,knor_p,knor_a,knor_l] = getMgenLevels( m, 'KNOR' );  %#ok<ASGLU>
    [strainret_i,strainret_p,strainret_a,strainret_l] = getMgenLevels( m, 'STRAINRET' );  %#ok<ASGLU>
    [arrest_i,arrest_p,arrest_a,arrest_l] = getMgenLevels( m, 'ARREST' );  %#ok<ASGLU>
    [v_leaf_i,v_leaf_p,v_leaf_a,v_leaf_l] = getMgenLevels( m, 'V_LEAF' );  %#ok<ASGLU>
    [id_proxorg_i,id_proxorg_p,id_proxorg_a,id_proxorg_l] = getMgenLevels( m, 'ID_PROXORG' );  %#ok<ASGLU>
    [id_inc_i,id_inc_p,id_inc_a,id_inc_l] = getMgenLevels( m, 'ID_INC' );  %#ok<ASGLU>
    [id_pgrad_i,id_pgrad_p,id_pgrad_a,id_pgrad_l] = getMgenLevels( m, 'ID_PGRAD' );  %#ok<ASGLU>
    [id_mid_i,id_mid_p,id_mid_a,id_mid_l] = getMgenLevels( m, 'ID_MID' );  %#ok<ASGLU>
    [id_lam_i,id_lam_p,id_lam_a,id_lam_l] = getMgenLevels( m, 'ID_LAM' );  %#ok<ASGLU>
    [karea_i,karea_p,karea_a,karea_l] = getMgenLevels( m, 'KAREA' );  %#ok<ASGLU>
    [id_late_i,id_late_p,id_late_a,id_late_l] = getMgenLevels( m, 'ID_LATE' );  %#ok<ASGLU>
    [id_distal_i,id_distal_p,id_distal_a,id_distal_l] = getMgenLevels( m, 'ID_DISTAL' );  %#ok<ASGLU>
    [id_cutedge_i,id_cutedge_p,id_cutedge_a,id_cutedge_l] = getMgenLevels( m, 'ID_CUTEDGE' );  %#ok<ASGLU>
    [id_boundary_i,id_boundary_p,id_boundary_a,id_boundary_l] = getMgenLevels( m, 'ID_BOUNDARY' );  %#ok<ASGLU>
    [id_highgro_i,id_highgro_p,id_highgro_a,id_highgro_l] = getMgenLevels( m, 'ID_HIGHGRO' );  %#ok<ASGLU>
    [id_highpar_i,id_highpar_p,id_highpar_a,id_highpar_l] = getMgenLevels( m, 'ID_HIGHPAR' );  %#ok<ASGLU>
    [id_highgrowlate_i,id_highgrowlate_p,id_highgrowlate_a,id_highgrowlate_l] = getMgenLevels( m, 'ID_HIGHGROWLATE' );  %#ok<ASGLU>
    [id_com_i,id_com_p,id_com_a,id_com_l] = getMgenLevels( m, 'ID_COM' );  %#ok<ASGLU>
    [id_commid_i,id_commid_p,id_commid_a,id_commid_l] = getMgenLevels( m, 'ID_COMMID' );  %#ok<ASGLU>
    [id_areathresh_i,id_areathresh_p,id_areathresh_a,id_areathresh_l] = getMgenLevels( m, 'ID_AREATHRESH' );  %#ok<ASGLU>
    [id_rim_i,id_rim_p,id_rim_a,id_rim_l] = getMgenLevels( m, 'ID_RIM' );  %#ok<ASGLU>
    [id_cdiv_i,id_cdiv_p,id_cdiv_a,id_cdiv_l] = getMgenLevels( m, 'ID_CDIV' );  %#ok<ASGLU>
    [id_early_i,id_early_p,id_early_a,id_early_l] = getMgenLevels( m, 'ID_EARLY' );  %#ok<ASGLU>
    [v_areathreshvis_i,v_areathreshvis_p,v_areathreshvis_a,v_areathreshvis_l] = getMgenLevels( m, 'V_AREATHRESHVIS' );  %#ok<ASGLU>
    [id_earlygrowth_i,id_earlygrowth_p,id_earlygrowth_a,id_earlygrowth_l] = getMgenLevels( m, 'ID_EARLYGROWTH' );  %#ok<ASGLU>
    [id_idiv_i,id_idiv_p,id_idiv_a,id_idiv_l] = getMgenLevels( m, 'ID_IDIV' );  %#ok<ASGLU>
    [v_justlam_i,v_justlam_p,v_justlam_a,v_justlam_l] = getMgenLevels( m, 'V_JUSTLAM' );  %#ok<ASGLU>
    [v_justpet_i,v_justpet_p,v_justpet_a,v_justpet_l] = getMgenLevels( m, 'V_JUSTPET' );  %#ok<ASGLU>
    [v_pgradthresh_i,v_pgradthresh_p,v_pgradthresh_a,v_pgradthresh_l] = getMgenLevels( m, 'V_PGRADTHRESH' );  %#ok<ASGLU>
    [id_mdgf_trunk_m3_i,id_mdgf_trunk_m3_p,id_mdgf_trunk_m3_a,id_mdgf_trunk_m3_l] = getMgenLevels( m, 'ID_MDGF_TRUNK_M3' );  %#ok<ASGLU>
    [s_mdgf_i,s_mdgf_p,s_mdgf_a,s_mdgf_l] = getMgenLevels( m, 'S_MDGF' );  %#ok<ASGLU>
    [id_mdgf_trunk_m2_i,id_mdgf_trunk_m2_p,id_mdgf_trunk_m2_a,id_mdgf_trunk_m2_l] = getMgenLevels( m, 'ID_MDGF_TRUNK_M2' );  %#ok<ASGLU>
    [id_mdgf_trunk_m1_i,id_mdgf_trunk_m1_p,id_mdgf_trunk_m1_a,id_mdgf_trunk_m1_l] = getMgenLevels( m, 'ID_MDGF_TRUNK_M1' );  %#ok<ASGLU>
    [v_justmid_i,v_justmid_p,v_justmid_a,v_justmid_l] = getMgenLevels( m, 'V_JUSTMID' );  %#ok<ASGLU>
    [r_anisotropy_i,r_anisotropy_p,r_anisotropy_a,r_anisotropy_l] = getMgenLevels( m, 'R_ANISOTROPY' );  %#ok<ASGLU>
    [v_midline_i,v_midline_p,v_midline_a,v_midline_l] = getMgenLevels( m, 'V_MIDLINE' );  %#ok<ASGLU>
    [v_cellareainhib_i,v_cellareainhib_p,v_cellareainhib_a,v_cellareainhib_l] = getMgenLevels( m, 'V_CELLAREAINHIB' );  %#ok<ASGLU>
    [s_area1_i,s_area1_p,s_area1_a,s_area1_l] = getMgenLevels( m, 'S_AREA1' );  %#ok<ASGLU>
    [s_area2_i,s_area2_p,s_area2_a,s_area2_l] = getMgenLevels( m, 'S_AREA2' );  %#ok<ASGLU>
    [s_distal_i,s_distal_p,s_distal_a,s_distal_l] = getMgenLevels( m, 'S_DISTAL' );  %#ok<ASGLU>
    [id_verydistal_i,id_verydistal_p,id_verydistal_a,id_verydistal_l] = getMgenLevels( m, 'ID_VERYDISTAL' );  %#ok<ASGLU>
    [id_moredistal_i,id_moredistal_p,id_moredistal_a,id_moredistal_l] = getMgenLevels( m, 'ID_MOREDISTAL' );  %#ok<ASGLU>
    [v_kx_i,v_kx_p,v_kx_a,v_kx_l] = getMgenLevels( m, 'V_KX' );  %#ok<ASGLU>
    [v_ky_i,v_ky_p,v_ky_a,v_ky_l] = getMgenLevels( m, 'V_KY' );  %#ok<ASGLU>
    [v_kz_i,v_kz_p,v_kz_a,v_kz_l] = getMgenLevels( m, 'V_KZ' );  %#ok<ASGLU>
    [c_split_i,c_split] = getCellFactorLevels( m, 'c_split' );
    [c_areathresh_i,c_areathresh] = getCellFactorLevels( m, 'c_areathresh' );
    [c_area_i,c_area] = getCellFactorLevels( m, 'c_area' );
    [c_splitrecord_i,c_splitrecord] = getCellFactorLevels( m, 'c_splitrecord' );
    [c_splitrecord2_i,c_splitrecord2] = getCellFactorLevels( m, 'c_splitrecord2' );
    [c_karea_i,c_karea] = getCellFactorLevels( m, 'c_karea' );
    [c_splitrecord3_i,c_splitrecord3] = getCellFactorLevels( m, 'c_splitrecord3' );
    [c_areavar_i,c_areavar] = getCellFactorLevels( m, 'c_areavar' );
    [c_justpet_i,c_justpet] = getCellFactorLevels( m, 'c_justpet' );
    [c_justlam_i,c_justlam] = getCellFactorLevels( m, 'c_justlam' );
    [c_justmid_i,c_justmid] = getCellFactorLevels( m, 'c_justmid' );
    [c_justnotpet_i,c_justnotpet] = getCellFactorLevels( m, 'c_justnotpet' );
    [c_justdist_i,c_justdist] = getCellFactorLevels( m, 'c_justdist' );
    [c_justverydist_i,c_justverydist] = getCellFactorLevels( m, 'c_justverydist' );
    [c_neighbours_i,c_neighbours] = getCellFactorLevels( m, 'c_neighbours' );

% Mesh type: lobes
%            base: 0
%        cylinder: 0
%          height: 0.053
%           lobes: 1
%          radius: 0.042
%      randomness: 0
%           rings: 15
%          strips: 35
%         version: 1

%            Morphogen    Diffusion   Decay   Dilution   Mutant
%            --------------------------------------------------
%                KAPAR         ----    ----       ----     ----
%                KAPER         ----    ----       ----     ----
%                KBPAR         ----    ----       ----     ----
%                KBPER         ----    ----       ----     ----
%                 KNOR         ----    ----       ----     ----
%            POLARISER         0.01     0.1       ----     ----
%            STRAINRET         ----    ----       ----     ----
%               ARREST         ----    ----       ----     ----
%               V_LEAF         ----    ----       ----     ----
%           ID_PROXORG         ----    ----       ----     ----
%               ID_INC         ----    ----       ----     ----
%             ID_PGRAD         ----    ----       ----     ----
%               ID_MID         ----    ----       ----     ----
%               ID_LAM         ----    ----       ----     ----
%                KAREA         ----    ----       ----     ----
%              ID_LATE         ----    ----       ----     ----
%            ID_DISTAL         ----    ----       ----     ----
%           ID_CUTEDGE         ----    ----       ----     ----
%          ID_BOUNDARY         ----    ----       ----     ----
%           ID_HIGHGRO         ----    ----       ----     ----
%           ID_HIGHPAR         ----    ----       ----     ----
%      ID_HIGHGROWLATE         ----    ----       ----     ----
%               ID_COM         ----    ----       ----     ----
%            ID_COMMID         ----    ----       ----     ----
%        ID_AREATHRESH         ----    ----       ----     ----
%               ID_RIM         ----    ----       ----     ----
%              ID_CDIV         ----    ----       ----     ----
%             ID_EARLY         ----    ----       ----     ----
%      V_AREATHRESHVIS         ----    ----       ----     ----
%       ID_EARLYGROWTH         ----    ----       ----     ----
%              ID_IDIV         ----    ----       ----     ----
%            V_JUSTLAM         ----    ----       ----     ----
%            V_JUSTPET         ----    ----       ----     ----
%        V_PGRADTHRESH         ----    ----       ----     ----
%     ID_MDGF_TRUNK_M3         ----    ----       ----     ----
%               S_MDGF         0.01     0.3       ----     ----
%     ID_MDGF_TRUNK_M2         ----    ----       ----     ----
%     ID_MDGF_TRUNK_M1         ----    ----       ----     ----
%            V_JUSTMID         ----    ----       ----     ----
%         R_ANISOTROPY         ----    ----       ----     ----
%            V_MIDLINE         ----    ----       ----     ----
%      V_CELLAREAINHIB         ----    ----       ----     ----
%              S_AREA1         ----    ----       ----     ----
%              S_AREA2         ----    ----       ----     ----
%             S_DISTAL         ----    ----       ----     ----
%        ID_VERYDISTAL         ----    ----       ----     ----
%        ID_MOREDISTAL         ----    ----       ----     ----
%                 V_KX         ----    ----       ----     ----
%                 V_KY         ----    ----       ----     ----
%                 V_KZ         ----    ----       ----     ----


%%% USER CODE: MORPHOGEN INTERACTIONS

% Get current area, in order to calculate relative growth rate in the
% preplotproc.
m.userdata.currentTissueArea = sum(m.cellareas);

%% model dynamics
if (Steps(m)==0) && m.globalDynamicProps.doinit  % Initialisation code.
    if PLOTTING_LINEAGE
        m = leaf_plotoptions( m, 'cellbodyvalue', 'c_splitrecord3' ); % The cell factor we use to record cells that split during the interval.
    end
        
    ZEROCOLOR = [1 1 1];  % Use [0 0 0] if you want them black.
    LIGHTGREEN = [0.31 0.45 0.13];
    DARKGREEN = [0.47 0.67 0.19];
    m = leaf_setcellcolorinfo( m, 'factor', 'c_splitrecord3', 'mode', 'indexed', 'colors', [ZEROCOLOR;DARKGREEN;LIGHTGREEN] );
%     m = leaf_setcellcolorinfo( m, 'factor', 'c_area', 'mode', 'rainbow', 'autorange', false, 'range', [0 OPTIONS.c_area_scaling/1e6] );
    
    % Set up the user data fields to handle recording of cell splits.    
    % The splits field is intended to be different at different times,
    % so it goes in the dynamic user data.
    if m.userdatastatic.collectsplitdata
        m = leaf_set_userdata( m, ...
            'splits', [] ); % Records information about cell splittings as they happen during the interval.
        % These items conceptually apply to the whole project rather than
        % separately to each point in time, so they are stored in the
        % static user data.
        m = leaf_set_userdatastatic( m, ...
            'firsttime', [115 132 140 147 156 164 178], ... % Times at which recording of cell splits begins.
            'lasttime',  [132 140 147 156 164 178 205], ... % Times at which recording of cell splits ends.
            'allsplits', cell(1,length(m.userdatastatic.firsttime)), ... % At the end of each interval, stores a list of all cells that split.
            'cellparent', [], ... % Records for every cell id the id of its parent, if any, otherwise 0.
            'celldaughters', [], ... % Records for every cell id, the ids of its two daughter cells, if it ever splits, otherwise 0.
            ... % (This is never used.)
            'celltimes', [], ... % Records for every cell id, its creation time and splitting time.  If it never splits, the splitting time is Inf.
            'splitareas', cell(1,length(m.userdatastatic.firsttime)), ... % For each interval, the areas of the cells that split
            ... % in that interval, at the moment they split.
            'unsplitareas', cell(1,length(m.userdatastatic.firsttime)), ... % For each interval, and each cell that did not split
            ... % in that interval, the areas of that cell at the
            ... % beginning and end of the interval.
            'cellareastartbyid', cell(1,length(m.userdatastatic.firsttime)), ... % For each interval, the area of every cell at the
            ... % beginning of the interval, indexed by cellid, the value
            ... % being zero for cells no longer existing.
            'cellareastart', cell(1,length(m.userdatastatic.firsttime)), ... % For each interval, the area of every cell at the
            ... % beginning of the interval.
            'cellareaend', cell(1,length(m.userdatastatic.firsttime)), ... % For each interval, the area of every cell at the
            ... % end of the interval.
            'cloneareas', cell(1,length(m.userdatastatic.firsttime)), ... % For each interval and each cell present at the start of that interval,
            ... % the area of its descendant clone at the end of the interval.
            'cellidstart', cell(1,length(m.userdatastatic.firsttime)), ... % For each interval, the id of every cell at the
            ... % beginning of the interval.
            'cellidend', cell(1,length(m.userdatastatic.firsttime)) ); % For each interval, the id of every cell at the
            ... % end of the interval.  (This is never used.)
            saveStaticPart( m );
    end
    
    % Initialise record of current relative growth rate.
    m.userdata.currentTissueRelGrowth = 0;
    
    % Initialise record of cell Size Feedback On Growth.
    % Turn on cell size feeback on growth for some models
    switch modelname
        case {'EPIMODEL-CELLS-DISTALGROWTH','SHIFT_LATE_CELLS_EARLY','SHIFT_LATE_CELLS_LATE', ...
                'CHANGE_GROWTH_RATE_CELLS_LOW','CHANGE_GROWTH_RATE_CELLS_HIGH', ...
                'CHANGE_CELL_DIVTHRESH_CELLS_LOW','CHANGE_CELL_DIVTHRESH_CELLS_HIGH'}
            m.userdata.cellSizeFeedbackOnGrowth = 1;
        otherwise
            m.userdata.cellSizeFeedbackOnGrowth = 0;
    end
    
    %set inter-cell space colour
    m = leaf_plotoptions(m, 'cellspacecolor', [1 1 1]);
    %set alpha of cellular morphogens
    m = leaf_plotoptions(m, 'bioAalpha', 1);
    
end
if isfield( OPTIONS,'c_area_scaling')
    m = leaf_setcellcolorinfo( m, 'factor', 'c_area', 'mode', 'rainbow', 'autorange', false, 'range', [0 OPTIONS.c_area_scaling/1e6] );
else
    m = leaf_setcellcolorinfo( m, 'factor', 'c_area', 'mode', 'rainbow', 'autorange', true );
end


if Steps(m) == 1 % a bit later to give the cluster the chance to catch up.
    
    % PROXORG and POLARISER (POL) setup
    proxorg_ind= m.nodes(:,2)<=min(m.nodes(:,2))+0.001;
    id_proxorg_p(proxorg_ind) = 1;
    P(proxorg_ind) = bpol;
    m = leaf_fix_mgen( m, polariser_i,'vertex',find(proxorg_ind),'fix', 1);
    m = leaf_mgen_absorption(m, polariser_i, 0.1);
    m = leaf_mgen_conductivity(m, polariser_i, 0.01);
    
    %% KRN FACTORS
    
    % Visulalisation Factors setup
    v_leaf_p(:) =1;
    
    %MIDLINE
    v_midline_p(:) = 0;
    v_midline_p((m.nodes(:,1) < 0.001) & (m.nodes(:,1) > -0.001)) = 1;
    
    %JUSTPET
    v_justpet_p(:) = 0;
    v_justpet_p(m.nodes(:,2)<(min(m.nodes(:,2))+0.011)) = 1;
    
    %PGRADTHRESH
    v_pgradthresh_p(:) = 0;
    v_pgradthresh_p(id_pgrad_p > PgradTHRESH ) = 1;
    
    % Identity Factors setup
    
    % PGRAD
    % set up a linear gradient
%     bpgrad = m.userdata.ranges.bpgrad.range(m.userdata.ranges.bpgrad.index);
    bpgrad = OPTIONS.bpgrad;
    id_pgrad_p(:) = max(m.nodes(:,2)) - m.nodes(:,2);
    id_pgrad_p(:) = id_pgrad_p/max(id_pgrad_p)*(1-bpgrad);
    id_pgrad_p(:) = id_pgrad_p+bpgrad;
    
    %BOUNDARY
    id_boundary_p(m.nodes(:,2)>(BOUNDRYPOS) & m.nodes(:,2)<(BOUNDRYPOS+0.002)) = 1;
    
    %MDGF
    s_mdgf_p(:) = 0;
    psig_region = (m.nodes(:,2)<(min(m.nodes(:,2))+0.013)) & ...
        (m.nodes(:,2)>(min(m.nodes(:,2))+0.009));
    m.userdata.psig_region = psig_region;
    s_mdgf_p(psig_region) = 1;
    m = leaf_fix_mgen( m, s_mdgf_i,'vertex',find(psig_region),'fix', 1);
    m = leaf_mgen_absorption(m, s_mdgf_i, 0.3);
    m = leaf_mgen_conductivity(m, s_mdgf_i, 0.01);
    
    %EARLYGROWTH - is on when LATE is increasing linear and off when it
    %goes exponential
    id_earlygrowth_p(:) = 1;
    
    %EARLY - this is the period of growth from the model start (66h) to the time
    %we have observations for TL02 (115h)
    id_early_p(:) = 1;
    
    %CDIV
    id_cdiv_p(:) = 1;
    
    %DISTAL
    s_distal_p(:) = 0;
    id_distal_p(:) = 0;
    id_verydistal_p(:) = 0;
    id_moredistal_p(:) = 0;
    
    %MDGF THRESHOLDS
    id_mdgf_trunk_m1_p(:) = 0;
    id_mdgf_trunk_m2_p(:) = 0;
    id_mdgf_trunk_m3_p(:) = 0;
        
    % LAMINA (LAM)
    id_lam_p(m.nodes(:,2)>(-0.038)& m.nodes(:,2)<0.005) = 1.5;
    m = leaf_mgen_conductivity(m, id_lam_i, 0.00082);
    
    % MIDVEIN (MID)
    id_mid_p(:) = 0;
    id_mid_p( (m.nodes(:,1)<0.025) & ...
        (m.nodes(:,1)>-0.025) & ...
        (m.nodes(:,2)< BOUNDRYPOS) ) = 1;
    
    %JUSTLAM
    v_justlam_p(:) = 1;
    v_justlam_p(id_mid_p > 0) = 0;
    v_justlam_p(m.nodes(:,2)<(min(m.nodes(:,2))+0.011)) = 0;
    
    %JUSTMID
    v_justmid_p(:) = id_mid_p(:);
    v_justmid_p(v_justpet_p > 0) = 0;
    
    % RK: Set up parameters for inhibition of growth by cell size.
    %in the lamina set the cell size thresh to 8000 and 10000
    s_area1_p(:) = 4000e-6;
    s_area2_p(:) = 8000e-6;
    %in the mid vein set the cell size thresh to 18000 20000
    s_area1_p(id_mid_p == 1) = 18000e-6;
    s_area2_p(id_mid_p == 1) = 20000e-6;
    
    m = leaf_set_userdatastatic( m, ...
        ... % 'cellareainhibparams', inhibparams, 'cellareaparam1', a, 'cellareaparam2', b, ...
        'growthsmallcells', 1, ...  % How much to scale the growth of cells smaller than s_area1_p.
        'growthlargecells', 0 );  % How much to scale the growth of cells larger than s_area2_p.
    
    % Storing these parameters in the static user data means that they can
    % be set once and will be effective always.  Even a stage file that has
    % never been recomputed will pick up the values set here when it is
    % loaded. This presumes that these values are never changed during a
    % run. If they do need to take different values during a single run,
    % then the dynamic user data (i.e. m.userdata) must be used.
    
elseif realtime > 67 && realtime < 68+dt
    
    % fixing LAM distribution and resetting LAM to a low value in
    % the petiole
    m = leaf_mgen_conductivity(m, id_lam_i, 0);
    id_lam_p(m.nodes(:,2)<(min(m.nodes(:,2))+0.011)) = 0.4;

    
    %% Biological Growth
elseif realtime >= GROWTH_START_TIME
    
    % read in the parameters and simplify the parameter names
%     ppgrad = m.userdata.ranges.ppgrad.range(m.userdata.ranges.ppgrad.index)
%     plam = m.userdata.ranges.plam.range(m.userdata.ranges.plam.index)
%     glate = m.userdata.ranges.glate.range(m.userdata.ranges.glate.index)
%     plate = m.userdata.ranges.plate.range(m.userdata.ranges.plate.index)
%     hlate = m.userdata.ranges.hlate.range(m.userdata.ranges.hlate.index)
%     hmid = m.userdata.ranges.hmid.range(m.userdata.ranges.hmid.index)
%     pmdgf = m.userdata.ranges.pmdgf.range(m.userdata.ranges.pmdgf.index)
    
    ppgrad = OPTIONS.ppgrad;
    plam = OPTIONS.plam;
    glate = OPTIONS.glate;
    plate = OPTIONS.plate;
    hlate = OPTIONS.hlate;
    hmid = OPTIONS.hmid;
    pmdgf = OPTIONS.pmdgf;

    
    %% GRN
    
    %DISTAL
    DISTAL_PGRAD_THRESH = 0.64;
    VERYDISTAL_PGRAD_THRESH = 0.38;
    MOREDISTAL_PGRAD_THRESH = 0.577;
    if atTime( physiotime, 112, physiotimestep )
        fprintf( 1, 'At time point %g, phys time %g, step %g.\n', 112, physiotime, physiotimestep );
        
        id_moredistal_p(id_pgrad_p < MOREDISTAL_PGRAD_THRESH) = 1;
        id_verydistal_p(id_pgrad_p < VERYDISTAL_PGRAD_THRESH) = 1;    
        
        %OVERWRITING DISTAL
        id_distal_p(id_pgrad_p < VERYDISTAL_PGRAD_THRESH) = 1;        
        
        s_distal_p = id_distal_p;
        m = leaf_fix_mgen( m, s_distal_i,'vertex',find(id_pgrad_p < VERYDISTAL_PGRAD_THRESH),'fix', 1);
        m = leaf_mgen_absorption(m, s_distal_i, 0.01);
        m = leaf_mgen_conductivity(m, s_distal_i, 0.0001);
    end
    if atTime( physiotime, 115, physiotimestep )
        fprintf( 1, 'At time point %g, phys time %g, step %g.\n', 115, physiotime, physiotimestep );
        
        m = leaf_mgen_absorption(m, s_distal_i, 0);
        m = leaf_mgen_conductivity(m, s_distal_i, 0);
    end
        
    switch modelname
        case { 'SPCH_IN_CHAMBER', 'SPCH_ON_PLATES', 'SPCH_ON_PLATES_NO_MID_COMP' }
            % The SPCH mutants are identical to EPIMODEL except for the
            % parameters defined here.
            
            SPCH_TWEAK = 1;  % To increase the inhibition of growth by id_late in the SPCH mutant.
            EXP_TURNON_SPCH_DELAY = 0;  % To increase the duration of the linear phase of id_late.
            EXP_TURNON_SLOWING = 0;
            EARLY_GROWTH_TAPER_TIME = 24;
        otherwise
             % Original model
            SPCH_TWEAK = 1;
            EXP_TURNON_SPCH_DELAY = 0;
            EXP_TURNON_SLOWING = 0;
            EARLY_GROWTH_TAPER_TIME = 24;
    end
    

    %LATE
    LINEAR_LATE_TURNON_TIME = 148;
    EXP_LATE_TURNON_TIME = 189 + EXP_TURNON_SPCH_DELAY;
    
    %if we are running the LATE mutant we should shift the onset of LATE
    switch modelname 
        case {'SHIFT_LATE_FACTOR_EARLY','SHIFT_LATE_CELLS_EARLY'}
            LATE_SHIFT_TIME = -6; %advance late for the SHIFT_LATE models
            LINEAR_LATE_TURNON_TIME = LINEAR_LATE_TURNON_TIME + LATE_SHIFT_TIME;
            EXP_LATE_TURNON_TIME = EXP_LATE_TURNON_TIME + LATE_SHIFT_TIME;
        case {'SHIFT_LATE_FACTOR_LATE','SHIFT_LATE_CELLS_LATE'}
            LATE_SHIFT_TIME = 6; %delay late for the SHIFT_LATE models
            LINEAR_LATE_TURNON_TIME = LINEAR_LATE_TURNON_TIME + LATE_SHIFT_TIME;
            EXP_LATE_TURNON_TIME = EXP_LATE_TURNON_TIME + LATE_SHIFT_TIME;
        otherwise
            LATE_SHIFT_TIME = 0;
    end
    
    %set late to increase linearly if physiotime > LINEAR_LATE_TURNON_TIME or
    %exponentially if physiotime >= EXP_LATE_TURNON_TIME
    if (physiotime >= LINEAR_LATE_TURNON_TIME) && (physiotime < EXP_LATE_TURNON_TIME)
        id_late_p = glate*(physiotime - LINEAR_LATE_TURNON_TIME);
    elseif physiotime >= EXP_LATE_TURNON_TIME
        %value of late at 205h is 0.2736 - value of late at 189h is 0.192
        DELTA_T = EXP_LATE_TURNON_TIME - LINEAR_LATE_TURNON_TIME + EXP_TURNON_SLOWING;
        B = 1/DELTA_T;
        A = DELTA_T*exp(-B*EXP_LATE_TURNON_TIME);
        id_late_p = glate*(A*exp(B*physiotime) - EXP_TURNON_SLOWING);
    end
    
    %EARLYGROWTH
    %EARLYGROWTH is on (1) during the stage of the model growth where LATE is
    %increasing linearly. When LATE starts increasing exponential then
    %EARLYGROWTH starts decreasing linearly to (0) over a time of
    %EARLY_GROWTH_TAPER_TIME (usual value is 24h).
    decrease_per_step = 1 / EARLY_GROWTH_TAPER_TIME;
    if physiotime >= EXP_LATE_TURNON_TIME;
        id_earlygrowth_p = id_earlygrowth_p - decrease_per_step;
        if sum(id_earlygrowth_p(:) <= 0)
            id_earlygrowth_p(:) = 0;
        end
    end      
    
    %Early off at 114h - early is the cell division early stage before
    %which all cells are dividing at a set threshold size.
    EARLY_OFF_TIME = 114;  
    if physiotime >= EARLY_OFF_TIME;
        id_early_p(:) = 0;
    end
    
    %% KRN
    %MDGF_TRUNK_M1
    mdgf_Thresh_M1 = 0.184; %This is M1=0.184=400um in the sub-epimodel
    id_mdgf_trunk_m1_p(:) = s_mdgf_p;
    ptr1 = id_mdgf_trunk_m1_p >= mdgf_Thresh_M1;
    id_mdgf_trunk_m1_p(ptr1) = mdgf_Thresh_M1;
    
    %MDGF_TRUNK_M2
    mdgf_Thresh_M2 = 0.295; %M2=0.295=300um in the epimodel
    id_mdgf_trunk_m2_p(:) = s_mdgf_p;
    ptr1 = id_mdgf_trunk_m2_p >= mdgf_Thresh_M2;
    id_mdgf_trunk_m2_p(ptr1) = mdgf_Thresh_M2;
    
    %MDGF_TRUNK_M3
    mdgf_Thresh_M3 = 0.51; %This is M3=0.51=150um in the epimodel
    id_mdgf_trunk_m3_p(:) = s_mdgf_p;
    ptr1 = id_mdgf_trunk_m3_p >= mdgf_Thresh_M3;
    id_mdgf_trunk_m3_p(ptr1) = mdgf_Thresh_M3;
    
    
    %AREATHRESH
    areaMod = id_mdgf_trunk_m3_p;
    ptr = areaMod >= mdgf_Thresh_M3;
    areaMod(ptr) = mdgf_Thresh_M3;
    switch modelname
        case{'SUBEPIMODEL'}
            ptr = areaMod <= mdgf_Thresh_M1;
            areaMod(ptr) = mdgf_Thresh_M1;
            mdgf_Thresh_M1M2 = mdgf_Thresh_M1;
            
        otherwise
            ptr = areaMod <= mdgf_Thresh_M2;
            areaMod(ptr) = mdgf_Thresh_M2;            
            mdgf_Thresh_M1M2 = mdgf_Thresh_M2;
    end
    
    %set the cell division thresholds here
    switch modelname
        case{'SUBEPIMODEL'}
            smallCellDivsionAreaThresh = 150;
            largeCellDivisionAreaThresh = 300;
            midVeinCellDivisionAreaThresh = 500;
        otherwise
            smallCellDivsionAreaThresh = 150; 
            largeCellDivisionAreaThresh = 300; 
            midVeinCellDivisionAreaThresh = 500; 
    end
    
    %if we are running the cell division threshold mutant then change the
    %division thresholds here
    switch modelname
        case {'CHANGE_CELL_DIVTHRESH_FACTOR_LOW','CHANGE_CELL_DIVTHRESH_CELLS_LOW'}
            CELL_DIV_THRESH_CHANGE = -85;
            smallCellDivsionAreaThresh = smallCellDivsionAreaThresh + CELL_DIV_THRESH_CHANGE;
            largeCellDivisionAreaThresh = largeCellDivisionAreaThresh + CELL_DIV_THRESH_CHANGE;
            midVeinCellDivisionAreaThresh = midVeinCellDivisionAreaThresh + CELL_DIV_THRESH_CHANGE;
        case {'CHANGE_CELL_DIVTHRESH_FACTOR_HIGH','CHANGE_CELL_DIVTHRESH_CELLS_HIGH'}
            CELL_DIV_THRESH_CHANGE = 85;
            smallCellDivsionAreaThresh = smallCellDivsionAreaThresh + CELL_DIV_THRESH_CHANGE;
            largeCellDivisionAreaThresh = largeCellDivisionAreaThresh + CELL_DIV_THRESH_CHANGE;
            midVeinCellDivisionAreaThresh = midVeinCellDivisionAreaThresh + CELL_DIV_THRESH_CHANGE;
    end
    
    %this sets the gradiated pattern of cell divisions in the lamina. All
    %cells within ~150um from lp-boundry divide at
    %smallCellDivsionAreaThresh. Cells at or above 300um(EPIMODEL) or
    %400um(SUBEPIMODEL) will divide at largeCellDivisionAreaThresh. And
    %cells inbetween will divide at a threshold set linearly inbetween
    v_areathreshvis_p(:) = ((largeCellDivisionAreaThresh-smallCellDivsionAreaThresh).*(areaMod - mdgf_Thresh_M3) ./ ...
        (mdgf_Thresh_M1M2 - mdgf_Thresh_M3)) + smallCellDivsionAreaThresh;
    
    %the midvein cell division threshold is set here. The subepimodel cell division threshold in the midvein is modulated
    %with pgrad so that the cell division threshold decreases
    %proximal-distal.
    switch modelname
        case{'SUBEPIMODEL'}
            v_areathreshvis_p(id_mid_p > 0) = midVeinCellDivisionAreaThresh .* id_pgrad_p(id_mid_p > 0) .* 1.5;            
        otherwise
            v_areathreshvis_p(id_mid_p > 0) = midVeinCellDivisionAreaThresh;
    end
    
    %if physiotime is < 114 then all cells everywhere will divide when they
    %reach the smallCellDivsionAreaThresh threshold size
    %if physiotime is >= 114 (at this point we now have cell division
    %observations) then 'turn on' the gradiated pattern of division as set above in v_areathreshvis_p    
    if physiotime >= EARLY_OFF_TIME    
        id_areathresh_p(:) = v_areathreshvis_p .* 1e-6;
    else
        id_areathresh_p(:) = smallCellDivsionAreaThresh * 1e-6;
    end
    
    %cells will stop dividing at CELL_DIVISION_OFF_TIME, midvein cells in
    %the epimodel will turn off earlier at EPI_MODEL_MID_CELL_DIVISION_OFF_TIME
    %CELL_DIVISION_OFF_TIME = 183;
    CELL_DIVISION_OFF_LATE_THRESHOLD = 0.1680; %the value of LATE at 183h
    %EPI_MODEL_MID_CELL_DIVISION_OFF_TIME = 164;
    EPI_MODEL_MID_CELL_DIVISION_OFF_LATE_THRESHOLD = 0.0768; %the value of LATE at 164h
    
    %define the CDIV region (in which cells may divide).
    switch modelname
        %in the subepimodel CDIV extends to ~400um based on an MDGF
        %threshold mdgf_Thresh_M1
        case{'SUBEPIMODEL'} 
            %CDIV
            id_cdiv_p(:) = 0;
            cdivZone =  (id_mdgf_trunk_m1_p(:) >= mdgf_Thresh_M1);
            id_cdiv_p( cdivZone) = 1;
            %turn off cell divisions at a specified time
            %physiotime -> late
            %if physiotime >= CELL_DIVISION_OFF_TIME + LATE_SHIFT_TIME
            if mean(id_late_p(:)) >= CELL_DIVISION_OFF_LATE_THRESHOLD
                id_cdiv_p(:) = 0;
            end
        %in the EPIMODEL CDIV is set with a pgrad threshold and extends with the growing tissue until it reaches
        %~300um (mdgf_Thresh_M2)
        otherwise
            %CDIV
            id_cdiv_p(:) = 0;
            cdivZone =  (id_mdgf_trunk_m2_p(:) >= mdgf_Thresh_M2);
            id_cdiv_p( cdivZone) = 1;            
            noMoreDivisionsZone = id_pgrad_p < PgradTHRESH .* (1-id_early_p(:));             
            id_cdiv_p(noMoreDivisionsZone == 1) = 0;            
            %turn off cell divisions at a specified time
            %physiotime -> late
            %if physiotime >= EPI_MODEL_MID_CELL_DIVISION_OFF_TIME + LATE_SHIFT_TIME
            if (mean(id_late_p(:)) >= EPI_MODEL_MID_CELL_DIVISION_OFF_LATE_THRESHOLD) ...
                    && ~strcmp( modelname, 'SPCH_ON_PLATES_NO_MID_COMP' ) ...
                    && ~strcmp( modelname, 'EPIMODEL_NO_MID_COMP' )
                % The excluded models do not turn off cell division on the midline.
                id_cdiv_p(id_mid_p == 1) = 0;
            end
            %physiotime -> late
            %if physiotime >= CELL_DIVISION_OFF_TIME + LATE_SHIFT_TIME
            if mean(id_late_p(:)) >= CELL_DIVISION_OFF_LATE_THRESHOLD
                id_cdiv_p(:) = 0;
            end
            
    end
    
    switch modelname
        case{'EPIMODEL-FACTOR-DISTALGROWTH', 'SHIFT_LATE_FACTOR_EARLY', 'SHIFT_LATE_FACTOR_LATE', ...
                'CHANGE_GROWTH_RATE_FACTOR_LOW','CHANGE_GROWTH_RATE_FACTOR_HIGH', ...
                'CHANGE_CELL_DIVTHRESH_FACTOR_LOW','CHANGE_CELL_DIVTHRESH_FACTOR_HIGH'}
            % KRN using PGRAD factor to enhance reduction of distal growth            
            kapar_p(:) = ppgrad.*id_pgrad_l... Kpar promotion by PGRAD (parallel to POL gradient)
                .*inh(hlate,id_late_l .* inh(1.5,(1-id_earlygrowth_l)))...
                .*inh(2,(1-id_lam_l) .* (1-id_earlygrowth_l))...
                .*inh(4,  (1-id_earlygrowth_l) .* (1-id_pgrad_p));
            kbpar_p(:) = kapar_p;
            kaper_p(:) = plam*id_lam_l...Kper promotion by LAM (perpendicular to POL gradient)
                .*pro(plate,id_late_l .* id_earlygrowth_l)...Kper promotion by LATE
                .*inh(1.2, id_late_l .* (1-id_earlygrowth_l))... %
                .*inh(4,(1-id_earlygrowth_l) .* (1-id_pgrad_p))...
                .*inh(hmid, id_mid_l)...
                .*pro(pmdgf,id_mdgf_trunk_m2_p);
            kbper_p(:) = kaper_p(:);

            
        otherwise
            % KRN using linear decreasing EARLYGROWTH to decrease distal
            % growth. Some models use Cell Size Feedback with this KRN to
            % limit growth of large cells.
            
            kapar_p(:) = ppgrad.*id_pgrad_l... Kpar promotion by PGRAD (parallel to POL gradient)
                .*inh(hlate,id_late_l)... %
                .*inh(0.24*SPCH_TWEAK,id_late_l .* (1-id_earlygrowth_l));
            kbpar_p(:) = kapar_p;
            kaper_p(:) = plam*id_lam_l...Kper promotion by LAM (perpendicular to POL gradient)
                .*pro(plate,id_late_l .* id_earlygrowth_l)...Kper promotion by LATE
                .*inh(2.8*SPCH_TWEAK,id_late_l .* (1-id_earlygrowth_l))... %
                .*inh(hmid, id_mid_l)...
                .*pro(pmdgf,id_mdgf_trunk_m2_p);
            kbper_p(:) = kaper_p(:);
    end
    
    % Correct growth for physiological slowing and growth slowing.
    % For most models these ratios are 1 and so have no effect.
    kapar_p = kapar_p * (PHYSIO_TIME_RATIO*GROWTH_RATIO);
    kbpar_p = kapar_p;
    kaper_p = kaper_p * (PHYSIO_TIME_RATIO*GROWTH_RATIO);
    kbper_p = kaper_p;
    
    %if we are running the growth rate mutant change the global growth
    %rate. Stop growth at 178h for the subepimodel or 412h for the EPIMODEL.
    switch modelname
        case {'CHANGE_GROWTH_RATE_FACTOR_LOW','CHANGE_GROWTH_RATE_CELLS_LOW'}
            growthRateChange = 0.95;
            kapar_p = kapar_p * growthRateChange;
            kbpar_p = kapar_p;
            kaper_p = kaper_p * growthRateChange;
            kbper_p = kaper_p;
        case {'CHANGE_GROWTH_RATE_FACTOR_HIGH','CHANGE_GROWTH_RATE_CELLS_HIGH'}
            growthRateChange = 1.05;
            kapar_p = kapar_p * growthRateChange;
            kbpar_p = kapar_p;
            kaper_p = kaper_p * growthRateChange;
            kbper_p = kaper_p;
        case { 'SPCH_IN_CHAMBER', 'SPCH_ON_PLATES', 'SPCH_ON_PLATES_NO_MID_COMP' }
            growthRateChange = 1;
            if physiotime >= 148
                kapar_p = kapar_p * growthRateChange;
                kbpar_p = kapar_p;
                kaper_p = kaper_p * growthRateChange;
                kbper_p = kaper_p;
            end
        case {'SUBEPIMODEL'}
            if physiotime > 178
                kapar_p(:) = 0;
                kbpar_p(:) = 0;
                kaper_p(:) = 0;
                kbper_p(:) = 0;
            end
        otherwise
            if physiotime > 412
                kapar_p(:) = 0;
                kbpar_p(:) = 0;
                kaper_p(:) = 0;
                kbper_p(:) = 0;
            end
    end    
    
end


%% Model Decorations - optional
switch modeldeco
    
    case 'CELLSFORDIV'
        if Steps(m) == 0
            if hasNonemptySecondLayer( m )
                fprintf( 1, 'Cellular layer already present, not rebuilding it.\n' );
            else
                % Import real cells
                fprintf( 1, 'Building cellular layer from data file.\n' );

                % Rotation and rescaling parameters to project cells onto the canvas
                ang = pi; scale = 1000;

                % Load in the position of each vertex and rescale (from meters to mm)
                path = [fileparts(which(m.globalProps.modelname)),filesep,'EarlyClonesAdjusted'];
                point = load([path,filesep,'allpoints.mat']);
                points = point.all_points(:,:,1)*scale;
                pointsC1 = points(:,1)-mean(points(:,1))-0.00075;
                pointsC2 = points(:,2)-mean(points(:,2))+0.004;

                % Rotate points
                R = [cos(ang) sin(ang); -sin(ang) cos(ang)];
                pointsC = [pointsC1 pointsC2]*R;
                points(:,3) = zeros(length(points),1); % to get z-axis component
                pointsC = [pointsC points(:,3)];

                % Load in the cell vertices
                cells = load([path,filesep,'cells.mat']);
                cells = cells.cells;
                cells = cells(1:116);

                l = 0;
                cellcount = 0:(length(cells)-1);
                for i=cellcount+1
                    l = l+1;
                    cellpts{l} = cells(i).pts;
                end

                % Assign green color to cells
                cols = repmat([0 1 0], length(cellcount), 1);
                cols = cols./255;
                cols = cols(1:length(cellcount),:);
            
                % Make the cells
                m = leaf_makesecondlayer( m, ...
                    'vertexdata',pointsC,...
                    'celldata',cellpts,...
                    'cellcolors',cols);
            end
            % Set up tables to track descendant relationships.
            m.userdata.cellid = (1:length(m.secondlayer.cells))';
            m.userdata.cellindex = m.userdata.cellid;
            m.userdatastatic.cellparent = zeros( length(m.secondlayer.cells), 1 );
            m.userdatastatic.celldaughters = zeros( length(m.secondlayer.cells), 2 );
            m.userdatastatic.celltimes = [ realtime + zeros( length(m.secondlayer.cells), 1 ), inf( length(m.secondlayer.cells), 1 ) ];
        end
        if (cellIntercellSpaces) %do we plot intercell spaces
            if(strcmp(modelname,'SUBEPIMODEL') && physiotime < 178)
                doIntercellularSpaces();
            end
        end
        
    otherwise
        % the model will not have a second layer.
        
end


%% cellular morphogens
% After the second layer is created, the c_**** variables are stale and
% need to be set to the current values.
[c_split_i,c_split] = getCellFactorLevels( m, 'c_split' );
[c_areathresh_i,c_areathresh] = getCellFactorLevels( m, 'c_areathresh' );
[c_areavar_i,c_areavar] = getCellFactorLevels( m, 'c_areavar' );
[c_area_i,c_area] = getCellFactorLevels( m, 'c_area' );
[c_splitrecord_i,c_splitrecord] = getCellFactorLevels( m, 'c_splitrecord' );
[c_splitrecord2_i,c_splitrecord2] = getCellFactorLevels( m, 'c_splitrecord2' );
[c_splitrecord3_i,c_splitrecord3] = getCellFactorLevels( m, 'c_splitrecord3' );
[c_karea_i,c_karea] = getCellFactorLevels( m, 'c_karea' );
if Steps(m)==0
    m.userdata.areavar = 0.2; %noise std value
    c_areavar = randn( size(c_areavar) ) * m.userdata.areavar;
end

%setup some indexing cellular morphogens
c_justpet = leaf_FEVertexToCellularValue( m, v_justpet_p, 'mode', 'cell');
c_justlam = leaf_FEVertexToCellularValue( m, v_justlam_p, 'mode', 'cell');
c_justmid = leaf_FEVertexToCellularValue( m, v_justmid_p, 'mode', 'cell');
c_justnotpet = leaf_FEVertexToCellularValue( m, ~v_justpet_p, 'mode', 'cell');
c_justdist = leaf_FEVertexToCellularValue( m, id_distal_p, 'mode', 'cell');
c_justverydist = leaf_FEVertexToCellularValue( m, id_verydistal_p, 'mode', 'cell');
c_justmoredist = leaf_FEVertexToCellularValue( m, id_moredistal_p, 'mode', 'cell');

%% Cell Layer setup and split rules
c_com = leaf_FEVertexToCellularValue( m, id_cdiv_p, 'mode', 'cell');
c_areathresh = leaf_FEVertexToCellularValue( m, id_areathresh_p, 'mode', 'cell');
c_area = m.secondlayer.cellarea;

if(celldivnoise) %if adding noise to areathresh
    c_split = (c_com > 0) & (m.secondlayer.cellarea > c_areathresh.*(1+c_areavar));
else
    c_split = (c_com > 0) & (m.secondlayer.cellarea > c_areathresh);
end

%% calculate area growthrate
% Specified karea
karea_p(:) = kapar_p+kaper_p;
c_karea = leaf_FEVertexToCellularValue( m, karea_p, 'mode', 'cell');

%add code to calc anisotropy here
actualgrowth = (m.outputs.actualstrain.A+m.outputs.actualstrain.B)/2;
[amounts,~] = tensorsToComponents( actualgrowth  );
Kmax = amounts(:,1);
Kmin = amounts(:,2);
anisotropyPerFE = 1 - Kmin./Kmax;
anisotropyPerVertex = perFEtoperVertex( m, anisotropyPerFE );
if(sum(isnan(anisotropyPerVertex) > 0))
    r_anisotropy_p(:) = 0;
else
    r_anisotropy_p(:) = anisotropyPerVertex;
end

%science model anisotropy bar plots
output = m.userdata.outputtimes(2:end);

if sum(ismember(output-24,realtime))
    m.userdata.oldpos = m.prismnodes;
    m.userdata.starttime = realtime;
    
elseif sum(ismember(output-1,realtime))
    displacements = m.prismnodes - m.userdata.oldpos;
    [growth,gf] = leaf_computeGrowthFromDisplacements( m, displacements, ...
        realtime - m.userdata.starttime,'axisorder', 'maxminnor', ...
        'anisotropythreshold', 0.05);
    
    % plot resultant areal growth rates over 24-h intervals.
    %m = leaf_plotoptions( m, 'pervertex',perFEtoperVertex(m,sum(growth(:,1:2),2)) ,'perelementaxes', gf(:,1,:), 'drawtensoraxes', true );
    
elseif sum(ismember(output,realtime)) && m.userdata.output ==1
    % Output - plot an image at high resolution
    path = fileparts(which(m.globalProps.modelname));
    [m,ok] = leaf_snapshot( m,[path,filesep,'snapshots',filesep,modelname,'_',modeltype,'.png'], 'resolution',[]); %
    
end

% RK
% m = leaf_plotoptions( m, 'morphogen', 'v_cellareainhib' );

%% Leaf Calculations

% data for csv file
% Leaf area
LeafArea = m.globalDynamicProps.currentArea;
% Leaf width
[minX, indminX] = min(m.nodes(:,1));
[maxX, indmaxX] = max(m.nodes(:,1));
LeafWidth = maxX + abs(minX);
% Leaf Length (tip to lam-pet boundry)
[minY, indminY] = min(m.nodes(:,2));
[maxY, indmaxY] = max(m.nodes(:,2));
petIdx = find(not(cellfun('isempty', strfind(m.mgenIndexToName,'V_JUSTPET'))));
somenodes = m.nodes(m.morphogens(:,petIdx)==1,:,:);
displace = max(somenodes(:,2)) - min(somenodes(:,2));
LeafLength = (maxY + abs(minY)) - abs(displace);
%plam
% plam = m.userdata.ranges.plam.range(m.userdata.ranges.plam.index);
plam = OPTIONS.plam;
%late
late = max(id_late_p);
%cells
numberOfCells = numel(m.secondlayer.cells);
meanCellArea = mean(m.secondlayer.cellarea);
numberOfCellInCDIV = numel(find(c_com > 0));
meanCellAreaInCDIV = mean(m.secondlayer.cellarea(c_com > 0));
cellAreaInCDIV = sum(m.secondlayer.cellarea(c_com > 0));
numberOfCellOutCDIV = numel(find(c_com == 0));
meanCellAreaOutCDIV = mean(m.secondlayer.cellarea(c_com == 0));
cellAreaOutCDIV = sum(m.secondlayer.cellarea(c_com == 0));

[centroids, areas, projections] = getCellAreasAndPositions( m,'morphogen', 'id_cdiv', 'threshold', 0.5 , 'mode', 'ave','axis',[1 0 0]);

dataToLog = [realtime LeafArea LeafWidth LeafLength late plam ...
    numberOfCells meanCellArea numberOfCellInCDIV meanCellAreaInCDIV ...
    cellAreaInCDIV numberOfCellOutCDIV meanCellAreaOutCDIV cellAreaOutCDIV];
pwd;
dlmwrite(['Extracted Growth', filesep, [modelname,'_LeafData.csv']], dataToLog, '-append')
if LeafWidth > 0.5
    xxxx = 1;
end

%get values of MDGF along midline
    midlineIdxs = v_midline_p(:) > 0;
[m.userdata.MDGFvaluesAlongMidline,m.userdata.positionAlongMidline] = find(s_mdgf_p(v_midline_p > 0));
%figure; plot(1:numel(EXTERNMESH.userdata.MDGFvaluesAlongMidline),EXTERNMESH.userdata.MDGFvaluesAlongMidline)

% Callback procedures must be installed on every iteration.
m = leaf_setproperty( m, ...
    'bioApostsplitproc', @bioApostsplitproc, ...
    'userpostiterateproc', @udstatictest_userpostiterateproc );
m = leaf_plotoptions( m, 'userpreplotproc', @udstatictest_userpreplotproc );

if m.userdatastatic.collectsplitdata
    % When recording of splits is about to begin, zero the relevant user data.
    atfirsttime = find( meshAtTime( m, m.userdatastatic.firsttime ), 1 );
    
    if ~isempty(atfirsttime)
        m.userdatastatic.allsplits{atfirsttime} = [];
        m.userdata.splits = [];
        saveStaticPart( m );
        fprintf( 1, 'Time %f, about to begin recording splits.\n', m.globalDynamicProps.currenttime );
    end
end

if(m.userdata.cellSizeFeedbackOnGrowth)
        %RK: Inhibition of growth by cell size.
        if isfield( m.userdatastatic, 'growthsmallcells' ) && ~isempty( m.userdatastatic.growthsmallcells )
            areaPerVertex = perCellToperFEVertex(m,c_area);
            szlo = s_area1_p;
            szhi = s_area2_p;
            glo = m.userdatastatic.growthsmallcells;
            ghi = m.userdatastatic.growthlargecells;
            v_cellareainhib_p = glo + (areaPerVertex - szlo) .* (ghi - glo) ./ (szhi - szlo);
            v_cellareainhib_p = max( ghi, min(glo,v_cellareainhib_p) );
            
            kapar_p = kapar_p .* v_cellareainhib_p;
            kaper_p = kaper_p .* v_cellareainhib_p;
            kbpar_p = kbpar_p .* v_cellareainhib_p;
            kbper_p = kbper_p .* v_cellareainhib_p;
        end
end

switch modelname
    case 'EPIMODEL_NO_POLARISER'
        P(:) = 0;
    otherwise
        % Nothing.
end
%%% END OF USER CODE: MORPHOGEN INTERACTIONS

%%% SECTION 3: INSTALLING MODIFIED VALUES BACK INTO MESH STRUCTURE
%%% AUTOMATICALLY GENERATED CODE: DO NOT EDIT.
    m.morphogens(:,polariser_i) = P;
    m.morphogens(:,kapar_i) = kapar_p;
    m.morphogens(:,kaper_i) = kaper_p;
    m.morphogens(:,kbpar_i) = kbpar_p;
    m.morphogens(:,kbper_i) = kbper_p;
    m.morphogens(:,knor_i) = knor_p;
    m.morphogens(:,strainret_i) = strainret_p;
    m.morphogens(:,arrest_i) = arrest_p;
    m.morphogens(:,v_leaf_i) = v_leaf_p;
    m.morphogens(:,id_proxorg_i) = id_proxorg_p;
    m.morphogens(:,id_inc_i) = id_inc_p;
    m.morphogens(:,id_pgrad_i) = id_pgrad_p;
    m.morphogens(:,id_mid_i) = id_mid_p;
    m.morphogens(:,id_lam_i) = id_lam_p;
    m.morphogens(:,karea_i) = karea_p;
    m.morphogens(:,id_late_i) = id_late_p;
    m.morphogens(:,id_distal_i) = id_distal_p;
    m.morphogens(:,id_cutedge_i) = id_cutedge_p;
    m.morphogens(:,id_boundary_i) = id_boundary_p;
    m.morphogens(:,id_highgro_i) = id_highgro_p;
    m.morphogens(:,id_highpar_i) = id_highpar_p;
    m.morphogens(:,id_highgrowlate_i) = id_highgrowlate_p;
    m.morphogens(:,id_com_i) = id_com_p;
    m.morphogens(:,id_commid_i) = id_commid_p;
    m.morphogens(:,id_areathresh_i) = id_areathresh_p;
    m.morphogens(:,id_rim_i) = id_rim_p;
    m.morphogens(:,id_cdiv_i) = id_cdiv_p;
    m.morphogens(:,id_early_i) = id_early_p;
    m.morphogens(:,v_areathreshvis_i) = v_areathreshvis_p;
    m.morphogens(:,id_earlygrowth_i) = id_earlygrowth_p;
    m.morphogens(:,id_idiv_i) = id_idiv_p;
    m.morphogens(:,v_justlam_i) = v_justlam_p;
    m.morphogens(:,v_justpet_i) = v_justpet_p;
    m.morphogens(:,v_pgradthresh_i) = v_pgradthresh_p;
    m.morphogens(:,id_mdgf_trunk_m3_i) = id_mdgf_trunk_m3_p;
    m.morphogens(:,s_mdgf_i) = s_mdgf_p;
    m.morphogens(:,id_mdgf_trunk_m2_i) = id_mdgf_trunk_m2_p;
    m.morphogens(:,id_mdgf_trunk_m1_i) = id_mdgf_trunk_m1_p;
    m.morphogens(:,v_justmid_i) = v_justmid_p;
    m.morphogens(:,r_anisotropy_i) = r_anisotropy_p;
    m.morphogens(:,v_midline_i) = v_midline_p;
    m.morphogens(:,v_cellareainhib_i) = v_cellareainhib_p;
    m.morphogens(:,s_area1_i) = s_area1_p;
    m.morphogens(:,s_area2_i) = s_area2_p;
    m.morphogens(:,s_distal_i) = s_distal_p;
    m.morphogens(:,id_verydistal_i) = id_verydistal_p;
    m.morphogens(:,id_moredistal_i) = id_moredistal_p;
    m.morphogens(:,v_kx_i) = v_kx_p;
    m.morphogens(:,v_ky_i) = v_ky_p;
    m.morphogens(:,v_kz_i) = v_kz_p;
    m.secondlayer.cellvalues(:,c_split_i) = c_split(:);
    m.secondlayer.cellvalues(:,c_areathresh_i) = c_areathresh(:);
    m.secondlayer.cellvalues(:,c_area_i) = c_area(:);
    m.secondlayer.cellvalues(:,c_splitrecord_i) = c_splitrecord(:);
    m.secondlayer.cellvalues(:,c_splitrecord2_i) = c_splitrecord2(:);
    m.secondlayer.cellvalues(:,c_karea_i) = c_karea(:);
    m.secondlayer.cellvalues(:,c_splitrecord3_i) = c_splitrecord3(:);
    m.secondlayer.cellvalues(:,c_areavar_i) = c_areavar(:);
    m.secondlayer.cellvalues(:,c_justpet_i) = c_justpet(:);
    m.secondlayer.cellvalues(:,c_justlam_i) = c_justlam(:);
    m.secondlayer.cellvalues(:,c_justmid_i) = c_justmid(:);
    m.secondlayer.cellvalues(:,c_justnotpet_i) = c_justnotpet(:);
    m.secondlayer.cellvalues(:,c_justdist_i) = c_justdist(:);
    m.secondlayer.cellvalues(:,c_justverydist_i) = c_justverydist(:);
    m.secondlayer.cellvalues(:,c_neighbours_i) = c_neighbours(:);

%%% USER CODE: FINALISATION

% In this section you may modify the mesh in any way whatsoever.


    function doIntercellularSpaces()
        if ~isempty(m.secondlayer.cells)
            % Choose a set of new vertexes.
            % Expand those vertexes and every other that already borders an
            % "air" space.
            m = leaf_setproperty( m, 'bioMinEdgeLength', 1e-9, 'bioSpacePullInRatio', 0.1 );
            %initiate = abs(physiotime - 124) < 0.1; % 87
            
            %at SUBEPIMODEL t=178 there are 1857 cells. at t = 120 there
            %are 146 cells.
            %146 * 0.2 = 29.2 holes
            %1857 * 0.35 = 649.95 holes
            %Therefore add 620 holes over 58h -> 10.68 holes per hour
            
            START_TIME = 120;
            INIT_SPACE_SIZE = 0.002;
            numOfStartingHoles = 30;
            
            if atTime( physiotime, START_TIME, physiotimestep )
                m = leaf_initiateICSpaces( m, 'number', numOfStartingHoles, 'abssize', INIT_SPACE_SIZE );
                m.userdata.numOfHoles = numOfStartingHoles;
            end
            
            times = [120 141 156 166];
            % Determine the amount of expansion of intercellular spaces, in
            % the form of the ratio of wall separation to distance of a
            % vertex from the centroid of its space.
            % This might be a single value to apply to all vertexes, or a
            % value per vertex.
            SPACE_EXPANSION_MODE = 'A';
            switch SPACE_EXPANSION_MODE
                case 'A'  % A single value for all vertexes, uniform in time
                    Gvs = 0.0125; %0.0325; %0.02;
                case 'B'  % A single value for all vertexes, but varying
                    % with time, being a constant fraction of overall
                    % leaf growth rate.
                    icFraction = 1;  % Wild guess, tweak this until you
                    % get a reasonable value for Gvs.
                    Gvs = icFraction * m.userdata.currentTissueRelGrowth;
                case 'C'  % Scale according to the areal growth rate of the
                    % patch we happened to measure.
                    patchRelGrowth = [ 0.03356 0.02195 0.01889 0.01964 ];  % Data from Sam's spreadsheet.
                    icFraction = 0.5;  % Wild guess, tweak this until you
                    % get a reasonable value for Gvs.
                    Gvs = icFraction * patchRelGrowth;
                case 'D'  % "Use vertex movement rates for each interval based on the patch relative to growth rate of the patch."
                    % I don't know what this means.
                    % Not implemented yet.
                    Gvs = 0;
                case 'E'  % "As C but relative to growth rate of the local growth rate."
                    % I don't know what this means.
                    % Not implemented yet.
                    Gvs = 0;
                case 'F'
                    % The method we used originally.
                    Gvs = [ 0.024229234 0.013264194 0.012025971 0.008136952 ];  % Data computed from Sam's spreadsheet.
                otherwise
                    % Shouldn't reach here.
                    Gvs = 0;
            end
            if numel(Gvs)==1
                Gvs = Gvs + zeros(size(times));
            end
            
            if physiotime >= START_TIME
                % Make new intercellular spaces.
                
                m = leaf_initiateICSpaces( m, 'number', 11, 'abssize', INIT_SPACE_SIZE );
                
                % Do wall separation.
                
                % Find which time interval we're in.
                [~,idx] = min(abs(times - realtime));
                
                % Select the current value of Gvs.
                % If positive, use it as the amount of wall separation.
                if Gvs(idx) > 0
                    m = leaf_growICSpaces( m, ...
                        'amount', Gvs(idx), ...
                        'amounttype', 'cellvertex', ...
                        'amountmode', 'reldist' );
                end
            end
        end
    end
%%% END OF USER CODE: FINALISATION

end


%%% USER CODE: SUBFUNCTIONS


function [m,splits,splitcentre,splitdirection] = bioApresplitproc( m, ci )
%   This is the pre-split procedure callback, which will be called when
%   GFtbox is considering which cells to split.  Its responsiblity is to
%   tell GFtbox which cells should split and in what way.
%
%   ci is an array of indexes of biological cells.
%
%   The results are;
%       splits: a boolean array the same size as ci, true if the
%               corresponding element of ci should be split.
%       splitcentre: an N*3 array, where N is the length of ci.  For each
%               cell that should be split, it contains the location of a
%               point that the new cell wall should pass through.  For a
%               cell that should not be split, set the corresponding row of
%               splitcentre to zero (it will be ignored anyway).
%       splitdirection: another N*3 array, this holds the direction
%               perpendicular to the new cell wall, represented by a unit
%               vector, and is zero for those cells that should not be
%               split.
%
% If it is desired to split a cell but to leave to GFtbox the decision
% about how to position the new cell wall, set the corresponding member of
% splits to true but set splitcentre and splitdirection to zero.

% [splits,asym] = doSplit( m, ci );
[c_split_i,splits] = getCellFactorLevels( m, 'c_split', ci );
splitcentre = zeros( length(ci), 3 );
splitdirection = zeros( length(ci), 3 );
for i=1:length(ci)
    if splits(i)
        cii = ci(i);
        % Split along shortest diameter through centre of cell.
        cvxs = m.secondlayer.cells(cii).vxs;
        cellcoords = m.secondlayer.cell3dcoords( cvxs, : );
        [splitcentre(i,:),splitdirection(i,:)] = splitShortestDiameter( cellcoords );
    end
end

% RK: update userdatastatic.splitareas and .unsplitareas.
% If we are in a time interval, record the area of every cell that has just
% been marked for splitting.  We record only the areas, not the cell
% indexes or any other information about the cell, since all we want to do
% with this data is plot a histogram of it.[No, actually we need to
% record the cell identities as well. E.g. use NaN for inapplicable cells.]
if m.userdatastatic.collectsplitdata
    [startInterval,~,atstart,atend] = findInterval( m );
    if atstart
        m.userdatastatic.splitareas{startInterval} = nan( length(m.secondlayer.cells), 1 );
        m.userdatastatic.unsplitareas{startInterval} = m.secondlayer.cellarea;
    end
    if ~isempty(startInterval)
        splitcells = ci(splits==1);
        foo = m.userdatastatic.splitareas{startInterval};
        foo(splitcells) = m.secondlayer.cellarea(splitcells);
        m.userdatastatic.splitareas{startInterval} = foo;
        if any(splits)==1
            xxxx = 1;
        end
        % The unsplitareas array must exclude cells that split.  [NOT! Need
        % to add up all daughter cells.]
        % Stupid Matlab syntax doesn't allow
        % m.userdatastatic.unsplitareas{startInterval}(splitcells) = 0.
        foo = m.userdatastatic.unsplitareas{startInterval};
        % Only interested in cells that were already present at the start
        % of the interval.
        splitcells = splitcells( splitcells <= size( foo, 1 ) );
        foo(splitcells) = 0;
        m.userdatastatic.unsplitareas{startInterval} = foo;
    end
    %     saveStaticPart( m );
end
end

function [splitcentre,splitdirection] = splitShortestDiameter( cellcoords )
% Split along shortest diameter through centre of cell.
% Only the X and Y components of cellcoords are used -- this is
% specialised for flat tissues.
[centroid,v,distance] = polysplitdiameter( cellcoords(:,[1 2]) );
splitcentre = [ (centroid + 0*randInCircle2( distance/4 )), 0 ];
splitdirection = [ -v(2), v(1), 0 ];
end

function m = bioApostsplitproc( m, ci, cei, newci, newcei, oe1, oe2, ne1, ne2, ne3 )
% Boilerplate section.
ne3check1 = m.secondlayer.cells(ci).edges(cei);
ne3check2 = m.secondlayer.cells(newci).edges(newcei);
ne = length( m.secondlayer.cells(ci).edges );
newne = length( m.secondlayer.cells(newci).edges );
ne2check = m.secondlayer.cells(ci).edges( mod(cei,ne)+1 );
ne1check = m.secondlayer.cells(newci).edges( mod(newcei,newne)+1 );
oe1check = cei-1; if oe1check==0, oe1check = ne; end
oe1check = m.secondlayer.cells(ci).edges( oe1check );
oe2check = newcei-1; if oe2check==0, oe2check = newne; end
oe2check = m.secondlayer.cells(newci).edges( oe2check );

if (ne3 ~= ne3check1) || (ne3 ~= ne3check2) ...
        || (ne1 ~= ne1check) ...
        || (ne2 ~= ne2check) ...
        || (oe1 ~= oe1check) ...
        || (oe2 ~= oe2check)
    fprintf( 1, '%s: oops\n', mfilename() );
end

% User code.

% Update the descendency tables.
numindexes = length(m.userdata.cellindex);
parentid = m.userdata.cellid(ci);
newindexes = [numindexes+1,numindexes+2];
m.userdata.cellid([ci,newci]) = newindexes;
m.userdata.cellindex([parentid,newindexes]) = [0,ci,newci];
m.userdatastatic.cellparent(newindexes) = [parentid,parentid];
m.userdatastatic.celldaughters(parentid,:) = newindexes;
m.userdatastatic.celltimes(parentid,2) = m.globalDynamicProps.currenttime;
m.userdatastatic.celltimes(newindexes,1) = m.globalDynamicProps.currenttime + m.globalProps.timestep;
m.userdatastatic.celltimes(newindexes,2) = Inf;

if m.userdatastatic.collectsplitdata
    % Determine if we are within the recording interval.
    [startInterval,endInterval] = findInterval( m );
    within = ~isempty(startInterval) || ~isempty(endInterval);
    if within
        % Record the indexes of the cells that split.
        m.userdata.splits([end+1, end+2]) = [ci; newci];
    end
end

% RK: set new area threshold perturbation for the daughter cells.
c_areavar_i = getCellFactorLevels( m, 'c_areavar' );
m.secondlayer.cellvalues([ci,newci],c_areavar_i) = ...
    randn( 2,1 ) * m.userdata.areavar;
end

function m = udstatictest_userpostiterateproc( m )
% This code runs in the post-iterate callback so that we can update
% cell areas after the cell splitting of the current iteration has
% happened, to keep it up to date.  Note that the time has been incremented
% by this point.

    c_neighbours_i = getCellFactorLevels( m, 'c_neighbours' );
    [~,c_neighbours] = numCellNeighbours( m, true );
    if ~isempty(c_neighbours)
        m.secondlayer.cellvalues(:,c_neighbours_i) = c_neighbours(:);
    end
    if m.userdatastatic.collectsplitdata
        [startInterval,endInterval,atstart,atend] = findInterval( m );
        if atstart
            foo = zeros( length(m.userdata.cellid), 1 );
            foo(m.userdata.cellid) = m.secondlayer.cellarea;
            m.userdatastatic.cellareastartbyid{startInterval} = foo;
            m.userdatastatic.cellareastart{startInterval} = m.secondlayer.cellarea;
            m.userdatastatic.cellidstart{startInterval} = m.userdata.cellid;
        end
        if atend
            numIDs = max(m.userdata.cellid);
            cellareaByID = zeros( numIDs, 1 );
            cellareaByID(m.userdata.cellid) = m.secondlayer.cellarea;
            m.userdatastatic.cellareaend{endInterval} = m.secondlayer.cellarea;
            m.userdatastatic.cellidend{endInterval} = m.userdata.cellid;

            cloneareas = cellareaByID;
            numcellidsatstart = length( m.userdatastatic.cellareastartbyid{endInterval} );
            %         wasAddedTo = false( length(cloneareas), 1 );
            for i=length(cloneareas):-1:(numcellidsatstart+1)
                pi = m.userdatastatic.cellparent(i);
                cloneareas(pi) = cloneareas(pi) + cloneareas(i);
                cloneareas(i) = 0;
                %             wasAddedTo(pi) = true;
            end
            cloneareas( (numcellidsatstart+1):end ) = [];
            %         wasAddedTo( (numcellidsatstart+1):end ) = [];
            m.userdatastatic.cloneareas{endInterval} = cloneareas( m.userdatastatic.cellidstart{endInterval} );

            % The ratio m.userdatastatic.cloneareas{n} ./ m.userdatastatic.cellareastart{n}
            % is the amount of relative growth for each cell extant at the start of interval n,
            % indexed by the cell indexes at that time.  All daughter cells at
            % the end of the interval are counted.
            xxxx = 1;
        end
        if atend
            % Copy the list of splits, removing duplicates (which will happen if a
            % cell splits more than once).
            m.userdatastatic.allsplits{endInterval} = unique( m.userdata.splits );
            saveStaticPart( m );
        end
    end
end


function m = udstatictest_userpreplotproc( m, theaxes )
if Steps(m)==0
    % Stuff isn't set up yet.
    return;
end

v_kx_i = getMgenLevels( m, 'V_KX' );
v_ky_i = getMgenLevels( m, 'V_KY' );
v_kz_i = getMgenLevels( m, 'V_KZ' );
m.morphogens(:,v_kx_i) = perFEtoperVertex( m, 0.5*(m.outputs.actualstrain.A(:,1) + m.outputs.actualstrain.B(:,1)) );
m.morphogens(:,v_ky_i) = perFEtoperVertex( m, 0.5*(m.outputs.actualstrain.A(:,2) + m.outputs.actualstrain.B(:,2)) );
m.morphogens(:,v_kz_i) = perFEtoperVertex( m, 0.5*(m.outputs.actualstrain.A(:,3) + m.outputs.actualstrain.B(:,3)) );

[startInterval,~,atstart,~] = findInterval( m );
if atstart
    starttime = m.userdatastatic.firsttime(startInterval);
    endtime = m.userdatastatic.lasttime(startInterval);
    finaltime = m.userdatastatic.lasttime(end);
    existsnow = (m.userdatastatic.celltimes(:,1) <= starttime) & (m.userdatastatic.celltimes(:,2) >= starttime);
    splitthisintervalmap = existsnow & (m.userdatastatic.celltimes(:,2) < endtime);
    splitthisorlaterintervalmap = existsnow & (m.userdatastatic.celltimes(:,2) < finaltime);
    splitthisinterval = m.userdata.cellindex( splitthisintervalmap(1:length(m.userdata.cellindex)) );
    splitthisorlaterinterval = m.userdata.cellindex( splitthisorlaterintervalmap(1:length(m.userdata.cellindex)) );
    
    c_splitrecord3 = zeros( size(m.secondlayer.cellvalues,1), 1 );
    c_splitrecord3(splitthisorlaterinterval) = 2;
    c_splitrecord3(splitthisinterval) = 1;
    c_splitrecord3_i = getCellFactorLevels( m, 'c_splitrecord3' );
    m.secondlayer.cellvalues(:,c_splitrecord3_i) = c_splitrecord3(:);
end


% Find the relative increase in area during this step.
m.userdata.currentTissueRelGrowth = sum(m.cellareas)/m.userdata.currentTissueArea - 1;
% This procedure does not change the static data, but previously called
% ones may have done, and this procedure is the last user code called
% before the simulation step ends, so this ensures that the save happens.
saveStaticPart( m );
end

function createLogFile( logFileHeaders,  filename )
%write exel file with project data to project directory
fid = fopen(['Extracted Growth', filesep, filename], 'wt');
csvFun = @(str)sprintf('%s,',str);
xchar = cellfun(csvFun, logFileHeaders, 'UniformOutput', false);
xchar = strcat(xchar{:});
xchar = strcat(xchar(1:end-1),'\n');
fprintf(fid,xchar);
fclose(fid);
end

function [startInterval,endInterval,atstart,atend] = findInterval( m )
atafterfirst = meshAtOrAfterTime( m, m.userdatastatic.firsttime );
atbeforelast = meshAtOrBeforeTime( m, m.userdatastatic.lasttime );
whichInterval = find( atafterfirst & atbeforelast );
switch length(whichInterval)
    case 0
        startInterval = [];
        endInterval = [];
    case 1
        startInterval = whichInterval;
        endInterval = whichInterval;
    otherwise
        endInterval = whichInterval(1);
        startInterval = whichInterval(2);
end
atstart = meshAtTime( m, m.userdatastatic.firsttime(startInterval) );
atend = meshAtTime( m, m.userdatastatic.lasttime(endInterval) );
end

function plotArabidopsisRun( m )
    projectdir = getModelDir( m );
%     modelname = m.userdata.ranges.modelname.range{m.userdata.ranges.modelname.index}
    modelname = getModelOption( m, 'modelname' );

    origdatafile = fullfile( projectdir, 'Extracted Growth', 'EPIMODEL_LeafData EPIMODEL ORIG.csv' );
    if ~exist( origdatafile, 'file' )
        return;
    end
    origdata = csvread( origdatafile, 1, 0 );
    origrealtime = origdata(:,1);
    origlogwidth = log( origdata(:,3) );

    origspchdatafile = fullfile( projectdir, 'Extracted Growth', 'SPCH_ORIG.csv' );
    if ~exist( origspchdatafile, 'file' )
        return;
    end
    origspchdata = csvread( origspchdatafile, 1, 0 );
    origspchrealtime = origspchdata(:,1);
    origspchlogwidth = origspchdata(:,2);

    origspchdatafile2 = fullfile( projectdir, 'Extracted Growth', 'spch_on_plates.csv' );
    if ~exist( origspchdatafile2, 'file' )
        return;
    end
    origspchdata2 = csvread( origspchdatafile2, 2, 0 );
    origspchrealtime2 = origspchdata2(:,1);
    origspchlogwidth2 = origspchdata2(:,2);

    origfulldatafile = fullfile( projectdir, 'Extracted Growth', 'origfulldata.csv' );
    if ~exist( origfulldatafile, 'file' )
        return;
    end
    fulldata = csvread( origfulldatafile, 2, 0 );
    fullrealtime = fulldata(:,1);
    lnlogisticwidth = fulldata(:,8);

    currentdatafile = fullfile( projectdir, 'Extracted Growth', [modelname,'_LeafData.csv'] );
    if ~exist( currentdatafile, 'file' )
        return;
    end
    try
        data = csvread( currentdatafile, 1, 0 );
    catch e
        % Because Matlab can't cope with there being no data in the file.
        return;
    end
    realtime = data(:,1);
    logwidth = log( data(:,3) );

    figure(1);
    plot( realtime, logwidth, 'ob', 'MarkerSize', 8 );
    hold on
    plot( origrealtime, origlogwidth, '.r' );
    plot( origspchrealtime, origspchlogwidth, '^g', 'MarkerFaceColor', 'g', 'MarkerSize', 12 );
    plot( origspchrealtime2, origspchlogwidth2, 'xk', 'MarkerFaceColor', 'k', 'MarkerSize', 12 );
    plot( fullrealtime, lnlogisticwidth, '.-k', 'MarkerFaceColor', 'k', 'MarkerSize', 12 );
    hold off
    axis([0 500 -3 3]);
end

function isat = atTimeX( currenttime, timepoint, timestep )
% Determine whether currenttime is close to timepoint, in such a way that
% if currentTime is successively increased by timestep, beginning when it
% is less than timepoint and ending at some time greater than timepoint,
% then this function will return true exactly once.

    timediff = currenttime - timepoint;
    isat = (timediff >= -timestep/2) && (timediff < timestep/2);
    if isat
        fprintf( 1, 'At time point %g, phys time %g, step %g.\n', timepoint, currenttime, timestep );
    end
end
