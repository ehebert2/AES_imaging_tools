function skel2_interp = GU_extractVesselness_Skel(fnrt, varargin)
% extract Vessels + medial axis + distance transform
% writes files in the same location with suffixess
% also writes interpolated volumes with isotropic voxel sizes
% Skeletonization code based off of:
% Kollmannsberger, Kerschnitzki et al., "The small world of osteocytes: connectomics of the lacuno-canalicular network in bone." New Journal of Physics 19:073019, 2017.
% T. Jerman, F. Pernus, B. Likar, Z. Spiclin, "Enhancement of Vascular Structures in 3D and 2D Angiographic Images", IEEE Transactions on Medical Imaging, 35(9), p. 2107-2118 (2016)
% T. Jerman, F. Pernus, B. Likar, Z. Spiclin, "Blob Enhancement and Visualization for Improved Intracranial Aneurysm Detection", IEEE Transactions on Visualization and Computer Graphics, 22(6), p. 1705-1717 (2016)
% with read/write libraries from:
% F. Aguet, S. Upadhyayula, R. Gaudin, et al., "Membrane dynamics of dividing cells imaged by lattice light-sheet microscopy", Molecular Biology of the Cell 2016 Nov 7;27(22):3418-3435.

% Gokul Upadhyayula, Oct 2019

%%%% input: fnrt = 'root+filename'; eg:'D:\Na\similar_gauss_region.tif'
% GU_extractVesselness_Skel('D:\Na\similar_gauss_region.tif');

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('fnrt', @isstr);
ip.addParameter('kernelSize', 3 , @isnumeric); % kernelSize for erosion and dilation
ip.addParameter('VoxelSizes', [0.8702; 0.8702; 1]); %0.8702 %1.3238
ip.addParameter('minBranchLength', 30 , @isnumeric); % for the skeletonization
ip.addParameter('tau', 0.95 , @isnumeric); % controls response uniformity 0.5-1. Smaller tau extracts more segments. 
ip.addParameter('bridgeSize', 1 , @isnumeric); % kernel size for connecting broken vessels 
ip.parse(fnrt, varargin{:});
p = ip.Results;

%%
im = readtiff(fnrt); % read image
se = strel('sphere',p.kernelSize); % dilation/erosion kernel
spacing = p.VoxelSizes; %voxel dimensions in microns  [y ,x ,z]
minBranchLength = p.minBranchLength; % prevent excessive branching on large vessels

%% segment vessels and clean
V = vesselness3D(im, 1:4, spacing, p.tau, true); % segment using Hessian + Eigen parameters

% Lower vesselness threshold to retain more weak signals
V_thresh = V;
V_thresh(V_thresh < 5e-3) = 0;  % Lower threshold

cm = imdilate(logical(V_thresh),se); %link connected components
cm = imerode(cm, se);

% Additional connection processing - use larger kernel to connect broken vessels
if p.bridgeSize > p.kernelSize
    bridge_se = strel('sphere', p.bridgeSize);
    cm_bridge = imdilate(cm, bridge_se);
    cm_bridge = imerode(cm_bridge, bridge_se);
    % Preserve original details but connect disconnected parts
    cm = cm | cm_bridge;
end

% discard small components
CC = bwconncomp(cm, 26);
csize = cellfun(@numel, CC.PixelIdxList);
minVox = prctile(csize, 90); % Changed to 90%, more lenient filtering
idx = csize>=minVox;
CC.NumObjects = sum(idx);
CC.PixelIdxList = CC.PixelIdxList(idx);
cm = labelmatrix(CC)~=0;

% cm = imclose(logical(V), se);
dtcm = bwdistsc(~cm,spacing);

%% skeletonization
pcm = padarray(cm,[3 3 3],0,'both');
tic
skel = skeleton3D(pcm);
toc
skel = skel(4:end-3,4:end-3,4:end-3);

% Add connection step after skeletonization
fprintf('Original skeleton generation completed, starting to connect disconnected parts...\n');
%skel_connected = connect_broken_vessels_improved(skel, spacing, 15, 2); % 15 micron intelligent connection, using 5x5x5 neighborhood search
skel_connected = skel; % Temporarily disable connection function to avoid false connections

tic
[A,node,link] = Skel2Graph3D(skel_connected,minBranchLength);
% generate a cleaned up graph
w = size(skel_connected,1);
l = size(skel_connected,2);
h = size(skel_connected,3);
skel2 = Graph2Skel3D(node,link,w,l,h);
toc
skeldtm = single(dtcm);
skeldtm(~skel2) = 0;

%% interpolate to make isotropic voxels --> interpolates to make all voxels 0.8702 um
[ny,nx,nz] = size(im);
[y,x,z] = ndgrid(1:ny,1:nx,1:nz);
%[Y,X,Z] = ndgrid(1:1/spacing(1):ny,1:1/spacing(2):nx,1:1/spacing(3):nz); % adjust this line if you need to adjust anisotropy
[Y,X,Z] = ndgrid(1:0.8702/spacing(1):ny,1:0.8702/spacing(2):nx,1:0.8702/spacing(3):nz); % adjust this line if you need to adjust anisotropy
skel2_interp = single(interp3(x,y,z,double(skel2),X,Y,Z,'nearest'));
clear x y z X Y Z;
%% Connectivity Analysis
fprintf('\n=== Connectivity Analysis Report ===\n');

% Analyze connectivity of original vessel segmentation
CC_vessels = bwconncomp(cm, 26);
fprintf('Number of connected components after vessel segmentation: %d\n', CC_vessels.NumObjects);

% Analyze connectivity of original skeleton
CC_skel_orig = bwconncomp(skel, 26);
fprintf('Number of connected components in original skeleton: %d\n', CC_skel_orig.NumObjects);

% Analyze connectivity of connected skeleton
CC_skel_connected = bwconncomp(skel_connected, 26);
fprintf('Number of connected components in connected skeleton: %d\n', CC_skel_connected.NumObjects);

% Analyze connectivity of final cleaned skeleton
CC_skel_final = bwconncomp(skel2, 26);
fprintf('Number of connected components in final skeleton: %d\n', CC_skel_final.NumObjects);

% Analyze sizes of each connected component
[component_stats, largest_component] = analyze_component_sizes(CC_skel_final);
fprintf('\nConnected Component Size Analysis:\n');
fprintf('- Largest component contains %d voxels (%.1f%%)\n', ...
    component_stats.max_size, component_stats.max_percentage);
fprintf('- Average component size: %.1f voxels\n', component_stats.mean_size);
fprintf('- Number of small components (<%d voxels): %d\n', component_stats.small_threshold, component_stats.num_small);

% Save labeled component image
labeled_components = labelmatrix(CC_skel_final);
writetiff(uint16(labeled_components), [fnrt(1:end-4) '_components.tif']);

% Create connectivity summary file
create_connectivity_report(fnrt, CC_vessels, CC_skel_orig, CC_skel_connected, CC_skel_final, component_stats);

%% Node Degree Analysis
fprintf('\n=== Starting Node Degree Analysis ===\n');
[degree_stats, node_degrees] = analyze_node_degrees(A, node, link);

% Save degree analysis results to file
degree_report_file = [fnrt(1:end-4) '_degree_analysis.txt'];
fid = fopen(degree_report_file, 'w');
if fid ~= -1
    try
        fprintf(fid, 'Vessel Network Node Degree Analysis Report\n');
        fprintf(fid, 'Generated at: %s\n', datestr(now));
        fprintf(fid, 'Data file: %s\n\n', fnrt);
        
        fprintf(fid, 'Node Degree Statistics:\n');
        fprintf(fid, '- Total number of nodes: %d\n', degree_stats.total_nodes);
        fprintf(fid, '- Minimum degree: %d\n', degree_stats.min_degree);
        fprintf(fid, '- Maximum degree: %d\n', degree_stats.max_degree);
        fprintf(fid, '- Average degree: %.2f\n', degree_stats.mean_degree);
        fprintf(fid, '- Median degree: %.1f\n', degree_stats.median_degree);
        
        fprintf(fid, '\nNode Type Distribution:\n');
        fprintf(fid, '- Isolated nodes (degree=0): %d (%.1f%%)\n', ...
            degree_stats.isolated_nodes, degree_stats.isolated_nodes/degree_stats.total_nodes*100);
        fprintf(fid, '- End nodes (degree=1): %d (%.1f%%)\n', ...
            degree_stats.end_nodes, degree_stats.end_nodes/degree_stats.total_nodes*100);
        fprintf(fid, '- Regular nodes (degree=2): %d (%.1f%%)\n', ...
            sum(node_degrees == 2), sum(node_degrees == 2)/degree_stats.total_nodes*100);
        fprintf(fid, '- Junction nodes (degree≥3): %d (%.1f%%)\n', ...
            degree_stats.junction_nodes, degree_stats.junction_nodes/degree_stats.total_nodes*100);
        fprintf(fid, '- Hub nodes (degree≥5): %d (%.1f%%)\n', ...
            degree_stats.hub_nodes, degree_stats.hub_nodes/degree_stats.total_nodes*100);
        
        fprintf(fid, '\nDegree Distribution Details:\n');
        for i = 1:length(degree_stats.degree_values)
            if degree_stats.degree_counts(i) > 0
                fprintf(fid, '- Degree %d: %d nodes (%.1f%%)\n', ...
                    degree_stats.degree_values(i), degree_stats.degree_counts(i), ...
                    degree_stats.degree_distribution(i));
            end
        end
        
        fprintf(fid, '\nNetwork Topology Assessment:\n');
        if degree_stats.junction_nodes > degree_stats.total_nodes * 0.1
            fprintf(fid, '- Network complexity: High (junction node ratio %.1f%%)\n', ...
                degree_stats.junction_nodes/degree_stats.total_nodes*100);
        elseif degree_stats.junction_nodes > degree_stats.total_nodes * 0.05
            fprintf(fid, '- Network complexity: Medium (junction node ratio %.1f%%)\n', ...
                degree_stats.junction_nodes/degree_stats.total_nodes*100);
        else
            fprintf(fid, '- Network complexity: Low (junction node ratio %.1f%%)\n', ...
                degree_stats.junction_nodes/degree_stats.total_nodes*100);
        end
        
        if degree_stats.end_nodes > degree_stats.total_nodes * 0.3
            fprintf(fid, '- Network integrity: Poor (end node ratio too high %.1f%%)\n', ...
                degree_stats.end_nodes/degree_stats.total_nodes*100);
        else
            fprintf(fid, '- Network integrity: Good (end node ratio %.1f%%)\n', ...
                degree_stats.end_nodes/degree_stats.total_nodes*100);
        end
        
        fclose(fid);
        fprintf('Degree analysis report saved to: %s\n', degree_report_file);
    catch ME
        fprintf('Failed to save degree analysis report: %s\n', ME.message);
        fclose(fid);
    end
end

%% write files
writetiff(uint8(cm), [fnrt(1:end-4) '_vessels.tif']);
writetiff(single(dtcm), [fnrt(1:end-4) '_dtm.tif']);
writetiff(uint8(skel), [fnrt(1:end-4) '_skel_original.tif']);  % Original skeleton
writetiff(uint8(skel_connected), [fnrt(1:end-4) '_skel_connected.tif']);  % Connected skeleton
writetiff(uint8(skel2), [fnrt(1:end-4) '_skel.tif']);
writetiff(uint8(skel2_interp), [fnrt(1:end-4) '_skelInterp.tif']);