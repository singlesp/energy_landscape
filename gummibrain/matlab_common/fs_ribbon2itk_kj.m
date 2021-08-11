function fs_ribbon2itk(subjID, outfile, fillWithCSF, alignTo)
% Read in a freesurfer ribbon.mgz file and convert to a nifti class file
% with our itkGray conventional labels.
%
% fs_ribbon2itk(subjID, [outfile], [fillWithCSF], [alignTo])
%
% Inputs:
%   subjID:      name of directory in freesurfer subejct directory (string)
%   outfile:     name of output nifti segmentation file (including path)
%                   [default = 't1_class.nii.gz']
%   fillWithCSF: if true, convert unlabeled voxels to CSF (0s => 1s)
%                   [default = false]        
%   alignTo:     optional nifti file which dertermines bounding box and
%                 alignment of output segmnentation
%
% The freesurfer automatic segmentation process produces many output files.
% From our preliminary experience, the gray-white segmentation is quite
% accurate, and creates white matter segmentations with minimal topological
% errors (handles, cavities). 
%
% The file called ribbon.mgz is similar to the class files we create when
% we segment in itkGray. It is aligned to the freesurfer file 't1.mgz', and
% it contains integer values: 
%
%    2: left white matter
%    3: left gray matter
%   41: right white matter
%   42: right gray matter
%    0: unlabeled
%
% If we want to use this freesurfer segmentation, then we need to (a)
% convert the mgz file to nifti, and (b) change the label values to 
%    3: left white matter
%    5: left gray matter
%    4: right white matter
%    6: right gray matter 
%
%    1: csf (if input argument fillWithCSF = true), or  
%    0: unlabeled
%
%
% For description of freesurfer automatic segmentation, see
%   http://surfer.nmr.mgh.harvard.edu/fswiki/ReconAllDevTable
% And for the ribbon file specifcally see:
%   http://surfer.nmr.mgh.harvard.edu/fswiki/cortribbon
%
% Example 1:
%   fs_ribbon2itk('bert')
%
% Example 2:
%   subjID      = 'andreas';
%   outfile     = 't1_freesurfer_class.nii.gz'
%   fillWithCSF = true; 
%   alignTo     = '/biac2/wandell2/data/anatomy/rauschecker/Anatomy081031/t1.nii.gz';
%   fs_ribbon2itk(subjID, outfile, fillWithCSF, alignTo);

%% Check Inputs
if ~exist('subjID', 'var')
    warning('Subject ID is required input'); %#ok<WNTAG>
    help fs_ribbon2itk;
    return
end

if notDefined('fillWithCSF'), fillWithCSF = false; end

%% Find paths
% subdir   = getenv('SUBJECTS_DIR');
% 
% if isempty(subdir), 
%     fshome = getenv('FREESURFER_HOME');
%     subdir = fullfile(fshome, 'subjects'); 
%     
% end
% 
%subdir   = element(unixvm('echo $SUBJECTS_DIR'),1);
subdir = freesurferdir;

if isempty(subdir), 
    fsdir   = element(unixvm('echo $FREESURFER_HOME'),1);
    subdir = fullfile(fshome, 'subjects'); 
    
end



ribbon = fullfile(subdir, subjID, 'mri', 'ribbon.mgz');

if ~exist(ribbon, 'file'),
    [fname pth] = uigetfile({'ribbon*.mgz', 'Ribbon files'; '*.mgz', '.mgz files'}, 'Cannot locate ribbon file. Please find it yourself.', pwd);
    ribbon = fullfile(pth, fname);
end    

if ~exist(ribbon, 'file'), error('Cannot locate ribbon.mgz file'); end

if notDefined('outfile'), 
    pth     = fileparts(ribbon);
    outfile = fullfile(pth, 't1_class.nii.gz'); 
end

%% Convert MGZ to NIFTI

% 
% if exist('alignTo', 'var'),
%     str = sprintf('!mri_convert  --out_orientation RAS --reslice_like %s -rt nearest %s  %s', alignTo, ribbon, outfile);
% else
%     str = sprintf('!mri_convert  --out_orientation RAS -rt nearest %s  %s', ribbon, outfile);
% end
% eval(str)
%     

if exist('alignTo', 'var'),
    str = sprintf('mri_convert  --out_orientation LIA --reslice_like %s -rt nearest %s  %s', alignTo, ribbon, outfile);
else
    str = sprintf('mri_convert  --out_orientation LIA -rt nearest %s  %s', ribbon, outfile);
end
unixvm(str)

%% Convert freesurfer label values to itkGray label values
% We want to convert
%   Left white:   2 => 3
%   Left gray:    3 => 5
%   Right white: 41 => 4
%   Right gray:  42 => 6
%   unlabeled:    0 => 0 (if fillWithCSF == 0) or 1 (if fillWithCSF == 1)          

% read in the nifti
ni = readFileNifti(outfile);

% check that we have the expected values in the ribbon file
vals = sort(unique(ni.data(:)));
if ~isequal(vals, [0 2 3 41 42]')
    warning('The values in the ribbon file - %s - do no match the expected values [2 3 41 42]. Proceeding anyway...') %#ok<WNTAG>
end

% map the replacement values
invals  = [3 2 41 42];
outvals = [5 3 4  6];
labels  = {'L Gray', 'L White', 'R White', 'R Gray'};

fprintf('\n\n****************\nConverting voxels....\n\n');
for ii = 1:4;
    inds = ni.data == invals(ii);
    ni.data(inds) = outvals(ii);
    fprintf('Number of %s voxels \t= %d\n', labels{ii}, sum(inds(:)));
end

if fillWithCSF, 
    ni.data(ni.data == 0) = 1;
end

% write out the nifti
writeFileNifti(ni)

% done.
return
