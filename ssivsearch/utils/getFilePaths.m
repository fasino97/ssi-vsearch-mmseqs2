function [ vsearchPath, usearchPath ] = getFilePaths( )
%GETFILEPATHS Return the paths to usearch and vsearch
%   Return the correct file path to VSEARCH and USEARCH depending on
%   what type of system you are on.
%   If this function is throwing errors:
%       - Run getenv('OS') on MATLAB and add it to the case statement below
%         along with the paths to usearch and vsearch save on your machine
%       - If your getenv('OS') is already listed, consider changing the
%         directory of your file paths to the ones here.

switch getenv('OS')
    case '' %HPC
        usearchPath = '/coe_data/spprl/old-uclust/old-uclust';
        vsearchPath = '/coe_data/spprl/vsearch/src/vsearch/bin/vsearch';
    case 'Windows_NT' %Windows Machine in 308
        usearchPath = [pwd, '/../../usearch5.2.236_win32'];
        vsearchPath = [pwd, '/../../vsearch-2.4.3-win-x86_64/vsearch'];
    otherwise
        error('VSEARCH and USEARCH paths not found')
end



end

