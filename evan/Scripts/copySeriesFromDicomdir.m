function copySeriesFromDicomdir(src, patt, dst)
%COPYSERIESFROMDICOMDIR    Copy images for series that match a pattern.
% 
% copySeriesFromDicomdir(src, patt, dst): read the DICOMDIR file at 'src'
% and copy all images from series described in this file whose descriptions
% match the specified regex pattern 'patt' into subdirectories of the
% destination directory 'dst'.
% 
% Arguments:
% - src (string): the root path of a DICOM media directory, containing a
%   DICOMDIR file.
% - patt (string): a regex pattern to compare series descriptions against.
%   If any part of a series' description matches 'patt' (case-insensitive,
%   using regexpi()), the images for that series will be copied.
% - dst (string): the path of an existing directory, in which a new
%   directory will be created for each matching series, and the images for
%   that series will be copied into that directory.
% 
% Returns: nothing.
% 
% See also: parseDicomdir

% Written 2010-05-25 by Jadrian Miles
	
	ffmt = {'DICOMDIR', 'dicomdir'};
	i = 0;
	while i < length(ffmt) && isempty(dir(fullfile(src, ffmt{i+1})))
		i = i+1;
	end
	if i == length(ffmt)
		error('copySeriesFromDicomdir:noDicomdirFound', ...
			'No DICOMDIR file found in source path "%s".', src);
	end
	
	fprintf('Reading %s ...\n', fullfile(src, ffmt{i+1}));
	patient = parseDicomdir(fullfile(src, ffmt{i+1}));
	
	fprintf('Searching for series that match "%s" ...\n', patt);
	for i=1:length(patient)
	for j=1:length(patient{i}.study)
	for k=1:length(patient{i}.study{j}.series)
		s = patient{i}.study{j}.series{k};
		if regexpi(s.info.SeriesDescription, patt)
			pname = strclean(patient{i}.info.PatientName.FamilyName);
			stname = strclean(patient{i}.study{j}.info.StudyDescription);
			srname = strclean(s.info.SeriesDescription);
			dstdir = fullfile(dst, strjoin({pname, stname, srname}, '.'));
			fprintf('Found %d files; copying into %s ...\n', ...
				length(s.image), dstdir);
			mkdir(dstdir);
			for l=1:length(s.image)
				relpath = strrep(s.image{l}.info.ReferencedFileID, '\', filesep);
				fsrc = fullfile(src, relpath);
				fdst = fullfile(dstdir, sprintf('img%03d.dcm', l));
				copyfile(fsrc, fdst);
			end
			fprintf('Done.\n');
		end
	end
	end
	end
end

function out = strclean(in)
	badchars = regexpi(in, '[^a-zA-Z0-9_]');
	out = in;
	out(badchars) = repmat('_', 1, length(badchars));
end

function out = strjoin(strs, sep)
	out = strs{1};
	for i=2:length(strs)
		out = [out, sep, strs{i}];
	end
end
