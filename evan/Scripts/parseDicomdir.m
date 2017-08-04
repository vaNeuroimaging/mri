function patient = parseDicomdir(fname)
%PARSEDICOMDIR    Parse a DICOMDIR file into a tree structure.
% 
% patient = parseDicomdir(fname): generate a parse tree of a DICOMDIR file.
% A DICOMDIR file describes a collection of (usually cryptically-named)
% DICOM image files and contains rich metadata for these images.
% 
% Arguments:
% - fname (string): the filename for a DICOMDIR file.
% 
% Returns:
% - patient (cell array of structs): the root of a metadata tree with the
%   following levels of nodes: patient, study, series, image.  Each node
%   is a struct, with an "info" field containing all the metadata as stored
%   in the DICOMDIR file and another field, named for the next deeper level
%   of nodes, that is a cell array of node structs.
% 
% Examples:
%   After calling the function,
%   - metadata for the i-th image in the j-th series in the k-th study
%     for the l-th patient described by this DICOMDIR is
%       patient{l}.study{k}.series{j}.image{i}.info;
%   - metadata for the k-th study for the l-th patient is
%       patient{l}.study{k}.info;
%   - the birth date for the l-th patient is
%       patient{l}.info.PatientBirthDate;
%   - the date of the k-th study for the l-th patient is
%       patient{l}.study{k}.info.StudyDate;
%   - the description of the j-th series in the k-th study for the l-th
%     patient is
%       patient{l}.study{k}.series{j}.info.SeriesDescription;
%   - the relative path from the DICOMDIR file to the actual image file for
%     the i-th image in the j-th series in the k-th study for the l-th
%     patient is
%       patient{l}.study{k}.series{j}.image{i}.info.ReferencedFileID;
%     Note that this relative path may be formatted using a path
%     separator different from the one used on your system.
% 
% Prerequisites:
%   This function relies on the dicominfo() function from the Matlab image
% processing toolbox.
% 

% Written 2010-05-25 by Jadrian Miles
	
	try
		dcmhdr = dicominfo(fname);
	catch me
		if strcmp(me.identifier, 'Images:dicominfo:notDICOM')
			error(me.identifier, '"%s" is not a DICOMDIR file.', fname);
		end
		error(me.identifier, me.message);
	end
	
	hdrfields = fieldnames(dcmhdr.DirectoryRecordSequence);
	nFields = length(hdrfields);
	nPatients = 0;
	nStudies = 0;
	nSeries = 0;
	nImages = 0;
	patient = cell(0);
	for i = 1:nFields
		record = dcmhdr.DirectoryRecordSequence.(hdrfields{i});
		switch lower(record.DirectoryRecordType)
			case 'patient'
				nPatients = nPatients + 1;
				nStudies = 0;
				patient{nPatients}.info = record;
				patient{nPatients}.study = cell(0);
			case 'study'
				nStudies = nStudies + 1;
				nSeries = 0;
				patient{nPatients}.study{nStudies}.info = record;
				patient{nPatients}.study{nStudies}.series = cell(0);
			case 'series'
				nSeries = nSeries + 1;
				nImages = 0;
				patient{nPatients}.study{nStudies}.series{nSeries}.info = record;
				patient{nPatients}.study{nStudies}.series{nSeries}.image = cell(0);
			case 'image'
				nImages = nImages + 1;
				patient{nPatients}.study{nStudies}.series{nSeries}.image{nImages}.info = record;
			otherwise
				warning('parseDicomdir:unknownDirtype', ...
					'Unknown directory record type "%s" at index %d', ...
					record.DirectoryRecordType, i);
		end
	end
end
