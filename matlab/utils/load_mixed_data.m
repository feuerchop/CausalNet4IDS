function [data, sym_tbl, hdls] = load_mixed_data( filename, varargin )
%IMPORTMIXEDDATA Summary of this function goes here
   % load data from text-numeric mixed file
   % data format can be numeric or norminal
   fid = fopen(filename);
   delimiter = strcmpi('delimiter', varargin);
   delimiter_idx = find(delimiter);
   if delimiter_idx > 0
   delimiter = varargin{delimiter_idx + 1};
   else
   delimiter = ',';
   end
   % set the headline
   headline = strcmpi('HeaderLines', varargin);
   headline_idx = find(headline);
   if headline_idx > 0
   headline = varargin{headline_idx + 1};
   hdls = cell(headline, 1);
   % skip head lines
   for j = 1:headline
      firstln = fgetl(fid);
      hdls{j, 1} = firstln;
   end
   else
   if (nargout > 2)
      fprintf('OUTPUT HEADLINE NOT ASSIGNED!\n');
      return;
   end
   headline = 0;
   firstln = fgetl(fid);
   end
   % read the line format
   aInst = strsplit(delimiter, firstln);
   fseek(fid, 0, 'bof');

   % construct the format string
   format = [];
   for j = 1:length(aInst)
   attr = aInst{j};
   if  isnan(str2double(attr))
      % it is a string value
      format = [format, '%s'];
   else
      % it is a numeric value
      format = [format, '%f'];
   end
   end

   % TODO: make it more robust
   data2 = textscan(fid, format, 'delimiter', delimiter, 'HeaderLines', headline);
   [data, sym_tbl] = data_merge(data2);
   fclose(fid);
end

