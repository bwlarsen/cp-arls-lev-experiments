function write_vector(fid,name,v,frmt,maxlen)
%WRITE_VECTOR Writes vector to specified output, padding with NaN's
%
%   WRITE_VECTOR(FID,NAME,X) writes the vector X to the output specified by
%   FID. The NAME is printed first, then a comma, then the entries of X
%   separated by commas. The last entry has no comma and is followed by a
%   new line.
%
%   WRITE_VECTOR(FID,NAME,X,FMT) specifies the format of the entries. The
%   default is FMT='%.6f'. Do not add commas or spaces.
%
%   WRITE_VECTOR(FID,NAME,X,FMT,LEN) specifies to pad the vector up to
%   length LEN with NaN's. This is useful when all vectors need to be the
%   same length.


if ~exist('maxlen','var')
    maxlen = 0;
end
    
if ~exist('frmt','var')
    frmt = '%.6f';
end


fprintf(fid,'%s, ', name);

len = length(v);

% Print vector
for i = 1:len
    fprintf(fid, frmt, v(i));
    if i < len
        fprintf(fid, ', ');
    end
end

% Pad, if needed
for i = len+1:maxlen
    fprintf(fid,', NaN');
end

fprintf(fid, '\n');


