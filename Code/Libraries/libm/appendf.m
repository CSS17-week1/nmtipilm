function appendf(filename,fmt,varargin)
  fid = fopen(filename,'a');
  fprintf(fid,fmt,varargin{:});
  fclose(fid);
