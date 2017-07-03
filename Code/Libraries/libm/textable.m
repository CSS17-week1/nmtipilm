% fmt can be cell array for each column, or just one format string
% '&' is inserted automatically
function tex_table(data,fmt,filehandle)
  [nr,nc] = size(data);
  for(i=1:nr)
    for(j=1:nc)
      if(iscell(fmt))
	fj = fmt{j};
      else
	fj = fmt;
      end
      fprintf(filehandle,fj,data(i,j));
      if(j<nc)
	fprintf(filehandle,' & ');
      else
	fprintf(filehandle,'\\\\\n');
      end
    end
  end
