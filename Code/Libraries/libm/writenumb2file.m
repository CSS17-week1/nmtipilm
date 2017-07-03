function writenumb2file(filename,x)
  FF = fopen(filename,'w');
  fprintf(FF,'%g  ',x);
  fprintf(FF,'\n');
  fclose(FF);