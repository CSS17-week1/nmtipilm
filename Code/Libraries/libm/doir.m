function IR = doir(G1,impact,VARNAMES,imp_sel,HORIZON,PERIOD,shocknames,scale_shocks,scale_y)
  ny = size(G1,1);
  nz = size(impact,2);

  if(nargin<9)
    scale_y = [];
  end;
  if(nargin<8)
    scale_shocks = [];
  end;
  if(nargin<7)
    shocknames = [];
    for(i=1:nz)
      ish = sprintf('shock%2d',i);
      shocknames = [shocknames;ish];
    end
  end;

  if(nargin<6)
    PERIOD = 4;
  end;

  if(nargin<5)
    HORIZON = 100;
  end;

  if(nargin<4)
    imp_sel = 1:ny;
  end;

  if(nargin<3)
    VARNAMES = [];
    for(i=1:ny)
      ivar = sprintf('var%2d',i);
      VARNAMES = [VARNAMES;ivar];
    end
  end;

  IR = impresp_sims(G1,impact,imp_sel,HORIZON,PERIOD,VARNAMES,shocknames,scale_shocks,scale_y);
