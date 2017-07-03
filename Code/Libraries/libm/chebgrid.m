% Inputs:
%   polord:  d-vector of polynomial degrees (3 for cubic, etc.)
%   lowb:    d-vector of lower bounds;
%   uppb:    d-vector of upper bounds;
% Outputs:
%   grid:  Cartesian grid of Chebyshev zeros, scaled to lowb,uppb
%   if nargout==2: same grid, but unscaled
%
function [grid,grid01] = chebgrid(polord,lowb,uppb)
  d = length(polord);

  nodes = cell(d,1);
  nodes_sc = nodes;

  for(i=1:d)
    nodesi = chebzero(polord(i)+1);
    gridi  = scalup(nodesi,lowb(i),uppb(i));
    nodes{i} = nodesi;
    nodes_sc{i} = gridi;
  end;

  if(nargout==2)
    grid01 = gridmake(nodes{:});
  end;
  grid = gridmake(nodes_sc{:});
