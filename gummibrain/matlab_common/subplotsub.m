function h = subplotsub(numrows,numcols,whichrow,whichcol,varargin)

n = (whichrow-1)*numcols + whichcol;
h=subplot(numrows,numcols,n,varargin{:});
