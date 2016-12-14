function fpl_vtk_write_field_quadmesh (basename, mesh, nodedata, celldata, endfile)

  ## Check input
  if nargin!=5
    error ("fpl_vtk_write_field: wrong number of input parameters");
  endif

  if (! ischar (basename))
    error ("fpl_vtk_write_field: basename should be a string");
  elseif (! isstruct (mesh))
    error ("fpl_vtk_write_field: mesh should be a struct");
  elseif (! (iscell (nodedata) && iscell (celldata)))
    error ("fpl_vtk_write_field: nodedata and celldata should be cell arrays");
  endif

  filename = [basename ".vtu"];

  if (! exist (filename, "file"))
    fid = fopen (filename, "w");

    ## Header
    fprintf (fid, "<?xml version=""1.0""?>\n");
    fprintf (fid, "<VTKFile type=""UnstructuredGrid"" version=""0.1"" byte_order=""LittleEndian"">\n");
    fprintf (fid, "<UnstructuredGrid>\n");
  else

    ## FIXME: the following should be performed in a cleaner way! Does a
    ## backward fgetl function exist?

    ## If file exist, check if it was already closed
    fid = fopen (filename, "r");
    fseek (fid, -10, SEEK_END);
    tst = fgetl (fid);
    if (strcmp (tst, "</VTKFile>"))
      error ("fpl_vtk_write_field: file %s exist and was already closed", filename);
    endif
    fclose (fid);
    fid = fopen (filename, "a");
  endif    

  p   = mesh.p;
  dim = rows (p); # 2D or 3D

  if dim == 2
    t = mesh.t (1:4, ! any(mesh.children));
  elseif dim == 3
    t = mesh.t (1:8, ! any(mesh.children));
  else
    error ("fpl_vtk_write_field: neither 2D quadrilaterals nor 3D hexahedral mesh");    
  endif
  
  t -= 1;
  
  nnodes = columns (p);
  nelems = columns (t);

  ## Header for <Piece>
  fprintf (fid, "<Piece NumberOfPoints=""%d"" NumberOfCells=""%d"">\n", nnodes, nelems);

  ## Print grid
  print_grid (fid, dim, p, nnodes, t, nelems);

  ## Print Data
  print_data_points (fid, nodedata, nnodes)
  print_cell_data  (fid, celldata, nelems)

  ## Footer for <Piece>
  fprintf (fid, "</Piece>\n");

  if (endfile)
    ## Footer
    fprintf (fid, "</UnstructuredGrid>\n");
    fprintf (fid, "</VTKFile>");
  endif

  fclose (fid);
	   
endfunction

## Print Points and Cells Data
function print_grid (fid, dim, p, nnodes, t, nelems)
  
  if dim == 2
    p      = [p; zeros(1,nnodes)];
    eltype = 9;
  else
    eltype = 12;
  endif
  
  ## VTK-Points (mesh nodes)
  fprintf (fid, "<Points>\n");
  fprintf (fid, "<DataArray type=""Float64"" Name=""Array"" NumberOfComponents=""3"" format=""ascii"">\n");  
  fprintf (fid, "%.17g %.17g %.17g\n", p);
  fprintf (fid, "</DataArray>\n");
  fprintf (fid, "</Points>\n");

  ## VTK-Cells (mesh elements)
  fprintf (fid, "<Cells>\n");
  fprintf (fid, "<DataArray type=""Int32"" Name=""connectivity"" format=""ascii"">\n");
  if dim == 2
    fprintf (fid, "%d %d %d %d \n", t);
  else
    fprintf (fid, "%d %d %d %d %d %d %d %d \n", t);
  endif
  fprintf (fid, "</DataArray>\n");
  fprintf (fid, "<DataArray type=""Int32"" Name=""offsets"" format=""ascii"">\n");
  fprintf (fid, "%d %d %d %d %d %d\n", (2^dim):(2^dim):((2^dim)*nelems));
  fprintf (fid, "</DataArray>\n");
  fprintf (fid, "<DataArray type=""Int32"" Name=""types"" format=""ascii"">\n");
  fprintf (fid, "%d %d %d %d %d %d\n", eltype*ones(nelems,1));
  fprintf (fid, "</DataArray>\n");
  fprintf (fid, "</Cells>\n");

endfunction

## Print DataPoints
function print_data_points (fid, nodedata, nnodes)
  
  ## # of data to print in 
  ## <PointData> field
  nvdata = size (nodedata, 1);  
  
  if (nvdata)
    fprintf (fid, "<PointData>\n");   
    for ii = 1:nvdata    
      data     = nodedata{ii,1};
      dataname = nodedata{ii,2};
      nsamples = rows (data);
      ncomp    = columns (data);
      if (nsamples != nnodes)
	error ("fpl_vtk_write_field: wrong number of samples in <PointData> ""%s""", dataname);
      endif
      fprintf (fid, "<DataArray type=""Float64"" Name=""%s"" ", dataname);
      fprintf (fid, "NumberOfComponents=""%d"" format=""ascii"">\n", ncomp);
	fprintf (fid, "%.17g ", data.');
	fprintf (fid, "\n");
      fprintf (fid, "</DataArray>\n"); 
    endfor
    fprintf (fid, "</PointData>\n");
  endif

endfunction

function print_cell_data (fid, celldata, nelems)
  
  ## # of data to print in 
  ## <CellData> field
  nvdata = size (celldata, 1); 

  if (nvdata)
    fprintf (fid, "<CellData>\n");
    for ii = 1:nvdata
      data     = celldata{ii,1};
      dataname = celldata{ii,2};
      nsamples = rows(data);
      ncomp    = columns(data);
      if nsamples != nelems
	error ("fpl_vtk_write_field: wrong number of samples in <CellData> ""%s""", dataname);
      endif
      fprintf (fid, "<DataArray type=""Float64"" Name=""%s"" ", dataname);
      fprintf (fid, "NumberOfComponents=""%d"" format=""ascii"">\n", ncomp);
	fprintf (fid, "%.17g ", data.');
	fprintf (fid, "\n");
      fprintf (fid, "</DataArray>\n");
    endfor
    fprintf (fid, "</CellData>\n"); 
  endif

endfunction
