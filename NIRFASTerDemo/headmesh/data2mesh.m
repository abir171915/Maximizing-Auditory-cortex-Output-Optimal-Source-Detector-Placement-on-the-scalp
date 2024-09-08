function [datamesh, mesh] = data2mesh(data, mesh)

if ~isfield(mesh, 'vol')
    warning('Field "vol" not found. Calculating integration matrices using the grid information in data');
    mesh = gen_intmat(mesh, data.xgrid, data.ygrid, data.zgrid);
end

allfields = fieldnames(data);
dim = max(length(data.xgrid),1) * max(length(data.ygrid),1) * max(length(data.zgrid),1);
for i=1:length(allfields)
    if ~isempty(strfind(allfields{i}, 'phi'))
        try
            eval(['datamesh.',allfields{i},'= mesh.vol.grid2mesh * reshape(data.',allfields{i},', dim, size(data.',allfields{i},', ndims(data.',allfields{i},')));']);
        catch
            % if the field is not calculated, set it to zero
            eval(['datamesh.',allfields{i},'=0;']);
        end
    elseif isempty(strfind(allfields{i}, 'grid'))
        % Copy over the other fields - they should not change
        eval(['datamesh.',allfields{i},'= data.',allfields{i},';']);
    end
end

