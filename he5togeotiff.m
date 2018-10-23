%The given script converts all *.he5 files in a given folder into
%geotiff in a gridded raster format. It was developed and written to
%convert a swath of Aura OMI data into gridded raster data and export the
%data into Geotiff.
%Place the script inside a folder with two additional folders on the same
%file level. One must be called "input_data", the other one "output_data".
%If the output folder does not exist, the script will create it for you.
clear

% Open the HDF5 File.

current_folder = pwd;
input_data_location = fullfile(pwd, 'input_data');
output_data_location = fullfile(pwd, 'output_data');

if 7==exist(input_data_location,'dir')
    list_files = dir(fullfile(input_data_location,'*.he5'));
    if 7==exist(output_data_location,'dir')
        convertHe5ToGeotiff(list_files);
    else 
        mkdir output_data;
        convertHe5ToGeotiff(list_files);
    end
    
else 
    error('Please create folder called input_data and place the *.he5 files in there')
end

%begin of the loop

function result = convertHe5ToGeotiff(data_list)
    %current_folder= pwd;
    input_data_location = fullfile(pwd, 'input_data');
    
    for j = 1:numel(data_list)
        file = data_list(j).name;
        %disp(isstruct(file));
        %disp(file);
        %file_in_ending = '.he5';
        file_name = fullfile(input_data_location, file);
        %disp(file_name);
        
        file_id = H5F.open(file_name);

        % Optain information about the file- uncomment and add breakpoint
        % to assess the he5 information
        
        %info = h5info(file_name);

        % Open the dataset.
        DATAFIELD_NAME = '/HDFEOS/SWATHS/ColumnAmountNO2/Data Fields/ColumnAmountNO2Trop';
        data_id = H5D.open (file_id, DATAFIELD_NAME);

        Lat_NAME='HDFEOS/SWATHS/ColumnAmountNO2/Geolocation Fields/Latitude';
        lat_id=H5D.open(file_id, Lat_NAME);

        Lon_NAME='HDFEOS/SWATHS/ColumnAmountNO2/Geolocation Fields/Longitude';
        lon_id=H5D.open(file_id, Lon_NAME);

        % Read the dataset.
        data=H5D.read (data_id,'H5T_NATIVE_DOUBLE', 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT');

        lat=H5D.read(lat_id,'H5T_NATIVE_DOUBLE', 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT');

        lon=H5D.read(lon_id,'H5T_NATIVE_DOUBLE', 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT');

        % Read the offset.
        ATTRIBUTE = 'Offset';
        attr_id = H5A.open_name (data_id, ATTRIBUTE);
        offset = H5A.read(attr_id, 'H5ML_DEFAULT');

        % Read the scale.
        ATTRIBUTE = 'ScaleFactor';
        attr_id = H5A.open_name (data_id, ATTRIBUTE);
        scale = H5A.read(attr_id, 'H5ML_DEFAULT');

        % Read the fill value.
        ATTRIBUTE = '_FillValue';
        attr_id = H5A.open_name (data_id, ATTRIBUTE);
        fillvalue=H5A.read (attr_id, 'H5T_NATIVE_DOUBLE');

        % Read the missing value.
        ATTRIBUTE = 'MissingValue';
        attr_id = H5A.open_name (data_id, ATTRIBUTE);
        missingvalue=H5A.read (attr_id, 'H5T_NATIVE_DOUBLE');

        % Read title attribute.
        ATTRIBUTE = 'Title';
        attr_id = H5A.open_name (data_id, ATTRIBUTE);
        long_name=H5A.read (attr_id, 'H5ML_DEFAULT');

        % Close and free up memory.
        H5A.close (attr_id)
        H5D.close (data_id);
        H5F.close (file_id);

        % Replace the fill value with NaN.
        data(data==fillvalue) = NaN;
        % Replace the missing value with NaN.
        data(data==missingvalue) = NaN;

        % Apply scale and offset.
        data = double(scale*(data-offset));

        % Attempt to convert data from molecules/cm^2 to micrograms/m^2
            %data = data*10^4*(1/(6.022*10^23))*46*10^6;


        % Following code partially taken from 
        % https://hdfeos.org/forums/showthread.php?t=694

        %Set the resolution in degrees for the gridded data
        cellsize = .25;

        %--> method 1
        %[Z, R] = vec2mtx(double(lat),double(lon), data, cellsize)

        [Z, refvec] = geoloc2grid(double(lat),double(lon), data, cellsize);

        %--> method 2

        % method = 'linear';
        % F = TriScatteredInterp(double(lon(:)),double(lat(:)),data(:),method);
        % halfcell = cellsize/2;
        % [lonmesh, latmesh] = meshgrid( ...
        %     ((lonlim(1)+halfcell):cellsize:(lonlim(2)-halfcell)),...
        %     ((latlim(1)+halfcell):cellsize:(latlim(2)-halfcell))');
        % Z = F(lonmesh, latmesh);
        % refvec = [1/cellsize, latlim(1) + cellsize*size(Z,1), lonlim(1)];

        %--> methode 3

        %[Z, refvec] = gridmap(double(lat),double(lon), data, cellsize);

        R = refvecToGeoRasterReference(refvec,size(Z));
        file_new = strrep(file, '.he5','.tif');
        file_out_location = fullfile(pwd, 'output_data', file_new);
        geotiffwrite(file_out_location,Z,R)
    end

    result = 'true';
end