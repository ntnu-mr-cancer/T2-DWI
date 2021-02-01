% ----------------------------------------------------------------------- %
% This function calculates and saves the mono-exponential ADC for a       %
% patient, using data from combined T2-DWI with two TEs and two b-values. %
% ----------------------------------------------------------------------- %

function calculateADC(patientNr, patientDataFilePath, dicomFilePath, coRegistrationFilePath, resultsFilePath, noiseThreshold, useCoRegisteredImages)

    
    % Find patient name string
    patientNrString = getPatientNrString(patientNr);


    % Write to output
    fprintf('\nMono-exponential ADC analysis of patient %s started %s\n', patientNrString, datetime(now,'ConvertFrom','datenum'));
    
    
    % Load images
    if useCoRegisteredImages
        
        [im_shortTE_lowB,im_shortTE_highB,im_longTE_lowB,im_longTE_highB] = loadCoRegisteredImages(patientNrString,coRegistrationFilePath);
        
    else
        
        [im_shortTE_lowB,im_shortTE_highB,im_longTE_lowB,im_longTE_highB] = loadDicoms(patientNrString,dicomFilePath);
        
    end
    
    
    % Extract image size
    [nRows, nCols, nSlices] = size(im_shortTE_lowB);


    % Take the logarithm of the images to prepare for linear regression
    ln_im_longTE_lowB = log(im_longTE_lowB);
    ln_im_longTE_highB = log(im_longTE_highB);
    
    
    % Find x and y limits of the box ROI containing the prostate
    [x_min,y_min,x_max,y_max] = findBoxROI(patientNr,patientDataFilePath);


    % Specify b-values
    b = [50; 700]; % s/mm²

    
    % Initialize ADC map
    map_ADC = NaN(nRows, nCols, nSlices);
    
    
    % Make variable that contains the number of excluded pixels
    nExcludedPixels = 0;
    
    
    % Specify that the fitting method is linear least squares
    fo = fitoptions('Method','LinearLeastSquares');
 
    
    % Start clock to measure time of model fit
    tic
    
    % Loops over slices, rows and columns to calculate mono-exponential ADC
    % for each voxel inside specified ROI ("box" containing the prostate)
     for slice = 1:nSlices
        for row = y_min:y_max % row = y
            for col = x_min:x_max % col = x
                
                % Check that the pixel values are above the threshold
                % value, and ensure that the ADC and T2 values will be
                % positive (increasing signal intensity with increasing
                % b-value or TE is considered "unphysical")
                if (im_longTE_highB(row,col,slice) > noiseThreshold) && ...
                        (im_shortTE_lowB(row,col,slice) > im_shortTE_highB(row,col,slice)) && ...
                        (im_longTE_lowB(row,col,slice) > im_longTE_highB(row,col,slice)) && ...
                        (im_shortTE_lowB(row,col,slice) > im_longTE_lowB(row,col,slice)) && ...
                        (im_shortTE_highB(row,col,slice) > im_longTE_highB(row,col,slice))
                    
                    % Linear regression to find ADC at the long TE
                    f = fit(b,[ln_im_longTE_lowB(row,col,slice); ln_im_longTE_highB(row,col,slice)],'poly1',fo); % fit curve
                    map_ADC(row,col,slice) = -f.p1*1E3; % ADC is the negative slope (unit: 1E-3 mm^2/s)
                    
                else
                    
                    % Exclude pixel
                    nExcludedPixels = nExcludedPixels + 1;
                    
                end % if

            end % for col
        end % for row
     end % for slice
    
    % Write to output
    fprintf('Time of model fit: %.f seconds.\n',toc);
    
    
    % Make ADC value array from the map, and remove NaN values
    values_ADC = rmmissing(reshape(map_ADC,[],1));

    
    % Save variables (and make directory if it does not already exist)
    if ~isfolder([resultsFilePath patientNrString])
        mkdir([resultsFilePath patientNrString])
    end
    
    if useCoRegisteredImages
        save([resultsFilePath patientNrString '/resultsADC_coRegistered.mat'], 'map_ADC');
        save([resultsFilePath patientNrString '/resultsADC_coRegistered_allVariables.mat']);
    else
        save([resultsFilePath patientNrString '/resultsADC.mat'], 'map_ADC');
        save([resultsFilePath patientNrString '/resultsADC_allVariables.mat']);
    end

    
end % function