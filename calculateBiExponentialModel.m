% ----------------------------------------------------------------------- %
% This function calculates and saves the values from the bi-exponential   %
% model (volume fractions) for a patient, using data from combined T2-DWI %
% with two TEs and two b-values.                                          %
% ----------------------------------------------------------------------- %

function calculateBiExponentialModel(patientNr, patientDataFilePath, dicomFilePath, coRegistrationFilePath, resultsFilePath, noiseThreshold, useCoRegisteredImages)


    % Find patient name string
    patientNrString = getPatientNrString(patientNr);


    % Write to output
    fprintf('\nBi-exponential model analysis of patient %s started %s\n', patientNrString, datetime(now,'ConvertFrom','datenum'));
    
    
    % Load images
    if useCoRegisteredImages
        
        [im_shortTE_lowB,im_shortTE_highB,im_longTE_lowB,im_longTE_highB] = loadCoRegisteredImages(patientNrString,coRegistrationFilePath);
        
    else
        
        [im_shortTE_lowB,im_shortTE_highB,im_longTE_lowB,im_longTE_highB] = loadDicoms(patientNrString,dicomFilePath);
        
    end
    
    
    % Extract image size
    [nRows, nCols, nSlices] = size(im_shortTE_lowB);
    
    
    % Find x and y limits of the box ROI containing the prostate
    [x_min,y_min,x_max,y_max] = findBoxROI(patientNr,patientDataFilePath);
    
    
    % Specify b-values
    b = [50e-3 700e-3]; % s/mm²
    
    
    % Initialize maps
    map_SF_fast = NaN(nRows, nCols, nSlices);
    map_SI_0 = NaN(nRows, nCols, nSlices);
    map_gof_sse = NaN(nRows, nCols, nSlices);
    map_gof_rsquare = NaN(nRows, nCols, nSlices);
    map_gof_dfe = NaN(nRows, nCols, nSlices);
    map_gof_adjrsquare = NaN(nRows, nCols, nSlices);
    map_gof_rmse = NaN(nRows, nCols, nSlices);
    
    
    % Make variable that contains the number of excluded pixels
    nExcludedPixels = 0;
    
    
    % Start clock to measure time of model fit
    tic
    
    % Loops over slices, rows and columns to fit the bi-exponential model
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
                    
                    % Perform model fit using the data at the longest TE
                    [fitresult, gof] = hybridfit_BEmodel(b, [im_longTE_lowB(row,col,slice) im_longTE_highB(row,col,slice)]);

                    % Fill maps with values from fit
                    map_SF_fast(row,col,slice) = fitresult.SF_fast;
                    map_SI_0(row,col,slice) = fitresult.SI_0;
                    map_gof_sse(row,col,slice) = gof.sse;
                    map_gof_rsquare(row,col,slice) = gof.rsquare;
                    map_gof_dfe(row,col,slice) = gof.dfe;
                    map_gof_adjrsquare(row,col,slice) = gof.adjrsquare;
                    map_gof_rmse(row,col,slice) = gof.rmse;
                    
                else
                    
                    % Exclude pixel
                    nExcludedPixels = nExcludedPixels + 1;
                    
                end % if

            end % for col
        end % for row
     end % for slice
    
    % Write to output
    fprintf('Time of model fit: %.f seconds.\n',toc);
    
    
    % Make SF_slow map (SF_slow + SF_fast = 1)
    map_SF_slow = 1 - map_SF_fast;
    
    % Make value arrays from the maps, and remove NaN values
    values_SF_slow = rmmissing(reshape(map_SF_slow,[],1));
    values_SF_fast = rmmissing(reshape(map_SF_fast,[],1));
    values_SI_0 = rmmissing(reshape(map_SI_0,[],1));
    values_gof_sse = rmmissing(reshape(map_gof_sse,[],1));
    values_gof_rsquare = rmmissing(reshape(map_gof_rsquare,[],1));
    values_gof_dfe = rmmissing(reshape(map_gof_dfe,[],1));
    values_gof_adjrsquare = rmmissing(reshape(map_gof_adjrsquare,[],1));
    values_gof_rmse = rmmissing(reshape(map_gof_rmse,[],1));
    
    
    % Save variables (and make directory if it does not already exist)
    if ~isfolder([resultsFilePath patientNrString])
        mkdir([resultsFilePath patientNrString])
    end
    
    if useCoRegisteredImages
        save([resultsFilePath patientNrString '/resultsBiExponentialModel_coRegistered.mat'], 'map_SF_slow', 'map_SF_fast', 'map_SI_0', 'map_gof_sse', 'map_gof_rsquare', 'map_gof_dfe', 'map_gof_adjrsquare', 'map_gof_rmse');
        save([resultsFilePath patientNrString '/resultsBiExponentialModel_coRegistered_allVariables.mat']);
    else
        save([resultsFilePath patientNrString '/resultsBiExponentialModel.mat'], 'map_SF_slow', 'map_SF_fast', 'map_SI_0', 'map_gof_sse', 'map_gof_rsquare', 'map_gof_dfe', 'map_gof_adjrsquare', 'map_gof_rmse');
        save([resultsFilePath patientNrString '/resultsBiExponentialModel_allVariables.mat']);
    end


end