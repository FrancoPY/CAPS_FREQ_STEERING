function sampleDef = getSampleAcousticParams(sampleName)

    % getSampleAcousticParams  Return acoustic properties for a sample
    %
    % Inputs
    %   sampleName : sample identifier
    %
    % Output
    %   sampleDef  : struct with fields ac, alphaPower, speedOfSound, and boa

    % Sample table
    switch lower(sampleName)

        case 'phantoma'
            sampleDef.ac = 0.34 / 8.686;
            sampleDef.alphaPower = 1.18;
            sampleDef.speedOfSound = 1534.2;
            sampleDef.boa = 6.75;
            sampleDef.profileType = 'homogeneous';

        case 'phantoma_6mhz'
            sampleDef.ac = 0.34 / 8.686;
            sampleDef.alphaPower = 1.18;
            sampleDef.speedOfSound = 1534.2;
            sampleDef.boa = 6.75;
            sampleDef.profileType = 'homogeneous';

        case 'phantoma_7mhz'
            sampleDef.ac = 0.34 / 8.686;
            sampleDef.alphaPower = 1.18;
            sampleDef.speedOfSound = 1534.2;
            sampleDef.boa = 6.75;
            sampleDef.profileType = 'homogeneous';

        case 'phantomb'
            sampleDef.ac = 0.70 / 8.686;
            sampleDef.alphaPower = 1;
            sampleDef.speedOfSound = 1538.3;
            sampleDef.boa = 6.78;
            sampleDef.profileType = 'homogeneous';

        case 'phantomb_6mhz'
            sampleDef.ac = 0.70 / 8.686;
            sampleDef.alphaPower = 1;
            sampleDef.speedOfSound = 1538.3;
            sampleDef.boa = 6.78;
            sampleDef.profileType = 'homogeneous';
        
        case 'phantomb_7mhz'
            sampleDef.ac = 0.70 / 8.686;
            sampleDef.alphaPower = 1;
            sampleDef.speedOfSound = 1538.3;
            sampleDef.boa = 6.78;
            sampleDef.profileType = 'homogeneous';

        case 'phantomagarv1'
            sampleDef.ac = 0.23 / 8.686;
            sampleDef.alphaPower = 1.13;
            sampleDef.speedOfSound = 1495.9;
            sampleDef.boa = 5.87;
            sampleDef.profileType = 'homogeneous';

        case 'phantomagarv1_6mhz'
            sampleDef.ac = 0.23 / 8.686;
            sampleDef.alphaPower = 1.13;
            sampleDef.speedOfSound = 1495.9;
            sampleDef.boa = 5.87;
            sampleDef.profileType = 'homogeneous';

        case 'phantomagarv1_7mhz'
            sampleDef.ac = 0.23 / 8.686;
            sampleDef.alphaPower = 1.13;
            sampleDef.speedOfSound = 1495.9;
            sampleDef.boa = 5.87;
            sampleDef.profileType = 'homogeneous';

        case 'phantomagarv1f2acq2'
            sampleDef.ac = 0.23 / 8.686;
            sampleDef.alphaPower = 1.13;
            sampleDef.speedOfSound = 1495.9;
            sampleDef.boa = 5.87;
            sampleDef.profileType = 'homogeneous';

        case 'phantomoilv1'
            sampleDef.ac = 0.83 / 8.686;
            sampleDef.alphaPower = 1.00;
            sampleDef.speedOfSound = 1470.4;
            sampleDef.boa = 9.21;
            sampleDef.profileType = 'homogeneous';

        case 'phantomoilv1_6mhz'
            sampleDef.ac = 0.83 / 8.686;
            sampleDef.alphaPower = 1.00;
            sampleDef.speedOfSound = 1470.4;
            sampleDef.boa = 9.21;
            sampleDef.profileType = 'homogeneous';

        case 'phantomoilv1_7mhz'
            sampleDef.ac = 0.83 / 8.686;
            sampleDef.alphaPower = 1.00;
            sampleDef.speedOfSound = 1470.4;
            sampleDef.boa = 9.21;
            sampleDef.profileType = 'homogeneous';

        case 'phantomoilv1f2'
            sampleDef.ac = 0.83 / 8.686;
            sampleDef.alphaPower = 1.00;
            sampleDef.speedOfSound = 1470.4;
            sampleDef.boa = 9.21;
            sampleDef.profileType = 'homogeneous';

        case 'layer69'
            sampleDef.ac = 0.1 / 8.686;
            sampleDef.alphaPower = 2.00;
            sampleDef.speedOfSound = 1500.0;
            sampleDef.boa = 9.00;
            sampleDef.profileType = 'homogeneous';

        case 'reference'
            sampleDef.ac = 0.1 / 8.686;
            sampleDef.alphaPower = 2.00;
            sampleDef.speedOfSound = 1500.0;
            sampleDef.boa = 6.00;
            sampleDef.profileType = 'homogeneous';
        
        case 'phantomincv1d1'
            sampleDef.profileType = 'inclusion';
        
            sampleDef.acInc  = 0.83 / 8.686;
            sampleDef.acBgnd = 0.23 / 8.686;
            sampleDef.alphaPowerInc  = 1.00;
            sampleDef.alphaPowerBgnd = 1.13;
            sampleDef.speedOfSoundInc  = 1470.4;
            sampleDef.speedOfSoundBgnd = 1495.9;
            sampleDef.boaInc  = 9.21;
            sampleDef.boaBgnd = 5.87;
        
            sampleDef.ac = 0.83 / 8.686;
            sampleDef.alphaPower = 1.00;
            sampleDef.speedOfSound = 1470.4;
            sampleDef.boa = 9.21;
        
            sampleDef.inclusionCenter = [19.2e-3, 22.50e-3];
            sampleDef.inclusionRadius = 7.5e-3;

        case 'phantomincv1d1_7mhz'
            sampleDef.profileType = 'inclusion';
        
            sampleDef.acInc  = 0.83 / 8.686;
            sampleDef.acBgnd = 0.23 / 8.686;
            sampleDef.alphaPowerInc  = 1.00;
            sampleDef.alphaPowerBgnd = 1.13;
            sampleDef.speedOfSoundInc  = 1470.4;
            sampleDef.speedOfSoundBgnd = 1495.9;
            sampleDef.boaInc  = 9.21;
            sampleDef.boaBgnd = 5.87;
        
            sampleDef.ac = 0.83 / 8.686;
            sampleDef.alphaPower = 1.00;
            sampleDef.speedOfSound = 1470.4;
            sampleDef.boa = 9.21;

            sampleDef.inclusionCenter = [19.46e-3, 22.24e-3];
            sampleDef.inclusionRadius = 8.33e-3;

        case 'phantomincv1d2'
            sampleDef.profileType = 'inclusion';
        
            sampleDef.acInc  = 0.83 / 8.686;
            sampleDef.acBgnd = 0.23 / 8.686;
            sampleDef.alphaPowerInc  = 1.00;
            sampleDef.alphaPowerBgnd = 1.13;
            sampleDef.speedOfSoundInc  = 1470.4;
            sampleDef.speedOfSoundBgnd = 1495.9;
            sampleDef.boaInc  = 9.21;
            sampleDef.boaBgnd = 5.87;
        
            sampleDef.ac = 0.83 / 8.686;
            sampleDef.alphaPower = 1.00;
            sampleDef.speedOfSound = 1470.4;
            sampleDef.boa = 9.21;
        
            sampleDef.inclusionCenter = [18.335e-3, 29.245e-3];
            sampleDef.inclusionRadius = 8.193e-3;

        case 'inc9'
            sampleDef.profileType = 'inclusion';
        
            sampleDef.acInc  = 0.16 / 8.686;
            sampleDef.acBgnd = 0.10 / 8.686;
            sampleDef.alphaPowerInc  = 2.00;
            sampleDef.alphaPowerBgnd = 2.00;
            sampleDef.speedOfSoundInc  = 1500.0;
            sampleDef.speedOfSoundBgnd = 1500.0;
            sampleDef.boaInc  = 9.00;
            sampleDef.boaBgnd = 6.00;
        
            sampleDef.ac = 0.16 / 8.686;
            sampleDef.alphaPower = 2.00;
            sampleDef.speedOfSound = 1500.0;
            sampleDef.boa = 9.00;
        
            sampleDef.inclusionCenter = [19.46e-3, 22.24e-3];
            sampleDef.inclusionRadius = 8.33e-3;

        case 'bgnd6'
            sampleDef.ac = 0.10 / 8.686;
            sampleDef.alphaPower = 2.00;
            sampleDef.speedOfSound = 1500.0;
            sampleDef.boa = 6.00;
            sampleDef.profileType = 'homogeneous';

        otherwise
            error('Unknown sample: %s', sampleName);

    end

end