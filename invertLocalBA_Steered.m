function BAL = invertLocalBA_Steered(BAeff, mu, invParam, thetaDeg)
% invertLocalBA_Steered  B/A local de UN ángulo de steering (paso previo al compounding)
%
% Idea (Sec. III-C del paper): con steering, la acumulacion de B/A ocurre a lo
% largo de rayos inclinados theta, no a lo largo de las columnas de la imagen.
% Hay que re-alinear la imagen para que cada COLUMNA sea un rayo, invertir
% columna a columna (invertBAeffTV) y volver al marco original.
%
% Aqui el re-alineamiento se hace muestreando la imagen directamente a lo largo
% de cada rayo (equivale a rotar -theta + relleno de ceros del paper, pero con
% la profundidad exacta: l = z/cos(theta), medida desde la cara del transductor).
%
%   rayo que entra en x0:   x(l) = x0 + l*sin(theta),   z(l) = l*cos(theta)
%   pixel (x,z)         ->  l = z/cos(theta),  x0 = x - z*tan(theta)
%
% Convencion de signo: theta > 0 inclina el haz hacia +x (la misma que usan
% bfPlaneWaveSimu: projDist = z*cos + x*sin).
%
% INPUTS
%   BAeff    : struct con .image (B/A acumulado de ESTE angulo, [nz x nx]),
%              .lateral [m] (1 x nx), .axial [m] (nz x 1, profundidad vertical)
%   mu       : peso TV
%   invParam : zCrop, zInv ([m], a lo largo del RAYO), gridSize [m], etc.
%              (los mismos campos que usa invertBAeffTV)
%   thetaDeg : angulo de steering en grados
%
% OUTPUT  BAL
%   .image   : B/A local en el marco ORIGINAL, NaN fuera de la zona valida
%   .lateral : eje x [m]  \  misma grilla para todos los angulos
%   .axial   : eje z [m]  /  (permite promediar directamente)
%   .BAloc   : salida cruda de invertBAeffTV (en el marco de rayos, diagnostico)

    if nargin < 3 || isempty(invParam); invParam = struct; end
    if nargin < 4; thetaDeg = 0; end

    gridSize = getOptionLocal(invParam, 'gridSize', 0.3e-3);
    zCrop    = getOptionLocal(invParam, 'zCrop', [BAeff.axial(1), BAeff.axial(end)]);
    zInv     = getOptionLocal(invParam, 'zInv', zCrop);

    x   = BAeff.lateral(:).';
    z   = BAeff.axial(:);
    img = BAeff.image;

    th = thetaDeg*pi/180;
    sT = sin(th);  cT = cos(th);  tT = tan(th);

    % ---- 1) Muestrear la imagen a lo largo de los rayos ---------------------
    % Filas = profundidad l a lo largo del rayo; columnas = punto de entrada x0
    % (solo rayos que nacen dentro de la apertura: x0 en [x(1), x(end)]).
    ellMax = min(zInv(2), z(end)/cT);
    ell    = (0:gridSize:ellMax).';
    x0     = x(1):gridSize:x(end);

    [X0, ELL] = meshgrid(x0, ell);
    Cray  = interp2(x, z, img, X0 + ELL*sT, ELL*cT, 'linear', NaN);
    valid = ~isnan(Cray);          % pixeles realmente dentro del FOV

    if ~any(valid(:))
        error('invertLocalBA_Steered: ningun rayo cae dentro de la imagen (theta=%g).', thetaDeg);
    end

    Cfill = Cray;  Cfill(isnan(Cfill)) = 0;

    % ---- 2) Inversion local con TV, columna a columna (rayo a rayo) --------
    BAeffRay.image   = Cfill;
    BAeffRay.lateral = x0;
    BAeffRay.axial   = ell;

    BAloc = invertBAeffTV(BAeffRay, mu, invParam);

    % Mascara de validez sobre la grilla de salida de la inversion
    [XL, LL]  = meshgrid(BAloc.lateral, BAloc.axial);
    validLoc  = interp2(x0, ell, double(valid), XL, LL, 'linear', 0);

    % ---- 3) Volver al marco original (x,z) ---------------------------------
    xOut = BAloc.lateral;
    zOut = BAloc.axial(:);
    [XO, ZO] = meshgrid(xOut, zOut);

    Lq  = ZO ./ cT;          % profundidad a lo largo del rayo
    X0q = XO - ZO .* tT;     % donde entro ese rayo

    img_out = interp2(BAloc.lateral, BAloc.axial, BAloc.image, X0q, Lq, 'linear', NaN);
    v_out   = interp2(BAloc.lateral, BAloc.axial, validLoc,    X0q, Lq, 'linear', 0);
    img_out(v_out < 0.999) = 0;      % fuera de la zona válida -> 0 (como el padding del paper)
    img_out(isnan(img_out)) = 0;     % lo que quede fuera de la grilla de rayos también 0

    BAL.image   = img_out;
    BAL.lateral = xOut;
    BAL.axial   = zOut;
    BAL.BAloc   = BAloc;
end

function val = getOptionLocal(s, fieldName, defaultVal)
    if isfield(s, fieldName) && ~isempty(s.(fieldName))
        val = s.(fieldName);
    else
        val = defaultVal;
    end
end
