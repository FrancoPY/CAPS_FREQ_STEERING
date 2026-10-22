function [h, refVal, db] = plotBmode(ax,x,z,env,viewBox,dynRange,refVal)

    % Draw a bone B-mode from a linear envelope and coordinates in meters
    % Use [xMin,xMax,zMin,zMax] for viewBox and [] for a local 0 dB reference
    if isvector(x) && isvector(z)
        [X,Z] = meshgrid(x,z);
    else
        X = x;
        Z = z;
    end
    viewMask = X >= viewBox(1) & X <= viewBox(2) & Z >= viewBox(3) & Z <= viewBox(4) & isfinite(env);
    if isempty(refVal)
        refVal = max(env(viewMask));
    end
    assert(isscalar(refVal) && isfinite(refVal) && refVal > 0,'The displayed B-mode region must contain a positive envelope.');
    db = 20 * log10(double(env) / double(refVal));
    db = min(0,max(-dynRange,db));
    db(~isfinite(db)) = -dynRange;
    rgb = ind2rgb(1 + round(255 * (db + dynRange) / dynRange),bone(256));

    % Preserve linear or curved acquisition geometry without resampling
    if isvector(x) && isvector(z)
        h = imagesc(ax,x * 100,z * 100,rgb);
    else
        h = surface(ax,X * 100,Z * 100,zeros(size(X)),rgb,'FaceColor','texturemap','EdgeColor','none');
        view(ax,2)
    end
    axis(ax,'image')
    set(ax,'YDir','reverse','FontName','Arial','FontSize',10,'LineWidth',0.9,'TickDir','out','Box','on','Layer','top')
    xlim(ax,viewBox(1:2) * 100)
    ylim(ax,viewBox(3:4) * 100)
    xlabel(ax,'Lateral (cm)','FontWeight','bold')
    ylabel(ax,'Depth (cm)','FontWeight','bold')
    colormap(ax,bone(256))
    clim(ax,[-dynRange,0])
    hold(ax,'on')

end
