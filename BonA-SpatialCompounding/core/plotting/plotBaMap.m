function h = plotBaMap(ax,x,z,ba,baLim,roiMask,opacity)

    % Draw a B/A map in centimeters from coordinates in meters
    % Use [] for the whole map, or a geometric metrics mask for an overlay
    valid = isfinite(ba);
    if ~isempty(roiMask)
        valid = valid & roiMask;
    end
    if isvector(x) && isvector(z)
        h = imagesc(ax,x * 100,z * 100,ba);
        h.AlphaData = opacity * double(valid);
    else
        h = surface(ax,x * 100,z * 100,ones(size(x)),ba,'FaceColor','flat','EdgeColor','none','FaceAlpha','flat','AlphaDataMapping','none');
        h.AlphaData = opacity * double(valid);
        view(ax,2)
    end
    daspect(ax,[1,1,1])
    set(ax,'YDir','reverse','Color',[0.8,0.8,0.8],'FontName','Arial','FontSize',10,'LineWidth',0.9,'TickDir','out','Box','on','Layer','top')
    colormap(ax,turbo(256))
    clim(ax,baLim)
    xlabel(ax,'Lateral (cm)','FontWeight','bold')
    ylabel(ax,'Depth (cm)','FontWeight','bold')
    hold(ax,'on')

end
