function plotScResult(maps,names,x,z,env,geom,outFile,overlay)

    % Display maintained B/A maps with optional B-mode and metrics overlays
    fig = figure('Visible','off','Color','w','Units','inches','Position',[1,1,3.1*numel(maps),4.1]);
    tiledlayout(fig,1,numel(maps),'TileSpacing','compact','Padding','compact')
    for iMap = 1:numel(maps)
        ax = nexttile;
        if overlay
            plotBmode(ax,x,z,env,geom.viewBox,60,[]);
            mask = false(size(maps{iMap}));
            for iRoi = 1:numel(geom.masks)
                mask = mask | geom.masks{iRoi};
            end
            plotBaMap(ax,x,z,maps{iMap},[4,13],mask,0.8);
        else
            plotBaMap(ax,x,z,maps{iMap},[4,13],[],1);
        end
        xlim(ax,geom.viewBox(1:2)*100)
        ylim(ax,geom.viewBox(3:4)*100)
        title(ax,names{iMap},'FontWeight','normal')
        cb = colorbar(ax);
        cb.Label.String = 'B/A';
        if overlay
            for iRoi = 1:numel(geom.masks)
                rect = geom.rects(iRoi,:);
                rectangle(ax,'Position',rect*100,'EdgeColor','k','LineWidth',1.5)
                val = mean(maps{iMap}(geom.masks{iRoi}),'omitnan');
                fontSize = 9;
                rotation = 0;
                if rect(3) < 3e-3
                    fontSize = 6;
                    rotation = 90;
                elseif rect(3) < 5e-3
                    fontSize = 7;
                    rotation = 90;
                end
                text(ax,(rect(1)+rect(3)/2)*100,(rect(2)+rect(4)/2)*100,sprintf('%.2f',val),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',fontSize,'Rotation',rotation,'FontWeight','bold')
            end
            for iInc = 1:size(geom.circles,1)
                circ = geom.circles(iInc,:);
                rectangle(ax,'Position',[circ(1:2)-circ(3),2*circ(3),2*circ(3)]*100,'Curvature',[1,1],'LineStyle','--','EdgeColor','k','LineWidth',1.5)
            end
        end
    end
    exportgraphics(fig,outFile,'Resolution',200)
    close(fig)

end
