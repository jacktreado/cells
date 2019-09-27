function visualizeMultipleCells(figureNum,NCELLS,L,xPos,yPos,l0,bodyColor,lineColor)
%% Function to plot a single frame of single cell
% with lab coordinates xPos and yPos, and colors defined as input RGB
% vectors

% size of repulsive zone
delMag = 1.0;

% plot box
figure(figureNum),clf,hold on, box on;

for nn = 1:NCELLS
    % get positions for this cell
    xp = xPos{:,nn};
    yp = yPos{:,nn};
    
    % get number of vertices
    nv = length(xp);
    
    % get width of rectangle
    del = delMag*l0(nn);
    del = 0.5*del;
    
    % plot pbc patch
    for xx = -3:3
        for yy = -3:3
            if nn > 0
                if (xx > 1 || xx < -1)
                    continue;
                elseif (yy > 1 || yy < -1)
                    continue;
                end
            end
            % create vertices
            V = [xp+L*xx, yp+L*yy];
            F = [1:size(xp,1) 1];
            
            % center of mass
            COM = sum(V)./size(V,1);
            
            % add a little bit to vertices to extend size
            relPos = V - COM;
            unitRelPos = relPos./sqrt(sum(relPos.^2,2));
            Vfull = V + 0.5*del.*unitRelPos; 

            % draw full patch
            patch('Faces',F,'Vertices',Vfull,'FaceColor',bodyColor(nn,:),'EdgeColor',lineColor(nn,:));
            
        end
    end
    
end

% Box
plot([0 L L 0 0],[0 0 L L 0],'k');
axis('equal'); box on; set(gca,'XTick',[],'YTick',[])
% axis([0 L 0 L]);
axis([-L/4 (5/4)*L -L/4 (5/4)*L]);
drawnow;

end



