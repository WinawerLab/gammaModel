
% Pretend contrast energy for the following hypothetical stimuli:
% rows: there are 9 locations (3x3)
% columns: there are 4 orientations

stim=zeros(7,9,4);

% Number of orientations:
% grating
% plaid
% squiggles

% grating:
stim(1,:,:) = repmat([1 0 0 0],9,1);
% plaid:
stim(2,:,:) = repmat([.5 0 .5 0],9,1);
% random:
temp_s = repmat([.5 0 .5 0],9,1);
stim(3,:,:) = reshape(temp_s(randperm(length(temp_s(:)))),9,4);
clear temp_s

% Variations across space
% size: small grating
stim(4,:,:) = stim(1,:,:);
stim(4,[1:4 6:9],:) = 0;

% pattern
stim(5,:,:) = stim(2,:,:);
stim(5,[1:4 6:9],:) = 0;

% squiggles
stim(6,:,:) = stim(3,:,:);
stim(6,[1:4 6:9],:) = 0;

% grating surrounded by other orientation
stim(7,:,:) = repmat([0 1 0 0],9,1);
stim(7,5,:) = [1 0 0 0];

%%
cm = gray(150); cm = cm(end:-1:1,:);
figure('Position',[0 0 300 600])
for kk = 1:size(stim,1)
    thisStim = reshape(squeeze(stim(kk,:,:)),3,3,4);
    subplot(size(stim,1),1,kk),hold on
    for rr = 1:4 % orientation   
        for ii = 1:size(thisStim(:,:,rr),1)
            for jj = 1:size(thisStim(:,:,rr),2)
                if thisStim(ii,jj,rr)>0
                    if rr==1
                        plot([ii ii],[jj-.4 jj+.4],'-',...
                            'Color',cm(thisStim(ii,jj,rr)*100+50,:),...
                            'LineWidth',thisStim(ii,jj,rr)*3)
                    elseif rr==2
                        plot([ii-.4 ii+.4],[jj-.4 jj+.4],'-',...
                            'Color',cm(thisStim(ii,jj,rr)*100+50,:),...
                            'LineWidth',thisStim(ii,jj,rr)*3)
                    elseif rr==3
                        plot([ii-.4 ii+.4],[jj jj],'-',...
                            'Color',cm(thisStim(ii,jj,rr)*100+50,:),...
                            'LineWidth',thisStim(ii,jj,rr)*3)
                    elseif rr==4
                        plot([ii+.4 ii-.4],[jj-.4 jj+.4],'-',...
                            'Color',cm(thisStim(ii,jj,rr)*100+50,:),...
                            'LineWidth',thisStim(ii,jj,rr)*3)
                    end
                end
            end
        end
        xlim([.5 3.5]),ylim([0.5 3.5])
        axis square, box on
        set(gca,'XTick',[],'YTick',[])
        set(gca,'Color',[1 1 1])
    end
    % sum over space 
    a = squeeze(sum(sum(thisStim,1),2));
%     title(num2str(std(a(:))))
end

set(gcf,'PaperPositionMode','auto')
print('-dpng','-r300',['./figures/testStim/testStim_contrastE2'])
print('-depsc','-r300',['./figures/testStim/testStim_contrastE2'])
