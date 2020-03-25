function [class]=confusion_matrix(lithos_num,lithos,lithos_sm,id)

[N_log,CMP]=size(lithos);

class=zeros(lithos_num,lithos_num);

for i=1:CMP
    for j=1:N_log
        switch lithos(j,i)
            case 1
                class(1,lithos_sm(j,i))=class(1,lithos_sm(j,i))+1;
            case 2
                class(2,lithos_sm(j,i))=class(2,lithos_sm(j,i))+1;
            case 3
                class(3,lithos_sm(j,i))=class(3,lithos_sm(j,i))+1;
            case 4
                class(4,lithos_sm(j,i))=class(4,lithos_sm(j,i))+1;
            case 5
                class(5,lithos_sm(j,i))=class(5,lithos_sm(j,i))+1;
            case 6
                class(6,lithos_sm(j,i))=class(6,lithos_sm(j,i))+1;
            case 7
                class(7,lithos_sm(j,i))=class(7,lithos_sm(j,i))+1;
            case 8
                class(8,lithos_sm(j,i))=class(8,lithos_sm(j,i))+1;
            case 9
                class(9,lithos_sm(j,i))=class(9,lithos_sm(j,i))+1;
            case 10
                class(10,lithos_sm(j,i))=class(10,lithos_sm(j,i))+1;
            case 11
                class(11,lithos_sm(j,i))=class(11,lithos_sm(j,i))+1;
            case 12
                class(12,lithos_sm(j,i))=class(12,lithos_sm(j,i))+1;
        end
    end
end

temp_clas=class;

if id

%plot the confusion matrix
mat = temp_clas;           %# A 5-by-5 matrix of random values from 0 to 1
figure;
% subplot(1,2,2)
imagesc(mat);            %# Create a colored plot of the matrix values
colormap(flipud(jet));  %# Change the colormap to gray (so higher values are black and lower values are white)

textStrings = num2str(mat(:),'%i');  %# Create strings from the matrix values %i for integer
textStrings = strtrim(cellstr(textStrings));  %# Remove any space padding


idx = find(strcmp(textStrings(:), '0'));
textStrings(idx) = {'   '};

[x,y] = meshgrid(1:lithos_num);   %# Create x and y coordinates for the strings
hStrings = text(x(:),y(:),textStrings(:),...      %# Plot the strings
                'HorizontalAlignment','center','VerticalAlignment','Bottom','fontsize',16);
            
midValue = mean(get(gca,'CLim'));  %# Get the middle value of the color range
textColors = repmat(mat(:) > midValue,1,3);  
set(hStrings,{'Color'},num2cell(textColors,2));  %# Change the text colors            
            

textStrings = num2str(mat(:)/sum(sum(temp_clas))*100,'%0.2f'); 
textStrings = strcat(textStrings,'%');
textStrings = strtrim(cellstr(textStrings));
idx = find(strcmp(textStrings(:), '0.00%'));
textStrings(idx) = {'   '};
hStrings = text(x(:),y(:),textStrings(:),'HorizontalAlignment','center','VerticalAlignment','top','fontsize',14);            
            
midValue = mean(get(gca,'CLim'));  %# Get the middle value of the color range
textColors = repmat(mat(:) > midValue,1,3);  %# Choose white or black for the text color of the strings so they can be easily seen over the background color
set(hStrings,{'Color'},num2cell(textColors,2));  %# Change the text colors

set(gca,'XTick',1:12,...                         %# Change the axes tick marks
        'XTickLabel',{'CS\_non','MS\_non','MS','FS\_non','FS','VFS\_non','VFS','SS\_non','SS','Clay\_non','Clay','Coal'},...  %#   and tick labels
        'YTick',1:12,...
        'YTickLabel',{'CS\_non','MS\_non','MS','FS\_non','FS','VFS\_non','VFS','SS\_non','SS','Clay\_non','Clay','Coal'},...
        'TickLength',[0 0],'fontsize',10);
    
% set(gca,'XTick',1:3,...                         %# Change the axes tick marks
%         'XTickLabel',{'Lithology A','Lithology B','Lithology C'},...  %#   and tick labels
%         'YTick',1:3,...
%         'YTickLabel',{'Lithology A','Lithology B','Lithology C'},...
%         'TickLength',[0 0],'fontsize',16);
title('Confusion Matrix','fontsize',20); 
% s=sprintf('1%c', char(176));
% title(s,'fontsize',24);
xlabel({'Predicted Lithologies'},'fontsize',18);
ylabel('True Lithologies ','fontsize',18);

end


end