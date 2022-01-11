

accept_dir = fullfile(output_dir,'accepted');
reject_dir = fullfile(output_dir,'rejected');
idk_dir = fullfile(output_dir,'idk');

if ~isfolder(accept_dir)
    mkdir(accept_dir);
end
if ~isfolder(reject_dir)
    mkdir(reject_dir);
end
if ~isfolder(idk_dir)
    mkdir(idk_dir);
end

for i = 1:numel(images)
    while 1
        fprintf('Input y/n\n');
        k = waitforbuttonpress;
        value = double(get(gcf,'CurrentCharacter')); hold on;
        switch value
            case 121         
                fprintf('Yes\n');
                text(0.5,0.5,'yes','Color','g','Fontsize',20);
                movefile(image,accept_dir)
                break
            case 110
                fprintf('No\n');
                text(0.5,0.5,'No','Color','r','Fontsize',20);
                movefile(image,reject_dir)
                break
            case 100
                fprintf("Don't know\n");
                text(0.5,0.5,"Don't know",'Color','r','Fontsize',20);
                movefile(image,idk_dir)
                break
            otherwise
                fprintf('Error: Press y/n\n');
        end
    end
    pause(1);
    close;
end