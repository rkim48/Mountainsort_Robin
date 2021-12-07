clu_last=0;

            tBefore=30;
            tAfter=45;
for layer_index=layernum
    %input_filename=strcat('test_190920_012044.rhdnoCMR_tracking-',num2str(layer_index));
%     input_filename=['128Layer' num2str(layer_index) '_tracking.mat'];
        input_filename=[recordingFolders(scheduleArray(fol,1)).folder '/' recordingFolders(scheduleArray(fol,1)).name 
            '_tracking-' num2str(layer_index*offset) '.mat'];


    load(input_filename);
    mkdir (input_filename(1:end-4))

    availCh = result.Dat_V_Map(:,2);
    
    list=1:numel(result.waveform);
    
    
    for clu=clu_last+1:clu_last+numel(list)
        tic
        proceed = 1;
        newList(clu)=1;
        
        
        if proceed == 1
            clu
            
%           row=[] 
%           col=[]
            pdfCellArray{clu}=[num2str(list(clu-clu_last)) '.pdf'];
            FiringTimeForThisUnit = result.time{clu-clu_last};
            
            ISI_in_MS_bins = diff(double(FiringTimeForThisUnit))/(Fs/1000); % Intervals in miliseconds.
            [counts] = histc(ISI_in_MS_bins,0:1:50);
            [counts_M] = histc(ISI_in_MS_bins,0:20:1000);
            [counts_L] = histc(ISI_in_MS_bins,0:1000:20000);
            
            
            Sec_bins = double(FiringTimeForThisUnit)/Fs;
            [counts_5s_1] = histc(Sec_bins,0:5:720);
            [counts_5s_2] = histc(Sec_bins,2.5:5:722.5);
            Max_Ins_5s_FR = max(max(counts_5s_1),max(counts_5s_2))/5;
            
            
            
            
            % ValleyAcrossTime = min(reshape(unit,[size(unit,1)*size(unit,2),size(unit,3)]));
             
              MeanPlot=result.waveform{clu-clu_last};
              
              StdPlot=result.waveformStd{clu-clu_last};


            % Mask_clu=Mask(cluster==list(clu),:);
            % if size(Mask_clu,1)~=1
            % Intensity=max(mean(Mask_clu),median(Mask_clu));
            % else
            % Intensity=Mask_clu;
            % end
            P2P = max(MeanPlot')-min(MeanPlot');
            ch=find(P2P==max(P2P));
            ch=ch(1);
            [sortP2P,sP2P]=sort(P2P,'descend');
            PeakAcrossTime = result.PeakAcrossTime{clu-clu_last};
            ValleyAcrossTime = result.ValleyAcrossTime{clu-clu_last};
            
            
            layer_index_r=5-layer_index;
            Valley = min(MeanPlot');
            Weighted_P2P = P2P.^2;
%             clear Weighted_P2P_Remap Valley_Remap P2P_Remap
            Ch_Map = zeros(8,16);Ch_Map_2 = nan(8,16);
            if HD ==0
                            Ch_Map(layer_index_r*2-1:layer_index_r*2,:)=[result.Ch_Map(1:2,:) result.Ch_Map(3:4,:)];
            StdPlot=zeros(128,tBefore+tAfter+1);
            MeanPlot=zeros(128,tBefore+tAfter+1);
            
            MeanPlot(availCh,:) = result.waveform_CM{clu-clu_last};
            StdPlot(availCh,:) = result.waveformStd_CM{clu-clu_last};
             
             Ch_Map_2(layer_index_r*2-1:layer_index_r*2,:)=[result.Ch_Map(1:2,:) result.Ch_Map(3:4,:)];

            elseif HD==1
                   Ch_Map=result.Ch_Map';
             Ch_Map_2=result.Ch_Map_2';
             Ch_Map = Ch_Map_2;
             
            elseif HD==2
                     Ch_Map = result.Ch_Map(1:16,1:4);
                     Ch_Map = [Ch_Map result.Ch_Map(17:32,1:4)]';
                     
                     Ch_Map_2 = result.Ch_Map_2(1:16,1:4);
                     Ch_Map_2 = [Ch_Map_2 result.Ch_Map_2(17:32,1:4)]';
            end
            
            WCList(clu,:)= result.Location{clu-clu_last};
            

            t=[-tBefore/Fs:1/Fs:tAfter/Fs]*1000;
            % 4-by-8 Channel Plot showing correlation by visual inspection
            for i=1:size(MeanPlot,1)
                Valley =Valley -MeanPlot(i,1);
                MeanPlot(i,:)=MeanPlot(i,:)-MeanPlot(i,1);
            end
            
            %% plot
            tic
            h=figure;
            set(h,'PaperUnits','centimeters','PaperPositionMode','Auto')
            set(h, 'Visible', 'off');
            
            Ox=1;
            Oy=4;
            gap=0.4;
            Gw=1.5;
            Gh=2.25;
            numGx = 16;
            numGy = 8;
            
            Gx=[Ox:(gap+Gw):Ox+(numGx-1)*(gap+Gw)];
            Gy=flip([Oy:(gap+Gh):Oy+(numGy-1)*(gap+Gh)]);
            [X,Y] = meshgrid(Gx,Gy);
            Dat_V_Map=result.Dat_V_Map;
            
            
            for row=1:size(Ch_Map,1)
                for col=1:size(Ch_Map,2)
                    
                    ax=axes('Units','centimeters','position', [X(row,col) Y(row,col) Gw Gh]);
                    set(ax,'FontSize',8)
                    %         set(ax,'FontWeight','bold')
                    
                    if Ch_Map(row,col)~=0
                        temp=Ch_Map(row,col);
                        errorbar(t,MeanPlot(temp,:),StdPlot(temp,:),'Color',[0.9 0.9 0.9])
                        
                        hold on;
                        if min(min(MeanPlot))<-100 || max(max(MeanPlot))>25
                            plot(t,MeanPlot(temp,:),'Color','r','LineWidth',2)
                            %                 axis([t(1) t(end) -150 50])
                            axis([t(1) t(end) min(min(MeanPlot)) max(max(MeanPlot))])
                            
                        else
                            plot(t,MeanPlot(temp,:),'Color','r','LineWidth',1)
                            %                 axis([t(1) t(end) -100 25])
                            axis([t(1) t(end) min(min(MeanPlot)) max(max(MeanPlot))])
                            
                        end
                        xlim([-0.75 1.25])
                        set(gca,'xticklabel',[])
                        
                        
                    end
                    
                    
                end
            end
            
%             Ox=32;
%             Oy=7;
%             
%             gapX=0.3;
%             gapY=1.65;
%             
%             Gw=0.6;
%             Gh=0.6;
%             
%             numGx = 16;
%             numGy = 8;
%             
%             Gx1=[Ox:(gapX+Gw):Ox+(numGx-1)*(gapX+Gw)];
%             Gy1=flip([Oy:(gapY+Gh):Oy+(numGy-1)*(gapY+Gh)]);
%             [X,Y] = meshgrid(Gx1,Gy1);
            
            %Ch_Map=flip(Ch_Map);
            %color=jet(500);
            %cl=randperm(500);
            %color=color(cl,:);
%             for row=1:size(Ch_Map,1)
%                 for col=1:size(Ch_Map,2)
%                     
%                     ax=axes('Units','centimeters','position', [X(row,col) Y(row,col) Gw Gh]);
%                     set(ax,'FontSize',8)
%                     %         set(ax,'FontWeight','bold')
%                     
%                     set(ax,'Box','on')
%                     
%                     set(ax,'LineWidth',2)
% %                     plot(0,0)
% %                     T = text(-0.2,0,num2str(floor(abs(P2P_Remap(row,col)))));
% %                     set(T,'FontSize',6)
% %                     
% %                     set(T,'Color',[0 0 1])
%                     
%                     set(ax,'xtick',[])
%                     set(ax,'ytick',[])
%                     set(ax,'XColor',[0.7 0.7 0.7])
%                     set(ax,'YColor',[0.7 0.7 0.7])
%                     
%                     
% %                     SelectedPriLoc = zeros(size(Ch_Map));
% %                     [rp,cp]= find(Ch_Map_2==result.clusterPriCh(clu-clu_last,2));
% %                     SelectedPriLoc(rp,cp)=1;
% %                     if SelectedPriLoc(row,col)==1
% %                         
% %                         set(ax,'XColor',[0 0 1])
% %                         set(ax,'YColor',[0 0 1])
% %                         
% %                     end
%                     
%                 end
%             end
            
            scaleCoeff = 0.9/60;
            
            
            
            
            ax=axes('Units','centimeters','position', [32 0.6 10 4.5]);
            set(ax,'FontSize',6)
            if ~isempty(counts)
                bar(0.5:1:50.5,counts)
                xlim([0 51])
            end
            
            ax=axes('Units','centimeters','position', [1 0.5 30 2.5]);
            set(ax,'FontSize',4)
            scatter(double(FiringTimeForThisUnit)/Fs/60,ValleyAcrossTime,5*ones(size(ValleyAcrossTime)),'+')
            
            
%             ax=axes('Units','centimeters','position', [38 0.6 5 4.5]);
%             set(ax,'FontSize',6)
%             
%             if ~isempty(counts_M)
%                 bar(10:20:1010,counts_M)
%                 xlim([0 1010])
%             end
            
            
            %     ax=axes('Units','centimeters','position', [12 0.6 5 4.5]);
            %     set(ax,'FontSize',6)
            % if ~isempty(counts_L)
            % bar(0.5:1:20.5,counts_L)
            % xlim([0 21])
            % end
            dataFile= result.date;
            ax=axes('Units','centimeters','position', [32 5.5 5 4.5]);
            set(ax,'FontSize',6)
%             T=text(.1,.1,[dataFile '-----' num2str(list(clu-clu_last))  ]); %
            T=text(.1,.2,[num2str(numel(result.time{clu-clu_last})) '-----']); %
            T=text(.1,.3,[num2str(numel(result.time{clu-clu_last})/(max(result.time{clu-clu_last})/Fs)) '  spk/sec']); %
            T=text(.1,.4,['2ms Violation % = ' sprintf('%0.1f%%',result.violation(clu-clu_last)*100) ]); %
            
            
            
            
            
            
            
            
            set(gcf, 'Position', [277         184        1892         952]);
            saveas(gcf,strcat([input_filename(1:end-4) '/'],num2str(clu)),'png')
            close(h)
% cdata = print('-RGBImage');
%    imwrite(cdata,[strcat('Waveform_Fig/',num2str(clu)) '.png'])
%                         print(gcf,[strcat('Waveform_Fig/',num2str(clu)) '.png'],'-dpng','-RGBImage')
% PlotAdd()
            %             F(clu)=getframe(h)
            
            sprintf('plot use time')
            toc
            
        end
        close all
    end
    
        clu_last=clu_last+numel(result.waveform);

    
end
% writerObj = VideoWriter('Waveform_summary','MPEG-4');
% writerObj.FrameRate=2;
% open(writerObj);
% pngFrames=dir('*.png');
% for K = 1 : numel(pngFrames)
%     K
%     filename = [num2str(numel(pngFrames)+1-K) '.png'];
%     thisimage = imread(filename);
%     writeVideo(writerObj, thisimage);
% end
% close(writerObj);
% cd ..
