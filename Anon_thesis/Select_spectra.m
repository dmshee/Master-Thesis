% Function to allow user to select spectra from an image and plot them
% Aoife Gowen 2023

function [Spec]=Select_spectra(Img,wavelength)
colrs=['r';'g';'b';'c';'m';'k'];
if isempty(wavelength)
    wavelength=1:1:size(Img,3);
end
clc
nclass= input('How many classes?: ', 's');
nclass=str2double(nclass);

reply = input('Do you want to select pixel (1) or ROI (2) spectra?: ', 's');

if reply == '1'
    [Spec,arg2]=Select_Spec(Img,nclass,wavelength);
else
    [Spec,arg2]=Select_ROI(Img,nclass,wavelength);
end



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %select pixel spectra
    function[Spec,Coord]=Select_Spec(Img,nclass,wavelength)
        [x,y,z]=size(Img);
        figure,subplot(1,nclass+2,1)
        for i=1:nclass
            sprintf('select spectra from class %d, then hit Enter',i)

            Img_mean=(mean(Img,3)-prctile(nonzeros(mean(Img,3)),5))./(prctile(nonzeros(mean(Img,3)),95)-prctile(nonzeros(mean(Img,3)),5));

            [r,c,P] = impixel(Img_mean);
            Coord{i}=[r,c];

        end

        for i=1:nclass
          
            for j=1:size(Coord{1,i},1)
                Spec{i}(j,:)=reshape((Img(Coord{1,i}(j,2),Coord{1,i}(j,1),:)),1,z);
                hold on,plot(Coord{1,i}(j,1),Coord{1,i}(j,2),'+', 'Color',colrs(i),'MarkerSize', 10,'LineWidth',2)

            end
        end




        for i=1:nclass
            subplot(1,nclass+2,i+1),plot(wavelength,(Spec{i}),'Color',colrs(i)),title(sprintf('Class %d',i))';

             xlim([wavelength(1),wavelength(end)])
             ylim([0,mean(Img(:))+3*std(Img(:))])
           
        end

       subplot(1,nclass+2,nclass+2)
        for i=1:nclass
            hold on,plot(wavelength,(Spec{i}),'Color',colrs(i)),title('All classes together');
            xlim([wavelength(1),wavelength(end)])
             ylim([0,mean(Img(:))+3*std(Img(:))])
           
        end

    end



        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function[Spec,BW]=Select_ROI(Img,nclass,wavelength)
            [x1,y1,z1]=size(Img);
            figure,
            for i=1:nclass
                Img_mean=(mean(Img,3)-prctile(nonzeros(mean(Img,3)),5))./(prctile(nonzeros(mean(Img,3)),95)-prctile(nonzeros(mean(Img,3)),5));

                sprintf('select spectra from class %d',i)
                [BW(:,:,i)] = roipoly(Img_mean);
            end
            Ind=[];
            Img_unfold=reshape(Img,x1*y1,z1);

        BW_sum=sum(BW,3)>0;

figure,
subplot(1,nclass+2,1),imshow(double(cat(3,BW_sum,BW_sum,BW_sum))+cat(3,Img(:,:,ceil(z1/2)),Img(:,:,ceil(z1/2)),Img(:,:,ceil(z1/2))))
for i=1:nclass
    %find centre coordinates of the ROI
   [~, rBW]=max(mean(BW(:,:,i)));
    [~,cBW]=max(mean(BW(:,:,i)'));

text(rBW,cBW,sprintf('Class %d',i),'EdgeColor',colrs(i),'FontSize',14)
end


title('Selected regions')

            for i=1:size(BW,3)
                Spec{i}=Img_unfold(BW(:,:,i)>0,:);
            end

           
            for i=1:nclass
                Mean_spec=mean(Spec{i});
                Std_spec=std(Spec{i});
                subplot(1,nclass+2,i+1),plot(wavelength,Mean_spec,'Color',colrs(i)),title(sprintf('Class %d',i)),hold on;
                plot(wavelength,Mean_spec+Std_spec,'--','Color',colrs(i))
                plot(wavelength,Mean_spec-Std_spec,'--','Color',colrs(i))
                xlim([wavelength(1),wavelength(end)])
                ylim([0,mean(Img(:))+3*std(Img(:))])
            end
             subplot(1,nclass+2,nclass+2),hold on,title('All Classes together')
               for i=1:nclass
                Mean_spec=mean(Spec{i});
                Std_spec=std(Spec{i});
                plot(wavelength,Mean_spec,'Color',colrs(i)),
                plot(wavelength,Mean_spec+Std_spec,'--','Color',colrs(i))
                plot(wavelength,Mean_spec-Std_spec,'--','Color',colrs(i))
                xlim([wavelength(1),wavelength(end)])
                ylim([0,mean(Img(:))+3*std(Img(:))])
            end

        end
end