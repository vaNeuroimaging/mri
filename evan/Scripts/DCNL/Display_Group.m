%Display_Group
%
%Allows fast loading of multiple subjects' images for display in SPM's
%checkreg function when given a common path pattern that the images share.
%
%Subject names and the image pattern are specified at the top of the
%script.
%
%The number of images to display in each window is also specified at the
%top of the script.  Note that more than 15 makes the images look really
%small. 6 is a pretty good number.
%
%The script will pause after displaying each set of images. Press any key
%in the main Matlab window to make the script continue.
%
%You must have SPM in your path for this to work.
%
%
%Created by E. Gordon 01/12



subject = 'Placeholder'; %DO NOT CHANGE


%USER INPUT
%--------------------------------------------------------------------------
%Subject names (as specified in the directory structure)
subjects = {'101','102','113','118','120','122','125','127','132','138','147'};%,'150','151','154','156','159','160','161','162','166','172','181','182','187','202','207','211','214','215','221','225','229','232','233','242','250','254','255','272','274','112','126','133','137','208','222','227','251','253','256','258','260','261','264','270','279','281','283','292','301','307','322','327','247','295','297','300','305','309','334','339','340','343','359','362','383','395','396','374','400','401','402','406','407','410','412','415','416','417','420'};

%Pattern of the image to find for each subject
basepattern = ['/fmri/data3/Evan/Gene-Rest-Nback/Data/' subject '/Struct/DC*.img'];

%Number of images to display in each graphics window. The script will pause
%after displaying each window.
%This number should not be more than 15 or the images get really small. I
%like using 6.
imagesperwindow = 6;

%END USER INPUT
%--------------------------------------------------------------------------

%Store the placeholder string
replacethisstring = subject;

%Calculate the number of windows that will be needed
numwindows = ceil(length(subjects)/imagesperwindow);

%Loop through the windows
for windownum = 1:numwindows;
    
    %Clear and reset various variables from the last loop
    clear images imagenames foundimageindex pattern thisimage
    foundimageindex = [];
    
    %Loop through the subjects to be shown in this window
    for subnum = 1:min(length(subjects) - imagesperwindow*(windownum-1) , imagesperwindow)
        
        %Get this subject's name
        subject = subjects{imagesperwindow*(windownum-1) + subnum};
        
        %Replace the placeholder string with this subject's name
        pattern = [basepattern(1:strfind(basepattern,replacethisstring)-1) subject basepattern(strfind(basepattern,replacethisstring)+length(replacethisstring):end)];
        
        %Get the directory the imnage is in
        directory = pattern(1:(find(pattern=='/',1,'last')));
        
        %Get the image name within the directory
        thisimage = dir(pattern);
        
        %If the image exists there
        if ~isempty(thisimage)
            
            %Store the image name
            imagenames(subnum,:) = sprintf('%-400s',[directory thisimage(1).name]);
            %And note that it was successfully found
            foundimageindex(subnum) = 1;
            
        else
            %Store a string reporting that the image wasn't found
            imagenames(subnum,:) = sprintf('%-400s',['Can''t find ' pattern]);
            %And note that fact
            foundimageindex(subnum) = 0;
        end
        
    end
    
    %Load the found images into SPM
    images = spm_vol(imagenames(find(foundimageindex),:));
    
    %Initialize the display window
    spm_figure('GetWin','Graphics');
    spm_figure('Clear','Graphics');
    spm_orthviews('Reset');
    
    %Display the images (this block stolen from SPM functions)
    mn = length(images);
    n  = round(mn^0.4);
    m  = ceil(mn/n);
    w  = 1/n;
    h  = 1/m;
    ds = (w+h)*0.02;
    for ij=1:mn
        i = 1-h*(floor((ij-1)/n)+1);
        j = w*rem(ij-1,n);
        handle = spm_orthviews('Image', images(ij),...
            [j+ds/2 i+ds/2 w-ds h-ds]);
        if ij==1, spm_orthviews('Space'); end
        spm_orthviews('AddContext',handle);
    end
    
    %List which images are being displayed
    disp('Showing Images:')
    disp(imagenames)
    
    %Wait for the user to press any button before moving on to the next
    %set of images to display
    if windownum < numwindows;
        disp('Press any key to see more images')
        pause
    end
    
    
end