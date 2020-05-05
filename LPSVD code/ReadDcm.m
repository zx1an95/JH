clc

FolderName = uigetdir;
FileList = dir('1.2.156.112605.75006881735495.20190207153015.3.26528.4');
uCSRImgSet=[];
InstNumSet=[];
Counter=0;

%Images are stored in uCSRImgSet
%For coding, each one of them will do the work
%For image effects testing, all of them should be tested
for i=1:length(FileList)
    CurName = FileList(i).name;
    if(length(strfind(CurName,'.dcm'))>0)
        Counter=Counter+1;
        Img=dicomread([FolderName '\' CurName]);
        Inf=dicominfo([FolderName '\' CurName]);        
        uCSRImgSet(:,:,Counter)=Img;
        InstNumSet(Counter)=Inf.InstanceNumber;      
    end
end
