function readTomodata(field_width)
global TPR SP OARXOPEN OARXLEAVES OARY; 

% Store the path of this function as a global. This is used by 
% readTomodata()
[start_path, ~, ~] = fileparts(mfilename('fullpath'));

%[fname pname] = uigetfile('*.txt');
SP.length=field_width; 
if field_width>=5.0
    fname='LA4data50.txt';
elseif field_width>=2.5
    fname='LA4data25.txt';
else
    fname='LA4data10.txt';
end
%pname='U:\Simon Thomas\comp_projects\TomoCheck\';
pname=start_path;
fid = fopen(fullfile(pname, fname)); % Modified by mwgeurts on 2017-04-19 to use fullfile to concatenate file name
while 1>0    
C=textscan(fid,'%[^[]'); 
C=textscan(fid,'%[^]]');
D=textscan(fid,'%c',1);
tagvalue=C{1};
if strcmp(tagvalue,'[TPRdata')
D=textscan(fid,'%d',1);
TPR.ndepths=D{1};
D=textscan(fid,'%f',TPR.ndepths);
TPR.depths=single(D{1});
D=textscan(fid,'%d',1);
TPR.nsizes=D{1};
D=textscan(fid,'%f',TPR.nsizes);
TPR.sizes=single(D{1});
form1='%f';
form2=' %f';
for i=2:TPR.nsizes
    form1=strcat(form1,form2);
end;
D=textscan(fid,form1,TPR.ndepths);
TPR.tpr=single(cell2mat(D));
end

if strcmp(tagvalue,'[SPdata')
            D=textscan(fid,'%d',1);
            SP.nsizes=D{1};
            D=textscan(fid,'%f',SP.nsizes);
            SP.sizes=single(D{1});
            D=textscan(fid,'%f',SP.nsizes);
            SP.values=single(D{1});
end 
if strcmp(tagvalue,'[OARXOPENdata')
    D=textscan(fid,'%d',1);
    OARXOPEN.ndepths=D{1};
    D=textscan(fid,'%f',OARXOPEN.ndepths);
    OARXOPEN.depths=single(D{1}); 
    D=textscan(fid,'%d',1);
    OARXOPEN.nvalues=D{1};
    OARXOPEN.indices=single(0:50);
    form1='%f';
    form2=' %f';
    for i=2:OARXOPEN.ndepths
        form1=strcat(form1,form2);
    end;
    D=textscan(fid,form1,OARXOPEN.nvalues);
    
    OARXOPEN.oar=single(cell2mat(D));
end  

if strcmp(tagvalue,'[OARXLEAVESdata')
    D=textscan(fid,'%d',1);    
    OARXLEAVES.nwidths=D{1};
    D=textscan(fid,'%d',OARXLEAVES.nwidths);
    OARXLEAVES.widths=single(D{1});
    D=textscan(fid,'%d',1);
    OARXLEAVES.ndepths=D{1};
    D=textscan(fid,'%f',OARXLEAVES.ndepths);
    OARXLEAVES.depths=single(D{1}); 
    D=textscan(fid,'%d',1);
    OARXLEAVES.nvalues=D{1};
    D=textscan(fid,'%f',OARXLEAVES.nvalues);
    OARXLEAVES.indices=single(D{1});
    OARXLEAVES.oar = single(zeros([OARXLEAVES.nvalues OARXLEAVES.ndepths OARXLEAVES.nwidths])); %initialise
    form1='%f';
    form2=' %f';
    for i=2:OARXLEAVES.ndepths
        form1=strcat(form1,form2);
    end;
    for w=1:OARXLEAVES.nwidths
       D=textscan(fid,form1,OARXLEAVES.nvalues); 
       OARXLEAVES.oar(:,:,w)=single(cell2mat(D)); 
    end
end  

if strcmp(tagvalue,'[OARYdata')
    D=textscan(fid,'%d',1);    
    OARY.nwidths=D{1};
    D=textscan(fid,'%f',OARY.nwidths);
    OARY.widths=single(D{1});
    
    D=textscan(fid,'%d',1);
    OARY.ndepths=D{1};
    D=textscan(fid,'%f',OARY.ndepths);
    OARY.depths=single(D{1}); 
    D=textscan(fid,'%d',1);
    OARY.nvalues=D{1};
    D=textscan(fid,'%f',OARY.nvalues);
    OARY.indices=single(D{1});
    OARY.oar = single(zeros([OARY.nvalues OARY.ndepths OARY.nwidths])); %initialise
    form1='%f';
    form2=' %f';
    for i=2:OARY.ndepths
        form1=strcat(form1,form2);
    end;
    for w=1:OARY.nwidths
        D=textscan(fid,form1,OARY.nvalues);
        OARY.oar(:,:,w)=single(cell2mat(D)); 
    end  
end

if isempty(tagvalue)
    fclose(fid);
    break
end

end