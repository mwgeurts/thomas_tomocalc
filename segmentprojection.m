function nproj=segmentprojection()
%read projection from text file
global projection segments;
%fname='projection.txt';
%pname='U:\Simon Thomas\comp_projects\TomoCheck\';
%fid = fopen([pname fname]);
%D=textscan(fid,'%f',64);
%projection=D{1};
nproj=split_projection(projection,64,0,0);
if (nproj<1)
    return
end
projval=1;
segments.nseg=0;
while (projval<=nproj)
    nproj=segmentsubproj(projval,nproj); % keeps increasing until done
    projval=projval+1;
end