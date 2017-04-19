function nproj=split_projection(proj,nleaf,nprojin,offsetval)
global subproj lensubproj;
%identify groups of adjacent leaves
nproj=nprojin;
inproj=false;
for leaf=1:nleaf
    if proj(leaf)>0
        if inproj
            lensubproj(nproj)=lensubproj(nproj)+1;
            subproj(nproj,lensubproj(nproj),1)=leaf+offsetval;
            subproj(nproj,lensubproj(nproj),2)=proj(leaf);           
        else
            nproj=nproj+1;
            subproj(nproj,1,1)=leaf+offsetval;
            subproj(nproj,1,2)=proj(leaf);
            lensubproj(nproj)=1;
            inproj=true;
        end
    else
        inproj=false;
    end
end
return 