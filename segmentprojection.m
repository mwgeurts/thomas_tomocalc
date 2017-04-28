function [nproj, segments] = segmentprojection(projection)

% Initialize nested function shared variables
segments = struct;
subproj = [];
lensubproj = [];

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


function nproj=split_projection(proj,nleaf,nprojin,offsetval)
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
end 


function newnproj=segmentsubproj(projval,nproj)
%extract the next segment from a projection.

    a=subproj(projval,:,2);
    mv=min(a(1:lensubproj(projval))); %smallest value - i.e. value for largest segment
    iseg=segments.nseg+1;
    segments.startval(iseg)=subproj(projval,1,1);
    segments.endval(iseg)=subproj(projval,lensubproj(projval),1);
    segments.value(iseg)=mv;
    subproj(projval,:,2)=subproj(projval,:,2)-mv;
    segments.nseg=iseg;
    if max(subproj(projval,:,2))<0.0001 
        newnproj=nproj; % finished with this subproj
        return  %end of segmentation
    end
    %you will now have reduced or split the group.  Therefore re-run splitting
    nleaves=1+segments.endval(iseg)-segments.startval(iseg);
    newnproj=split_projection(subproj(projval,:,2),nleaves,nproj,subproj(projval,1,1)-1);
 
end



end