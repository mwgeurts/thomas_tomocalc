function newnproj=segmentsubproj(projval,nproj)
%extract the next segment from a projection.
global segments subproj lensubproj;

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
 
return