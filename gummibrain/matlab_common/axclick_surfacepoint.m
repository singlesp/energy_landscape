function axclick_surfacepoint(gcbo,eventdata,handles)

    global calc_idx_tmp;
    
    curpoint = get(gca,'CurrentPoint');
    q1 = curpoint(1,:);
    q2 = curpoint(2,:);
    
    d = q2-q1;
    d = d/norm(d);
    
    verts = get(gcbo,'Vertices');
    norms = get(gcbo,'VertexNormals');
    faces = get(gcbo,'Faces');
    
    [intersect,indx,dist,u,v] = ray_intersect(faces,verts,q1,d,'b');
    
    %find the closest face to the starting point
    [dmin dindx] = min(dist);
    closest_face_idx = indx(dindx(1));
    closest_face_verts = verts(faces(closest_face_idx,:),:);
    
    %find the closest of the face's 3 verts to the line
    d = DistanceFromPointToLine(closest_face_verts, q1,q2);
    [dmin dindx] = min(d);

    vert_idx = faces(closest_face_idx,dindx);
    closest_vert = verts(vert_idx,:);
    
    calc_idx_tmp(end+1) = vert_idx;
    
    nextplot = get(gca,'nextplot');
    set(gca,'nextplot','add');
    %plot3(curpoint(:,1),curpoint(:,2),curpoint(:,3),'kx-');
    %plot3(closest_face_verts(:,1),closest_face_verts(:,2),closest_face_verts(:,3),'-ok','MarkerSize',5,'MarkerFaceColor','r');
    plot3(closest_vert(:,1),closest_vert(:,2),closest_vert(:,3),'-ok','MarkerSize',10,'MarkerFaceColor','w');
    set(gca,'nextplot',nextplot);
    
    fprintf('%d (%.2f,%.2f,%.2f)\n',vert_idx, closest_vert);
    