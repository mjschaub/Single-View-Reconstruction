close all;
im = imread('west-quad-pool.jpg');
  
imshow(im);
%%
[im_y, im_x, c] = size(im);
[vp_x,vp_y,irx,iry,orx,ory] = TIP_GUI(im);
[bim,bim_alpha,vp_x,vp_y,ceilrx,ceilry,floorrx,floorry,...
    leftrx,leftry,rightrx,rightry,backrx,backry] = ...
    TIP_get5rects(im,vp_x,vp_y,irx,iry,orx,ory);

figure(2);
imshow(bim);
figure(2);
hold on;
plot(vp_x,vp_y,'w*');
plot([ceilrx ceilrx(1)], [ceilry ceilry(1)], 'y-');
plot([floorrx floorrx(1)], [floorry floorry(1)], 'm-');
plot([leftrx leftrx(1)], [leftry leftry(1)], 'c-');
plot([rightrx rightrx(1)], [rightry rightry(1)], 'g-');
hold off;

%% Calculate focal length

%%
vp_x = getVanishingPoint_shell(im);
%%
vp_y = getVanishingPoint_shell(im);
%%
vp_z = getVanishingPoint_shell(im);
%%



temp_1_x = vp_x(1)+vp_y(1);
temp_2_x = vp_y(1)+vp_z(1);
temp_3_x = vp_z(1)+vp_x(1);
temp_1_y = vp_x(2)+vp_y(2);
temp_2_y = vp_y(2)+vp_z(2);
temp_3_y = vp_z(2)+vp_x(2);
temp_xy = vp_x(1)*vp_y(1)+vp_x(2)*vp_y(2);
temp_yz = vp_z(1)*vp_y(1)+vp_z(2)*vp_y(2);
temp_xz = vp_x(1)*vp_z(1)+vp_x(2)*vp_z(2);

%solve for image center ax+by+c=0
syms u_0 v_0;
[final_u, final_v] = solve(-(temp_1_x)*u_0 - (temp_1_y)*v_0 + temp_xy == -(temp_2_x)*u_0 -(temp_2_y)*v_0 + temp_yz,...
-(temp_3_x)*u_0 - (temp_3_y)*v_0 + temp_xz == -(temp_2_x)*u_0 - (temp_2_y)*v_0 + temp_yz);

u = double(final_u);
v = double(final_v);
%focal length solving
eq_f = (u - vp_x(1))*(u - vp_y(1)) + (v - vp_x(2))*(v - vp_y(2));
f = sqrt(-eq_f);
f = double(f);


K = [f      0       u;
     0      f       v;
     0      0       1];   




%% calculate depth
focal_length = 100;
%focal_length = f;

ratio = (iry(4) - iry(1))/(ory(4) - ory(1));
depth_left = (focal_length - ratio*focal_length) / ratio; 

ratio = (iry(3) - iry(2))/(ory(3) - ory(2));
depth_right = (focal_length - ratio*focal_length) / ratio; 

ratio = (irx(2) - irx(1))/(orx(2) - orx(1));
depth_up = (focal_length - ratio*focal_length) / ratio; 

ratio = (irx(3) - irx(4))/(orx(3) - orx(4));
depth_down = (focal_length - ratio*focal_length) / ratio;

%% Crop the texture faces

h = size(im,1);
w = size(im,2);

img_middle = zeros(h, w, 3);
img_left = zeros(h, w, 3);
img_right = zeros(h, w, 3);
img_up = zeros(h, w, 3);
img_down = zeros(h, w, 3);

% middle mask
mask = poly2mask([irx(1) irx(2) irx(3) irx(4)], [iry(1) iry(2) iry(3) iry(4)], h, w);
for i=1:h
    for j=1:w
        if(mask(i,j)==1)
            img_middle(i,j,1) = double(im(i,j,1))/255;
            img_middle(i,j,2) = double(im(i,j,2))/255;
            img_middle(i,j,3) = double(im(i,j,3))/255;
        end
    end
end
%figure(2), imshow(img_middle);
%left mask
mask = poly2mask([orx(1) irx(1) irx(4) orx(4)], [ory(1) iry(1) iry(4) ory(4)], h, w);
for i=1:h
    for j=1:w
        if(mask(i,j)==1)
            img_left(i,j,1) = double(im(i,j,1))/255;
            img_left(i,j,2) = double(im(i,j,2))/255;
            img_left(i,j,3) = double(im(i,j,3))/255;
        end
    end
end
%figure(3), imshow(img_left);

%right mask
mask = poly2mask([orx(2) irx(3) irx(3) orx(2)], [ory(2) iry(2) iry(3) ory(4)], h, w);
for i=1:h
    for j=1:w
        if(mask(i,j)==1)
            img_right(i,j,1) = double(im(i,j,1))/255;
            img_right(i,j,2) = double(im(i,j,2))/255;
            img_right(i,j,3) = double(im(i,j,3))/255;
        end
    end
end
%figure(4), imshow(img_right);

%ceiling mask
mask = poly2mask([orx(1) orx(2) irx(2) irx(1)], [ory(1) ory(2) iry(2) iry(1)], h, w);
for i=1:h
    for j=1:w
        if(mask(i,j)==1)
            img_up(i,j,1) = double(im(i,j,1))/255;
            img_up(i,j,2) = double(im(i,j,2))/255;
            img_up(i,j,3) = double(im(i,j,3))/255;
        end
    end
end

%figure(5), imshow(img_up);


%floor mask
mask = poly2mask([orx(3) orx(4) irx(4) irx(3)], [ory(3) ory(4) iry(4) iry(3)], h, w);
for i=1:h
    for j=1:w
        if(mask(i,j)==1)
            img_down(i,j,1) = double(im(i,j,1))/255;
            img_down(i,j,2) = double(im(i,j,2))/255;
            img_down(i,j,3) = double(im(i,j,3))/255;
        end
    end
end
%figure(6), imshow(img_down);
 
%% calculate homographies

points_one = [orx(1) orx(2) irx(2) irx(1) ; ory(1) ory(2) iry(2) iry(1) ; 1 1 1 1];
points_two = [irx(1) irx(2) irx(2) irx(1) ; 1 1 depth_up depth_up ; 1 1 1 1];
H_Up = calc_homography(points_one, points_two);
T = maketform('projective', H_Up.');
img_up_transform = imtransform(img_up, T, 'XData',[irx(1) irx(2)],'YData',[0 depth_up]);
%figure(3), imshow(img_up_transform);

points_one = [orx(3) orx(4) irx(4) irx(3) ; ory(3) ory(4) iry(4) iry(3) ; 1 1 1 1];
points_two = [irx(3) irx(4) irx(4) irx(3) ; depth_down depth_down 1 1 ; 1 1 1 1];
H_down = calc_homography(points_one, points_two);
T = maketform('projective', H_down.');
img_down_transform = imtransform(img_down, T, 'XData',[irx(4) irx(3)],'YData',[depth_down 0]);
%figure(4), imshow(img_down_transform);

points_one = [orx(1) irx(1) irx(4) orx(1) ; ory(1) iry(1) iry(4) ory(4) ; 1 1 1 1];
points_two = [1 depth_left depth_left 1 ; iry(1) iry(1) iry(4) iry(4) ; 1 1 1 1];
H_Left = calc_homography(points_one, points_two);
T = maketform('projective', H_Left.');
img_left_transform = imtransform(img_left, T, 'XData',[0 depth_left],'YData',[iry(1) iry(4)]);
%figure(5), imshow(img_left_transform);

points_one = [orx(3) irx(3) irx(2) orx(2) ; ory(3) iry(3) iry(2) ory(2) ; 1 1 1 1];
points_two = [depth_right 1 1 depth_right ; iry(3) iry(3) iry(2) iry(2) ; 1 1 1 1];
H_Right = calc_homography(points_one, points_two);
T = maketform('projective', H_Right.');
img_right_transform = imtransform(img_right, T, 'XData',[0 depth_right],'YData',[iry(2) iry(3)]);
%figure(6), imshow(img_right_transform);

img_middle_transform = imcrop(img_middle, [irx(1) iry(1) (irx(2)-irx(1)) (iry(4)-iry(1))]);
%figure(7), imshow(img_middle_transform);

%% warp each plane
x_max = irx(2)-irx(1);
y_max = iry(4)-iry(1);

up_x = [depth_up depth_up depth_up; 0 0 0];
up_y = [x_max x_max/2 0; x_max x_max/2 0];
up_z = [y_max y_max y_max; y_max y_max y_max];

bottom_x = [depth_down depth_down depth_down; 0 0 0];
bottom_y = [x_max x_max/2 0; x_max x_max/2 0];
bottom_z = [0 0 0; 0 0 0];

middle_x = [0 0 0; 0 0 0];
middle_y = [x_max x_max/2 0; x_max x_max/2 0];
middle_z = [y_max y_max y_max; 0 0 0];

left_x = [0 depth_left/2 depth_left; 0 depth_left/2 depth_left];
left_y = [0 0 0; 0 0 0];
left_z = [y_max y_max y_max; 0 0 0];

right_x = [depth_right depth_right/2 0; depth_right depth_right/2 0];
right_y = [x_max x_max x_max; x_max x_max x_max];
right_z = [y_max y_max y_max; 0 0 0];


view = figure('name','3DViewer: Movement[WASD] Zoom in[Q] Zoom out[E] Exit[ESC]');
set(view,'windowkeypressfcn','set(gcbf,''Userdata'',get(gcbf,''CurrentCharacter''))') ;
set(view,'windowkeyreleasefcn','set(gcbf,''Userdata'','''')') ;
set(view,'Color','black')
hold on

warp(up_x, up_y, up_z, img_up_transform);
warp(bottom_x, bottom_y, bottom_z, img_down_transform);
warp(right_x, right_y, right_z, img_left_transform);
warp(left_x, left_y, left_z,img_right_transform);
warp(middle_x, middle_y, middle_z, img_middle_transform);

%% 3D model

axis equal; 
axis vis3d; 
axis off;
camproj('perspective');

step = 0.1;
camup([0,0,1]);

key = 0;
while (~key),
    waitforbuttonpress;
    key = get(view, 'currentch');
    
    switch key
        case 'd'
            camdolly(-step,0,0,'fixtarget');
        case 'a'
            camdolly(step,0,0,'fixtarget');
        case 's'
            camdolly(0,step,0,'fixtarget');
        case 'w'
            camdolly(0,-step,0,'fixtarget');
        case 'q'
            camdolly(0,0,step,'fixtarget');
        case 'e'
            camdolly(0,0,-step,'fixtarget');
        case 'p'
            break;
    end
    
    key = 0;
    pause(.001);

    
end;


