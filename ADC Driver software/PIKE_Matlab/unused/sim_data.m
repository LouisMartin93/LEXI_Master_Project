rows = 480;
columns =640;
numbit = 2^8;
Lumper = 0.7;
bg = 40;
c1=[220,150];
R=100;
c2 = [290,440];

data = zeros(rows,columns);

[x,y] = meshgrid(1:R*2,1:R*2);

rc1 = sqrt((x-R).^2+(y-R).^2);
rc2 = sqrt((x-R).^2+(y-R).^2);

pup = (rc1<R);

data(c1(1)-R:c1(1)+R-1,c1(2)-R:c1(2)+R-1)=pup;
data(c2(1)-R:c2(1)+R-1,c2(2)-R:c2(2)+R-1)=pup;


imshow(data,[]);
