%% 
Sh(1)=regions.rect('name','original');
Sh(2)=Sh(1)+[1.01,1.01];
Sh(2).Name='plus';

Sh(3)=Sh(1)*[0.5,2]+[-1,-1];
Sh(3).Name='times+plus';

Sh(4)=Sh(1).rotate(45)+[-1,1];
Sh(3).Name='rotated+plus';
figure; 
Sh.draw();
axis equal
legend('Location','southeast')
%% 
clear Sh;
figure;
R=regions.rect('name','R');
C=regions.circle([0.5, 0.5],0.5);
Sh(1)=R+[-1,1];

Sh(2)=(R+C)+[1,1];
Sh(2).Name='R+C';
Sh(3)=(R-C)+[-1,-1];
Sh(3).Name='R-C';
Sh(4)=(R&C)+[1,-1];
Sh(4).Name='R&C';
Sh.draw();
legend('location','NorthEastOutside');
%% 
figure;
R=regions.rect('name','Emmental');
for r=0:4
    for c=0:4
        R=R-regions.circle([r/4-.5,c/4-.5],[0.1, 0.1]);
    end
end
R.draw();
