%=========================================================================%
%===========      Program for generating composite network      ==========%
%===========         (actin filament and microtubule)           ==========%
%===========             Jane 20, 2019.  Bo Gong                ==========%
%=========================================================================%

%================================2023.10.31===============================%
%--------------------------------修改建模顺序------------------------------%
%-----------------------微管纤维生成顺序互换，变量名称互换-------------------%
%=========================================================================%
clc
clear

%== Generate the network by placing filament and microtuble one by one. ==%

% define parameters of the simulation box, filaments and microtubules
Wx=1.0e3;  % the box size
Wy=1.0e3;
Wz=1.0e3; 

d_fil=4; % diameter of the filament bead
d_mt=6; % diameter of the MT bead
len_fil=6.16e2;   % length of the filament 600
len_MT=6.00e2;  % length of the microtubule 606  用mt长度换算
beads_fil=len_fil/d_fil;  % number of beads for each filament  100个珠子
beads_MT=len_MT/d_mt; % 101个珠子

del_fil=30; % the minimum distance between beads for two different filaments   来自不同丝的两个珠子的最小距离
del_filMT=25;  % 13 the minimum distance between beads for filaments and microtubules  来自微管和细丝的两个珠子的最小距离
del_MT=6;   %  23 the minimum distance between beads for two different microtubules   来自不同微管的两个珠子的最小距离

N_fil=93;  % number of filaments  细丝的数量 500
N_MT=525;  % number of microtubules  微管的数量 25
N_max=1e8;  % the maximum loops for generate filaments  创建细丝的最大循环数




%============= generate filaments and microtubue one by one ==============%

%--------------------- generate filaments one by one----------------------%
Position_atoms=zeros(1,6);
MT_vector=zeros(1,7);         % record the direction vectors and intial endpoint of the microtubule
atomtype_fil=2;  % beads type of filaments  1-细丝
atomtype_MT=1;   % beads type of MT         2-mt
moletype=0;
n_atoms=0;
num_fil=0;  % count the filament number 长丝的数量

for i=1:N_max
    point_1x=rand(1,1)*Wx; % 取两个点
    point_1y=rand(1,1)*Wy;
    point_1z=rand(1,1)*Wz;

    point_2x=rand(1,1)*Wx;
    point_2y=rand(1,1)*Wy;
    point_2z=rand(1,1)*Wz;

    direct_vectx=point_2x-point_1x;  % x方向的方向向量
    direct_vecty=point_2y-point_1y;
    direct_vectz=point_2z-point_1z;

    % length of the direction vector
    len_vect=sqrt(direct_vectx^2+direct_vecty^2+direct_vectz^2); % 计算平方根 两个点构成向量的长度

    % unit direction vector
    vect_ix=direct_vectx/len_vect;      % x方向上的单位向量
    vect_iy=direct_vecty/len_vect;
    vect_iz=direct_vectz/len_vect;


    % generate the first endpoint of each filament
    endpoint_1x=rand(1,1)*Wx;   % 随机一个端点的x坐标 
    endpoint_1y=rand(1,1)*Wy;
    endpoint_1z=rand(1,1)*Wz;

    % generate the coordinates of the filament  1根坐标
    Temp_atoms=zeros(1,6);
    for j=1:beads_fil                       % 100个细丝珠子
        temp_x=endpoint_1x+j*vect_ix*6;     % 第j个珠子的x坐标
        temp_y=endpoint_1y+j*vect_iy*6;
        temp_z=endpoint_1z+j*vect_iz*6;
        if temp_x<=0                        % 珠子坐标的周期性边界条件
           temp_x=temp_x+Wx;
        elseif temp_x>=Wx
               temp_x=temp_x-Wx;
        end
        if temp_y<=0
            temp_y=temp_y+Wy;
        elseif temp_y>=Wy
            temp_y=temp_y-Wy;
        end
        if temp_z<=0
            temp_z=temp_z+Wz;
        elseif temp_z>=Wz
            temp_z=temp_z-Wz;
        end
        Temp_atoms(j,4)=temp_x;             % 4、5、6列分别为x，y，z坐标
        Temp_atoms(j,5)=temp_y;
        Temp_atoms(j,6)=temp_z;
    end
        
    % delete the filament if the atoms comes from the two filament are      ？
    % closed enough 
    Temp_natoms=length(Position_atoms(:,1)); % 求数组的维数行数和列数的最大值 1
    delete=zeros(1,3);
    n_delete=0;
    if Temp_natoms>=beads_fil 
        for ii=1:Temp_natoms
            for jj=1:beads_fil
                del_x=Temp_atoms(jj,4)-Position_atoms(ii,4);
                del_y=Temp_atoms(jj,5)-Position_atoms(ii,5);
                del_z=Temp_atoms(jj,6)-Position_atoms(ii,6);
                del_filfil=sqrt(del_x^2+del_y^2+del_z^2);
                if del_filfil<=del_fil
                    n_delete=n_delete+1;
                    delete(n_delete,:)=[ii,jj,1];
                end
            end
        end
    end
    
    temp_delete=sum(delete(:,3));
    if temp_delete==0
        num_fil=num_fil+1;
        MT_vector(num_fil,:)=[num_fil, vect_ix, vect_iy, vect_iz, Temp_atoms(1,4), Temp_atoms(1,5), Temp_atoms(1,6)];
        Temp_atoms(:,1)=beads_fil*(num_fil-1)+1:beads_fil*num_fil;
        Temp_atoms(:,2)=num_fil;
        Temp_atoms(:,3)=atomtype_fil;
        if Temp_natoms==1
            Position_atoms(1:beads_fil,:)=Temp_atoms;
        else
            Position_atoms(Temp_natoms+1:Temp_natoms+beads_fil,:)=Temp_atoms; 
        end
    end

    if num_fil>=N_fil
        break
    end
end

% for k=1:length(Position_atoms(:,1))
%     plot3(Position_atoms(k,4), Position_atoms(k,5), Position_atoms(k,6),"o",'MarkerSize',3);
%     hold on
% end

save MT_direction_vector.txt MT_vector -ascii
save MT_Position_beads.txt Position_atoms -ascii

save ver_all.mat


%---------------------------------delete bonds and angles of boundary---------------------------------%

% pos_x=Position_atoms(:,4);    % fil原子的xyz坐标
% pos_y=Position_atoms(:,5);
% pos_z=Position_atoms(:,6);
% 
% temp_boundx_atoms=find(pos_x>=(Wx-5.2)&pos_x<=Wx);
% temp_boundy_atoms=find(pos_y>=(Wy-5.2)&pos_y<=Wy);
% temp_boundz_atoms=find(pos_z>=(Wz-5.2)&pos_z<=Wz);

%---------------------------------矩阵索引超出删除范围-------------------------------------------------%

% for i=1:N_MAX
%     temp_boundx_atoms=find(pos_x>=(Wx-5.2)&pos_x<=Wx);
%     if temp_boundx_atoms(end,1)<=size(bonds_fil,1)
%        temp_boundx_atoms=temp_boundxx_atoms;
%     end
%     i=i+1;
% end
% 
% for i=1:N_MAX
%     temp_boundx_atoms=find(pos_x>=(Wx-5.2)&pos_x<=Wx);
%     if temp_boundx_atoms(end,1)<=size(bonds_fil,1)
%        temp_boundx_atoms=temp_boundxx_atoms;
%     end
%     i=i+1;
% end
% 
% for i=1:N_MAX
%     temp_boundx_atoms=find(pos_x>=(Wx-5.2)&pos_x<=Wx);
%     if temp_boundx_atoms(end,1)<=size(bonds_fil,1)
%        temp_boundx_atoms=temp_boundxx_atoms;
%     end
%     i=i+1;
% end

%-------------------------------------------------------------------------%

% bonds_fil(temp_boundx_atoms,:)=[];
% bonds_fil(temp_boundy_atoms,:)=[];
% bonds_fil(temp_boundz_atoms,:)=[];




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%----------- generate the microtubule position one by one ---------------%
Temp_MTbeads=zeros(1,6);     
Temp_nMTbead=0;  
Position_beads=zeros(1,6);   % the matrix for microtubule beads space
num_MT=0;
for i=1:N_max
    point_1x=rand(1,1)*Wx;
    point_1y=rand(1,1)*Wy;
    point_1z=rand(1,1)*Wz;

    point_2x=rand(1,1)*Wx;
    point_2y=rand(1,1)*Wy;
    point_2z=rand(1,1)*Wz;

    direct_vectx=point_2x-point_1x;
    direct_vecty=point_2y-point_1y;
    direct_vectz=point_2z-point_1z;

    % length of the direction vector
    len_vect=sqrt(direct_vectx^2+direct_vecty^2+direct_vectz^2);

    % unit direction vector
    vect_ix=direct_vectx/len_vect;
    vect_iy=direct_vecty/len_vect;
    vect_iz=direct_vectz/len_vect;


    % generate the first endpoint of each filament
    endpoint_1x=rand(1,1)*Wx;
    endpoint_1y=rand(1,1)*Wy;
    endpoint_1z=rand(1,1)*Wz;

    % generate the coordinates of the filament
    for j=1:beads_MT
        temp_x=endpoint_1x+j*vect_ix*6;
        temp_y=endpoint_1y+j*vect_iy*6;
        temp_z=endpoint_1z+j*vect_iz*6;
        if temp_x<=0
           temp_x=temp_x+Wx;
        elseif temp_x>=Wx
               temp_x=temp_x-Wx;
        end
        if temp_y<=0
            temp_y=temp_y+Wy;
        elseif temp_y>=Wy
            temp_y=temp_y-Wy;
        end
        if temp_z<=0
            temp_z=temp_z+Wz;
        elseif temp_z>=Wz
            temp_z=temp_z-Wz;
        end
        Temp_MTbeads(j,4)=temp_x;
        Temp_MTbeads(j,5)=temp_y;
        Temp_MTbeads(j,6)=temp_z; 
    end 
    
    % delete the MT if the atoms comes from the two filament are
    % closed enough 
    Temp_nMTbeads=length(Position_beads(:,1));
    Natoms_fil=length(Position_atoms(:,1));  % The number of atoms in all filaments  5000
    delete_MT=zeros(1,3);
    n_deleteMT=0;
    for ii=1:Natoms_fil
        for jj=1:beads_MT       % 101  
            del_xMTfil=Temp_MTbeads(jj,4)-Position_atoms(ii,4);
            del_yMTfil=Temp_MTbeads(jj,5)-Position_atoms(ii,5);
            del_zMTfil=Temp_MTbeads(jj,6)-Position_atoms(ii,6);
            del_mtfil=sqrt(del_xMTfil^2+del_yMTfil^2+del_zMTfil^2);
            if del_mtfil<=del_filMT
                n_deleteMT=n_deleteMT+1;
                delete_MT(n_deleteMT,:)=[ii,jj,1];
            end                
        end
    end
    
    if Temp_nMTbeads>=beads_MT
        for ii=1:Temp_nMTbeads
            for jj=1:beads_MT
                del_xMTMT=Temp_MTbeads(jj,4)-Position_beads(ii,4);
                del_yMTMT=Temp_MTbeads(jj,5)-Position_beads(ii,5);
                del_zMTMT=Temp_MTbeads(jj,6)-Position_beads(ii,6);
                del_MTMT=sqrt(del_xMTMT^2+del_yMTMT^2+del_zMTMT^2);
                if del_MTMT<=del_MT
                    n_deleteMT=n_deleteMT+1;
                    delete_MT(n_deleteMT,:)=[ii,jj,1];
                end
            end
        end
    end

    output=[n_deleteMT,num_MT];
    disp('已经删除纤维数量：')
    disp(output)
    
    temp_deleteMT=sum(delete_MT(:,3));
    if temp_deleteMT==0
        num_MT=num_MT+1;
        Temp_MTbeads(:,1)=beads_MT*(num_MT-1)+1:beads_MT*num_MT;
        Temp_MTbeads(:,2)=num_MT;
        Temp_MTbeads(:,3)=atomtype_MT;
        if Temp_nMTbeads==1
            Position_beads(1:beads_MT,:)=Temp_MTbeads;
        else
            Position_beads(Temp_nMTbeads+1:Temp_nMTbeads+beads_MT,:)=Temp_MTbeads; 
        end
    end

    if num_MT>=N_MT
        break
    end
end

%---------------- generate the bonds and angles of the filaments---------------------%
Temp_atoms_fil=zeros(1,6);
bonds_fil=zeros(1,4);    % record the bonds information of the filaments
angles_fil=zeros(1,5);
bonds_type=1;
angles_type=1;
n_bonds=0;
n_angles=0;
for ii=1:N_MT  % 500根细丝
    Temp_index=find(Position_beads(:,2)==ii);       % 第ii根细丝
    Temp_atoms_fil=Position_beads(Temp_index,:);    % 第ii根细丝的坐标信息
    % generate bonds 
    for jj=1:beads_MT-1        % 100个珠子 99个键
        atom_bonds_01=Temp_atoms_fil(jj,1);
        atom_bonds_02=Temp_atoms_fil(jj+1,1);
        n_bonds=n_bonds+1;
        bonds_fil(n_bonds,:)=[n_bonds, bonds_type, atom_bonds_01, atom_bonds_02];
    end
    % generate angles
    for jj=1:beads_MT-2
        atom_angles_01=Temp_atoms_fil(jj,1);
        atom_angles_02=Temp_atoms_fil(jj+1,1);
        atom_angles_03=Temp_atoms_fil(jj+2,1);
        n_angles=n_angles+1;
        angles_fil(n_angles,:)=[n_angles, angles_type, atom_angles_01, atom_angles_02, atom_angles_03];
    end   
end

save Atoms_fil.txt Position_beads -ascii
save Bonds_fil.txt bonds_fil -ascii
save Angles_fil.txt angles_fil -ascii



