%=========================================================================%
%==========              Generate composite network             ==========%
%==========   generate filaments atoms and microtubules atoms   ==========%
%==========               June 21th, 2019, Bo Gong              ==========%
%=========================================================================%

% the program using the data generate from the program "Filaments_atoms_and MT_vector_position.m"
% to generate the whole composite networks atoms, bonds, angles and
% dihedrals, but not including the crosslinking between filaments or
% filaments and microtubules or microtubules and microtubules

clc
clear 

% load the data generate by program "Filaments_atoms_and MT_vector_position.m"
atoms_fil=load('Atoms_fil.txt'); %（1，6） 原子序号 细丝序号 原子类型 xyz
bonds_fil=load('Bonds_fil.txt'); % 键序号 类型 原子1 原子2
angles_fil=load('Angles_fil.txt'); % 角序号 类型 原子1 原子2 原子3 

MT_Position_vector=load('MT_direction_vector.txt');  % positon vector of microtubules  微管的位置向量
% MT数量 xyz方向单位向量  MT珠子的xyz坐标

Wx=1.0e3;  % the box size
Wy=1.0e3;
Wz=1.0e3; 
sigma=4.0;  % the units length for dimensionless    无量纲单位

% load the MT atoms, bond, angles and dihedrals (units nano)
One_MT_atoms=load('Atoms_oneMT_nano.txt');
One_MT_bonds=load('Bonds_oneMT.txt');
One_MT_angles=load('Angles_oneMT.txt');
One_MT_dihedrals=load('Dihedrals_oneMT.txt');

% define variables
MT_X=One_MT_atoms(:,4);     % MT x坐标
MT_Y=One_MT_atoms(:,5);     % MT y坐标
MT_Z=One_MT_atoms(:,6);     % MT z坐标
MT_X_trans=zeros(1,1);      %创建1*1的0矩阵
MT_Y_trans=zeros(1,1);
MT_Z_trans=zeros(1,1);

%---------------  parameters of one microtubule  -------------------------% 
Num_MT=length(MT_Position_vector(:,1));        % 微管位置向量维度长度 微管数量
Nfil_fil=atoms_fil(end,2);                     % 最后一行，第二列 蛋白丝数量
Natoms_fil=atoms_fil(end,1);                   % 蛋白丝原子数量
Nbonds_fil=bonds_fil(end,1);                   % 键的数量 100014 设置处理好边界交联的纤维网络的键角数量
Nangles_fil=angles_fil(end,1);                 % 角的数量 97093
Ndihedrals_fil=0;

Natomtype_fil=max(atoms_fil(:,3));             % 返回第三列的最大值 蛋白丝原子类型数量
Nbondtype_fil=max(bonds_fil(:,2));             % 键类型的数量
Nangletype_fil=max(angles_fil(:,2));           % 角类型的数量
Ndihedraltype_fil=0;

Natoms_oneMT=length(One_MT_atoms(:,1));         % MT 原子数量
Nfil_oneMT=One_MT_atoms(end,2);                 % MT 中细丝的数量
Nbonds_oneMT=length(One_MT_bonds(:,1));         % MT 中键的数量
Nangles_oneMT=length(One_MT_angles(:,1));       % MT 中角的数量
Ndihedrals_oneMT=length(One_MT_dihedrals(:,1)); % MT 中二面角的数量

MT_atoms=[];
MT_bonds=[];
MT_angles=[];
MT_dihedrals=[];

for n=1:Num_MT % 微管数量
    Temp_atoms_one_MT=zeros(Natoms_oneMT,6); % 微管原子 预留矩阵
    Temp_bonds_one_MT=zeros(Nbonds_oneMT,4);
    Temp_angles_one_MT=zeros(Nangles_oneMT,5);
    Temp_dihedrals_one_MT=zeros(Ndihedrals_oneMT,6);
    
    vect_x=MT_Position_vector(n,2);  % x单位向量
    vect_y=MT_Position_vector(n,3);
    vect_z=MT_Position_vector(n,4);
    x_1=MT_Position_vector(n,5);     % 珠子x坐标
    y_1=MT_Position_vector(n,6);
    z_1=MT_Position_vector(n,7);

    %------------- Compute the Euler angles of the vector -------------------%
    % the yaw angle  (pay attention to the direction of the angle) (Z=0)
    % the projection of the vector in xy plane   在xy面上的投影
    vect_xy_x=vect_x; % 单位向量x 在xy面上的投影
    vect_xy_y=vect_y;
    vect_xy_z=0;

    vect_ix=1;   % the units vector in x axis
    vect_iy=0;
    vect_iz=0;

    if vect_xy_y>=0
        if sqrt(vect_xy_x^2+vect_xy_y^2)==0 % 垂直于xy平面的
            cos_yangle=1; % y方向角度为0度
        else
            cos_yangle=(vect_xy_x*vect_ix+vect_xy_y*vect_iy)/sqrt(vect_xy_x^2+vect_xy_y^2); % xy面上两个向量夹角
        end
        yaw_angle=acos(cos_yangle); % 与y轴的角度
    elseif vect_xy_y<0
        if sqrt(vect_xy_x^2+vect_xy_y^2)==0
            cos_yangle=1;
        else
            cos_yangle=(vect_xy_x*vect_ix+vect_xy_y*vect_iy)/sqrt(vect_xy_x^2+vect_xy_y^2);  %(theta's direction should be minus)
        end
        yaw_angle=-acos(cos_yangle);
    end

    % yaw_angle=0;

    % the pitch angle 
    % angle between the vector and the projection of the vector in xy plane

    if vect_z>=0
        if (sqrt(vect_xy_x^2+vect_xy_y^2+vect_xy_z^2)*sqrt(vect_x^2+vect_y^2+vect_z^2))==0
            cos_pangle=0;
        else
            cos_pangle=(vect_xy_x*vect_x+vect_xy_y*vect_y+vect_xy_z*vect_z)/(sqrt(vect_xy_x^2+vect_xy_y^2+vect_xy_z^2)*sqrt(vect_x^2+vect_y^2+vect_z^2));
        end
        pitch_angle=-acos(cos_pangle);
    elseif vect_z<0
        if (sqrt(vect_xy_x^2+vect_xy_y^2+vect_xy_z^2)*sqrt(vect_x^2+vect_y^2+vect_z^2))==0
            cos_pangle=0;
        else
            cos_pangle=(vect_xy_x*vect_x+vect_xy_y*vect_y+vect_xy_z*vect_z)/(sqrt(vect_xy_x^2+vect_xy_y^2+vect_xy_z^2)*sqrt(vect_x^2+vect_y^2+vect_z^2));
        end
        pitch_angle=acos(cos_pangle);  
    end

    % pitch_angle=0;

    % the roll angle
    % random generate the roll angle
    roll_angle=rand(1,1)*2*pi;

    %------------------- transformation coordinates ---------------%

    N_trans=length(MT_X); % 1300个原子转换坐标
    for ii=1:N_trans
        x_ii=MT_X(ii,1); % 第ii个原子的x坐标
        y_ii=MT_Y(ii,1);
        z_ii=MT_Z(ii,1);

        x_ii_trans=cos(yaw_angle)*cos(pitch_angle)*x_ii+(cos(yaw_angle)*sin(pitch_angle)*sin(roll_angle)-sin(yaw_angle)*cos(roll_angle))*y_ii+(cos(yaw_angle)*sin(pitch_angle)*cos(roll_angle)+sin(yaw_angle)*sin(roll_angle))*z_ii+x_1;
        y_ii_trans=sin(yaw_angle)*cos(pitch_angle)*x_ii+(sin(yaw_angle)*sin(pitch_angle)*sin(roll_angle)+cos(yaw_angle)*cos(roll_angle))*y_ii+(sin(yaw_angle)*sin(pitch_angle)*cos(roll_angle)-cos(yaw_angle)*sin(roll_angle))*z_ii+y_1;
        z_ii_trans=-sin(pitch_angle)*x_ii+cos(pitch_angle)*sin(roll_angle)*y_ii+cos(pitch_angle)*cos(roll_angle)*z_ii+z_1;

        MT_X_trans(ii,1)=x_ii_trans;
        MT_Y_trans(ii,1)=y_ii_trans;
        MT_Z_trans(ii,1)=z_ii_trans;
    %     plot3(x_ii_trans, y_ii_trans, z_ii_trans,"o");
    %     hold on
    end
    
    %------ periodic boundarty condition for transformed coordinates -----% 周期性边界条件
    for jj=1:N_trans
        temp_x=MT_X_trans(jj,1);    % 第jj个原子的转化后坐标
        temp_y=MT_Y_trans(jj,1);
        temp_z=MT_Z_trans(jj,1);
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
        Temp_atoms_one_MT(jj,4)=temp_x;
        Temp_atoms_one_MT(jj,5)=temp_y;
        Temp_atoms_one_MT(jj,6)=temp_z;
    end
    
    Temp_atoms_one_MT(:,1)=Natoms_fil+(n-1)*Natoms_oneMT+1:Natoms_fil+n*Natoms_oneMT;
    Temp_atoms_one_MT(:,2)=One_MT_atoms(:,2)+(Nfil_fil+(n-1)*Nfil_oneMT);
    Temp_atoms_one_MT(:,3)=One_MT_atoms(:,3)+Natomtype_fil;
    
    Temp_bonds_one_MT(:,1)=Nbonds_fil+(n-1)*Nbonds_oneMT+1:Nbonds_fil+n*Nbonds_oneMT;
    Temp_bonds_one_MT(:,2)=One_MT_bonds(:,2)+Nbondtype_fil;
    Temp_bonds_one_MT(:,3:4)=One_MT_bonds(:,3:4)+(Natoms_fil+(n-1)*Natoms_oneMT);
    
    Temp_angles_one_MT(:,1)=Nangles_fil+(n-1)*Nangles_oneMT+1:Nangles_fil+n*Nangles_oneMT;
    Temp_angles_one_MT(:,2)=One_MT_angles(:,2)+Nangletype_fil;
    Temp_angles_one_MT(:,3:5)=One_MT_angles(:,3:5)+(Natoms_fil+(n-1)*Natoms_oneMT);
    
    Temp_dihedrals_one_MT(:,1)=Ndihedrals_fil+(n-1)*Ndihedrals_oneMT+1:Ndihedrals_fil+n*Ndihedrals_oneMT;
    Temp_dihedrals_one_MT(:,2)=One_MT_dihedrals(:,2)+Ndihedraltype_fil;
    Temp_dihedrals_one_MT(:,3:6)=One_MT_dihedrals(:,3:6)+(Natoms_fil+(n-1)*Natoms_oneMT);
    
    MT_atoms=[MT_atoms;Temp_atoms_one_MT]; 
    MT_bonds=[MT_bonds;Temp_bonds_one_MT];
    MT_angles=[MT_angles;Temp_angles_one_MT];
    MT_dihedrals=[MT_dihedrals;Temp_dihedrals_one_MT];
end
%---------------dimensionless the coordinates of atoms -------------------%

MT_atoms_lj(:,1:3)=MT_atoms(:,1:3);
MT_atoms_lj(:,4:6)=MT_atoms(:,4:6)/sigma;
Fil_atoms_lj(:,1:3)=atoms_fil(:,1:3);
Fil_atoms_lj(:,4:6)=atoms_fil(:,4:6)/sigma;

save MT_atoms_lj.txt MT_atoms_lj -ascii
save Fil_atoms_lj.txt Fil_atoms_lj -ascii

save MT_atoms.txt MT_atoms -ascii
save MT_bonds.txt MT_bonds -ascii
save MT_angles.txt MT_angles -ascii
save MT_dihedrals.txt MT_dihedrals -ascii

