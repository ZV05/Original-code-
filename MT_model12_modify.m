%%=======================================================================%%
%%================= Program for generate Microtubule ====================%%
%%=================      Bo Gong, May, 30, 2019      ====================%%    
%%=======================================================================%%

clc
clear

%========================== define variables ============================== 
N_profil=13;      % the number of protofilaments in each microtubule
N_tubulin=150;   % the number of tubulins in each profilament
D_tub=21.52866;   % diameter of the microtubule
Atoms=zeros(N_profil,6);

l_long=4.0;    % lattice parameter of longitudinal  纵向晶格参数
l_ring=5.2;    % lattice parameter of laterial      横向
delta=3*l_long/13;  % screw hight difference        螺旋高度

index_atom=0;  % atom index 

% 9.12 合并原子类型 1 和键的类型 1,2,角类型改为1，2，3

%=========================== Generate Atoms ==============================

atom_type=1; 
for i=1:N_profil               % 生成原子的坐标xyz
        x=0;
    y=D_tub/2*cos(3/2*pi-(i-1)*(2*pi/N_profil));   % 微管最底层的13个原子坐标
    z=D_tub/2*sin(3/2*pi-(i-1)*(2*pi/N_profil));
    
    for j=1:N_tubulin
%         if mod(j,2)==0         % 取模运算 j/2的余数
%             atom_type=2;       % belta tubulin 2的倍数
%         else 
%             atom_type=1;       % alpha tubulin
%         end
        X=x+(j-1)*l_long+delta*(i-1);       % 微管轴向坐标
        Y=y;
        Z=z;
        index_atom=index_atom+1;
        molecular_type=i;
        Atoms(index_atom,1)=index_atom;         %原子序号
        Atoms(index_atom,2)=molecular_type;     %分子类型（细丝编号）
        Atoms(index_atom,3)=atom_type;          %原子类型
        Atoms(index_atom,4)=X;
        Atoms(index_atom,5)=Y;
        Atoms(index_atom,6)=Z;
            
      plot3(X, Y, Z,"o")
      hold on    
    end       
end

save Atoms_oneMT_nano.txt Atoms -ascii

% ========================== Generate Bonds ===============================

N_atoms=length(Atoms(:,1));         % 1300个原子
% the bonds in each protofilament
n_bonds01=0; % 轴向键数量
bond_type01=1;
bond_type02=1;
bonds_part01=zeros(1,4);

for ii=1:N_profil
    Temp=find(Atoms(:,2)==ii);      % 找到等于ii的细丝编号
    for jj=1:N_tubulin-1            % 99个纵向键
        atom_01=Temp(jj,1);         % 一列一列提出来
        atom_02=Temp(jj+1,1);
        n_bonds01=n_bonds01+1;      % 键数+1
        if mod(jj,2)==1
            bond_type=bond_type01;  % 如果是奇数 01类型
        else
            bond_type=bond_type02;  % 偶数 02类型
        end
        bonds_part01(n_bonds01,:)=[n_bonds01, bond_type, atom_01, atom_02];
    end  
end

% the bonds for circumference       % 横向键
bond_type03=2;                      % 双 belta 双 alpha 键 
bond_type04=2;
bonds_part02=zeros(1,4);
bonds_part03=zeros(1,4);
bond_type05=2;                      % 第 13 根丝和第 1 根丝的键
bond_type06=2;
n_bonds02=0;
n_bonds03=0;

for ii=1:N_profil
    if ii<=N_profil-1
        Temp_01=find(Atoms(:,2)==ii);       % 找到ii编号的细丝
        Temp_02=find(Atoms(:,2)==ii+1);
        for jj=1:N_tubulin
            n_bonds02=n_bonds02+1;
            if mod(n_bonds02,2)==1
                bond_type=bond_type03;      % 键编号是奇数03类型
            else
                bond_type=bond_type04;      % 键编号是偶数04类型
            end
            atom_01=Temp_01(jj,1);          % ii丝的第jj个原子
            atom_02=Temp_02(jj,1);
            bonds_part02(n_bonds02,:)=[n_bonds02, bond_type, atom_01, atom_02];
        end
    else
        Temp_01=find(Atoms(:,2)==N_profil); % 第13 根丝
        Temp_02=find(Atoms(:,2)==1);        % 第 1 根丝
        for jj=1:N_tubulin-3
            n_bonds03=n_bonds03+1;
            if mod(n_bonds03,2)==1
                bond_type=bond_type05;
            else
                bond_type=bond_type06;
            end
            atom_01=Temp_01(jj,1);
            atom_02=Temp_02(jj+3,1);
            bonds_part03(n_bonds03,:)=[n_bonds03, bond_type, atom_01, atom_02];
        end
    end   
end

Bonds=[bonds_part01;bonds_part02;bonds_part03]; % 键数量编号 键类型 成键原子1 成键原子2 
N_bonds=length(Bonds(:,1));
Bonds(:,1)=1:N_bonds;
save Bonds_oneMT.txt Bonds -ascii

%========================== Generate Angles ===============================

% generate angles in each protofilament
angle_type01=1;                     % 纵向键角
n_angle01=0;
angles01=zeros(1,5);

for i=1:N_profil
    Temp=find(Atoms(:,2)==i);       % 第i根细丝
    for j=1:N_tubulin-2
        atom_01=Temp(j,1);          % 成角原子
        atom_02=Temp(j+1,1);
        atom_03=Temp(j+2,1);
        n_angle01=n_angle01+1;
        angles01(n_angle01,:)=[n_angle01, angle_type01, atom_01, atom_02, atom_03];
    end
end

% generate angles of the circumference orientation to prevent slip between
% two protofilaments

angle_type02=2;                     % 圆周当向上的键角
n_angle02=0;
angles02=zeros(1,5);                % 预留矩阵

angle_type03=3;
n_angle03=0;
angles03=zeros(1,5);

for ii=1:N_profil
    if ii==1
        Temp_01=find(Atoms(:,2)==ii+1);         % 成键丝第2根   
        Temp_02=find(Atoms(:,2)==ii);           % 第1根
        Temp_03=find(Atoms(:,2)==N_profil);     % 第13根
        for jj=1:N_tubulin
            if jj==1
                atom_01=Temp_01(jj,1);          % 第2根丝上第1个原子
                atom_02=Temp_02(jj,1);          % 第1根丝上第1个原子
                atom_03=Temp_02(jj+1,1);        % 第13根丝上第2个原子    （第13根丝上的第4个原子？） 
                n_angle02=n_angle02+1;
                angles02(n_angle02,:)=[n_angle02, angle_type02, atom_01, atom_02, atom_03];
            elseif jj<=3
                atom_01=Temp_01(jj,1);          % 第2根丝上的第2个原子
                atom_02=Temp_02(jj,1);          % 第1根丝上的第2个原子  
                atom_03=Temp_02(jj+1,1);        % 第1根丝上的第3个原子
                atom_04=Temp_02(jj-1,1);        % 第1根丝上的第1个原子
                n_angle02=n_angle02+1;
                angles02(n_angle02,:)=[n_angle02, angle_type02, atom_01, atom_02, atom_03];
                n_angle03=n_angle03+1;
                angles03(n_angle03,:)=[n_angle03, angle_type03, atom_01, atom_02, atom_04];
            elseif jj<=N_tubulin-1
                atom_01=Temp_01(jj,1);
                atom_02=Temp_02(jj,1);
                atom_03=Temp_02(jj+1,1);
                atom_04=Temp_02(jj-1,1);
                atom_05=Temp_03(jj-3,1);
                n_angle02=n_angle02+1;
                angles02(n_angle02,:)=[n_angle02, angle_type02, atom_01, atom_02, atom_03];
                n_angle02=n_angle02+1;
                angles02(n_angle02,:)=[n_angle02, angle_type02, atom_04, atom_02, atom_05];
                n_angle03=n_angle03+1;
                angles03(n_angle03,:)=[n_angle03, angle_type03, atom_01, atom_02, atom_04];
                n_angle03=n_angle03+1;
                angles03(n_angle03,:)=[n_angle03, angle_type03, atom_03, atom_02, atom_05];
            elseif jj==N_tubulin
                atom_01=Temp_01(jj,1);
                atom_02=Temp_02(jj,1);
                atom_04=Temp_02(jj-1,1);
                atom_05=Temp_03(jj-3,1);
                n_angle02=n_angle02+1;
                angles02(n_angle02,:)=[n_angle02, angle_type02, atom_04, atom_02, atom_05];
                n_angle03=n_angle03+1;
                angles03(n_angle03,:)=[n_angle03, angle_type03, atom_01, atom_02, atom_04];                              
            end
        end
    elseif ii<=N_profil-1                       
        Temp_01=find(Atoms(:,2)==ii+1);     
        Temp_02=find(Atoms(:,2)==ii);
        Temp_03=find(Atoms(:,2)==ii-1);
        for jj=1:N_tubulin
            if jj==1
                atom_01=Temp_01(jj,1);
                atom_02=Temp_02(jj,1);
                atom_03=Temp_02(jj+1,1);
                atom_05=Temp_03(jj,1);
                n_angle02=n_angle02+1;
                angles02(n_angle02,:)=[n_angle02, angle_type02, atom_01, atom_02, atom_03];
                n_angle03=n_angle03+1;
                angles03(n_angle03,:)=[n_angle03, angle_type03, atom_03, atom_02, atom_05];
            elseif jj<=N_tubulin-1
                atom_01=Temp_01(jj,1);
                atom_02=Temp_02(jj,1);
                atom_03=Temp_02(jj+1,1);
                atom_04=Temp_02(jj-1,1);
                atom_05=Temp_03(jj,1);
                n_angle02=n_angle02+1;
                angles02(n_angle02,:)=[n_angle02, angle_type02, atom_01, atom_02, atom_03];
                n_angle02=n_angle02+1;
                angles02(n_angle02,:)=[n_angle02, angle_type02, atom_04, atom_02, atom_05];
                n_angle03=n_angle03+1;
                angles03(n_angle03,:)=[n_angle03, angle_type03, atom_01, atom_02, atom_04];
                n_angle03=n_angle03+1;
                angles03(n_angle03,:)=[n_angle03, angle_type03, atom_03, atom_02, atom_05];
            elseif jj==N_tubulin
                atom_01=Temp_01(jj,1);
                atom_02=Temp_02(jj,1);
                atom_04=Temp_02(jj-1,1);
                atom_05=Temp_03(jj,1);
                n_angle02=n_angle02+1;
                angles02(n_angle02,:)=[n_angle02, angle_type02, atom_04, atom_02, atom_05];
                n_angle03=n_angle03+1;
                angles03(n_angle03,:)=[n_angle03, angle_type03, atom_01, atom_02, atom_04];
            end
        end
    elseif ii==N_profil
        Temp_01=find(Atoms(:,2)==1);
        Temp_02=find(Atoms(:,2)==ii);
        Temp_03=find(Atoms(:,2)==ii-1);
        for jj=1:N_tubulin
            if jj==1
                atom_01=Temp_01(jj+3,1);
                atom_02=Temp_02(jj,1);
                atom_03=Temp_02(jj+1,1);
                atom_05=Temp_03(jj,1);
                n_angle02=n_angle02+1;
                angles02(n_angle02,:)=[n_angle02, angle_type02, atom_01, atom_02, atom_03];
                n_angle03=n_angle03+1;
                angles03(n_angle03,:)=[n_angle03, angle_type03, atom_03, atom_02, atom_05];
            elseif jj<=N_tubulin-3
                atom_01=Temp_01(jj+3,1);
                atom_02=Temp_02(jj,1);
                atom_03=Temp_02(jj+1,1);
                atom_04=Temp_02(jj-1,1);
                atom_05=Temp_03(jj,1);
                n_angle02=n_angle02+1;
                angles02(n_angle02,:)=[n_angle02, angle_type02, atom_01, atom_02, atom_03];
                n_angle02=n_angle02+1;
                angles02(n_angle02,:)=[n_angle02, angle_type02, atom_04, atom_02, atom_05];
                n_angle03=n_angle03+1;
                angles03(n_angle03,:)=[n_angle03, angle_type03, atom_01, atom_02, atom_04];
                n_angle03=n_angle03+1;
                angles03(n_angle03,:)=[n_angle03, angle_type03, atom_03, atom_02, atom_05];
            elseif jj<=N_tubulin-1
                atom_02=Temp_02(jj,1);
                atom_03=Temp_02(jj+1);
                atom_04=Temp_02(jj-1,1);
                atom_05=Temp_03(jj,1);
                n_angle02=n_angle02+1;
                angles02(n_angle02,:)=[n_angle02, angle_type02, atom_04, atom_02, atom_05];
                n_angle03=n_angle03+1;
                angles03(n_angle03,:)=[n_angle03, angle_type03, atom_03, atom_02, atom_05];
            elseif jj==N_tubulin
                atom_02=Temp_02(jj,1);
                atom_04=Temp_02(jj-1,1);
                atom_05=Temp_03(jj,1);
                n_angle02=n_angle02+1;
                angles02(n_angle02,:)=[n_angle02, angle_type02, atom_04, atom_02, atom_05];
            end
        end       
    end
end


Angles=[angles01; angles02; angles03];
N_angles=length(Angles(:,1));
Angles(:,1)=1:N_angles;
save Angles_oneMT.txt Angles -ascii

%========================= Generate Dihedrals =============================

% generate dihedrals in each three protofilaments

dihedral_type=1;
n_dihedrals=0;
dihedrals=zeros(1,6);

for ii=1:N_profil
    if ii==1
        Temp_01=find(Atoms(:,2)==ii+1);         % 第2根丝
        Temp_02=find(Atoms(:,2)==ii);           % 第1根丝
        Temp_03=find(Atoms(:,2)==N_profil);     % 第13根丝
        for jj=3:N_tubulin-1
            if jj==3
                atom_01=Temp_01(jj,1);          % 第2根丝第3个原子
                atom_02=Temp_02(jj,1);          % 第1根丝第3个原子
                atom_03=Temp_02(jj+1,1);        % 第1根丝第4个原子
                atom_04=Temp_03(jj-2,1);        % 第13根丝第1个原子
                n_dihedrals=n_dihedrals+1;
                dihedrals(n_dihedrals,:)=[n_dihedrals, dihedral_type, atom_01, atom_02, atom_03, atom_04];
            else
                atom_01=Temp_01(jj,1);          % 2，5   （以5为例）
                atom_02=Temp_02(jj,1);          % 1，5
                atom_03=Temp_02(jj+1,1);        % 1，6
                atom_04=Temp_03(jj-2,1);        % 13，3
                atom_05=Temp_01(jj+1,1);        % 2，6
                atom_06=Temp_03(jj-3,1);        % 13，2
                n_dihedrals=n_dihedrals+1;
                dihedrals(n_dihedrals,:)=[n_dihedrals, dihedral_type, atom_01, atom_02, atom_03, atom_04];
                n_dihedrals=n_dihedrals+1;
                dihedrals(n_dihedrals,:)=[n_dihedrals, dihedral_type, atom_05, atom_03, atom_02, atom_06];
            end                
        end
    elseif ii<=N_profil-1                       % 以5根丝为例
        Temp_01=find(Atoms(:,2)==ii+1);
        Temp_02=find(Atoms(:,2)==ii);
        Temp_03=find(Atoms(:,2)==ii-1);
        for jj=1:N_tubulin-1
            atom_01=Temp_01(jj,1);              % 6，5
            atom_02=Temp_02(jj,1);              % 5，5
            atom_03=Temp_02(jj+1,1);            % 5，6
            atom_04=Temp_03(jj+1,1);            % 4，6
            atom_05=Temp_01(jj+1,1);            % 6，6
            atom_06=Temp_03(jj,1);              % 4，5
            n_dihedrals=n_dihedrals+1;
            dihedrals(n_dihedrals,:)=[n_dihedrals, dihedral_type, atom_01, atom_02, atom_03, atom_04];
            n_dihedrals=n_dihedrals+1;
            dihedrals(n_dihedrals,:)=[n_dihedrals, dihedral_type, atom_05, atom_03, atom_02, atom_06];
        end
    elseif ii==N_profil                         %第13根丝
        Temp_01=find(Atoms(:,2)==1);
        Temp_02=find(Atoms(:,2)==ii);
        Temp_03=find(Atoms(:,2)==ii-1);
        for jj=1:N_tubulin-3
            if jj==N_tubulin-3
                atom_01=Temp_01(jj+3,1);        % 1，8
                atom_02=Temp_02(jj,1);          % 13，5
                atom_03=Temp_02(jj+1,1);        % 13，6
                atom_04=Temp_03(jj+1,1);        % 12，6
                n_dihedrals=n_dihedrals+1;
                dihedrals(n_dihedrals,:)=[n_dihedrals, dihedral_type, atom_01, atom_02, atom_03, atom_04];
            else
                atom_01=Temp_01(jj+3,1);        % 1，4
                atom_02=Temp_02(jj,1);          % 13，1
                atom_03=Temp_02(jj+1,1);        % 13，2
                atom_04=Temp_03(jj+1,1);        % 12，2
                atom_05=Temp_01(jj+4,1);        % 1，5
                atom_06=Temp_03(jj,1);          % 12，1
                n_dihedrals=n_dihedrals+1;
                dihedrals(n_dihedrals,:)=[n_dihedrals, dihedral_type, atom_01, atom_02, atom_03, atom_04];
                n_dihedrals=n_dihedrals+1;
                dihedrals(n_dihedrals,:)=[n_dihedrals, dihedral_type, atom_05, atom_03, atom_02, atom_06];
            end
        end
    end
end

Dihedrals=dihedrals;
N_dihedrals=length(Dihedrals(:,1));
Dihedrals(:,1)=1:N_dihedrals;
save Dihedrals_oneMT.txt Dihedrals -ascii

%=================== Generate Boundary Atom Indexs ========================

% Lower=zeros(1,N_profil);
% Upper=zeros(1,N_profil);
% n_low=0;
% n_up=0;
% 
% for ii=1:N_profil
%     Temp=find(Atoms(:,2)==ii);
%     atom_lower=Temp(1,1);
%     atom_upper=Temp(end,1);
%     n_low=n_low+1;
%     n_up=n_up+1;
%     Lower(1,n_low)=atom_lower;
%     Upper(1,n_up)=atom_upper;
% end

lower_index=find(Atoms(:,6)<=14);
Lower=lower_index';

high_pz=Atoms(end,6);
upper_index=find(Atoms(:,6)>=high_pz-14);
Upper=upper_index';

save Lower.txt Lower -ascii
save Upper.txt Upper -ascii

%%%%%%%%%%%%%%%%%%%%%%%%%%%生成微管外挂二聚体%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%=========================生成单根微管的反应位点===========================%

% num_group_fil_atoms=25;         % 原丝上每组的原子数量
% num_mt_fil_beads=125;           % 微管上每根原丝的原子数量
% num_group=5;
% change_atoms=[];
% % 随机生成所在纤维数量
% s1 = randi([3,11],1,num_group);          % 取4-9随机数避免微管接缝处的反应位点
% % 5段原子反应位点所在纤维的序号
% s2 = randi([2,23],1,num_group);          % 取2-4防止出现反应位点相邻
% for k=1:num_group   % 共5组数据
%     c1=s1(1,k);    % 纤维序号
%     c2=s2(1,k);    % 原子顺序号码
%     % 筛选第c1根纤维数据
%     start_fil_idx = (c1 - 1) * num_mt_fil_beads + 1;
%     over_fil_idx=start_fil_idx+num_mt_fil_beads-1;
%     atoms_fil_mt=Atoms(start_fil_idx:over_fil_idx,:);
%     % 在c1根纤维上筛选第k组数据，找到第c2个原子
%     start_group_idx=(k - 1) * num_group_fil_atoms + 1;
%     over_group_idx=start_group_idx+num_group_fil_atoms-1;
%     atoms_group=atoms_fil_mt(start_group_idx:over_group_idx,:);
%     
%     % 保存所筛选原子的原子序号
%     temp_atoms=atoms_group(c2,1);
%     change_atoms=[change_atoms;temp_atoms];     % 每根微管上的反应位点
%  end
% 
% % 修改反应位点的原子类型
% 
% a1=size(change_atoms,1);
% a2=size(Atoms,1);
% 
% for m=1:a1
%     for n=1:a2
%         temp_data=Atoms(n,1);
%         temp_change=change_atoms(m,1);
%         if temp_data==temp_change
%             Atoms(n,3)=2;   % 原子类型改为2
%         end
%     end
% end
% 
% %=========================反应位点挂带二聚体===============================%
% 
% % react_point=load('reactpoint.txt');
% 
% react_point=change_atoms;
% 
% num_fil=0;
% 
% D_tub_dimer1=29.52866;      % 二聚体第一个原子的微管直径
% D_tub_dimer2=37.52866;      % 二聚体第二个原子的微管直径
% % idxdimer1=0;
% idxdimer1=size(Atoms,1);    % 初始微管原子序号
% idxdimer1_fil=13;
% idxdimer2=0;
% idxdimer2_fil=13;
% 
% dimer_type=2;               % 二聚体原子类型
% 
% idx_bonds=size(Bonds,1);
% dimer_bond_type=3;
% 
% len_atoms=size(Atoms,1);
% len_react=size(react_point,1);
% dimer1_atoms=[];
% dimer2_atoms=[];
% bonds_dimer_mt=[];
% 
% % 在反应为点生成二聚体坐标
% 
% for n=1:len_react
%     react_idx=react_point(n,1);
%     for m=1:len_atoms
%         atoms_idx=Atoms(m,1);
%         if react_idx==atoms_idx
%             dimer_profil=mod((Atoms(m,2)-num_fil),N_profil);
% 
%             dimer1_x=Atoms(m,4);
%             dimer1_y=D_tub_dimer1/2*cos(3/2*pi-(dimer_profil-1)*(2*pi/N_profil));
%             dimer1_z=D_tub_dimer1/2*sin(3/2*pi-(dimer_profil-1)*(2*pi/N_profil));
% 
%             dimer2_x=Atoms(m,4);
%             dimer2_y=D_tub_dimer2/2*cos(3/2*pi-(dimer_profil-1)*(2*pi/N_profil));
%             dimer2_z=D_tub_dimer2/2*sin(3/2*pi-(dimer_profil-1)*(2*pi/N_profil));
%             
%             idxdimer1=idxdimer1+1;
%             idxdimer1_fil=idxdimer1_fil+1;
%             dimer1=[idxdimer1,idxdimer1_fil,dimer_type,dimer1_x,dimer1_y,dimer1_z];
%             % 微管反应位点与二聚体成键
%             bonds_dimer_mt1=[idx_bonds,dimer_bond_type,react_idx,idxdimer1];
%             bonds_dimer_mt=[bonds_dimer_mt;bonds_dimer_mt1];
% 
%             idxdimer1=idxdimer1+1;
%             idxdimer2_fil=idxdimer1_fil;
%             dimer2=[idxdimer1,idxdimer2_fil,dimer_type,dimer2_x,dimer2_y,dimer2_z];
% 
%             dimer1_atoms=[dimer1_atoms;dimer1;dimer2];  % 生成二聚体的原子信息
% 
%         end
%     end
% end
% 
% % 将二聚体的原子信息加入之前原子信息
% Atoms=[Atoms;dimer1_atoms];
% Bonds=[Bonds;bonds_dimer_mt];
% % save Atoms_oneMT_nano.txt Atoms -ascii
% 
% %============================反应位点和二聚体成键==========================%
% 
% 
% 
% %================================二聚体成键===============================%
% len_dim_atoms=size(dimer1_atoms,1);
% % idx_bonds=0;
% idx_bonds=size(Bonds,1);
% bonds_dimer=[];
% for nn=1:2:len_dim_atoms
%     bonds_dim_atoms=dimer1_atoms(nn,1);
%     idx_bonds=idx_bonds+1;
%     dimer_bond_type=3;
%     bonds_dimer1=[idx_bonds,dimer_bond_type,bonds_dim_atoms,bonds_dim_atoms+1];
%     bonds_dimer=[bonds_dimer;bonds_dimer1];
% end
% 
% % 将二聚体键的信息加入之前键的信息
% Bonds=[Bonds;bonds_dimer];
% save Bonds_oneMT.txt Bonds -ascii

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%====================== Dimensionless Atoms =============================%

Atoms_lj_part2=Atoms(:,4:6)/l_long;
Atoms_lj_part1=Atoms(:,1:3);
Atoms_lj=[Atoms_lj_part1, Atoms_lj_part2];
save Atoms_lj.txt Atoms_lj -ascii