function [g,dg] = opf_Gmis_fcn(x,vv,mpc)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
[Va,Vm,Pg,Qg,Prs,PGs,LCg] = deal(x{:});
Pg = 100 * Pg; %100MVA
%exclude the LCe (because of LCe, the order of original gen is mixed
LCeIndex = find(mpc.gen(:,22)==2);
Pg(LCeIndex,:)=[];
mpc.gen(LCeIndex,:) = [];
g = calculate_g(Pg,Prs,PGs,LCg,vv,mpc);
dg = calculate_dg(Pg,Prs,PGs,LCg,vv,mpc);
%%
    function g = calculate_g(Pgen,Prs,PGs,LCg,vv,mpc)
        % numbers
        Gnb=size(mpc.Gbus,1);Gnl=size(mpc.Gline,1);
        %% 首先计算不涉及到线路注入潮流的几个分量
        %计算节点本身净产气量，注意要用变量x里面的值
        Pgd=zeros(Gnb,1); % nodal gas injection
        gtpEfficiency = 200;
        for i=1:Gnb
            Pgd(i)=-mpc.Gbus(i,3);
            if ismember(i,mpc.Gsou(:,1))
               Pgd(i)=Pgd(i)+sum(PGs(find(i==mpc.Gsou(:,1))));%-负荷+气源
            end
        end
        % we need to consider the gas consumption of gfu
        for i = 1:size(mpc.gfuIndex,1)
            electricityBus = mpc.gen(mpc.gfuIndex(i),1);
            gasBus = mpc.GEcon(find(mpc.GEcon(:,2)==electricityBus),1);
            Pgd(gasBus)=Pgd(gasBus)-Pgen(mpc.gfuIndex(i)) / gtpEfficiency;%找到相应的Ggtp编号,再-GTP
        end
        %% 然后计算线路流量（流出为正）
        GY = zeros(vv.N.Prs,vv.N.Prs);
        for k = 1:Gnl
            I = mpc.Gline(k,1);
            J = mpc.Gline(k,2);
            GY(I,J)= mpc.Gline(k,3);
            GY(J,I)= mpc.Gline(k,3);
        end
        dGP = zeros(Gnb,1);

        for i=1:Gnb%维数是节点
            for j=1:Gnb
                  Sub = abs(Prs(i)^2-Prs(j)^2);
                  dGP(i,1)=dGP(i,1)+sign(Prs(i)-Prs(j))* GY(i,j)*sqrt(Sub);%流出
            end
        end
        dPgd=Pgd-dGP;  
        addg=dPgd;
        %% 然后修改对应的LCg
        LCg_bus = find(mpc.Gbus(:,3)~=0);
        addg(LCg_bus) = addg(LCg_bus) + LCg;

        g = addg;
    end
    function dg = calculate_dg(Pg,Prs,PGs,LCg,vv,mpc)
        % numbers
        Gnb=size(mpc.Gbus,1);Gnl=size(mpc.Gline,1);
        % GY
        GY = zeros(vv.N.Prs,vv.N.Prs);
        for k = 1:Gnl
            I = mpc.Gline(k,1);
            J = mpc.Gline(k,2);
            GY(I,J)= mpc.Gline(k,3);
            GY(J,I)= mpc.Gline(k,3);
        end
        %% 
        % adddg=zeros(size(x,1),Gnb);
        %% g=PGs-PGd-Ggtp+Gptg-sigma(Gflow)
        %% 首先对+PGs-PGd-Ggtp+Gptg偏导（只对有PGs的节点值为1）
        dg_dPGs=zeros(vv.N.PGs,Gnb);dg_dGgtp=zeros(vv.N.Pg,Gnb);
        gfuEfficiency = 200;
        for i=1:Gnb
            %PGs
            if ismember(i,mpc.Gsou(:,1))
               dg_dPGs(find(mpc.Gsou(:,1)==i),i)=1;
            end
        end
        % set the same demenssion as mpc.gen, than renew dg_dgfu by adding to
        % dg_dgen
        for i = 1:size(mpc.gfuIndex,1)
            electricityBus = mpc.gen(mpc.gfuIndex(i),1);
            gasBus = mpc.GEcon(find(mpc.GEcon(:,2)==electricityBus),1);
            dg_dGgtp(gasBus)=dg_dGgtp(gasBus) - 1 / gfuEfficiency;%找到相应的Ggtp编号,再-GTP
        end
        %% 接着对sigma(Gflow)偏导：可以先求各个f对变量偏导，最后再累加起来
        %但是在prs(i)=prs(j)的时候导数无穷大
        df_dPrs=zeros(Gnb,Gnl);%最后不用这个形式的
        dPgd_dPrs=zeros(Gnb,Gnb);%最后用这个形式的
        for i=1:Gnl
            fb=mpc.Gline(i,1);tb=mpc.Gline(i,2);%from bus and to bus
            if Prs(fb)>=Prs(tb)
            df_dPrs(fb,i)=GY(fb,tb)*Prs(fb)/sqrt(Prs(fb)^2-Prs(tb)^2);
            df_dPrs(tb,i)=-GY(fb,tb)*Prs(tb)/sqrt(Prs(fb)^2-Prs(tb)^2);
            end
            if Prs(fb)<Prs(tb)
            df_dPrs(fb,i)=GY(fb,tb)*Prs(fb)/sqrt(Prs(tb)^2-Prs(fb)^2);
            df_dPrs(tb,i)=-GY(fb,tb)*Prs(tb)/sqrt(Prs(tb)^2-Prs(fb)^2);
            end
        %     if Prs(fb)==Prs(tb)%把无穷大变成有限值，待测试
        %        df_dPrs(i,fb)=99;
        %        df_dPrs(i,tb)=-99;
        %     end
        end
        %累加得到节点流出功率对Prs偏导

        for iF=1:Gnb
                connect_fb=find(mpc.Gline(:,1)==iF);connect_tb=find(mpc.Gline(:,2)==iF);
                dPgd_dPrs(:,iF)=sum(df_dPrs(:,connect_fb'),2)-sum(df_dPrs(:,connect_tb'),2);
        %     dPgd_dPrs(:,i)=sum(df_dPrs(:,(find(mpc.Gline(:,1)==i))'),2)-sum(df_dPrs(:,(find(mpc.Gline(:,2)==i))'),2);%横着相加，流出的减去流入的
        end
        %-sigma
        dPgd_dPrs=-dPgd_dPrs;
        %% 将各个矩阵方块拼接
        adddg=[zeros((vv.N.Va),Gnb);zeros((vv.N.Vm),Gnb);dg_dGgtp;zeros((vv.N.Qg),Gnb);dPgd_dPrs;dg_dPGs];
        % 修正LCg
        LCg_bus=find(mpc.Gbus(:,3)~=0);
        for i=1:vv.N.LCg
            adddg((vv.i1.LCg-1+i),LCg_bus(i))=1;
        end
        dg = adddg'; 
    end
end

