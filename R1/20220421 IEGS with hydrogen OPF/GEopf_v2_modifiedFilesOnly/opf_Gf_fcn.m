function [h,dh] = opf_Gf_fcn(x,vv,mpc)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
[Va,Vm,Pg,Qg,Prs,PGs,LCg] = deal(x{:});
h = add_h(Prs,vv,mpc);
dh =add_dh(Prs,vv,mpc);

    function h = add_h(Prs,vv,mpc)
    %UNTITLED2 Summary of this function goes here
    %   Detailed explanation goes here
    %% ����mpc�Ĳ�����x��״̬������·�ϵ�����������һ����ʽ�����漰������ϵͳ�ĳ�������͵���
    %% unpack x
%     Prs= x(vv.i1.Prs:vv.iN.Prs);
    % numbers
    Gnb=size(mpc.Gbus,1);
    Gnl=size(mpc.Gline,1);
    %flow max
    Gasflow_max = (mpc.Gline(:,5)).^2;
    Gasflow_max(Gasflow_max == 0) = Inf;
    %% ������·����������Ϊ����
    GY = zeros(vv.N.Prs,vv.N.Prs);
    for k = 1:Gnl
        I = mpc.Gline(k,1);
        J = mpc.Gline(k,2);
        GY(I,J)= mpc.Gline(k,3);
        GY(J,I)= mpc.Gline(k,3);
    end
    Gasflow=zeros(Gnl,1);
    for i=1:Gnl
        fb=mpc.Gline(i,1);tb=mpc.Gline(i,2);
        Sub = abs(Prs(fb)^2-Prs(tb)^2);%�ڵ��һ������ڶ���Ϊ��
        Gasflow(i)=sign(Prs(fb)-Prs(tb))* GY(fb,tb)*sqrt(Sub);
    end
    %% ����addh
    addh=Gasflow.^2-Gasflow_max;%Gasflow_max�Ѿ�ƽ������

    h = addh;
    end
    function dh =add_dh(Prs,vv,mpc)
    %UNTITLED Summary of this function goes here
    %   Detailed explanation goes here
    %% unpack x

%     Prs= x(vv.i1.Prs:vv.iN.Prs);
    % numbers
    Gnb=size(mpc.Gbus,1);Gnl=size(mpc.Gline,1);
    GY = zeros(vv.N.Prs,vv.N.Prs);
    for k = 1:Gnl
        I = mpc.Gline(k,1);
        J = mpc.Gline(k,2);
        GY(I,J)= mpc.Gline(k,3);
        GY(J,I)= mpc.Gline(k,3);
    end
    dGf_dPrs=zeros(vv.N.Prs,Gnl);
    for i=1:Gnl
        fb=mpc.Gline(i,1);tb=mpc.Gline(i,2);%frombus and tobus
        if Prs(fb)>Prs(tb)
            dGf_dPrs(fb,i)=2*GY(fb,tb)^2*Prs(fb);
            dGf_dPrs(tb,i)=-2*GY(fb,tb)^2*Prs(tb);
        end
        if Prs(fb)<Prs(tb)
            dGf_dPrs(fb,i)=-2*GY(fb,tb)^2*Prs(fb);
            dGf_dPrs(tb,i)=2*GY(fb,tb)^2*Prs(tb);
        end
    %     if Prs(fb)==Prs(tb)
    %         dGf_dPrs(fb,i)=0;
    %         dGf_dPrs(tb,i)=0;%���ձ��ʽ�ڸõ�û�е���ֵ����ȡ���ߵ�ƽ��=0��������
    %     end
    end
    %% ƴ��
    % adddh=zeros(nxyz,Gnl);
    adddh=[zeros(vv.iN.Qg,Gnl);dGf_dPrs;zeros(vv.N.PGs+vv.N.LCg,Gnl)];
    dh=adddh';  
    end
end

