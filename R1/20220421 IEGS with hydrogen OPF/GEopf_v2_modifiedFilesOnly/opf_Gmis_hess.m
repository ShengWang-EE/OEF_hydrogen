function [d2G] = opf_Gmis_hess(x,lam,vv,mpc)
%hess function for nodal gas power balance
%   Detailed explanation goes here
[Va,Vm,Pg,Qg,Prs,PGs,LCg] = deal(x{:});
nx = vv.iN.LCg;%last variable
d2G = zeros(nx,nx);
[d2GPrs] = d2PGb_dPrs2(Prs,mpc,lam);
d2G(vv.i1.Prs:vv.iN.Prs,vv.i1.Prs:vv.iN.Prs) = d2GPrs;
%%
function [ d2GPrs] = d2PGb_dPrs2(Prs,mpc,lamGb)
%% numbers
Gnb=size(mpc.Gbus,1);Gnl=size(mpc.Gline,1);
% GY
GY = zeros(Gnl,Gnl);
for k = 1:Gnl
    I = mpc.Gline(k,1);
    J = mpc.Gline(k,2);
    GY(I,J)= mpc.Gline(k,3);
    GY(J,I)= mpc.Gline(k,3);
end
%% ���Ŷ�sigma(Gflow)ƫ���������������f�Ա���ƫ����������ۼ�����
d2Gf_dPrs2=zeros(Gnl,Gnb,Gnb);%ÿ����·һ������
d2GPrs=zeros(Gnb,Gnb);
for i=1:Gnl
    fb=mpc.Gline(i,1);tb=mpc.Gline(i,2);
    if Prs(fb)>Prs(tb)
       d2Gf_dPrs2(i,fb,fb)=GY(fb,tb)*(-Prs(tb)^2)/(Prs(fb)^2-Prs(tb)^2)^1.5;%d2Gf_dpi2
       d2Gf_dPrs2(i,fb,tb)=GY(fb,tb)*Prs(fb)*Prs(tb)/(Prs(fb)^2-Prs(tb)^2)^1.5;%d2Gf_dpidpj
       d2Gf_dPrs2(i,tb,fb)=d2Gf_dPrs2(i,fb,tb);%%d2Gf_dpjdpi
       d2Gf_dPrs2(i,tb,tb)=GY(fb,tb)*(-Prs(fb)^2)/(Prs(fb)^2-Prs(tb)^2)^1.5;%d2Gf_dpi2
    end
    if Prs(fb)<Prs(tb)
       d2Gf_dPrs2(i,fb,fb)=GY(fb,tb)*(Prs(tb)^2)/(Prs(tb)^2-Prs(fb)^2)^1.5;%d2Gf_dpi2
       d2Gf_dPrs2(i,fb,tb)=GY(fb,tb)*(-Prs(fb)*Prs(tb))/(Prs(tb)^2-Prs(fb)^2)^1.5;%d2Gf_dpidpj
       d2Gf_dPrs2(i,tb,fb)=d2Gf_dPrs2(i,fb,tb);%%d2Gf_dpjdpi
       d2Gf_dPrs2(i,tb,tb)=GY(fb,tb)*(Prs(fb)^2)/(Prs(tb)^2-Prs(fb)^2)^1.5;%d2Gf_dpi2
    end
%     if Prs(fb)==Prs(tb)
%        d2Gf_dPrs2(i,fb,fb)=9999999999;d2Gf_dPrs2(i,fb,tb)=9999999999;
%        d2Gf_dPrs2(i,tb,fb)=9999999999;d2Gf_dPrs2(i,tb,tb)=9999999999;
%     end
% ������ȵ�ʱ����Ե���0
end
%% ���սڵ��ۼ�
for i=1:Gnb%Լ��������
    %squeeze���ã�����ά�����ά
    d2GPrs=d2GPrs+squeeze(lamGb(i)*(sum(d2Gf_dPrs2(find(i==mpc.Gline(:,1)),:,:),1)-sum(d2Gf_dPrs2(find(i==mpc.Gline(:,2)),:,:),1)));
end

end
end

