function [d2H] = opf_Gf_hess(x,lam,vv,mpc)
%hess of limit of gas pipeline
%   Detailed explanation goes here
[Va,Vm,Pg,Qg,Prs,PGs,LCg] = deal(x{:});
Pg = 100 * Pg;
nx = vv.iN.LCg;%last variable
d2H = zeros(nx,nx);
[ d2HPrs ] = d2Gasf_dPrs2( Prs,mpc,lam );
d2H(vv.i1.Prs:vv.iN.Prs,vv.i1.Prs:vv.iN.Prs) = d2HPrs;
function [ d2HPrs ] = d2Gasf_dPrs2( Prs,mpc,muGl )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
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
%% 对Gf^2-Gfmax^2求偏导
d2HPrs=zeros(Gnb,Gnb);%最后输出的
d2Gf_dPrs2=zeros(Gnl,Gnb,Gnb);%中间过程用的
%由于h不等式约束平方了，所以没有方向上的正负
for i=1:Gnl
    fb=mpc.Gline(i,1);tb=mpc.Gline(i,2);
    if Prs(fb)>Prs(tb)
       d2Gf_dPrs2(i,fb,fb)=2*GY(fb,tb)^2;
       d2Gf_dPrs2(i,tb,tb)=-2*GY(fb,tb)^2;
       %偏pi偏pj=0
    end
    if Prs(fb)<Prs(tb)
        d2Gf_dPrs2(i,fb,fb)=-2*GY(fb,fb)^2;
        d2Gf_dPrs2(i,tb,tb)=2*GY(tb,tb)^2;
    end
    %气压相等的时候等于0，但是给点微小的值？让其不能限于局部最优？
%      if Prs(fb)==Prs(tb)
%         d2Gf_dPrs2(i,fb,fb)=0;
%         d2Gf_dPrs2(i,tb,tb)=0;
%      end
end
%% 整合
for i=1:Gnl
    d2HPrs=d2HPrs+squeeze(muGl(i)*d2Gf_dPrs2(i,:,:));
end
end

end

