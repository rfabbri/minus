%Análise dos experimentos realizados com P problemas gerados de forma aleatoria
%MAI25

clear all
close all

addpath(genpath('/home/juliana/Documents/doutorado/pesquisa/pesquisa2025/Experimentos'));
addpath(genpath('frpt'));
P=140; %numero de problemas
aux=0;
for np=1:P %número de problemas
    i=string(np);
    str1=join(['gt_sol-P' i 'txt'], ["", "."]); %Quando não encontra a solucão escreve 313
    str2=join(['Nrealsol-P' i 'txt'],["", "."]);
    str3=join(['sol_step-P' i 'txt'],["", "."]); %O arquivo sol_step-Pi.txt tem 1 coluna com o número de iteracões em cada raiz
    str4=join(['vec_gamma-P' i 'txt'],["", "."]);
    gtsol(500*(np-1)+1:500*np,:)=load(str1);
    real=load(str2);
    raizrealN(aux+1:aux+length(real),:)=real; %raízes reais que não são solucão
    aux=aux+length(real);
    steps(312*500*(np-1)+1:312*500*np,:)=load(str3);
    %randoms(248*500*(np-1)+1:248*500*np,:)=load(str4);  % COM ERRO

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NÚMERO DE ITERAÇOES

stepNsol=steps; % stepNsol tem o número de iteracões nas raízes que Não são solucão
 for i=1:500*np
      gtsol(i,2)=gtsol(i,1)==313; 
      if gtsol(i,1)~=313
          stepNsol((i-1)*312+gtsol(i,1)+1,1)=NaN;
          stepsol(i,1)=steps((i-1)*312+gtsol(i,1)+1,1); %# de iteracoes na raiz que é solucao no teste i
      else
          stepsol(i,1)=NaN;
      end
 end

 % Boxplot do número de iteracões nas raízes que SÃO e nas que NÃO SÃO solucão 
 
 for i=1:np
        stepsolP(:,i)=stepsol((i-1)*500+1:(i-1)*500+500,1);
        stepNsolP(:,i)=stepNsol((i-1)*312*500+1:312*500*i,1);
 end
 
 [a,b]=min(stepsolP);  % a=min value, b=position
 [c,d]=max(stepsolP);
 
 % LEGENDA DOS PROBLEMAS NO EIXO HORIZONTAL
 for i=1:np
     s=string(i);
     j=join(['P' s],[""]);
     strp(i)=j;   % P+número
 end
 
% for i=1:np
%     MeP(i)=median(stepsolP(:,i),'omitnan'); % median value for each problem
% end

MeS=median(stepsol,'omitnan');  % median of all prblems - gt solution
xm=linspace(0,140,280);
ym=MeS+0*xm;

MeNS=median(stepNsol,'omitnan');  % median of all prblems - not gt solutions
ymN=MeNS+0*xm;

 % To print all #min interation
 for i=1:np/2
     xt(2*i)=2*i;
     xt(2*i-1)=2*i-1;
     yt(2*i)=a(2*i)-7;
     yt(2*i-1)=a(2*i-1)-7;
     strt(2*i)=string(a(2*i));
     strt(2*i-1)=string(a(2*i-1));
 end

 % To print only the #min iteration<35
% for i=1:np
%      if a(i)<35
%          amin(i)=a(i);   % minimum number of iterations less than 35
%      else
%          amin(i)=NaN;
%      end
%  end
% 
%  for i=1:np/2
%      xt(2*i)=2*i-0.6;
%      xt(2*i-1)=2*i-1-0.6;
%      yt(2*i)=a(2*i)-7;
%      yt(2*i-1)=a(2*i-1)-7;
%      strt(2*i)=string(amin(2*i));
%      strt(2*i-1)=string(amin(2*i-1));
%  end
   
 %-------------------------------------------------------------------------------------------------------------%
% PRINT ALL 140 PROBLEMS TOGETHER


figure
boxplot(stepsolP, 'Symbol', '+k', PlotStyle='compact')
m = findobj(gcf, 'Tag', 'Outliers');
set(m, 'color', [0.9 0.9 0.9], 'markerfacecolor',[0.9 0.9 0.9], 'markeredgecolor', [0.7 0.7 0.7], 'markersize',1);
hold on
plot(a,'o','MarkerSize', 4, 'MarkerFaceColor', 'r');
plot(xm,ym,'Color','g', 'LineWidth',2);
hold off
text(xt,yt,strt)
ylabel('\fontsize{18}Number of iterations per root')
xlabel('\fontsize{18}Problem')
title('\fontsize{20}Time spent in paths leading to ground-truth solution')

figure
boxplot(stepNsolP, 'Symbol', '+k', PlotStyle="compact");
m = findobj(gcf, 'Tag', 'Outliers');
set(m, 'color', [0.9 0.9 0.9], 'markerfacecolor',[0.9 0.9 0.9], 'markeredgecolor', [0.7 0.7 0.7], 'markersize',1);
hold on
plot(xm,ymN,'Color','g', 'LineWidth',2);
hold off
ylabel('\fontsize{18}Number of iterations per root')
xlabel('\fontsize{18}Problem')
title('\fontsize{20}Time spent in paths leading to undesired roots')

 %-----------------------------------------------------------------------------------------------------------%

% figure
% boxplot(stepsolP(:,1:70),strp(1:70));
% hold on
% plot(a(1:70), 'o','MarkerSize', 7, 'MarkerFaceColor', 'r');
% hold off
% text(xt(1:70),yt(1:70),strt(1:70))
% ylabel('\fontsize{18}Number of iterations per root')
% xlabel('\fontsize{18}Problem')
% title('\fontsize{20}Time spent in paths leading to ground-truth solution')
% 
% figure
% boxplot(stepsolP(:,71:140),strp(71:140));
% hold on
% plot(a(71:140), 'o','MarkerSize', 7, 'MarkerFaceColor', 'r');
% hold off
% text(xt(71:140)-70,yt(71:140),strt(71:140))
% ylabel('\fontsize{18}Number of iterations per root')
% xlabel('\fontsize{18}Problem')
% title('\fontsize{20}Time spent in paths leading to ground-truth solution')
% 
% figure
% boxplot(stepNsolP(:,1:70),strp(1:70));
% ylabel('\fontsize{18}Number of iterations per root')
% xlabel('\fontsize{18}Problem')
% title('\fontsize{20}Time spent in paths leading to undesired roots')
% 
% figure
% boxplot(stepNsolP(:,71:140),strp(71:140));
% ylabel('\fontsize{18}Number of iterations per root')
% xlabel('\fontsize{18}Problem')
% title('\fontsize{20}Time spent in paths leading to undesired roots')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TAXA DE ACERTO
%Contar quantas vezes acha a solução verdadeira

% ---------------------------------------------------------------------------------------%
% POR PROBLEMA
numgtsol=zeros(np,1); %matriz com # de vezes em que acha a solução verdadeira em cada problema
for i=1:np
    aux=sum(gtsol((i-1)*500+1:i*500,2)); %quantas vezes não acha solução verdadeira
    numgtsol(i,1)=500-aux;
end

figure
plot(0.2*numgtsol,'o');
title('Taxa de acerto');

txm=0.2*sum(numgtsol)/103; %tx média de encontro da solução verdadeira

%---------------------------------------------------------------------------------------------%
%TODOS OS PROBLEMAS JUNTOS de acordo com # max iterações
nt=["500" "400" "300" "200" "190" "180" "170" "160" "150" "130" "100" "90" "70" "50"]; %#max iterações
S=zeros(length(nt),10);  %S: matriz que conta em quantas rodadas encontra solução verdadeira pelo menos 1 vez, de acordo com o número de testes por rodada(coluna) e #máx de iterações (linha)
TA=zeros(length(nt),10); %taxa de acerto de S
for k=1:length(nt)
    rn=str2num(nt(k)); %#max iterações
    gtsolk(:,k)=gtsol(:,2); %mtz binária: 0:#iterações<=#max iterações, 1:caso contrário; cada coluna corresponde a um #max
    for r=1:length(gtsol)
        if stepsol(r,1)>rn %#iterações na raiz > #max iterações
            gtsolk(r,k)=1;
        end
    end
    for i=1:10%# testes por rodada
        for j=0:-1+length(gtsol)/i %número de rodadas
            for a=1:i %testes da rodada
                if gtsolk(j*i+a,k)==0
                    S(k,i)=S(k,i)+1;
                break
                end
            end
        end
        TA(k,i)=100*S(k,i)/(length(gtsol)/i);
    end
end

figure
hold on
for i=1:length(nt)
    plot(TA(i,:),'-s');
end
hold off
ylabel('\fontsize{18}Success rate(%)')
xlabel('\fontsize{18}Number of repeated runs')
title('\fontsize{20}Tradeoff of success rate vs. number of iterations per root')
lgd=legend(nt);
title(lgd,'Maximum number of iterations per path')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FEITO INICIALMENTE PARA DETECTAR POSSÍVEIS PROBLEMAS DEGENERADOS
% s=sum(numgtsol==0);
% pNgtsol=zeros(s,1); % Problemas para os quais não acha soluçao verdadeira
% aux=1;
% for j=1:np
%     if numgtsol(j,1)==0
%     pNgtsol(aux,1)=j;
%     aux=aux+1;
%     end 
% end
% 
% writematrix(pNgtsol,'problemas_degenerados.txt');
% 
% addpath(genpath('Dfrpt'));

% A ANÁLISE DOS PROBLEMAS LEVOU A UMA CORREÇÃO NO MINUS, NO TESTE DA ÁREA
% DO TRIÂNGULO EM CADA IMAGEM. EPSILON ALTERADO DE 10-4 PARA 10-10.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% REAL ROOTS

aux1=1;
aux2=1;
for i=1:length(raizrealN)
    if raizrealN(i,1)==313
        idreal(aux1,1)=i-1; %posicao da última raiz real, not gt, em cd teste
        idreal(aux1,2)=aux1; % id do teste, considerando todos os problemas juntos
        idreal(aux1,3)=i-aux2; % # de raízes reais que nao sao solucao em cada teste
        aux1=aux1+1;
        aux2=i+1;
    end
end

% # of real roots not gt solution in each step for each problem
for i=1:np
    for j=1:500
        NrealP(j,i)=idreal(500*(i-1)+j,3);
    end
end


% # of iterations of each real root for each problem
aux=1;
for r=1:500*np % # do teste, considerando todos os problemas juntos
    for s=aux:idreal(r,1) %intervalo das raízes reais do mesmo teste
        raizrealN(s,2)=steps(raizrealN(s,1)+1+312*(idreal(r,2)-1)); % # iteracões na raiz
    end
    aux=idreal(r,1)+2;
end


% Boxplot real roots
figure
boxplot(NrealP, PlotStyle="compact")
m = findobj(gcf, 'Tag', 'Outliers');
set(m, 'color', [0.9 0.9 0.9], 'markerfacecolor',[0.9 0.9 0.9], 'markeredgecolor', [0.7 0.7 0.7], 'markersize',1);
ylabel('\fontsize{18}Number of real roots')
xlabel('\fontsize{18}Problem')
yticks(0:10:80)
title('\fontsize{20}Number of real roots')

%-----------------------------------------------------------------------------------------------%

% ALL REAL ROOTS
%Matrix of all real roots
aux1=idreal(500,1);
aux2=1;
for i=1:np-1
   raizrealP(1:500,2*i-1)=gtsol((i-1)*500+1:500*i,1); % lines with index of gt solution (that are real roots)
   raizrealP(1:500,2*i)=stepsolP(:,i); % # of steps for each gt solution
   raizrealP(501:aux1+500,(i-1)*2+1:2*i)=raizrealN(aux2:idreal(500*i,1),:); % real roots (not gt) and # of steps
   aux1=idreal(500*(i+1),1)-idreal(500*i);
   aux2=idreal(500*i,1)+1;
end
raizrealP(1:500,2*np-1)=gtsol((np-1)*500+1:500*np,1);
raizrealP(1:500,2*np)=stepsolP(:,np);
raizrealP(501:idreal(500*np,1)-idreal(500*(np-1),1)+500,(np-1)*2+1:2*np)=raizrealN(idreal(500*(np-1),1)+1:idreal(500*np,1),:);

for i=2:2:2*np
    m=find(raizrealP(:,i));
    n=m(length(m),1)+1;
    for j=n:length(raizrealP)
    raizrealP(j,i-1:i)=[NaN NaN];
    end
end
for i=1:length(raizrealP)
    for j=2:2:2*np
        if raizrealP(i,j)==0
            raizrealP(i,j-1:j)=[NaN NaN];
        end
    end
end
%stepreal: matrix of # iterations foreach real root for all problems
stepreal=raizrealN(:,1);
stepreal(length(raizrealN)+1:length(raizrealN)+70000,1)=stepsol(:,1);

MeR=median(stepreal, "omitnan");
ymR=MeR'+0*xm;
figure
boxplot(raizrealP(:,2:2:2*np), 'Symbol', '+k', PlotStyle='compact')
m = findobj(gcf, 'Tag', 'Outliers');
set(m, 'color', [0.9 0.9 0.9], 'markerfacecolor',[0.9 0.9 0.9], 'markeredgecolor', [0.7 0.7 0.7], 'markersize',1);
hold on
plot(xm,ymR,'Color','g', 'LineWidth',2);
hold off
ylabel('\fontsize{18}Number of iterations per root')
xlabel('\fontsize{18}Problem')
title('\fontsize{20}Time spent in paths leading to all real roots')

%----------------------------------------------------------------------------------------------------------%
