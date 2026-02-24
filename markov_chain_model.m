%probabilities of transitions  between states
% see output table Tp: Tp(i,j) - transition from state i to state j

% select healers or non-healers in line 24-25
% select main or validation cohort in lines 8-9

clear
Clusters_analysis
C1=ClustVal(2,:);
C3=ClustVal(1,:);
cl=h;

for i=1:numel(p)
    if C1(i)>bnd1      % vertical boundary in cluster plane  
        ST(i)=2;
    else if C3(i)<bnd2  % horizontal boundary in cluster plane 
            ST(i)=1;
        else
            ST(i)=3;
        end
    end
end
clear i C1 C3
%% Number and fraction of healers and non-healers in each state
ST0=ST(cl==0); n0=numel(ST0);
ST1=ST(cl==1); n1=numel(ST1);
for i=1:3
   STT(i,1)=numel(find(ST0==i));
   STT(i,2)=numel(find(ST1==i));
   STP(i,1)=STT(i,1)/n0;
   STP(i,2)=STT(i,2)/n1;
end
%%
%select healers only / non-healers only
% s=0 for non healers;   s=1 for healers
%%
f5=figure;
fsz=7;
for s=[0 1]
%s=0;
ws=w(cl==s);
ps=p(cl==s);
STs=ST(cl==s);

T=zeros(3,3); %matrix of transitions 1st index - from, second index - to

psu=unique(ps);
for i=1:numel(psu)
    %select samples from one patient
    k=find(ismember(ps,psu(i)));
    STu=STs(k);
    wsu=ws(k);
    %rearrange samples in order of time of collection
    [I,i1]=mink(wsu,numel(wsu));
    STu=STu(i1); wsu=wsu(i1);
    for j=1:numel(wsu)-1
        if wsu(j+1)-wsu(j)==1
            T(STu(j),STu(j+1))=T(STu(j),STu(j+1))+1;
        end
    end   
end

tmp=sum(T');
Tp=T./tmp';

subplot(1,2,2-s)
circle(0,0,1,[0.9647    0.7765    0.9529]);
circle(3,0,1,[1 1 1]);
circle(1.5,3,1,[0.9686    0.9961    0.6863]);
xlim([-2 5])
ylim([-2 5])
text(-0.7, 0, 'Impairment', 'FontSize', fsz);
text(2.3, 0, 'Inflammation', 'FontSize', fsz);
text(0.8, 3, 'Proliferation', 'FontSize', fsz);
hold on, quiver(1,0.2,1.1,0,'LineWidth',1,'Color','black')
hold on, quiver(2,-0.2,-1.1,0,'LineWidth',1,'Color','black')
hold on, quiver(0.6,1,0.5,1,'LineWidth',1,'Color','black')
hold on, quiver(0.8,2,-0.5,-1,'LineWidth',1,'Color','black')
hold on, quiver(2.4,1,-0.5,1,'LineWidth',1,'Color','black')
hold on, quiver(2.2,2,0.5,-1,'LineWidth',1,'Color','black')

text(-1, -1,string(Tp(1,1)), 'FontSize', fsz);
text(0.8, 1.1,string(Tp(1,2)), 'FontSize', fsz);
text(1.3, 0.4,string(Tp(1,3)), 'FontSize', fsz);

text(0, 1.6,string(round(Tp(2,1),2)), 'FontSize', fsz);
text(1.2, 4.2,string(round(Tp(2,2),2)), 'FontSize', fsz);
text(2.6, 1.6,string(round(Tp(2,3),2)), 'FontSize', fsz);

text(1.3, -0.5,string(round(Tp(3,1),2)), 'FontSize', fsz);
text(1.7, 1.3,string(round(Tp(3,2),2)), 'FontSize', fsz);
text(3.6, -1,string(round(Tp(3,3),2)), 'FontSize', fsz);

xticklabels([]);
yticklabels([]);

if s==1
    title('Healers')
    %picturefilename='Healers.pdf';
else
    title('Non-healers')
    %picturefilename='Nonealers.pdf';
end

axis square
ax = gca;
ax.YColor = 'w';
ax.XColor = 'w';
end
exportgraphics(f5,'markov_chain.pdf','ContentType','vector');

function h = circle(x,y,r,color)
hold on
th = 0:pi/50:2*pi;
xunit = r * cos(th) + x;
yunit = r * sin(th) + y;
h = plot(xunit, yunit,'k');
fill(xunit,yunit,color)
hold off
end