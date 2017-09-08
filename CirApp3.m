clear all; clc;close all;
%%
%DIMITROU VASILEIOS 5380
My=5000000;                           %Nmm Bending Moment
f=100000;                               %axial force
numellipse=8;                          %data 1 - No circles/ellipses
step=pi/40;                            %data 2 - angle
d=26.6;                                %data 3 - b\o std ID
%%
%-------------------first part---------------------------------------------
%---------------------meshing----------------------------------------------
%reg1
h1=50; xbound=50; y1=d/2; h=47.5; maxcoeffellipse=h/y1;
h2=85; ybound=50;  check=0; check2=0; thcr=pi/4; scale2=1/100;
yin=linspace(1,maxcoeffellipse,numellipse)*y1; scale1=1/100;
i=1;  [c r]=size(yin); th=0:step:pi/2;Elem1=[]; l=0; ll=[]; 
for theta=th
    phi=pi/2-theta;
    b=yin(:);
    x(:,i)=b*cos(phi);
    y(:,i)=b*sin(phi);
    i=i+1;
end
[cn1,rn1]=size(x);xx=x(1:(cn1-1),:); yy=y(1:(cn1-1),:); 
NodeID=[[1:length(reshape(xx,1,(cn1-1)*(rn1))')]' reshape(xx,1,(cn1-1)*(rn1))' reshape(yy,1,(cn1-1)*(rn1))'];
[cx rx]=size(x); itemp=1; 
x2=[];y2=[]; nth=[]; j=1;  nodout=[];Elem2=[];
%παρένθεση περιοχής 2
for i=2:length(th)-1
    difx=x(cx,i)-x(cx-1,i);
    dify=y(cx,i)-y(cx-1,i);
    
    if th(i)==thcr
        x2=[x2 x(cx,i):difx:xbound  xbound];
        y2=[y2 y(cx,i):dify:ybound  ybound];
        nth(i)=length(x2);         holdthi=i;
        check2=1;
        if y2(length(x2))-y2(length(x2)-1)<0.5
          saveempos(j)=length(x2)-1;             %hold empty posiition
          j=j+1;
          nodout=[nodout i];
          
        end
    end 
     
     if th(i)<thcr
        x2=[x2 x(cx,i):difx:ybound*tan(th(i))  ybound*tan(th(i))];
        y2=[y2 y(cx,i):dify:ybound  ybound];
        nth(i)=length(x2);
        if y2(length(x2))-y2(length(x2)-1)<1
          saveempos(j)=length(x2)-1;             %η προηγούμενη θέση θέλω να απαλείφεται 
          j=j+1;
          nodout=[nodout i];  holdthi2=i;
        end
     end
     
     if th(i)>thcr
        x2=[x2 x(cx,i):difx:xbound  xbound];
        y2=[y2 y(cx,i):dify:xbound*tan((pi/2)-th(i))  xbound*tan((pi/2)-th(i))];
        nth(i)=length(x2);
        if ( x2(length(x2))-x2(length(x2)-1)<0.7) %&& th(i)<(pi/2)*0.9 y2(length(x2))-y2(length(x2)-1)<1 %||
          saveempos(j)=length(x2)-1;
          j=j+1;
          nodout=[nodout i];  
        end
      end
end
for ii=1:length(nodout)
    nth=[nth(1:nodout(ii)-1) [nth(nodout(ii):length(nth))]-1 ] ;  
    %αφαιρώ για κάθε ακτίνα και αποτις επόμενες της τον αριθμό των αφαιρούμενων κόμβων
end
for j=1:length(th)-1
 for i=1:cx-1
     if j==1 && i~=cx-1
Elem1=[Elem1;i i+(cx-1) i+1;i+(cx-1) i+1+(cx-1) i+1 ];     
     end
     if j~=1 && i~=cx-1
Elem1=[Elem1;i+(cx-1)*(j-1) i+(cx-1)*j i+1+(cx-1)*(j-1) ; i+(cx-1)*j i+1+(cx-1)*j i+1+(cx-1)*(j-1)  ];
     end
     if j==1 && i==cx-1
Elem1=[Elem1];       
     end  
     if j~=1 && i==cx-1
Elem1=[Elem1;i+(cx-1)*(j-1) i+(cx-1)*j ((cx-1)*rx)+1+nth(j-1) ; i+(cx-1)*j ((cx-1)*rx)+1+nth(j) ((cx-1)*rx)+1+nth(j-1)  ];
     end
xoddelemreg1(itemp,:)=[x(i,j) x(i,j+1) x(i+1,j)];
xevenelemreg1(itemp,:)=[ x(i,j+1) x(i+1,j+1) x(i+1,j) ];
yoddelemreg1(itemp,:)=[y(i,j) y(i,j+1) y(i+1,j)];
yevenelemreg1(itemp,:)=[ y(i,j+1) y(i+1,j+1) y(i+1,j) ]; 
itemp=itemp+1;
 end
end
%%
%reg2
x2=[];y2=[]; nth=[]; j=1;  nodout=[];Elem2=[];
%node determination
for i=2:length(th)-1
    difx=x(cx,i)-x(cx-1,i);
    dify=y(cx,i)-y(cx-1,i);
    
    if th(i)==thcr
        x2=[x2 x(cx,i):difx:xbound  xbound];
        y2=[y2 y(cx,i):dify:ybound  ybound];
        nth(i)=length(x2);         holdthi=i;
        check2=1;
        if y2(length(x2))-y2(length(x2)-1)<0.5
          saveempos(j)=length(x2)-1;             %hold empty posiition
          j=j+1;
          nodout=[nodout i];
          
        end
    end 
     
     if th(i)<thcr
        x2=[x2 x(cx,i):difx:ybound*tan(th(i))  ybound*tan(th(i))];
        y2=[y2 y(cx,i):dify:ybound  ybound];
        nth(i)=length(x2);
        if y2(length(x2))-y2(length(x2)-1)<1
          saveempos(j)=length(x2)-1;             %η προηγούμενη θέση θέλω να απαλείφεται 
          j=j+1;
          nodout=[nodout i];  holdthi2=i;
        end
     end
     
     if th(i)>thcr
        x2=[x2 x(cx,i):difx:xbound  xbound];
        y2=[y2 y(cx,i):dify:xbound*tan((pi/2)-th(i))  xbound*tan((pi/2)-th(i))];
        nth(i)=length(x2);
        if ( x2(length(x2))-x2(length(x2)-1)<0.7) %&& th(i)<(pi/2)*0.9 y2(length(x2))-y2(length(x2)-1)<1 %||
          saveempos(j)=length(x2)-1;
          j=j+1;
          nodout=[nodout i];  
        end
      end
end
for ii=1:length(nodout)
    nth=[nth(1:nodout(ii)-1) [nth(nodout(ii):length(nth))]-1 ] ;  
    %αφαιρώ για κάθε ακτίνα και αποτις επόμενες της τον αριθμό των αφαιρούμενων κόμβων
end
x2(saveempos)=[];  y2(saveempos)=[];          %empty positions
x2me=x2;y2me=y2; diff=100;diff1=200;
for i=2:length(th)-1
    diff=sqrt((x2(nth(i))-85)^2+(y2(nth(i))-50)^2);
    if diff1>diff
       diff1=diff; chpos=i;
    end
end
x2(nth(chpos))=xbound; y2(nth(chpos))=ybound;
%plot(x2,y2,'ko');
x2me=x2;y2me=y2;
%plot(x2me,y2me,'go');
NodeID=[NodeID ;[(cn1-1)*(rn1)+(1:nth(length(nth)))]' x2me((1:length(x2me)))' y2me((1:length(x2me)))'];
itemp=0;  ie=1; j=2; sth=0;
%element determination
for j=2:length(th)-2
    if th(j)<thcr

            if nth(j)-nth(j-1)==nth(j+1)-nth(j)
            for ii=1:(nth(j)-nth(j-1))-1
              i=ii+nth(j-1);
              Elem2=[Elem2; (cx-1)*rx+i (cx-1)*rx+i+(nth(j)-nth(j-1)) (cx-1)*rx+i+1 ; (cx-1)*rx+i+(nth(j)-nth(j-1))  (cx-1)*rx+i+(nth(j)-nth(j-1))+1  (cx-1)*rx+i+1  ];
              xoddelemreg2(ie,:)=[x2(i) x2(i+(nth(j)-nth(j-1))) x2(i+1)];
              xevenelemreg2(ie,:)=[x2(i+(nth(j)-nth(j-1))) x2(i+(nth(j)-nth(j-1))+1) x2(i+1) ];
              yoddelemreg2(ie,:)=[y2(i) y2(i+(nth(j)-nth(j-1))) y2(i+1)];
              yevenelemreg2(ie,:)=[y2(i+(nth(j)-nth(j-1))) y2(i+(nth(j)-nth(j-1))+1) y2(i+1) ]; 
              ie=ie+1;
            end 

            end
            if nth(j)-nth(j-1)<nth(j+1)-nth(j)
            for ii=1:(nth(j)-nth(j-1))-1
              i=ii+nth(j-1);
             Elem2=[Elem2; (cx-1)*rx+i (cx-1)*rx+i+(nth(j)-nth(j-1)) (cx-1)*rx+i+1 ; (cx-1)*rx+i+(nth(j)-nth(j-1))  (cx-1)*rx+i+(nth(j)-nth(j-1))+1  (cx-1)*rx+i+1  ];
              xoddelemreg2(ie,:)=[x2(i) x2(i+(nth(j)-nth(j-1))) x2(i+1)];
              xevenelemreg2(ie,:)=[x2(i+(nth(j)-nth(j-1))) x2(i+(nth(j)-nth(j-1))+1) x2(i+1) ];
              yoddelemreg2(ie,:)=[y2(i) y2(i+(nth(j)-nth(j-1))) y2(i+1)];
              yevenelemreg2(ie,:)=[y2(i+(nth(j)-nth(j-1))) y2(i+(nth(j)-nth(j-1))+1) y2(i+1) ]; 
              ie=ie+1;
            end
              if nth(j+1)-nth(j)-nth(j)+nth(j-1)==1
              xevenelemreg2(ie,:)=[x2(nth(j+1)-1) x2(nth(j+1)) x2(nth(j))];
              yevenelemreg2(ie,:)=[y2(nth(j+1)-1) y2(nth(j+1)) y2(nth(j))];
              ie=ie+1;
              Elem2=[Elem2;(cx-1)*rx+nth(j+1)-1   (cx-1)*rx+nth(j+1) (cx-1)*rx+nth(j) ];
              end
                if nth(j+1)-nth(j)-nth(j)+nth(j-1)>1
                         for a=1:nth(j+1)-nth(j)-nth(j)+nth(j-1)
                             xevenelemreg2(ie,:)=[ x2(nth(j)+nth(j)-nth(j-1)+a-1)  x2(nth(j)+nth(j)-nth(j-1)+a) x2(nth(j))];
                             yevenelemreg2(ie,:)=[ y2(nth(j)+nth(j)-nth(j-1)+a-1) y2(nth(j)+nth(j)-nth(j-1)+a) y2(nth(j))];
                             ie=ie+1; 
                             Elem2=[Elem2;((cn1)*(rn1))+n(j)+n(j)-n(j-1)+a-1   ((cn1)*(rn1))+n(j)+n(j)-n(j-1)+a ((cn1)*(rn1))+n(j)];
                         end
                end 
             end
     end

             
            
      %μέχρι εδώ για  διάστημα θ=[0,45]      
            
            
            
            
            
      %απο εδώ και κάτω για διάστημα θ=[45,90]     
            
           if th(j)==thcr 
           if nth(j)-nth(j-1)==nth(j+1)-nth(j)
            for ii=1:(nth(j)-nth(j-1))-1
              i=ii+nth(j-1);
              Elem2=[Elem2; (cx-1)*rx+i (cx-1)*rx+i+(nth(j)-nth(j-1)) (cx-1)*rx+i+1 ; (cx-1)*rx+i+(nth(j)-nth(j-1))  (cx-1)*rx+i+(nth(j)-nth(j-1))+1  (cx-1)*rx+i+1  ];
              xoddelemreg2(ie,:)=[x2(i) x2(i+(nth(j)-nth(j-1))) x2(i+1)];
              xevenelemreg2(ie,:)=[x2(i+(nth(j)-nth(j-1))) x2(i+(nth(j)-nth(j-1))+1) x2(i+1) ];
              yoddelemreg2(ie,:)=[y2(i) y2(i+(nth(j)-nth(j-1))) y2(i+1)];
              yevenelemreg2(ie,:)=[y2(i+(nth(j)-nth(j-1))) y2(i+(nth(j)-nth(j-1))+1) y2(i+1) ]; 
              ie=ie+1;
            end 
          end
       else if nth(j)-nth(j-1)>nth(j+1)-nth(j)
            for ii=1:(nth(j+1)-nth(j))-1
              i=ii+nth(j-1);
             Elem2=[Elem2; (cx-1)*rx+i (cx-1)*rx+i+(nth(j)-nth(j-1)) (cx-1)*rx+i+1 ; (cx-1)*rx+i+(nth(j)-nth(j-1))  (cx-1)*rx+i+(nth(j)-nth(j-1))+1  (cx-1)*rx+i+1  ];
              xoddelemreg2(ie,:)=[x2(i) x2(i+(nth(j)-nth(j-1))) x2(i+1)];
              xevenelemreg2(ie,:)=[x2(i+(nth(j)-nth(j-1))) x2(i+(nth(j)-nth(j-1))+1) x2(i+1) ];
              yoddelemreg2(ie,:)=[y2(i) y2(i+(nth(j)-nth(j-1))) y2(i+1)];
              yevenelemreg2(ie,:)=[y2(i+(nth(j)-nth(j-1))) y2(i+(nth(j)-nth(j-1))+1) y2(i+1) ]; 
              ie=ie+1;
            end
              if nth(j)-nth(j-1)-nth(j+1)+nth(j)==1
              xevenelemreg2(ie,:)=[x2(nth(j+1)-1) x2(nth(j+1)) x2(nth(j))];
              yevenelemreg2(ie,:)=[y2(nth(j+1)-1) y2(nth(j+1)) y2(nth(j))];
              ie=ie+1;
              Elem2=[Elem2;(cx-1)*rx+nth(j+1)-1   (cx-1)*rx+nth(j+1) (cx-1)*rx+nth(j) ];
              else if nth(j+1)-nth(j)-nth(j)+nth(j-1)<1
                         for a=1:nth(j+1)-nth(j)-nth(j)+nth(j-1)
                             xevenelemreg2(ie,:)=[ x2(nth(j-1)+nth(j+1)-nth(j)+a) x2(nth(j-1)+nth(j+1)-nth(j)+a-1) x2(nth(j+1))];
                             yevenelemreg2(ie,:)=[ y2(nth(j-1)+nth(j+1)-nth(j)+a) y2(nth(j-1)+nth(j+1)-nth(j)+a-1) y2(nth(j+1))];
                             ie=ie+1; 
                             Elem2=[Elem2;((cn1)*(rn1))+n(j-1)+n(j+1)-n(j)+a   ((cn1)*(rn1))+n(j-1)+n(j+1)-n(j)+a-1 ((cn1)*(rn1))+n(j+1)];
                         end
                   end 
              end
           end
           end
       
     
           
          if th(j)>thcr
          if nth(j)-nth(j-1)==nth(j+1)-nth(j)
            for ii=1:(nth(j)-nth(j-1))-1
              i=ii+nth(j-1);
              Elem2=[Elem2; (cx-1)*rx+i (cx-1)*rx+i+(nth(j)-nth(j-1)) (cx-1)*rx+i+1 ; (cx-1)*rx+i+(nth(j)-nth(j-1))  (cx-1)*rx+i+(nth(j)-nth(j-1))+1  (cx-1)*rx+i+1  ];
              xoddelemreg2(ie,:)=[x2(i) x2(i+(nth(j)-nth(j-1))) x2(i+1)];
              xevenelemreg2(ie,:)=[x2(i+(nth(j)-nth(j-1))) x2(i+(nth(j)-nth(j-1))+1) x2(i+1) ];
              yoddelemreg2(ie,:)=[y2(i) y2(i+(nth(j)-nth(j-1))) y2(i+1)];
              yevenelemreg2(ie,:)=[y2(i+(nth(j)-nth(j-1))) y2(i+(nth(j)-nth(j-1))+1) y2(i+1) ]; 
              ie=ie+1;
            end 
          end
       else if nth(j)-nth(j-1)>nth(j+1)-nth(j)
            for ii=1:(nth(j+1)-nth(j))-1
              i=ii+nth(j-1);
              Elem2=[Elem2; (cx-1)*rx+i (cx-1)*rx+i+(nth(j)-nth(j-1)) (cx-1)*rx+i+1 ; (cx-1)*rx+i+(nth(j)-nth(j-1))  (cx-1)*rx+i+(nth(j)-nth(j-1))+1  (cx-1)*rx+i+1  ];
              xoddelemreg2(ie,:)=[x2(i) x2(i+(nth(j)-nth(j-1))) x2(i+1)];
              xevenelemreg2(ie,:)=[x2(i+(nth(j)-nth(j-1))) x2(i+(nth(j)-nth(j-1))+1) x2(i+1) ];
              yoddelemreg2(ie,:)=[y2(i) y2(i+(nth(j)-nth(j-1))) y2(i+1)];
              yevenelemreg2(ie,:)=[y2(i+(nth(j)-nth(j-1))) y2(i+(nth(j)-nth(j-1))+1) y2(i+1) ]; 
              ie=ie+1;
            end
              if nth(j)-nth(j-1)-nth(j+1)+nth(j)==1
              xevenelemreg2(ie,:)=[x2(nth(j+1)-1) x2(nth(j+1)) x2(nth(j))];
              yevenelemreg2(ie,:)=[y2(nth(j+1)-1) y2(nth(j+1)) y2(nth(j))];
              ie=ie+1;
              Elem2=[Elem2;(cx-1)*rx+nth(j+1)-1   (cx-1)*rx+nth(j+1) (cx-1)*rx+nth(j) ];
              else if nth(j+1)-nth(j)-nth(j)+nth(j-1)<1
                         for a=1:nth(j+1)-nth(j)-nth(j)+nth(j-1)
                             xevenelemreg2(ie,:)=[ x2(nth(j-1)+nth(j+1)-nth(j)+a) x2(nth(j-1)+nth(j+1)-nth(j)+a-1) x2(nth(j+1))];
                             yevenelemreg2(ie,:)=[ y2(nth(j-1)+nth(j+1)-nth(j)+a) y2(nth(j-1)+nth(j+1)-nth(j)+a-1) y2(nth(j+1))];
                             ie=ie+1; 
                             Elem2=[Elem2;((cn1)*(rn1))+n(j-1)+n(j+1)-n(j)+a   ((cn1)*(rn1))+n(j-1)+n(j+1)-n(j)+a-1 ((cn1)*(rn1))+n(j+1)];
                         end
                   end 
              end
           end
           end
end
for i=1:length(x2) %ελεγχος για το σημείο στην κορυφή αν έχει ορισθεί 
    if x2(i)==xbound && y2(i)==ybound
       check=1;
    end   
end
%ορίζω τα πρώτα 2 και τα τελευταία 2
x22=[ x(cx,1) 0 ]; y22=[y(cx,1) ybound];
%plot(x22,y22,'ko')          %για θ = 0
x22e=[ x(cx,length(th)) xbound ]; y22e=[y(cx,length(th)) 0];
%plot(x22e,y22e,'ko')      %για θ = 90
xoddelemreg2=[x22(1) x2(1) x22(2);xoddelemreg2; x2(length(x2)-1) x22e(1) x2(length(x2))];
yoddelemreg2=[y22(1) y2(1) y22(2);yoddelemreg2; y2(length(x2)-1) y22e(1) y2(length(x2))];
xevenelemreg2=[x22(2) x2(1) x2(2);xevenelemreg2 ; x22e(1) x22e(2) x2(length(x2))];
yevenelemreg2=[y22(2) y2(1) y2(2);yevenelemreg2 ; y22e(1) y22e(2) y2(length(x2))];  %όριζω τα πρώτα 2 στοιχέια και τα 2 τελευταία
[ieo2,ro2]=size(xevenelemreg2);[iee2,re2]=size(xoddelemreg2); [cn2,rn2]=size([x22 x2 x22e]);

XY=[NodeID(Elem2(:,1),[2 3]) NodeID(Elem2(:,2),[2 3]) NodeID(Elem2(:,3),[2 3])];
XY3=NodeID(Elem2(:,3),[2 3]); XY2=NodeID(Elem2(:,2),[2 3]); XY1=NodeID(Elem2(:,1),[2 3]);

Elem1=Elem1(1:length(Elem1)-1,:);
Elem=[Elem1 ; Elem2];
NodeID=[NodeID ; NodeID(length(NodeID),1)+1 0 h; NodeID(length(NodeID),1)+2 0 ybound ; NodeID(length(NodeID),1)+3  h 0 ;NodeID(length(NodeID),1)+4  ybound 0 ];

Elem=[Elem ; [  NodeID(length(NodeID),1)-3  ((cx-1)*rx)+1  NodeID(length(NodeID),1)-2 ];...
    [((cx-1)*rx)+1  ((cx-1)*rx)+2 NodeID(length(NodeID),1)-2 ];...
    [NodeID(length(NodeID),1)-5  NodeID(length(NodeID),1)-1  NodeID(length(NodeID),1)-4];...
    [NodeID(length(NodeID),1)-1  NodeID(length(NodeID),1)  NodeID(length(NodeID),1)-4 ];...
    [cx-1 2*(cx-1)  NodeID(length(NodeID),1)-3];[2*(cx-1)  ((cx-1)*rx)+1  NodeID(length(NodeID),1)-3];...
    [((cx-1)*rx) NodeID(length(NodeID),1)-1 NodeID(length(NodeID),1)-4]]; %οκ τα στοιχεία


%%
%3rd region

%nodes
j=1; 
for i=holdthi:length(th)-1 
y3_temp(j)=y2(nth(i));
j=j+1;
end
y3_temp=[y3_temp'; 0]; y3=[]; x3=[]; c=10;
for i=1:length(linspace(50,85,c))
y3=[y3 y3_temp];
end
for i=1:length(y3_temp)
x3=[x3;linspace(50,85,c)];
end

[c r]=size(x3);numnode=c*r; Elem3=[];x33=[];y33=[];
x33=x3(:,2:r);y33=y3(:,2:r);
%elements
ie=1; Elem3=[]; Tempval=find(NodeID(:,2)==xbound);  NodeID3=[max(max(Elem))+(1:length(reshape(x33',1,c*r-c))) ; reshape(x33',1,c*r-c) ; reshape(y33',1,c*r-c)]';
for j=1:c-1
     for i=1:r-1 
     if i==1
Elem3=[Elem3;Tempval(j+1,1) i+(j-1)*(r-1)+max(max(Elem)) Tempval(j,1) ;...
    Tempval(j+1,1)  i+(j)*(r-1)+max(max(Elem)) i+(j-1)*(r-1)+max(max(Elem))  ];
     end
     if i~=1
Elem3=[Elem3; i-1+(r-1)*j+max(max(Elem))  i+(j-1)*(r-1)+max(max(Elem))  (i-1)+(j-1)*(r-1)+max(max(Elem))  ;...
      (i-1)+(r-1)*j+max(max(Elem))     i+(r-1)*j+max(max(Elem))  i+(j-1)*(r-1)+max(max(Elem)) ];  
     end
     end
end
%αυτά με ενδιαφέρουν απο την παραπάνω διαδικασία. είναι οι κομβοι και τα
%στοιχεία ενός τεραρτημορίου 
x2h=x2;y2h=y2;
if check==0
   for i=1:length(th)-1
       if x2h(nth(i))<xbound
       holddati=i;
       end
   end
x2=[x2 85]; y2=[y2 50];
xoddelemreg2=[xoddelemreg2 ; x2h(nth(holddati)) x2h(nth(holddati+1)) xbound];
yoddelemreg2=[yoddelemreg2 ; y2h(nth(holddati)) y2h(nth(holddati+1)) ybound];
end


Elem=[Elem; Elem3]; NodeID=[NodeID;NodeID3];

XYtet=[NodeID(Elem(:,1),[2 3]) NodeID(Elem(:,2),[2 3]) NodeID(Elem(:,3),[2 3])];


%figure(1000)
%for i=1:length(XYtet)
%    plot([XYtet(i,1) XYtet(i,3) XYtet(i,5)],[XYtet(i,2),XYtet(i,4),XYtet(i,6)],'b');hold on
%end

%%
%Έτοιμα τα στοιχεια, τώρα θα κάνω mirror, μέχρι εδώ σωστά
%Στοιχεία για ολόκληρη την πλάκα

NodeIDmirr=NodeID([find(NodeID(:,3))],:); Noholdval=NodeID([find(~NodeID(:,3))],:);%NodeIDmirr osa den einai 0
NodeIDtot=[NodeID;NodeIDmirr(:,1)+NodeID(length(NodeID),1) NodeIDmirr(:,2) -NodeIDmirr(:,3)];
Elemm=Elem+NodeID(length(NodeID),1);
tofind=Noholdval(:,1)+NodeID(length(NodeID),1);        %τους κομβους για τους οποίους έχουμε 0 , άρα είναι κοινοί
savefind=[]; len=[];
for i=1:length(tofind)
    for j=1:3
    savefind1=find(Elemm(:,j)==tofind(i));
    Elemm(savefind1,j)=Elemm(savefind1,j)-NodeID(length(NodeID),1);  
    end
end

Elemtot=[Elem;Elemm(:,3) Elemm(:,2) Elemm(:,1)];
%element determination
for i=1:length(Elemtot(:,1))
XYtot(i,[1 2])=[NodeIDtot(find(NodeIDtot(:,1)==Elemtot(i,1)),[2 3])];
XYtot(i,[3 4])=[NodeIDtot(find(NodeIDtot(:,1)==Elemtot(i,2)),[2 3])];
XYtot(i,[5 6])=[NodeIDtot(find(NodeIDtot(:,1)==Elemtot(i,3)),[2 3])];
end


for i=1:length(NodeID)
if NodeID(i,2)<0.001
NodeID(i,2)=0;
end
end

for i=1:length(NodeIDtot)
if NodeIDtot(i,2)<0.001
NodeIDtot(i,2)=0;
end
end

NodeIDmirr2=NodeIDtot([find(NodeIDtot(:,2))],:); Noholdval2=NodeIDtot([find(~NodeIDtot(:,2))],:);
NodeIDtot2=[NodeIDtot;NodeIDmirr2(:,1)+NodeIDtot(length(NodeIDtot),1) -NodeIDmirr2(:,2) NodeIDmirr2(:,3)];
Elemtot2=Elemtot+NodeIDtot(length(NodeIDtot),1);
tofind=Noholdval2(:,1)+NodeIDtot(length(NodeIDtot),1);        %τους κομβους για τους οποίους έχουμε 0 , άρα είναι κοινοί
savefind=[]; len=[];
for i=1:length(tofind)
    for j=1:3
    savefind1=find(Elemtot2(:,j)==tofind(i));
    Elemtot2(savefind1,j)=Elemtot2(savefind1,j)-NodeIDtot(length(NodeIDtot),1);
    end
end
Elemtot=[Elemtot; Elemtot2(:,3) Elemtot2(:,2) Elemtot2(:,1)];
XYtotf=zeros(length(Elemtot(:,1)),6);
%element determination
for i=1:length(Elemtot(:,1))
XYtotf(i,[1 2])=[NodeIDtot2(find(NodeIDtot2(:,1)==Elemtot(i,1)),[2 3])];
XYtotf(i,[3 4])=[NodeIDtot2(find(NodeIDtot2(:,1)==Elemtot(i,2)),[2 3])];
XYtotf(i,[5 6])=[NodeIDtot2(find(NodeIDtot2(:,1)==Elemtot(i,3)),[2 3])];
end


Nodes=NodeIDtot2;Elements=Elemtot;
v=0.3; E=210000; t=4;Ktot=zeros(length(Nodes)*2);
D=E/(1-v^2)*[1 v 0;v 1 0;0 0 (1-v)/2]; As=[];
for i=1:length(Elements)
    b= [ Nodes(find(Nodes(:,1)==Elements(i,2)),3)-Nodes(find(Nodes(:,1)==Elements(i,3)),3);...
        Nodes(find(Nodes(:,1)==Elements(i,3)),3)-Nodes(find(Nodes(:,1)==Elements(i,1)),3);...
        Nodes(find(Nodes(:,1)==Elements(i,1)),3)-Nodes(find(Nodes(:,1)==Elements(i,2)),3)];
    c=[Nodes(find(Nodes(:,1)==Elements(i,3)),2)-Nodes(find(Nodes(:,1)==Elements(i,2)),2);...
        Nodes(find(Nodes(:,1)==Elements(i,1)),2)-Nodes(find(Nodes(:,1)==Elements(i,3)),2);...
        Nodes(find(Nodes(:,1)==Elements(i,2)),2)-Nodes(find(Nodes(:,1)==Elements(i,1)),2)];
   
    A=1/2*det([1 Nodes(find(Nodes(:,1)==Elements(i,1)),2)  Nodes(find(Nodes(:,1)==Elements(i,1)),3);...
        1 Nodes(find(Nodes(:,1)==Elements(i,2)),2) Nodes(find(Nodes(:,1)==Elements(i,2)),3);...
        1 Nodes(find(Nodes(:,1)==Elements(i,3)),2) Nodes(find(Nodes(:,1)==Elements(i,3)),3)]);
    B=(1/(2*A))*[b(1,1) 0 b(2,1) 0 b(3,1) 0;0 c(1,1) 0 c(2,1) 0 c(3,1);c(1,1) b(1,1) c(2,1) b(2,1) c(3,1) b(3,1)];
    ke=(B'*D*B)*t*A;
    DOF_temp=[2*find(Nodes(:,1)==Elements(i,1))-1,2*find(Nodes(:,1)==Elements(i,1)),...
        2*find(Nodes(:,1)==Elements(i,2))-1,2*find(Nodes(:,1)==Elements(i,2)),...
        2*find(Nodes(:,1)==Elements(i,3))-1,2*find(Nodes(:,1)==Elements(i,3))];
    Ktot(DOF_temp,DOF_temp)=Ktot(DOF_temp,DOF_temp)+ke; As=[As;A];
end




%%
%--------------------Static Load Modeling----------------------------------

%determine the boundary condition - stable-movable
HINodes=find(Nodes(:,2)==0);HINodes(find(Nodes(find(Nodes(:,2)==0),3)==-(d/2)));
HoldBoundaryNodes=[HINodes(find(Nodes(find(Nodes(:,2)==0),3)==(d/2))),HINodes(find(Nodes(find(Nodes(:,2)==0),3)==-(d/2)))]; 
HoldRestiDOFs=[2*HoldBoundaryNodes(1)-1 ; 2*HoldBoundaryNodes(1) ; 2*HoldBoundaryNodes(2)-1];
actDOFs=(1:2*length(Nodes))'; actDOFs(HoldRestiDOFs)=[];
K=Ktot; maxd=(max(diag(Ktot)))*10^11;



%determine the Loaded condition

HoldLoadedNodes=find(Nodes(:,2)==85);HoldMeanDistperElem=zeros(1,length(HoldLoadedNodes)-1);
for i=1:(length(HoldLoadedNodes)-1)/2
    HoldMeanDistperElem(i)=(Nodes(HoldLoadedNodes(i),3)+Nodes(HoldLoadedNodes(i+1),3))/2 ;
end
for i=1:(length(HoldLoadedNodes)-1)/2
    HoldMeanDistperElem(i+(length(HoldLoadedNodes)-1)/2)=-HoldMeanDistperElem(1+(length(HoldLoadedNodes)-1)/2-i);
end
alph=(My)*(1/(sum(HoldMeanDistperElem.^2))); Fbtw=zeros(length(HoldMeanDistperElem),1);Fh=zeros(length(HoldLoadedNodes),1);
Fbtw=(alph*HoldMeanDistperElem)';
LNodes=2*HoldLoadedNodes;
for i=1:length(HoldLoadedNodes)
  if i==1
      Fh(i)=Fbtw(i)/2;
  end
  if i==length(HoldLoadedNodes)
      Fh(i)=Fbtw(i-1)/2;
  end
  if i~=1 && i~=length(HoldLoadedNodes)
      Fh(i)=(Fbtw(i)/2)+(Fbtw(i-1)/2);
  end
end
for i=1:(length(Fh)-1)/2
   Fh(i+1+(length(HoldLoadedNodes)-1)/2)=-Fh(i);
end
M=Nodes(HoldLoadedNodes,3)'*Fh;
F=zeros(2*length(Nodes),1);F(2*HoldLoadedNodes-1,1)=Fh;

HoldLoadedNodes=find(Nodes(:,2)==-85);HoldMeanDistperElem=zeros(1,length(HoldLoadedNodes)-1);
for i=1:(length(HoldLoadedNodes)-1)/2
    HoldMeanDistperElem(i)=(Nodes(HoldLoadedNodes(i),3)+Nodes(HoldLoadedNodes(i+1),3))/2 ;
end
for i=1:(length(HoldLoadedNodes)-1)/2
    HoldMeanDistperElem(i+(length(HoldLoadedNodes)-1)/2)=-HoldMeanDistperElem(1+(length(HoldLoadedNodes)-1)/2-i);
end
alph=(My)*(1/(sum(HoldMeanDistperElem.^2))); Fbtw=zeros(length(HoldMeanDistperElem),1);Fh=zeros(length(HoldLoadedNodes),1);
Fbtw=-(alph*HoldMeanDistperElem)';
LNodes=2*HoldLoadedNodes;
for i=1:length(HoldLoadedNodes)
  if i==1
      Fh(i)=Fbtw(i)/2;
  end
  if i==length(HoldLoadedNodes)
      Fh(i)=Fbtw(i-1)/2;
  end
  if i~=1 && i~=length(HoldLoadedNodes)
      Fh(i)=(Fbtw(i)/2)+(Fbtw(i-1)/2);
  end
end
for i=1:(length(Fh)-1)/2
   Fh(i+1+(length(HoldLoadedNodes)-1)/2)=-Fh(i);
end
M=Nodes(HoldLoadedNodes,3)'*Fh;
F(2*HoldLoadedNodes-1,1)=Fh;

ua=K(actDOFs,actDOFs)\F(actDOFs);
u=zeros(2*length(Nodes),1);u(actDOFs,1)=ua;   %eliminatiion method
XYtotfnew=zeros(length(Elements(:,1)),6);
Nodesst=zeros(length(Nodes),3);Nodesst(:,[2 3])=Nodes(:,[2 3])+reshape(u,2,length(Nodes))';Nodesst(:,1)=Nodes(:,1);
for i=1:length(Elements(:,1))
XYtotfnew(i,[1 2])=[Nodesst(find(Nodes(:,1)==Elements(i,1)),[2 3])];
XYtotfnew(i,[3 4])=[Nodesst(find(Nodes(:,1)==Elements(i,2)),[2 3])];
XYtotfnew(i,[5 6])=[Nodesst(find(Nodes(:,1)==Elements(i,3)),[2 3])];
end
figure(500)
for i=1:length(XYtotfnew)
    plot([XYtotfnew(i,1) XYtotfnew(i,3) XYtotfnew(i,5)],[XYtotfnew(i,2),XYtotfnew(i,4),XYtotfnew(i,6)],'b');hold on
end
title('kampsh katharh-prasino prin mple meta')
for i=1:length(HoldBoundaryNodes)
plot(Nodesst(HoldBoundaryNodes(i),2),Nodesst(HoldBoundaryNodes(i),3),'ro')
end
for i=1:length(XYtotf)
    plot([XYtotf(i,1) XYtotf(i,3) XYtotf(i,5)],[XYtotf(i,2),XYtotf(i,4),XYtotf(i,6)],'g');hold on
end
eps=zeros(3,size(Elements,1)); tension=zeros(3,size(Elements,1)); u2=u;
for i=1:length(Elements)
    A=1/2*det([1 Nodes(find(Nodes(:,1)==Elements(i,1)),2)  Nodes(find(Nodes(:,1)==Elements(i,1)),3);...
        1 Nodes(find(Nodes(:,1)==Elements(i,2)),2) Nodes(find(Nodes(:,1)==Elements(i,2)),3); 1 Nodes(find(Nodes(:,1)==Elements(i,3)),2) Nodes(find(Nodes(:,1)==Elements(i,3)),3)]);
    b= [ Nodes(find(Nodes(:,1)==Elements(i,2)),3)-Nodes(find(Nodes(:,1)==Elements(i,3)),3);...
        Nodes(find(Nodes(:,1)==Elements(i,3)),3)-Nodes(find(Nodes(:,1)==Elements(i,1)),3);...
        Nodes(find(Nodes(:,1)==Elements(i,1)),3)-Nodes(find(Nodes(:,1)==Elements(i,2)),3)];
    c=[Nodes(find(Nodes(:,1)==Elements(i,3)),2)-Nodes(find(Nodes(:,1)==Elements(i,2)),2);...
        Nodes(find(Nodes(:,1)==Elements(i,1)),2)-Nodes(find(Nodes(:,1)==Elements(i,3)),2);...
        Nodes(find(Nodes(:,1)==Elements(i,2)),2)-Nodes(find(Nodes(:,1)==Elements(i,1)),2)]; 
    B=(1/(2*A))*[b(1,1) 0 b(2,1) 0 b(3,1) 0;0 c(1,1) 0 c(2,1) 0 c(3,1);c(1,1) b(1,1) c(2,1) b(2,1) c(3,1) b(3,1)];ke=(B'*D*B)*t*A;
    eps(:,i)=B*[u2(2*(find(Nodes(:,1)==Elements(i,1)))-1);...
        u2(2*(find(Nodes(:,1)==Elements(i,1))));...
        u2(2*(find(Nodes(:,1)==Elements(i,2)))-1);...
        u2(2*(find(Nodes(:,1)==Elements(i,2))));...
        u2(2*(find(Nodes(:,1)==Elements(i,3)))-1);...
        u2(2*(find(Nodes(:,1)==Elements(i,3))))]; tension(:,i)=D*eps(:,i);   %per element
end
sigx=zeros(length(Nodes),1);sigy=zeros(length(Nodes),1);dutot=zeros(length(Nodes),1);tau=zeros(length(Nodes),1);
ex=zeros(length(Nodes),1);ey=zeros(length(Nodes),1);gama=zeros(length(Nodes),1);vonmises=zeros(length(Nodes),1);
u1=reshape(u,2,length(Nodes))';ihold=zeros(1,length(Nodes(:,1)));
for i=1:length(Nodes)
    ihold(i)=0;
  for j=1:size(Elements,1)
   if (Elements(j,1)==Nodes(i,1) || Elements(j,2)==Nodes(i,1) || Elements(j,3)==Nodes(i,1))
       sigx(i)=sigx(i)+tension(1,j);sigy(i)=sigy(i)+tension(2,j);tau(i)=tau(i)+tension(3,j);
       ex(i)=ex(i)+eps(1,j);ey(i)=ey(i)+eps(2,j);gama(i)=gama(i)+eps(3,j);
       ihold(i)=ihold(i)+1;      %#of elems. sourrounding one node 
       dutot(i)=sqrt((u1(i,1)^2+u1(i,2)^2));
       vonmises(i)=sqrt(sigx(i)^2-(sigx(i)*sigy(i))+(sigy(i)^2)+(3*(tau(i)^2)));
   end
  end
   sigx(i)=sigx(i)/ihold(i);sigy(i)=sigy(i)/ihold(i);tau(i)=tau(i)/ihold(i);  %tension mean val
   ex(i)=ex(i)/ihold(i);ey(i)=ey(i)/ihold(i);gama(i)=gama(i)/ihold(i);        %deformation mean val
   vonmises(i)=vonmises(i)/ihold(i);
end
sigma1=zeros(length(Nodes),1);sigma2=zeros(length(Nodes),1);
for i=1:size(Nodes,1)
sigma1(i,1)=((sigx(i)+sigy(i))/2)+sqrt(((sigx(i)-sigy(i))/2)^2+tau(i)^2);
sigma2(i,1)=((sigx(i)+sigy(i))/2)-sqrt(((sigx(i)-sigy(i))/2)^2+tau(i)^2);
end
colormap;Sigx=zeros(length(Elements(:,1)),3);
figure(2);title('normal stress {\sigma}_x[N/{mm}^2]');
for i=1:size(Elements,1)
Sigx(i,:)=[sigx(find(Nodes(:,1)==Elements(i,1))) sigx(find(Nodes(:,1)==Elements(i,2)))...
     sigx(find(Nodes(:,1)==Elements(i,3)))];
end
for i=1:size(Elements,1)
  patch([XYtotfnew(i,1) XYtotfnew(i,3) XYtotfnew(i,5)],[XYtotfnew(i,2),XYtotfnew(i,4),XYtotfnew(i,6)],Sigx(i,:)');
  hold on
end
axis equal;colorbar;
colormap;Sigy=zeros(length(Elements(:,1)),3);
figure(3);title('normal stress {\sigma}_y[N/{mm}^2]');
for i=1:size(Elements,1)
Sigy(i,:)=[sigy(find(Nodes(:,1)==Elements(i,1))) sigy(find(Nodes(:,1)==Elements(i,2)))...
     sigy(find(Nodes(:,1)==Elements(i,3))) ];
end
for i=1:size(Elements,1)
   patch([XYtotfnew(i,1) XYtotfnew(i,3) XYtotfnew(i,5)],[XYtotfnew(i,2),XYtotfnew(i,4),XYtotfnew(i,6)],Sigy(i,:)');
   hold on
end
axis equal;colorbar;
colormap;Tau=zeros(length(Elements(:,1)),3);
figure(4);title('shear stress {\tau}_{xy}[N/{mm}^2]');
for i=1:size(Elements,1)
Tau(i,:)=[tau(find(Nodes(:,1)==Elements(i,1))) tau(find(Nodes(:,1)==Elements(i,2)))...
    tau(find(Nodes(:,1)==Elements(i,3))) ];
end
 for i=1:size(Elements,1)
   patch([XYtotfnew(i,1) XYtotfnew(i,3) XYtotfnew(i,5)],[XYtotfnew(i,2),XYtotfnew(i,4),XYtotfnew(i,6)],Tau(i,:)');
   hold on
end
axis equal;colorbar;
colormap;Du=zeros(length(Elements(:,1)),3); 
 figure(5);title('total du[mm]')
for i=1:size(Elements,1)
Du(i,:)=[dutot(find(Nodes(:,1)==Elements(i,1))) dutot(find(Nodes(:,1)==Elements(i,2)))...
     dutot(find(Nodes(:,1)==Elements(i,3))) ];
end
 for i=1:size(Elements,1)
   patch([XYtotfnew(i,1) XYtotfnew(i,3) XYtotfnew(i,5)],[XYtotfnew(i,2),XYtotfnew(i,4),XYtotfnew(i,6)],Du(i,:)');
   hold on
 end
axis equal;colorbar;
 colormap;Ex=zeros(length(Elements(:,1)),3);
figure(6);title('strain displacement ex'); 
for i=1:size(Elements,1)
Ex(i,:)=[ex(find(Nodes(:,1)==Elements(i,1))) ex(find(Nodes(:,1)==Elements(i,2)))...
     ex(find(Nodes(:,1)==Elements(i,3))) ];
end
for i=1:size(Elements,1)
   patch([XYtotfnew(i,1) XYtotfnew(i,3) XYtotfnew(i,5)],[XYtotfnew(i,2),XYtotfnew(i,4),XYtotfnew(i,6)],Ex(i,:)');
   hold on
end
axis equal;colorbar;
colormap;Ey=zeros(length(Elements(:,1)),3);
figure(7);title('strain displacement ey');
for i=1:size(Elements,1)
Ey(i,:)=[ey(find(Nodes(:,1)==Elements(i,1))) ey(find(Nodes(:,1)==Elements(i,2)))...
     ey(find(Nodes(:,1)==Elements(i,3)))];
end
for i=1:size(Elements,1)
   patch([XYtotfnew(i,1) XYtotfnew(i,3) XYtotfnew(i,5)],[XYtotfnew(i,2),XYtotfnew(i,4),XYtotfnew(i,6)],Ey(i,:)');
   hold on
end
axis equal;colorbar;
colormap;Gam=zeros(length(Elements(:,1)),3);
figure(8);title('strain displacement {\gamma}_{xy}');
for i=1:size(Elements,1)
Gam(i,:)=[gama(find(Nodes(:,1)==Elements(i,1))) gama(find(Nodes(:,1)==Elements(i,2)))...
    gama(find(Nodes(:,1)==Elements(i,3)))];
end
for i=1:size(Elements,1)
   patch([XYtotfnew(i,1) XYtotfnew(i,3) XYtotfnew(i,5)],[XYtotfnew(i,2),XYtotfnew(i,4),XYtotfnew(i,6)],Gam(i,:)');
   hold on
end
axis equal;colorbar;
colormap;Sig1=zeros(length(Elements(:,1)),3);
figure(9);title(' {\sigma}_{1}[N/mm^2]');
for i=1:size(Elements,1)
Sig1(i,:)=[sigma1(find(Nodes(:,1)==Elements(i,1))) sigma1(find(Nodes(:,1)==Elements(i,2)))...
     sigma1(find(Nodes(:,1)==Elements(i,3))) ];
end
for i=1:size(Elements,1)
   patch([XYtotfnew(i,1) XYtotfnew(i,3) XYtotfnew(i,5)],[XYtotfnew(i,2),XYtotfnew(i,4),XYtotfnew(i,6)],Sig1(i,:)');
   hold on
end
axis equal;colorbar;
colormap;Sig2=zeros(length(Elements(:,1)),3);
figure(10);title(' {\sigma}_{2}[N/mm^2]');
for i=1:size(Elements,1)
Sig2(i,:)=[sigma2(find(Nodes(:,1)==Elements(i,1))) sigma2(find(Nodes(:,1)==Elements(i,2)))...
     sigma2(find(Nodes(:,1)==Elements(i,3))) ];
end
for i=1:size(Elements,1)
   patch([XYtotfnew(i,1) XYtotfnew(i,3) XYtotfnew(i,5)],[XYtotfnew(i,2),XYtotfnew(i,4),XYtotfnew(i,6)],Sig2(i,:)');
   hold on
end
axis equal;colorbar;
colormap;VM=zeros(length(Elements(:,1)),3);
figure(25);title('vonmises[N/mm^2]');
for i=1:size(Elements,1)
VM(i,:)=[vonmises(find(Nodes(:,1)==Elements(i,1))) vonmises(find(Nodes(:,1)==Elements(i,2)))...
     vonmises(find(Nodes(:,1)==Elements(i,3)))];
end
for i=1:size(Elements,1)
   patch([XYtotfnew(i,1) XYtotfnew(i,3) XYtotfnew(i,5)],[XYtotfnew(i,2),XYtotfnew(i,4),XYtotfnew(i,6)],VM(i,:)');
   hold on
end
axis equal;colorbar;
figure(22);
for i=1:length(Nodes(:,1))
    dx=sigma1(i).*cos((1/2)*atan((2*tau(i))./(sigx(i)-sigy(i))));
    dy=sigma1(i).*sin((1/2)*atan((2*tau(i))./(sigx(i)-sigy(i))));
   quiver(Nodes(i,2),Nodes(i,3),dx*scale1,dy*scale1,'b')
   hold on;
end
title(['print arrows for {\sigma}_{1} scale for arrows:' num2str(scale1)])
figure(23);
for i=1:length(Nodes(:,1))
    dx=sigma2(i).*cos(pi/2-((1/2)*atan((2*tau(i))./(sigx(i)-sigy(i)))));
    dy=sigma2(i).*sin(pi/2-((1/2)*atan((2*tau(i))./(sigx(i)-sigy(i)))));
   quiver(Nodes(i,2),Nodes(i,3),dx*scale1,dy*scale1,'b')
   hold on;
end
title(['print arrows for {\sigma}_{2} scale for arrows:' num2str(scale1)])
%%
%in txt file
[Ysx,Isx]=max(sigx);[Ysy,Isy]=max(sigy);[Yt,It]=max(tau);[Ys1,Is1]=max(sigma1);[Yvm,Ivm]=max(vonmises);
[Ydu,Idu]=max(dutot);[Yex,Iex]=max(ex);[Yey,Iey]=max(ey);[Yg,Ig]=max(gama);
fpr=[Ysx,Nodes(Isx,:);Ysy,Nodes(Isy,:);Yt,Nodes(It,:);Ys1,Nodes(Is1,:);Yvm,Nodes(Ivm,:);...
    Ydu,Nodes(Idu,:);Yex,Nodes(Iex,:);Yey,Nodes(Iey,:);Yg,Nodes(Ig,:)];


nod=fopen('AEM5380_2nd_A_ResultsBending_MeshingCir.txt','w+t');
fprintf(nod,  'Maximum_Value #Code_of_NODE      X[mm]    Y[mm]  \n');fprintf('\n');
fprintf(nod,'  sigma_{x}[N/mm2] \n');fprintf('\n');
fprintf(nod,'  %f      %.0f        %f     %f  \n',[ fpr(1,:) ] );fprintf('\n');
fprintf(nod,'  sigma_{y}[N/mm2] \n');fprintf('\n');
fprintf(nod,'  %f      %.0f        %f     %f  \n',[ fpr(2,:) ] );fprintf('\n');
fprintf(nod,'  tau_{xy}[N/mm2] \n');fprintf('\n');
fprintf(nod,'  %f      %.0f        %f     %f  \n',[ fpr(3,:) ] );fprintf('\n');
fprintf(nod,'  principal stress sigma_{1}[N/mm2] \n');fprintf('\n');
fprintf(nod,'  %f      %.0f        %f     %f  \n',[ fpr(4,:) ] );fprintf('\n');
fprintf(nod,'  Max_Von_Mises_stress[N/mm2] \n');fprintf('\n');
fprintf(nod,'  %f      %.0f        %f     %f  \n',[ fpr(5,:) ] );fprintf('\n');
fprintf(nod,'  du \n');fprintf('\n');
fprintf(nod,'  %f         %.0f        %f     %f  \n',[ fpr(6,:) ] );fprintf('\n');
fprintf(nod,'  epsilon_{x} \n');fprintf('\n');
fprintf(nod,'  %f         %.0f        %f     %f  \n',[ fpr(7,:) ] );fprintf('\n');
fprintf(nod,'  epsilon_{y} \n');fprintf('\n');
fprintf(nod,'  %f         %.0f        %f     %f  \n',[ fpr(8,:) ] );fprintf('\n');
fprintf(nod,'  gamma_{xy} \n');fprintf('\n');
fprintf(nod,'  %f         %.0f        %f     %f  \n',[ fpr(9,:) ] );fprintf('\n');
alpha_kappa=Ysx/mean(abs(sigx(find(Nodes(:,2)==Nodes(Isx,2)))));
fprintf(nod,'  alpha_{kappa} \n');fprintf('\n');
fprintf(nod,'      %f  \n',[ alpha_kappa ] );fprintf('\n');
fclose(nod);

%figure(1)
%for i=1:length(Elements)
%  patch([XYtotfnew(i,1) XYtotfnew(i,3) XYtotfnew(i,5)],[XYtotfnew(i,2),XYtotfnew(i,4),XYtotfnew(i,6)],tension(1,i));
%  hold on
%end

%%
%anamenomenos syntelesths egkophs kt gia kampsh

Kttttkampsh=2

sigmanominal1=6*My*d/(t*(100^3-d^3));
holdsimgamax11=vonmises(find(Nodes(find(Nodes(:,2)==0),3)==13.3));

Ktttkampsh=holdsimgamax11/sigmanominal1
%%

HoldLoadedNodes=find(Nodes(:,2)==85);
f=f/(length(HoldLoadedNodes)-1);
F=zeros(2*length(Nodes),1);F(2*HoldLoadedNodes-1,1)=f;
F(2*HoldLoadedNodes(1)-1,1)=f/2;F(2*HoldLoadedNodes(length(HoldLoadedNodes))-1,1)=f/2;
%u=K\F;                                       %penalty method

HoldLoadedNodes=find(Nodes(:,2)==-85);
F(2*HoldLoadedNodes-1,1)=-f;
F(2*HoldLoadedNodes(1)-1,1)=-f/2;F(2*HoldLoadedNodes(length(HoldLoadedNodes))-1,1)=-f/2;

ua=K(actDOFs,actDOFs)\F(actDOFs);
u=zeros(2*length(Nodes),1);u(actDOFs,1)=ua;   %eliminatiion method
XYtotfnew=zeros(length(Elements(:,1)),6);
Nodesst=zeros(length(Nodes),3);Nodesst(:,[2 3])=Nodes(:,[2 3])+reshape(u,2,length(Nodes))';Nodesst(:,1)=Nodes(:,1);
for i=1:length(Elements(:,1))
XYtotfnew(i,[1 2])=[Nodesst(find(Nodes(:,1)==Elements(i,1)),[2 3])];
XYtotfnew(i,[3 4])=[Nodesst(find(Nodes(:,1)==Elements(i,2)),[2 3])];
XYtotfnew(i,[5 6])=[Nodesst(find(Nodes(:,1)==Elements(i,3)),[2 3])];
end
figure(600)
for i=1:length(XYtotfnew)
    plot([XYtotfnew(i,1) XYtotfnew(i,3) XYtotfnew(i,5)],[XYtotfnew(i,2),XYtotfnew(i,4),XYtotfnew(i,6)],'b');hold on
end
title('efelkusmos katharos-prasino prin mple meta')
for i=1:length(HoldBoundaryNodes)
plot(Nodesst(HoldBoundaryNodes(i),2),Nodesst(HoldBoundaryNodes(i),3),'ro')
end
for i=1:length(XYtotf)
    plot([XYtotf(i,1) XYtotf(i,3) XYtotf(i,5)],[XYtotf(i,2),XYtotf(i,4),XYtotf(i,6)],'g');hold on
end

eps=zeros(3,size(Elements,1)); tension=zeros(3,size(Elements,1)); u2=u;
for i=1:length(Elements)
    A=1/2*det([1 Nodes(find(Nodes(:,1)==Elements(i,1)),2)  Nodes(find(Nodes(:,1)==Elements(i,1)),3);...
        1 Nodes(find(Nodes(:,1)==Elements(i,2)),2) Nodes(find(Nodes(:,1)==Elements(i,2)),3); 1 Nodes(find(Nodes(:,1)==Elements(i,3)),2) Nodes(find(Nodes(:,1)==Elements(i,3)),3)]);
    b= [ Nodes(find(Nodes(:,1)==Elements(i,2)),3)-Nodes(find(Nodes(:,1)==Elements(i,3)),3);...
        Nodes(find(Nodes(:,1)==Elements(i,3)),3)-Nodes(find(Nodes(:,1)==Elements(i,1)),3);...
        Nodes(find(Nodes(:,1)==Elements(i,1)),3)-Nodes(find(Nodes(:,1)==Elements(i,2)),3)];
    c=[Nodes(find(Nodes(:,1)==Elements(i,3)),2)-Nodes(find(Nodes(:,1)==Elements(i,2)),2);...
        Nodes(find(Nodes(:,1)==Elements(i,1)),2)-Nodes(find(Nodes(:,1)==Elements(i,3)),2);...
        Nodes(find(Nodes(:,1)==Elements(i,2)),2)-Nodes(find(Nodes(:,1)==Elements(i,1)),2)]; 
    B=(1/(2*A))*[b(1,1) 0 b(2,1) 0 b(3,1) 0;0 c(1,1) 0 c(2,1) 0 c(3,1);c(1,1) b(1,1) c(2,1) b(2,1) c(3,1) b(3,1)];ke=(B'*D*B)*t*A;
   eps(:,i)=B*[u2(2*(find(Nodes(:,1)==Elements(i,1)))-1);...
        u2(2*(find(Nodes(:,1)==Elements(i,1))));...
        u2(2*(find(Nodes(:,1)==Elements(i,2)))-1);...
        u2(2*(find(Nodes(:,1)==Elements(i,2))));...
        u2(2*(find(Nodes(:,1)==Elements(i,3)))-1);...
        u2(2*(find(Nodes(:,1)==Elements(i,3))))]; tension(:,i)=D*eps(:,i);   %per element
end
sigx=zeros(length(Nodes),1);sigy=zeros(length(Nodes),1);dutot=zeros(length(Nodes),1);tau=zeros(length(Nodes),1);
ex=zeros(length(Nodes),1);ey=zeros(length(Nodes),1);gama=zeros(length(Nodes),1);vonmises=zeros(length(Nodes),1);
u1=reshape(u,2,length(Nodes))';ihold=zeros(1,length(Nodes(:,1)));
for i=1:length(Nodes)
    ihold(i)=0;
  for j=1:size(Elements,1)
   if (Elements(j,1)==Nodes(i,1) || Elements(j,2)==Nodes(i,1) || Elements(j,3)==Nodes(i,1))
       sigx(i)=sigx(i)+tension(1,j);sigy(i)=sigy(i)+tension(2,j);tau(i)=tau(i)+tension(3,j);
       ex(i)=ex(i)+eps(1,j);ey(i)=ey(i)+eps(2,j);gama(i)=gama(i)+eps(3,j);
       ihold(i)=ihold(i)+1;      %#of elems. sourrounding one node 
       dutot(i)=sqrt((u1(i,1)^2+u1(i,2)^2));
       vonmises(i)=sqrt(sigx(i)^2-(sigx(i)*sigy(i))+(sigy(i)^2)+(3*(tau(i)^2)));
   end
  end
   sigx(i)=sigx(i)/ihold(i);sigy(i)=sigy(i)/ihold(i);tau(i)=tau(i)/ihold(i);  %tension mean val
   ex(i)=ex(i)/ihold(i);ey(i)=ey(i)/ihold(i);gama(i)=gama(i)/ihold(i);        %deformation mean val
   vonmises(i)=vonmises(i)/ihold(i);
end
sigma1=zeros(length(Nodes),1);sigma2=zeros(length(Nodes),1);
for i=1:size(Nodes,1)
sigma1(i,1)=((sigx(i)+sigy(i))/2)+sqrt(((sigx(i)-sigy(i))/2)^2+tau(i)^2);
sigma2(i,1)=((sigx(i)+sigy(i))/2)-sqrt(((sigx(i)-sigy(i))/2)^2+tau(i)^2);
end
colormap;Sigx=zeros(length(Elements(:,1)),3);
figure(11);title('normal stress {\sigma}_x[N/{mm}^2]');
for i=1:size(Elements,1)
Sigx(i,:)=[sigx(find(Nodes(:,1)==Elements(i,1))) sigx(find(Nodes(:,1)==Elements(i,2)))...
     sigx(find(Nodes(:,1)==Elements(i,3))) ];
end
for i=1:size(Elements,1)
  patch([XYtotfnew(i,1) XYtotfnew(i,3) XYtotfnew(i,5)],[XYtotfnew(i,2),XYtotfnew(i,4),XYtotfnew(i,6)],Sigx(i,:)');
  hold on
end
axis equal;colorbar;
colormap;Sigy=zeros(length(Elements(:,1)),3);
figure(12);title('normal stress {\sigma}_y[N/{mm}^2]');
for i=1:size(Elements,1)
Sigy(i,:)=[sigy(find(Nodes(:,1)==Elements(i,1))) sigy(find(Nodes(:,1)==Elements(i,2)))...
    sigy(find(Nodes(:,1)==Elements(i,3))) ];
end
for i=1:size(Elements,1)
   patch([XYtotfnew(i,1) XYtotfnew(i,3) XYtotfnew(i,5)],[XYtotfnew(i,2),XYtotfnew(i,4),XYtotfnew(i,6)],Sigy(i,:)');
   hold on
end
axis equal;colorbar;
colormap;Tau=zeros(length(Elements(:,1)),3);
figure(13);title('shear stress {\tau}_{xy}[N/{mm}^2]');
for i=1:size(Elements,1)
Tau(i,:)=[tau(find(Nodes(:,1)==Elements(i,1))) tau(find(Nodes(:,1)==Elements(i,2)))...
     tau(find(Nodes(:,1)==Elements(i,3))) ];
end
 for i=1:size(Elements,1)
   patch([XYtotfnew(i,1) XYtotfnew(i,3) XYtotfnew(i,5)],[XYtotfnew(i,2),XYtotfnew(i,4),XYtotfnew(i,6)],Tau(i,:)');
   hold on
end
axis equal;colorbar;
colormap;Du=zeros(length(Elements(:,1)),3); 
figure(14);title('total du[mm]')
for i=1:size(Elements,1)
Du(i,:)=[dutot(find(Nodes(:,1)==Elements(i,1))) dutot(find(Nodes(:,1)==Elements(i,2)))...
    dutot(find(Nodes(:,1)==Elements(i,3))) ];
end
 for i=1:size(Elements,1)
   patch([XYtotfnew(i,1) XYtotfnew(i,3) XYtotfnew(i,5)],[XYtotfnew(i,2),XYtotfnew(i,4),XYtotfnew(i,6)],Du(i,:)');
   hold on
 end
axis equal;colorbar;
 colormap;Ex=zeros(length(Elements(:,1)),3);
figure(15);title('strain displacement ex'); 
for i=1:size(Elements,1)
Ex(i,:)=[ex(find(Nodes(:,1)==Elements(i,1))) ex(find(Nodes(:,1)==Elements(i,2)))...
     ex(find(Nodes(:,1)==Elements(i,3))) ];
end
for i=1:size(Elements,1)
   patch([XYtotfnew(i,1) XYtotfnew(i,3) XYtotfnew(i,5)],[XYtotfnew(i,2),XYtotfnew(i,4),XYtotfnew(i,6)],Ex(i,:)');
   hold on
end
axis equal;colorbar;
colormap;Ey=zeros(length(Elements(:,1)),3);
figure(16);title('strain displacement ey');
for i=1:size(Elements,1)
Ey(i,:)=[ey(find(Nodes(:,1)==Elements(i,1))) ey(find(Nodes(:,1)==Elements(i,2)))...
     ey(find(Nodes(:,1)==Elements(i,3))) ];
end
for i=1:size(Elements,1)
   patch([XYtotfnew(i,1) XYtotfnew(i,3) XYtotfnew(i,5)],[XYtotfnew(i,2),XYtotfnew(i,4),XYtotfnew(i,6)],Ey(i,:)');
   hold on
end
axis equal;colorbar;
colormap;Gam=zeros(length(Elements(:,1)),3);
figure(17);title('strain displacement {\gamma}_{xy}');
for i=1:size(Elements,1)
Gam(i,:)=[gama(find(Nodes(:,1)==Elements(i,1))) gama(find(Nodes(:,1)==Elements(i,2)))...
     gama(find(Nodes(:,1)==Elements(i,3)))];
end
for i=1:size(Elements,1)
   patch([XYtotfnew(i,1) XYtotfnew(i,3) XYtotfnew(i,5)],[XYtotfnew(i,2),XYtotfnew(i,4),XYtotfnew(i,6)],Gam(i,:)');
   hold on
end
axis equal;colorbar;
colormap;Sig1=zeros(length(Elements(:,1)),3);
figure(18);title(' {\sigma}_{1}[N/mm^2]');
for i=1:size(Elements,1)
Sig1(i,:)=[sigma1(find(Nodes(:,1)==Elements(i,1))) sigma1(find(Nodes(:,1)==Elements(i,2)))...
     sigma1(find(Nodes(:,1)==Elements(i,3))) ];
end
for i=1:size(Elements,1)
   patch([XYtotfnew(i,1) XYtotfnew(i,3) XYtotfnew(i,5)],[XYtotfnew(i,2),XYtotfnew(i,4),XYtotfnew(i,6)],Sig1(i,:)');
   hold on
end
axis equal;colorbar;
colormap;Sig2=zeros(length(Elements(:,1)),3);
figure(19);title(' {\sigma}_{2}[N/mm^2]');
for i=1:size(Elements,1)
Sig2(i,:)=[sigma2(find(Nodes(:,1)==Elements(i,1))) sigma2(find(Nodes(:,1)==Elements(i,2)))...
     sigma2(find(Nodes(:,1)==Elements(i,3))) ];
end
for i=1:size(Elements,1)
   patch([XYtotfnew(i,1) XYtotfnew(i,3) XYtotfnew(i,5)],[XYtotfnew(i,2),XYtotfnew(i,4),XYtotfnew(i,6)],Sig2(i,:)');
   hold on
end
axis equal;colorbar;
colormap;VM=zeros(length(Elements(:,1)),3);
figure(24);title('vonmises[N/mm^2]');
for i=1:size(Elements,1)
VM(i,:)=[vonmises(find(Nodes(:,1)==Elements(i,1))) vonmises(find(Nodes(:,1)==Elements(i,2)))...
    vonmises(find(Nodes(:,1)==Elements(i,3)))];
end
for i=1:size(Elements,1)
   patch([XYtotfnew(i,1) XYtotfnew(i,3) XYtotfnew(i,5)],[XYtotfnew(i,2),XYtotfnew(i,4),XYtotfnew(i,6)],VM(i,:)');
   hold on
end
axis equal;colorbar;
figure(20);
for i=1:length(Nodes(:,1))
    dx=sigma1(i).*cos((1/2)*atan((2*tau(i))./(sigx(i)-sigy(i))));
    dy=sigma1(i).*sin((1/2)*atan((2*tau(i))./(sigx(i)-sigy(i))));
   quiver(Nodes(i,2),Nodes(i,3),dx*scale2,dy*scale2,'b')
   hold on;
end
title(['print arrows for {\sigma}_{1} scale for arrows:' num2str(scale2)])
figure(21);
for i=1:length(Nodes(:,1))
    dx=sigma2(i).*cos(pi/2-((1/2)*atan((2*tau(i))./(sigx(i)-sigy(i)))));
    dy=sigma2(i).*sin(pi/2-((1/2)*atan((2*tau(i))./(sigx(i)-sigy(i)))));
   quiver(Nodes(i,2),Nodes(i,3),dx*scale2,dy*scale2,'b')
   hold on;
end
title(['print arrows for {\sigma}_{2} scale for arrows:' num2str(scale2)])
%%
%in txt file
[Ysx,Isx]=max(sigx);[Ysy,Isy]=max(sigy);[Yt,It]=max(tau);[Ys1,Is1]=max(sigma1);[Yvm,Ivm]=max(vonmises);
[Ydu,Idu]=max(dutot);[Yex,Iex]=max(ex);[Yey,Iey]=max(ey);[Yg,Ig]=max(gama);
fpr=[Ysx,Nodes(Isx,:);Ysy,Nodes(Isy,:);Yt,Nodes(It,:);Ys1,Nodes(Is1,:);Yvm,Nodes(Ivm,:);...
    Ydu,Nodes(Idu,:);Yex,Nodes(Iex,:);Yey,Nodes(Iey,:);Yg,Nodes(Ig,:)];


nod=fopen('AEM5380_2nd_A_ResultsTensile_MeshingCir.txt','w+t');
fprintf(nod,  'Maximum_Value #Code_of_NODE      X[mm]    Y[mm]  \n');fprintf('\n');
fprintf(nod,'  sigma_{x}[N/mm2] \n');fprintf('\n');
fprintf(nod,'  %f         %.0f        %f     %f  \n',[ fpr(1,:) ] );fprintf('\n');
fprintf(nod,'  sigma_{y}[N/mm2] \n');fprintf('\n');
fprintf(nod,'  %f         %.0f        %f     %f  \n',[ fpr(2,:) ] );fprintf('\n');
fprintf(nod,'  tau_{xy}[N/mm2] \n');fprintf('\n');
fprintf(nod,'  %f         %.0f        %f     %f  \n',[ fpr(3,:) ] );fprintf('\n');
fprintf(nod,'  principal stress sigma_{1}[N/mm2] \n');fprintf('\n');
fprintf(nod,'  %f         %.0f        %f     %f  \n',[ fpr(4,:) ] );fprintf('\n');
fprintf(nod,'  Max_Von_Mises_stress[N/mm2] \n');fprintf('\n');
fprintf(nod,'  %f         %.0f        %f     %f  \n',[ fpr(5,:) ] );fprintf('\n');
fprintf(nod,'  du \n');fprintf('\n');
fprintf(nod,'  %f         %.0f        %f     %f  \n',[ fpr(6,:) ] );fprintf('\n');
fprintf(nod,'  epsilon_{x} \n');fprintf('\n');
fprintf(nod,'  %f         %.0f        %f     %f  \n',[ fpr(7,:) ] );fprintf('\n');
fprintf(nod,'  epsilon_{y} \n');fprintf('\n');
fprintf(nod,'  %f         %.0f        %f     %f  \n',[ fpr(8,:) ] );fprintf('\n');
fprintf(nod,'  gamma_{xy} \n');fprintf('\n');
fprintf(nod,'  %f         %.0f        %f     %f  \n',[ fpr(9,:) ] );fprintf('\n');
alpha_kappa=Ysx/mean(abs(sigx(find(Nodes(:,2)==Nodes(Isx,2)))));
fprintf(nod,'  alpha_{kappa} \n');fprintf('\n');
fprintf(nod,'      %f  \n',[ alpha_kappa ] );fprintf('\n');
fclose(nod);clc;


%figure(1000)
%for i=1:length(Elements)
%  patch([XYtotfnew(i,1) XYtotfnew(i,3) XYtotfnew(i,5)],[XYtotfnew(i,2),XYtotfnew(i,4),XYtotfnew(i,6)],tension(1,i));
%  hold on
%end

%%
%anamenomenos syntelesths egkophs kt gia efelkusmo

Kttttefelkusmo=3-3.14*(d/100)+3.667*(d/100)^2-1.527*(d/100)^3

sigmanominal=(100/(100-d))*sum(vonmises(find(Nodes(:,2)==85)))/length(find(Nodes(:,2)==85));
holdsimgamax1=vonmises(find(Nodes(find(Nodes(:,2)==0),3)==13.3));

Ktttefelkusmo=holdsimgamax1/sigmanominal