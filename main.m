%--------------------------------PROJECT 1---------------------------------------------%
% Written by Sohanjit 
clc,clear 
    filename = 'Uniaxial_linear.txt';
%   filename = 'Uniaxial_quadratic.txt';
% filename = 'demo.txt';

syms e


[node, element, elemType, nel, nen, nip, nnd, Eparam, Force_Node, qparam, disp_node, disp_BC, trac_el, trac_BC]=Read_input(filename);


if elemType == 'Q2'

[N,B] = linear_sf(e);

elseif elemType == 'Q3'

[N,B] = quadratic_sf(e);

end

%       nip =1 ; %Uncomment this line if interested in 1-point Gauss integration

[w,IC] = GQ(nip);

K = zeros(nnd,nnd);
ke = zeros(nen,nen);

h = 1/4;
J= h/2;

%------Element stiffness matrix

for i=1:nen
    for j=1:nen
        sum=0;
        for k=1:nip
            sum=sum+subs(B(i)*B(j)/J,e,IC(k))*w(k);
        end
            ke(i,j)= ke(i,j)+sum;
     end
end

%------Global stiffness matrix

for el =1:nel

    for i=1:nen
        for j=1:nen

K(element(el,i),element(el,j)) = K(element(el,i),element(el,j))+ke(i,j);

        end
    end
end

%------Element force vector

F = zeros(nnd,1);

for el=1:nel

    el_nodes=zeros(nen,1);

for t =1:nen
    el_nodes(t)= node(element(el,t));
end

fe = zeros(nen,1);

for i=1:nen
    sigma=0;
    
    for j=1:nen
        sigma = sigma+(N(j)*el_nodes(j));
    end

        sum=0;
        for k=1:nip
            sum=sum+subs(N(i)*qparam*sigma*J,e,IC(k))*w(k);
        end
            fe(i,1)= fe(i,1)+sum;
   
end


for i=1:nen
F(element(el,i),1) = F(element(el,i),1)+fe(i,1);
end

end



F_traction = zeros(nnd,1);
id = element(trac_el(1),trac_el(2));

for i =1:nnd

    if i==id

        F_traction(i)=trac_BC;
    end

end

F = F+F_traction;

%------Solving for nodal displacements-------%
 

CK = K;CF = F; %Copied matrices

d = zeros(nnd,1);

% for i=1:nnd
%     for j=1:nnd
% if i == disp_node || j==disp_node
%     CK(i,j)=0;
% end
% if i == disp_node && j==disp_node
%     CK(i,j)=1;
% end
% end
% end    

CK(:,disp_node)=[];
CK(disp_node,:)=[];

for i=1:nnd

if i ~= disp_node

    CF(i,1)=F(i,1)-K(i,disp_node)*disp_BC;
end


end

CF(disp_node,:)=[];

d=inv(CK)*CF;


d(disp_node)=disp_BC;


%------Element mass matrix----%

me=zeros(nen,nen);
M=zeros(nnd,nnd);

for i = 1:nen
    for j=1:nen

        sum=0;
        for k=1:nip
            sum=sum+subs(N(i)*N(j)*J,e,IC(k))*w(k);
        end
            me(i,j)= me(i,j)+sum;

    end
end


for el =1:nel

    for i=1:nen
        for j=1:nen

M(element(el,i),element(el,j)) = M(element(el,i),element(el,j))+me(i,j);

        end
    end
end


%--------Projection vector

P = zeros(nnd,1);

for el=1:nel

    el_nodes=zeros(nen,1);

for t =1:nen
    el_nodes(t)= element(el,t);
end


pe = zeros(nen,1);

for i=1:nen
    sigma=0;
    
    for j=1:nen
        sigma = sigma+(B(j)*d(el_nodes(j)));
    end

        sum=0;
        for k=1:nip
            sum=sum+subs(N(i)*qparam*sigma,e,IC(k))*w(k);
        end
            pe(i,1)= pe(i,1)+sum;
   
end


for i=1:nen
P(element(el,i),1) = P(element(el,i),1)+pe(i,1);
end

end


%---------Plotting ---- %

syms u(x) 
Dy = diff(u);

ode = diff(u,x,2) == -qparam*x;
cond1 = u(node(disp_node)) == disp_BC;
cond2 = Dy(node(element(trac_el(1),trac_el(2)))) == -1*trac_BC;

conds = [cond1 cond2];
uSol(x) = dsolve(ode,conds);

%---Nodal Displacement plot---%
figure(1)
exact_disp=zeros(nnd,1);

for i =1:nnd
exact_disp(i)=subs(uSol(x),x,node(i));

end

plot(node,d,'b-*','LineWidth',3)
hold on
plot(node,exact_disp,'r:*','LineWidth',3)
xlabel('Domain','FontSize',20, 'FontName','Arial','FontWeight','bold')
ylabel('Nodal displacement','FontSize',20, 'FontName','Arial','FontWeight','bold')

ax = gca;
ax.FontWeight= 'bold';
ax.FontSize= 16;
legend('FEM','analytical')

%---Stress@nodes plot---%
figure(2)
exact_stress_nodes=zeros(nnd,1);

for i =1:nnd
exact_stress_nodes(i)=subs(Eparam*diff(uSol(x)),x,node(i));

end

fem_stress_nodes=Eparam*inv(M)*P;

plot(node,exact_stress_nodes,'b-*','LineWidth',3)
hold on
plot(node,fem_stress_nodes,'r:*','LineWidth',3)

xlabel('Domain','FontSize',20, 'FontName','Arial','FontWeight','bold')
ylabel('Stress @ Nodes','FontSize',20, 'FontName','Arial','FontWeight','bold')

ax = gca;
ax.FontWeight= 'bold';
ax.FontSize= 16;
legend('FEM','analytical')


%---Stress@IPs plot---%

exact_stress_ips = zeros(nip*nel,1);

ip_coords= zeros(nip*nel,1);
count=1;

for el=1:nel

    el_nodes=zeros(nen,1);
    sum=0;
for t =1:nen
    el_nodes(t)= node(element(el,t));
    sum=sum+N(t)*el_nodes(t);
end

for i=1:nip
ip_coords(count)= subs(sum,e,IC(i));
count=count+1;

end

end


for i =1:length(ip_coords)
   exact_stress_ips(i)=subs(Eparam*diff(uSol(x)),x,ip_coords(i));
        
end

fem_stress_ips=zeros(nip*nel,1);
count=1;

for el=1:nel

    el_nodes=zeros(nip,1);
    sum=0;
for t =1:nen
    el_nodes(t)= element(el,t);
    sum=sum+B(t)*d(el_nodes(t));
end

for i=1:nip
fem_stress_ips(count)=Eparam*subs(sum/J,e,IC(i));
count=count+1;
end
end

% [fem_stress_ips exact_stress_ips]
% 
figure(3)
plot(ip_coords,exact_stress_ips,'b-*','LineWidth',3)
hold on
plot(ip_coords,fem_stress_ips,'r:*','LineWidth',3)

xlabel('Domain','FontSize',20, 'FontName','Arial','FontWeight','bold')
ylabel('Stress @ IPs','FontSize',20, 'FontName','Arial','FontWeight','bold')

ax = gca;
ax.FontWeight= 'bold';
ax.FontSize= 16;
legend('FEM','analytical')

%---------Reaction forces------------% 
% 
% 

R=K*d-F;
disp(R(Force_Node))







function [w,IC] = GQ(nip)

if nip == 1

    w(1)=2;
    IC(1)=0;

elseif nip ==2
w(1)=1;w(2)=1;
IC(1)= - 1/sqrt(3);IC(2)= 1/sqrt(3);

elseif nip == 3
w(1)=5/9;w(2)=8/9;w(3)=5/9;
IC(1)= - sqrt(3/5); IC(2)= 0; IC(3)= sqrt(3/5);

else 
    disp('Exceeded number of IPs')

end

end

function [N,B] = linear_sf(e)

N(1) = (1-e)/2;
N(2) = (1+e)/2;
B(1) = diff(N(1));
B(2) = diff(N(2));


end

function [N,B] = quadratic_sf(e)

N(1) = e*(e-1)/2;
N(2) = 1-e*e;
N(3) = e*(e+1)/2;
B(1) = diff(N(1));
B(2) = diff(N(2));
B(3) = diff(N(3));

end

function [node, element, elemType, nel, nen, nip, nnd, Eparam, Force_Node, qparam, disp_node, disp_BC, trac_el, trac_BC]=Read_input(filename)

fid=fopen(filename, 'r');

C=textscan(fid,'%s',1);
C=textscan(fid,'%n',1);
nnd=C{1,1};
C=textscan(fid,'%n','MultipleDelimsAsOne', 1);
node=C{1};

C=textscan(fid,'%s',2);
C=textscan(fid,'%s',3);
C=textscan(fid,'%n %n','MultipleDelimsAsOne', 1);
disp_node = C{1};
disp_BC = C{2};

C=textscan(fid,'%s',4);
C=textscan(fid,'%n','MultipleDelimsAsOne', 1);
Force_Node=C{1};

C=textscan(fid,'%s',1);
C=textscan(fid,'%n',1);
nel=C{1,1};
C=textscan(fid,'%n',1);
nen=C{1,1};
if (nen==2)
    elemType='Q2';
    C=textscan(fid,'%n %n','MultipleDelimsAsOne', 1);
    element=cell2mat(C(1:2));
elseif (nen==3)
    elemType='Q3';
    C=textscan(fid,'%n %n %n','MultipleDelimsAsOne', 1);
    element=cell2mat(C(1:3));
else
    disp('nen must be 2 or 3');
    return;
end

C=textscan(fid,'%s',2);
C=textscan(fid,'%n',1);
nip=C{1,1};

C=textscan(fid,'%s',2);
C=textscan(fid,'%s',5);
C=textscan(fid,'%n %n %n','MultipleDelimsAsOne', 1);
trac_el = [C{1}, C{2}];
trac_BC = C{3};

C=textscan(fid,'%s',4);
C=textscan(fid,'%n',1);
qparam=C{1,1};

C=textscan(fid,'%s',2);
C=textscan(fid,'%n',1);
Eparam=C{1,1};

fclose(fid);
end