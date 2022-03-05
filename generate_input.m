function [fileName,element,node]=generate_input(l,SF,IP)

% make sure to change the file name
fileName    = 'input.txt';
L           = 1;        % length of the domain
nel         = L/l;        % total number of elements
elemType    = SF;     % type of element: linear and quadratic
nip         = IP;        % no of integration points
Eparam      = 1;        % Young's modulus
qparam      = 1;        % constant of distributed loading
disp_BC     = 0.1;      % value of essential BC
disp_loc    = 1;        % x-coordinate/location of displacement BC applied
trac_BC     = 0.1;      % value of traction BC
trac_loc    = 0;        % location of traction BC to be applied
react_loc   = 1;        % location of reaction force to be calculated

    
fid=fopen(fileName, 'w');

if strcmp(elemType,'Q2')
    nen     = 2;
elseif strcmp(elemType,'Q3')
    nen     = 3;
else
    disp('Invalid element type');
    return;
end

nnd         = (nen-1)*nel+1;    % total no of nodes

delL        = L/(nnd-1);        % node to node distance
node        = (0:delL:L)';      % x-coordinate of the nodes

element     = zeros(nel,nen);   % element connectivity matrix    

% generate element connectivity
for i = 1:nel
    if nen == 2
        element(i,:)  = [i i+1];
    elseif nen==3
        element(i,:)  = [2*i-1 2*i 2*i+1];
    end
end

disp_node   = find(node==disp_loc);
react_node  = find(node==react_loc);
trac_node   = find(node==trac_loc);
[trac_elem, local_node]     = find(element == trac_node);

% writes # node and coordinate
fprintf(fid,'#nodes \n');
fprintf(fid,'%d\n',nnd);
fprintf(fid,'%.4f\n',node);

% writes essential BC and location
fprintf(fid,'#displacement BC \n');
fprintf(fid,'#node  -  value \n');
fprintf(fid,'%d %10.4f\n',disp_node,disp_BC);

% writes node to calculate reaction force
fprintf(fid,'#reaction force at node \n');
fprintf(fid,'%d\n',react_node);

% writes element connectivity
fprintf(fid,'#elements \n');
fprintf(fid,'%d\n',nel);
fprintf(fid,'%d\n',nen);
if nen == 2
    fprintf(fid,'%d %6d\n',element.');
elseif nen == 3
    fprintf(fid,'%d %6d %6d\n',element.');
else
    fprintf('Invalid element type. Input file not correct');
end

% writes # integration pts
fprintf(fid,'#integration points \n');
fprintf(fid,'%d\n',nip);

% writes the element, it's locatl node and traction value
fprintf(fid,'#traction BC \n');
fprintf(fid,'#element - local_node - value \n');
fprintf(fid,'%d %6d %10.4f\n',trac_elem,local_node,trac_BC);

% writes body force constant
fprintf(fid,'#body force parameter q \n');
fprintf(fid,'%.2f\n',qparam);

% writes young's modulus
fprintf(fid,'Young''s Modulus\n');
fprintf(fid,'%.2f',Eparam);
end