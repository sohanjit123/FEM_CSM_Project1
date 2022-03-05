

function [node, element, elemType, nel, nen, nip, nnd, Eparam, Force_Node, qparam, disp_node, disp_BC, trac_el, trac_BC]=Read_input(filename)

% you can remove the two cd ../ line if your files are in the same
% directories
%cd ../file
%filename
fid=fopen(filename, 'r');
%cd ../Solution


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

