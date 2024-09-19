# Beam_variable
clc; clear; lengthwo_sup=input('length of beam w/o support=  ' );
a=lengthwo_sup/0.01+1;
sup_vector=zeros(a,1);
fprintf('Number of supports')
supp=input('no.' );
sup_matrix=zeros(supp,2);

for new=1:supp
    sup_matrix(new,1)=input('1.Restrained (disp +rotation ) \n 2.Restrained(disp)');
    sup_matrix(new,2)=input('location present');
end

start_length=0;

for i=1:a
    for k=1:supp
        if sup_matrix(k,1)==1 && abs(sup_matrix(k,2)-start_length)<1e-06
            sup_vector(i,1)=2;
        elseif sup_matrix(k,1)==2 && abs(sup_matrix(k,2)-start_length)<1e-06
            sup_vector(i,1)=1;
        end
    end
    
        start_length=start_length+0.01;
end

supd_vector=zeros(2*a,1)

start=1;
modified_indices = [];

for j=1:2:2*a
         if sup_vector(start,1)==2
           
               modified_indices = [modified_indices; j; j+1];
               
         elseif sup_vector(start,1)==1
             
                modified_indices = [modified_indices; j];
                
                
         elseif sup_vector(start,1)==0
             
         end
   start=start+1;
end
% 
% supd_vector(modified_indices,:)=[];
% x1=[1;2]; x2=[2,3]; y=x1/x2;
double_n=2*a;
dx=0.01;
E = 210e9; % Pascal
I = 1e-6; 
K_elem = (E * I / dx^3) * [12 6 * dx -12 6 * dx;
                           6 * dx 4 * dx^2 -6 * dx 2 * dx^2;
                           -12 -6 * dx 12 -6 * dx;
                           6 * dx 2 * dx^2 -6 * dx 4 * dx^2];
                       
                       
K_global = zeros(2 * a);

for i = 1:a-1
    K_global(2*i-1:2*i+2, 2*i-1:2*i+2) = K_global(2*i-1:2*i+2, 2*i-1:2*i+2) + K_elem;
end
Kred=K_global;
Kred(modified_indices,:)=[];
Kred(:,modified_indices)=[];

d = input('Enter number of loads in beam specimen: ');
start_node = 0.01;
input_matrix = zeros(d, 2);


for i = 1:d
    input_matrix(i, 1) = input('Enter load in Newton: ');
    input_matrix(i, 2) = input('Enter position from left (in meters): ');
end

F_global = zeros(2*a, 1);

  for i=1:2*a
        F_global(i,1)=0;
     
  end
    
j=0;

for i=1:2:2*a
     
    for k=1:d
%            disp('entering');
      if (abs(input_matrix(k,2)-j)<(1e-06))
       
            F_global(i,1)=input_matrix(k,1);
            
        end
    end
      j=j+0.01;
end

disp(F_global);

Fred=F_global;
Fred(modified_indices,:)=[];

displacements = Kred \ Fred;

disp_full = zeros(2*a, 1);


free_dofs = setdiff(1:2*a, modified_indices); 


if length(free_dofs) == length(displacements)
    disp_full(free_dofs) = displacements;
end

vertical_disp=disp_full(1:2:end);


