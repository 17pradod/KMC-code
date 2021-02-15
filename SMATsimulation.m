%Diego Prado 
% A kintetic monte carlo simulation for mechanical alloying with a randomly
% jumping vacany. 

%FCCLattice
% "a" = Lattice Parameter 
% "V" = number of units on each side
a = .5;
V = 3;

%Experimental Values 
T = 300;                       %temp   
steps = 100;                    %number of steps    
sigma = -0.05533;              %function of interaction energy 
k = 8.62*10^-5;                %Boltzmann 
w_shear = 10^3;                %Shearing Frequency 
TotalTime = 0;
chi = 2.4*10^-22;
delta = .1;
%Vacancy Values
w_attempt = 10^15;             %Frequency of vacancy jump attempts 
Eo = .8;                       %Config. Indep. energy

%w_VacJumps = (w_attempt)*exp(-(Eo-Exv)/(K*T))

%Atomic positions of FCC cells in the bulk and at every edge
%The last valur will be used to identify atom type 

FCCatoms1 = [0 0 0 0;(a/2) (a/2) 0 0;(a/2) 0 (a/2) 0;0 (a/2) (a/2) 0;a 0 0 0; 0 a 0 0; 0 0 a 0;a a 0 0; a 0 a 0; 0 a a 0;a a a 0;a (a/2) (a/2) 0;(a/2) a (a/2) 0;(a/2) (a/2) a 0];
FCCatoms2 = [(a/2) 0 (a/2) 0;0 (a/2) (a/2) 0; 0 0 a 0; a 0 a 0; 0 a a 0;a a a 0;a (a/2) (a/2) 0;(a/2) a (a/2) 0;(a/2) (a/2) a 0];
FCCatoms3 = [(a/2) (a/2) 0 0;0 (a/2) (a/2) 0; 0 a 0 0;a a 0 0; 0 a a 0;a a a 0;a (a/2) (a/2) 0;(a/2) a (a/2) 0;(a/2) (a/2) a 0];
FCCatoms4 = [0 (a/2) (a/2) 0; 0 a a 0;a a a 0;a (a/2) (a/2) 0;(a/2) a (a/2) 0;(a/2) (a/2) a 0];
FCCatoms5 = [(a/2) (a/2) 0 0;(a/2) 0 (a/2) 0;a 0 0 0; a a 0 0; a 0 a 0; a a a 0;a (a/2) (a/2) 0;(a/2) a (a/2) 0;(a/2) (a/2) a 0];
FCCatoms6 = [(a/2) 0 (a/2) 0; a 0 a 0;a a a 0;a (a/2) (a/2) 0;(a/2) a (a/2) 0;(a/2) (a/2) a 0];
FCCatoms7 = [(a/2) (a/2) 0 0;a a 0 0;a a a 0;a (a/2) (a/2) 0;(a/2) a (a/2) 0;(a/2) (a/2) a 0];
FCCatoms8 = [a a a 0;a (a/2) (a/2) 0;(a/2) a (a/2) 0;(a/2) (a/2) a 0];



FCCLattice = [];

n = 0;
%Create 3 For loops to mix the positions between xyz 
for x = 1:V 
    for y = 1:V 
        for z = 1:V 
            if x==1 && y==1 && z==1 
                for i=1:length(FCCatoms1) 
                    coordinatestranslation = a*[x y z 0]; 
                    n = n+1; 
                    FCCLattice(n,:) =  FCCatoms1(i,:); 
                end 
            elseif x==1 && y==1 && z>1 
                for i=1:length(FCCatoms2) 
                    coordinatestranslation = a*[x y z 0]-a; 
                    n = n+1; 
                    FCCLattice(n,:) = coordinatestranslation +FCCatoms2(i,:);
                end
            elseif x==1 && y>1 && z==1 
                for i=1:length(FCCatoms3)
                    coordinatestranslation = a*[x y z 0]-a;
                    n = n+1;
                    FCCLattice(n,:) = coordinatestranslation +FCCatoms3(i,:);
                end
            elseif x==1 && y>1 && z>1
                for i=1:length(FCCatoms4)
                    coordinatestranslation = a*[x y z 0]-a;
                    n = n+1;
                    FCCLattice(n,:) = coordinatestranslation +FCCatoms4(i,:);
                end    
            elseif x>1 && y==1 && z==1
                for i=1:length(FCCatoms5)
                    coordinatestranslation = a*[x y z 0]-a;
                    n = n+1;
                    FCCLattice(n,:) = coordinatestranslation +FCCatoms5(i,:);
                end   
            elseif x>1 && y==1 && z>1
                for i=1:length(FCCatoms6)
                    coordinatestranslation = a*[x y z 0]-a;
                    n = n+1;
                    FCCLattice(n,:) = coordinatestranslation +FCCatoms6(i,:);
                end 
            elseif x>1 && y>1 && z==1
                for i=1:length(FCCatoms7)
                    coordinatestranslation = a*[x y z 0]-a;
                    n = n+1;
                    FCCLattice(n,:) = coordinatestranslation +FCCatoms7(i,:);
                end 
            else
                for i=1:length(FCCatoms8)
                    coordinatestranslation = a*[x y z 0]-a;
                    n = n+1;
                    FCCLattice(n,:) = coordinatestranslation +FCCatoms8(i,:);
                end
            end
        end
    end
end

%Divide into two different elements 

for i = 1:length(FCCLattice) 
    if FCCLattice(i,3) >= (a*V*.5)
        %plot3(FCCLattice2(i,1),FCCLattice2(i,2),FCCLattice2(i,3), 'o', 'MarkerSize',20, 'MarkerEdgeColor', 'red', 'MarkerFaceColor', 'red')
        FCCLattice(i,4) =  2;
        hold on
    else
       % plot3(FCCLattice2(i,1),FCCLattice2(i,2),FCCLattice2(i,3), 'o', 'MarkerSize',20, 'MarkerEdgeColor', 'blue', 'MarkerFaceColor', 'blue')
        FCCLattice(i,4) = 1;
        hold on
    end 
end

%introduce one random vacancy to the system 
index = randi([1 length(FCCLattice)]);
Vac = FCCLattice(index,:);              %matrix row of vacancy atom 
Vac1 = FCCLattice(index,:);
VacColor = Vac(4);
FCCLattice(index,:) = [FCCLattice(index,1:3) 0];
hold off

%Plot new lattice with vacancy 
for i = 1:length(FCCLattice)
    if FCCLattice(i,4) == 1
        plot3(FCCLattice(i,1),FCCLattice(i,2),FCCLattice(i,3), 'o', 'MarkerSize',20, 'MarkerEdgeColor', 'red', 'MarkerFaceColor', 'red')
        hold on
    elseif FCCLattice(i,4) == 2 
        plot3(FCCLattice(i,1),FCCLattice(i,2),FCCLattice(i,3), 'o', 'MarkerSize',20, 'MarkerEdgeColor', 'blue', 'MarkerFaceColor', 'blue')
        hold on
    else
        plot3(FCCLattice(i,1),FCCLattice(i,2),FCCLattice(i,3), 'o', 'MarkerSize',20, 'MarkerEdgeColor', 'yellow', 'MarkerFaceColor', 'yellow')
        hold on
    end 
end
axis equal

title('Initial State with One Vacancy')  
pause;

%%%%%%%%%%%%%%% Run Experiment %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%array of possible jumps for vacancy
Possible_Jumps = []; indices = []; nn = 0;
for n = 1:length(FCCLattice) %Goes through the array of atoms
    x = Vac(1);  %X, Y, and Z values of the vacancy 
    y = Vac(2);
    z = Vac(3);
    xmin = Vac(1) - a/2;  %Finds the max and min X, Y, and Z values that all nearest neighbors could have around the vacancy 
    xmax = Vac(1) + a/2;
    ymin = Vac(2) - a/2;
    ymax = Vac(2) + a/2;
    zmin = Vac(3) - a/2;
    zmax = Vac(3) + a/2;
    if FCCLattice(n,1)>=xmin & FCCLattice(n,1)<=xmax  %If the atom is close enough that it has to be a nearest neighbor
        if FCCLattice(n,2)>=ymin & FCCLattice(n,2)<=ymax
            if FCCLattice(n,3)>=zmin & FCCLattice(n,3)<=zmax
                if FCCLattice(n,1)==x & FCCLattice(n,2)==y & FCCLattice(n,3)==z %If the chosen atoms is not the vacancy 
                else
                    nn = nn + 1;
                    Possible_Jumps(nn,:) = FCCLattice(n,:); %add the information about the atom to the array of possible jumps
                    indices(nn) = n; %Gets the indices of the atoms in the main array 
                end
            end
        end
    end
end
    
  
%number of steps in the experiment 
for j = 1:steps    
%%%%%%%% Time algorithm to decide if this step is a jump or a shearing event
    
    %Jump Frequency 
    numopp = 0; len = size(Possible_Jumps);
    for q = 1:len(1)                      %Go through the array of possible jumps 
        x = Possible_Jumps(q,1);          %for each possible jump, these are the coordinates of all it's nearest neighbors
        y = Possible_Jumps(q,2);
        z = Possible_Jumps(q,3);
        xmin = Possible_Jumps(q,1) - a/2;
        xmax = Possible_Jumps(q,1) + a/2;
        ymin = Possible_Jumps(q,2) - a/2;
        ymax = Possible_Jumps(q,2) + a/2;
        zmin = Possible_Jumps(q,3) - a/2;
        zmax = Possible_Jumps(q,3) + a/2;
        for qq = 1:length(FCCLattice)                               %Go through the array of atoms and 
            if FCCLattice(qq,1)>=xmin & FCCLattice(qq,1)<=xmax      %find all the nearest neighbors for
                if FCCLattice(qq,2)>=ymin & FCCLattice(qq,2)<=ymax  % each atom in possible jumps
                    if FCCLattice(qq,3)>=zmin & FCCLattice(qq,3)<=zmax
                        if FCCLattice(qq,1)==x & FCCLattice(qq,2)==y & FCCLattice(qq,3)==z
                        else
                            if Possible_Jumps(q,4) == FCCLattice(qq,4)
                            else
                                numopp = numopp +1;          %finds n, the number of nearest neighbors that are of the oppsoits kind
                            end
                        end
                    end
                end
            end
        end
    end

Exv = (sigma*numopp)/2;
w_VacJumps = (w_attempt)*exp(-(Eo+Exv)/(k*T));
t_jump = 1/(w_VacJumps);
    
%Residence time algorithm 
                                            
% Rj = w_VacJumps + w_shear;
% ran = randi([1 round(Rj)]);  
% % if ran <= w_shear
% %     mech = 1;            This determines if a jump or shearing event
% % else                     occurs, commented for right now
% %     mech = 2;
% % end 

numevent = 0; 
%%%%%%% Shearing Events %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%select random slip system 
if mech ==1
   FCCLattice = Shearing(FCCLattice, numevent, a, V);
end



hold off

%%%%%%VACANCY MOVEMENT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Switch Vacany with selected atom 

if mech == 2
    numevent = numevent +1;
    Possible_Jumps = []; indices = []; nn = 0;
for n = 1:length(FCCLattice)                      %again finds n, the number of nearst neighbors of the opposite type
    x = Vac(1);                                   %but for the new vacancy after it jumps
    y = Vac(2);
    z = Vac(3);
    xmin = Vac(1) - a/2;
    xmax = Vac(1) + a/2;
    ymin = Vac(2) - a/2;
    ymax = Vac(2) + a/2;
    zmin = Vac(3) - a/2;
    zmax = Vac(3) + a/2;
    if FCCLattice(n,1)>=xmin & FCCLattice(n,1)<=xmax 
        if FCCLattice(n,2)>=ymin & FCCLattice(n,2)<=ymax
            if FCCLattice(n,3)>=zmin & FCCLattice(n,3)<=zmax
                if FCCLattice(n,1)==x & FCCLattice(n,2)==y & FCCLattice(n,3)==z
                else
                    nn = nn + 1;
                    Possible_Jumps(nn,:) = FCCLattice(n,:);
                    indices(nn) = n;
                end
            end
        end
    end
end

len = size(Possible_Jumps);       %Gets the length of the array of possible jumps
if len == 0                       %If the array is zero for some reason then nothing happens
else
    RandAtom = randi([1 len(1)]);          %a random number between 1 and the number of possible jumps  
    VacColor = Possible_Jumps(RandAtom,4); %New vacancy type
    FCCLattice(index,4) = VacColor;        %replace current vacany with color of new vacancy
    index = indices(RandAtom); 
    Vac = FCCLattice(index,:);
    FCCLattice(index,4) = 0;
end
    for i = 1:length(FCCLattice)
        if FCCLattice(i,4) == 1
            plot3(FCCLattice(i,1),FCCLattice(i,2),FCCLattice(i,3), 'o', 'MarkerSize',20, 'MarkerEdgeColor', 'red', 'MarkerFaceColor', 'red')
            hold on
        elseif FCCLattice(i,4) == 2 
            plot3(FCCLattice(i,1),FCCLattice(i,2),FCCLattice(i,3), 'o', 'MarkerSize',20, 'MarkerEdgeColor', 'blue', 'MarkerFaceColor', 'blue')
            hold on
        else
            plot3(FCCLattice(i,1),FCCLattice(i,2),FCCLattice(i,3), 'o', 'MarkerSize',20, 'MarkerEdgeColor', 'yellow', 'MarkerFaceColor', 'yellow')
            hold on
        end 
    end
axis equal
end
end