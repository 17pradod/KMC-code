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

for i = 1:length(FCCLattice)
    plot3(FCCLattice(i,1),FCCLattice(i,2),FCCLattice(i,3), 'o', 'MarkerSize',20, 'MarkerEdgeColor', 'green', 'MarkerFaceColor', 'green')
    hold on
end
axis equal 
title('Square FCC Lattice')
pause

hold off
%Divide into two different elements 

for i = 1:length(FCCLattice) 
    if FCCLattice(i,3) >= (a*V*.5)
        plot3(FCCLattice(i,1),FCCLattice(i,2),FCCLattice(i,3), 'o', 'MarkerSize',20, 'MarkerEdgeColor', 'red', 'MarkerFaceColor', 'red')
        FCCLattice(i,4) =  2;
        hold on
    else
        plot3(FCCLattice(i,1),FCCLattice(i,2),FCCLattice(i,3), 'o', 'MarkerSize',20, 'MarkerEdgeColor', 'blue', 'MarkerFaceColor', 'blue')
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
title('Square FCC Lattice Divided into Two Materials') 
axis equal 

hold off