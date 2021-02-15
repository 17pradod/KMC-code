function [NewLattice] = Shearing(FCCLattice, numevent, a, V)
numevent = numevent +1;
sysnum = randi([1 12]);
%An array of all possible slip planes
SlipSys = [1 1 1 1 -1 0; 1 1 1 0 1 -1; 1 1 1 1 0 -1; -1 1 1 1 1 0; -1 1 1 0 1 -1; -1 1 1 1 0 1; 1 -1 1 1 1 0; 1 -1 1 0 1 1; 1 -1 1 1 0 -1; 1 1 -1 1 -1 0; 1 1 -1 0 1 1; 1 1 -1 1 0 1];
Sys = SlipSys(sysnum,:); %select the proper system 

%After selecting a system, create arrays of all of the points in each
%parallel plane 
p = 0;
w = 1;
SlipPlanes = [];  %Array with 4 columns for each plane, X, Y, Z, and atom type 
    if sysnum == 1     
        x = [a 0 0];
        y = [0 a 0];
        z = [0 0 a];
        Vec1 = y-x;
        Vec2 = z-x;
        C = cross(Vec1,Vec2);
        Plane = [C(1) C(2) C(3) C(1)*-x(1)+C(2)*-x(2)+C(3)*-x(3)];
        %Slip Planes is an array where every 4 columns are the X Y Z and
        %color value of the atoms in each plane
        StuffToAdd = [];       %All of the constants that determine the different parallel planes 
        for m = -Plane(4):Plane(4):-V*a*2*2  %scan through all of the planes parallel to the chosen one
            for ppp = 1:length(FCCLattice)   %Run through the main lattice array
                if FCCLattice(ppp,1)*Plane(1)+ FCCLattice(ppp,2)*Plane(2) + FCCLattice(ppp,3)*Plane(3) + (Plane(4) + m) == 0 %see if that atom is on the current plane
                    p = p+1;    %Counter to place the next atom in the plane in the next row of the array
                    SlipPlanes(p,w:w+3) = FCCLattice(ppp,:);  %Put the X, Y, Z, and atom type in the array
                    StuffToAdd = [StuffToAdd (Plane(4) + m)]; %Puts the constants for each plane in an array 
                    StuffToAdd = unique(StuffToAdd); %Gets rid of all of the repeated values so there are the 
                end                                  %same number of constants as planes with atoms in them
            end
            w = w+4;
        end
        
        Dims = size(SlipPlanes);           %dimenstions of the total array with the slip planes
        NumsPerPlane = zeros(1,3*Dims(2)/4); %This array will have the Naa, Nbb, and Nab for each plane
        w3 = 0;     %Various counters 
        w1 = 1;
        w2 = 1;
        w4 = 4;
        for ii= 1:Dims(2)/4  %Loop form 1 to the number of planes
            w3 = w3 + 4;
            CurrentPlane = []; %This will have just the atoms in each plane
            for o = 1:Dims(1)
                if SlipPlanes(o,w4)== 1 | SlipPlanes(o,w4)== 2
                    CurrentPlane = [CurrentPlane; SlipPlanes(o,w4-3:w4)];
                end
            end
            w4 = w4+4;
            numOpp = 0;   % Nab per plane
            numSam1 = 0;  % Naa per plane
            numSam2 = 0;  % Nbb per plane
            DimsN = size(CurrentPlane); %dimensions of the array with just the atoms in a plane
            for iii = 1:DimsN(1)           %This is going through each atom in the plane
                for iiii = 1:DimsN(1)      %For each atom, it is compared to every other atom in the plane
                    dist = sqrt((0-0)^2+(a/2 - 0)^2+(a/2 - a)^2)+ a/10;
                    if sqrt((CurrentPlane(iii,1)-CurrentPlane(iiii,1))^2+(CurrentPlane(iii,2)-CurrentPlane(iiii,2))^2+(CurrentPlane(iii,3)-CurrentPlane(iiii,3))^2) <= dist  %if distace between the atoms is equal to or less than the nearest neighbor distance 
                        if sqrt((CurrentPlane(iii,1)-CurrentPlane(iiii,1))^2+(CurrentPlane(iii,2)-CurrentPlane(iiii,2))^2+(CurrentPlane(iii,3)-CurrentPlane(iiii,3))^2) ~= 0 %If the atom is not itself
                            if CurrentPlane(iii,4) == 1 & CurrentPlane(iiii,4) == 1
                                numSam1 = numSam1 + 1;
                            elseif CurrentPlane(iii,4) == 2 & CurrentPlane(iiii,4) == 2
                                numSam2 = numSam2 + 1;
                            else 
                                numOpp = numOpp + 1;
                            end
                        end
                    end
                end
            end
                NumsPerPlane(w1,w2) = NumsPerPlane(w1,w2) + numSam1;    %Assigns the number of aa, bb, and ab bonds to the array set up above
                NumsPerPlane(w1,w2+1) = NumsPerPlane(w2+1) +  numSam2;    
                NumsPerPlane(w1,w2+2) = NumsPerPlane(w2+2) + numOpp;
                w2 = w2 + 3;
        end
        R = [];
        jj = 1;
            for j = 1:Dims(2)/4         %Plugs the values for each plane into the equation in the paper
                Naa = NumsPerPlane(jj);
                Nbb = NumsPerPlane(jj+1);
                Nab = NumsPerPlane(jj+2);
                R(j) = exp(-.0001*(Naa+Nbb*(1-.001)+Nab*(1-.001/2))/1);
                jj = jj+3;
            end
            
            Sorted = sort(R);          %Sorts all the R values base on value
            Maxes = Sorted(end-1:end); %Selects the largest two (In my code that will always the corners :(
            Plane1 = 1;
            Plane2 = 1;
            num = 1;
            for jjj = 1:Dims(2)/4            %Finds the plane that the two largest values belong to
                if  R(jjj) == Maxes(1)
                    Plane1 = num;
                elseif R(jjj) == Maxes(2)
                    Plane2 = num;
                end
                num = num + 1;
            end
            
            NewLattice = FCCLattice;
            ww = 2;                            % Since it the corners never move it is easier to start at 2 and skip it
            Random = randi([1 3]);             % The reason this is between one and 3 is because randi([1 2]) for some reason kept giving me 1 a lot more than 2
            for jjjj = 1:length(StuffToAdd)-2  % Starts after the first plane and ends before the last, hence the -2
                for kkkk = 1:Dims(1)           % Runs though the main array of atoms
                    if FCCLattice(kkkk,1)*Plane(1) + FCCLattice(kkkk,2)*Plane(2) + FCCLattice(kkkk,3)*Plane(3) + StuffToAdd(ww) == 0 %if the atoms is on the planes between the selected two 
                        if Random == 1         % This determines weather the atoms will be sheared in the burgers vector direction or in the negative burgers vector direction 
                            NewLattice(kkkk,1) = NewLattice(kkkk,1) + a/2;
                            NewLattice(kkkk,2) = NewLattice(kkkk,2) - a/2;   %Sum of a/2 in one direction and a/2 in a perpendiculr direction gives the burgers vector 
                            NewLattice(kkkk,3) = NewLattice(kkkk,3) + 0;
                        else
                            NewLattice(kkkk,1) = NewLattice(kkkk,1) - a/2;
                            NewLattice(kkkk,2) = NewLattice(kkkk,2) + a/2;
                            NewLattice(kkkk,3) = NewLattice(kkkk,3) + 0;
                        end
                    end
                end
                ww = ww +1;
            end
            clf
            for i = 1:length(NewLattice)
                if FCCLattice(i,4) == 1
                    plot3(NewLattice(i,1),NewLattice(i,2),NewLattice(i,3), 'o', 'MarkerSize',20, 'MarkerEdgeColor', 'red', 'MarkerFaceColor', 'red')
                    hold on
                elseif FCCLattice(i,4) == 2 
                    plot3(NewLattice(i,1),NewLattice(i,2),NewLattice(i,3), 'o', 'MarkerSize',20, 'MarkerEdgeColor', 'blue', 'MarkerFaceColor', 'blue')
                    hold on
                else
                    plot3(NewLattice(i,1),NewLattice(i,2),NewLattice(i,3), 'o', 'MarkerSize',20, 'MarkerEdgeColor', 'yellow', 'MarkerFaceColor', 'yellow')
                    hold on
                end 
            end
            axis equal
        
    elseif sysnum == 2
        x = [a 0 0];
        y = [0 a 0];
        z = [0 0 a];
        Vec1 = y-x;
        Vec2 = z-x;
        C = cross(Vec1,Vec2);
        Plane = [C(1) C(2) C(3) C(1)*-x(1)+C(2)*-x(2)+C(3)*-x(3)]; %(Plane1 + Plane2 +Plane3 = -Plane(4))
        StuffToAdd = [];
        for m = -Plane(4):Plane(4):-V*a*2*2
            for ppp = 1:length(FCCLattice)
                if FCCLattice(ppp,1)*Plane(1)+ FCCLattice(ppp,2)*Plane(2) + FCCLattice(ppp,3)*Plane(3) + (Plane(4) + m) == 0
                    p = p+1;
                    SlipPlanes(p,w:w+3) = FCCLattice(ppp,:);
                    StuffToAdd = [StuffToAdd (Plane(4) + m)];
                    StuffToAdd = unique(StuffToAdd);
                end
            end
            w = w+4;
        end
        
        Dims = size(SlipPlanes);
        NumsPerPlane = zeros(1,3*Dims(2)/4);
        w3 = 0;
        w1 = 1;
        w2 = 1;
        w4 = 4;
        for ii= 1:Dims(2)/4 
            w3 = w3 + 4;
            CurrentPlane = [];
            for o = 1:Dims(1)
                if SlipPlanes(o,w4)== 1 | SlipPlanes(o,w4)== 2
                    CurrentPlane = [CurrentPlane; SlipPlanes(o,w4-3:w4)];
                end
            end
            w4 = w4+4;
            numOpp = 0;
            numSam1 = 0;
            numSam2 = 0;
            DimsN = size(CurrentPlane);
            for iii = 1:DimsN(1)
                for iiii = 1:DimsN(1)
                    dist = sqrt((0-0)^2+(a/2 - 0)^2+(a/2 - a)^2)+ a/10;
                    if sqrt((CurrentPlane(iii,1)-CurrentPlane(iiii,1))^2+(CurrentPlane(iii,2)-CurrentPlane(iiii,2))^2+(CurrentPlane(iii,3)-CurrentPlane(iiii,3))^2) <= dist
                        if sqrt((CurrentPlane(iii,1)-CurrentPlane(iiii,1))^2+(CurrentPlane(iii,2)-CurrentPlane(iiii,2))^2+(CurrentPlane(iii,3)-CurrentPlane(iiii,3))^2) ~= 0
                            if CurrentPlane(iii,4) == 1 & CurrentPlane(iiii,4) == 1
                                numSam1 = numSam1 + 1;
                            elseif CurrentPlane(iii,4) == 2 & CurrentPlane(iiii,4) == 2
                                numSam2 = numSam2 + 1;
                            else 
                                numOpp = numOpp + 1;
                            end
                        end
                    end
                end
            end
                NumsPerPlane(w1,w2) = NumsPerPlane(w1,w2) + numSam1;
                NumsPerPlane(w1,w2+1) = NumsPerPlane(w2+1) +  numSam2;
                NumsPerPlane(w1,w2+2) = NumsPerPlane(w2+2) + numOpp;
                w2 = w2 + 3;
        end
        R = [];
        jj = 1;
            for j = 1:Dims(2)/4
                Naa = NumsPerPlane(jj);
                Nbb = NumsPerPlane(jj+1);
                Nab = NumsPerPlane(jj+2);
                R(j) = exp(-.0001*(Naa+Nbb*(1-.001)+Nab*(1-.001/2))/1);
                jj = jj+3;
            end
            
            Sorted = sort(R);
            Maxes = Sorted(end-1:end);
            Plane1 = 1;
            Plane2 = 1;
            num = 1;
            for jjj = 1:Dims(2)/4
                if  R(jjj) == Maxes(1)
                    Plane1 = num;
                elseif R(jjj) == Maxes(2)
                    Plane2 = num;
                end
                num = num + 1;
            end
            
            NewLattice = FCCLattice;
            ww = 2;
            Random = randi([1 3]);
            for jjjj = 1:length(StuffToAdd)-2
                for kkkk = 1:Dims(1)
                    if FCCLattice(kkkk,1)*Plane(1) + FCCLattice(kkkk,2)*Plane(2) + FCCLattice(kkkk,3)*Plane(3) + StuffToAdd(ww) == 0
                        if Random == 1
                            NewLattice(kkkk,1) = NewLattice(kkkk,1) + 0;
                            NewLattice(kkkk,2) = NewLattice(kkkk,2) + a/2;
                            NewLattice(kkkk,3) = NewLattice(kkkk,3) - a/2;
                        else
                            NewLattice(kkkk,1) = NewLattice(kkkk,1) + 0;
                            NewLattice(kkkk,2) = NewLattice(kkkk,2) - a/2;
                            NewLattice(kkkk,3) = NewLattice(kkkk,3) + a/2;
                        end
                        
                    end
                end
                ww = ww +1;
            end
            clf
            for i = 1:length(NewLattice)
                if FCCLattice(i,4) == 1
                    plot3(NewLattice(i,1),NewLattice(i,2),NewLattice(i,3), 'o', 'MarkerSize',20, 'MarkerEdgeColor', 'red', 'MarkerFaceColor', 'red')
                    hold on
                elseif FCCLattice(i,4) == 2 
                    plot3(NewLattice(i,1),NewLattice(i,2),NewLattice(i,3), 'o', 'MarkerSize',20, 'MarkerEdgeColor', 'blue', 'MarkerFaceColor', 'blue')
                    hold on
                else
                    plot3(NewLattice(i,1),NewLattice(i,2),NewLattice(i,3), 'o', 'MarkerSize',20, 'MarkerEdgeColor', 'yellow', 'MarkerFaceColor', 'yellow')
                    hold on
                end 
            end
            axis equal
    elseif sysnum == 3
        x = [a 0 0];
        y = [0 a 0];
        z = [0 0 a];
        Vec1 = y-x;
        Vec2 = z-x;
        C = cross(Vec1,Vec2);
        Plane = [C(1) C(2) C(3) C(1)*-x(1)+C(2)*-x(2)+C(3)*-x(3)];
        StuffToAdd = [];
        for m = -Plane(4):Plane(4):-V*a*2*2
            for ppp = 1:length(FCCLattice)
                if FCCLattice(ppp,1)*Plane(1)+ FCCLattice(ppp,2)*Plane(2) + FCCLattice(ppp,3)*Plane(3) + (Plane(4) + m) == 0
                    p = p+1;
                    SlipPlanes(p,w:w+3) = FCCLattice(ppp,:);
                    StuffToAdd = [StuffToAdd (Plane(4) + m)];
                    StuffToAdd = unique(StuffToAdd);
                end
            end
            w = w+4;
        end
        
        Dims = size(SlipPlanes);
        NumsPerPlane = zeros(1,3*Dims(2)/4);
        w3 = 0;
        w1 = 1;
        w2 = 1;
        w4 = 4;
        for ii= 1:Dims(2)/4 
            w3 = w3 + 4;
            CurrentPlane = [];
            for o = 1:Dims(1)
                if SlipPlanes(o,w4)== 1 | SlipPlanes(o,w4)== 2
                    CurrentPlane = [CurrentPlane; SlipPlanes(o,w4-3:w4)];
                end
            end
            w4 = w4+4;
            numOpp = 0;
            numSam1 = 0;
            numSam2 = 0;
            DimsN = size(CurrentPlane);
            for iii = 1:DimsN(1)
                for iiii = 1:DimsN(1)
                    dist = sqrt((0-0)^2+(a/2 - 0)^2+(a/2 - a)^2)+ a/10;
                    if sqrt((CurrentPlane(iii,1)-CurrentPlane(iiii,1))^2+(CurrentPlane(iii,2)-CurrentPlane(iiii,2))^2+(CurrentPlane(iii,3)-CurrentPlane(iiii,3))^2) <= dist
                        if sqrt((CurrentPlane(iii,1)-CurrentPlane(iiii,1))^2+(CurrentPlane(iii,2)-CurrentPlane(iiii,2))^2+(CurrentPlane(iii,3)-CurrentPlane(iiii,3))^2) ~= 0
                            if CurrentPlane(iii,4) == 1 & CurrentPlane(iiii,4) == 1
                                numSam1 = numSam1 + 1;
                            elseif CurrentPlane(iii,4) == 2 & CurrentPlane(iiii,4) == 2
                                numSam2 = numSam2 + 1;
                            else 
                                numOpp = numOpp + 1;
                            end
                        end
                    end
                end
            end
                NumsPerPlane(w1,w2) = NumsPerPlane(w1,w2) + numSam1;
                NumsPerPlane(w1,w2+1) = NumsPerPlane(w2+1) +  numSam2;
                NumsPerPlane(w1,w2+2) = NumsPerPlane(w2+2) + numOpp;
                w2 = w2 + 3;
        end
        R = [];
        jj = 1;
            for j = 1:Dims(2)/4
                Naa = NumsPerPlane(jj);
                Nbb = NumsPerPlane(jj+1);
                Nab = NumsPerPlane(jj+2);
                R(j) = exp(-.0001*(Naa+Nbb*(1-.001)+Nab*(1-.001/2))/1);
                jj = jj+3;
            end
            
            Sorted = sort(R);
            Maxes = Sorted(end-1:end);
            Plane1 = 1;
            Plane2 = 1;
            num = 1;
            for jjj = 1:Dims(2)/4
                if  R(jjj) == Maxes(1)
                    Plane1 = num;
                elseif R(jjj) == Maxes(2)
                    Plane2 = num;
                end
                num = num + 1;
            end
            
            NewLattice = FCCLattice;
            ww = 2;
            Random = randi([1 3]);
            for jjjj = 1:length(StuffToAdd)-2
                for kkkk = 1:Dims(1)
                    if FCCLattice(kkkk,1)*Plane(1) + FCCLattice(kkkk,2)*Plane(2) + FCCLattice(kkkk,3)*Plane(3) + StuffToAdd(ww) == 0
                        if Random == 1
                            NewLattice(kkkk,1) = NewLattice(kkkk,1) + a/2;
                            NewLattice(kkkk,2) = NewLattice(kkkk,2) + 0;
                            NewLattice(kkkk,3) = NewLattice(kkkk,3) - a/2;
                        else
                            NewLattice(kkkk,1) = NewLattice(kkkk,1) - a/2;
                            NewLattice(kkkk,2) = NewLattice(kkkk,2) + 0;
                            NewLattice(kkkk,3) = NewLattice(kkkk,3) + a/2;
                        end
                    end
                end
                ww = ww +1;
            end
            clf
            for i = 1:length(NewLattice)
                if FCCLattice(i,4) == 1
                    plot3(NewLattice(i,1),NewLattice(i,2),NewLattice(i,3), 'o', 'MarkerSize',20, 'MarkerEdgeColor', 'red', 'MarkerFaceColor', 'red')
                    hold on
                elseif FCCLattice(i,4) == 2 
                    plot3(NewLattice(i,1),NewLattice(i,2),NewLattice(i,3), 'o', 'MarkerSize',20, 'MarkerEdgeColor', 'blue', 'MarkerFaceColor', 'blue')
                    hold on
                else
                    plot3(NewLattice(i,1),NewLattice(i,2),NewLattice(i,3), 'o', 'MarkerSize',20, 'MarkerEdgeColor', 'yellow', 'MarkerFaceColor', 'yellow')
                    hold on
                end 
            end
            axis equal
    elseif sysnum == 4
        x = [-a 0 0];
        y = [0 a 0];
        z = [0 0 a];
        Vec1 = y-x;
        Vec2 = z-x;
        C = cross(Vec1,Vec2);
        Plane = [C(1) C(2) C(3) C(1)*-x(1)+C(2)*-x(2)+C(3)*-x(3)];
        StuffToAdd = [];
        for m = -V*a*2:Plane(4):V*a*2 
            for ppp = 1:length(FCCLattice)
                if FCCLattice(ppp,1)*Plane(1)+ FCCLattice(ppp,2)*Plane(2) + FCCLattice(ppp,3)*Plane(3) + (Plane(4) + m) == 0
                    p = p+1;
                    SlipPlanes(p,w:w+3) = FCCLattice(ppp,:);
                    StuffToAdd = [StuffToAdd (Plane(4) + m)];
                    StuffToAdd = unique(StuffToAdd);
                end
            end
            w = w+4;
        end
        another = 1;
        PlaneFix = [];
        Dims = size(SlipPlanes);
        for l = 1:Dims(2)/4
            if max(SlipPlanes(:,4*l)) > 0
                PlaneFix(:,another*4-3:another*4) = SlipPlanes(:,4*l-3:4*l);
                another = another + 1;
            end
        end
        SlipPlanes = PlaneFix;
        Dims = size(SlipPlanes);
        NumsPerPlane = zeros(1,3*Dims(2)/4);
        w3 = 0;
        w1 = 1;
        w2 = 1;
        w4 = 4;
        for ii= 1:Dims(2)/4 
            w3 = w3 + 4;
            CurrentPlane = [];
            for o = 1:Dims(1)
                if SlipPlanes(o,w4)== 1 | SlipPlanes(o,w4)== 2
                    CurrentPlane = [CurrentPlane; SlipPlanes(o,w4-3:w4)];
                end
            end
            w4 = w4+4;
            numOpp = 0;
            numSam1 = 0;
            numSam2 = 0;
            DimsN = size(CurrentPlane);
            for iii = 1:DimsN(1)
                for iiii = 1:DimsN(1)
                    dist = sqrt((0-0)^2+(a/2 - 0)^2+(a/2 - a)^2)+ a/10;
                    if sqrt((CurrentPlane(iii,1)-CurrentPlane(iiii,1))^2+(CurrentPlane(iii,2)-CurrentPlane(iiii,2))^2+(CurrentPlane(iii,3)-CurrentPlane(iiii,3))^2) <= dist
                        if sqrt((CurrentPlane(iii,1)-CurrentPlane(iiii,1))^2+(CurrentPlane(iii,2)-CurrentPlane(iiii,2))^2+(CurrentPlane(iii,3)-CurrentPlane(iiii,3))^2) ~= 0
                            if CurrentPlane(iii,4) == 1 & CurrentPlane(iiii,4) == 1
                                numSam1 = numSam1 + 1;
                            elseif CurrentPlane(iii,4) == 2 & CurrentPlane(iiii,4) == 2
                                numSam2 = numSam2 + 1;
                            else 
                                numOpp = numOpp + 1;
                            end
                        end
                    end
                end
            end
                NumsPerPlane(w1,w2) = NumsPerPlane(w1,w2) + numSam1;
                NumsPerPlane(w1,w2+1) = NumsPerPlane(w2+1) +  numSam2;
                NumsPerPlane(w1,w2+2) = NumsPerPlane(w2+2) + numOpp;
                w2 = w2 + 3;
        end
        R = [];
        jj = 1;
            for j = 1:Dims(2)/4
                Naa = NumsPerPlane(jj);
                Nbb = NumsPerPlane(jj+1);
                Nab = NumsPerPlane(jj+2);
                R(j) = exp(-.0001*(Naa+Nbb*(1-.001)+Nab*(1-.001/2))/1);
                jj = jj+3;
            end
            
            Sorted = sort(R);
            Maxes = Sorted(end-1:end);
            Plane1 = 1;
            Plane2 = 1;
            num = 1;
            for jjj = 1:Dims(2)/4
                if  R(jjj) == Maxes(1)
                    Plane1 = num;
                elseif R(jjj) == Maxes(2)
                    Plane2 = num;
                end
                num = num + 1;
            end
            
            NewLattice = FCCLattice;
            ww = 2;
            Random = randi([1 3]);
            for jjjj = 1:length(StuffToAdd)-2
                for kkkk = 1:Dims(1)
                    if FCCLattice(kkkk,1)*Plane(1) + FCCLattice(kkkk,2)*Plane(2) + FCCLattice(kkkk,3)*Plane(3) + StuffToAdd(ww) == 0
                        if Random == 1
                            NewLattice(kkkk,1) = NewLattice(kkkk,1) + a/2;
                            NewLattice(kkkk,2) = NewLattice(kkkk,2) + a/2;
                            NewLattice(kkkk,3) = NewLattice(kkkk,3) + 0;
                        else
                            NewLattice(kkkk,1) = NewLattice(kkkk,1) - a/2;
                            NewLattice(kkkk,2) = NewLattice(kkkk,2) - a/2;
                            NewLattice(kkkk,3) = NewLattice(kkkk,3) + 0;
                        end
                        
                    end
                end
                ww = ww +1;
            end
            clf
            for i = 1:length(NewLattice)
                if FCCLattice(i,4) == 1
                    plot3(NewLattice(i,1),NewLattice(i,2),NewLattice(i,3), 'o', 'MarkerSize',20, 'MarkerEdgeColor', 'red', 'MarkerFaceColor', 'red')
                    hold on
                elseif FCCLattice(i,4) == 2 
                    plot3(NewLattice(i,1),NewLattice(i,2),NewLattice(i,3), 'o', 'MarkerSize',20, 'MarkerEdgeColor', 'blue', 'MarkerFaceColor', 'blue')
                    hold on
                else
                    plot3(NewLattice(i,1),NewLattice(i,2),NewLattice(i,3), 'o', 'MarkerSize',20, 'MarkerEdgeColor', 'yellow', 'MarkerFaceColor', 'yellow')
                    hold on
                end 
            end
            axis equal
    elseif sysnum == 5
        x = [-a 0 0];
        y = [0 a 0];
        z = [0 0 a];
        Vec1 = y-x;
        Vec2 = z-x;
        C = cross(Vec1,Vec2);
        Plane = [C(1) C(2) C(3) C(1)*-x(1)+C(2)*-x(2)+C(3)*-x(3)];
        StuffToAdd = [];
        for m = -V*a*2:Plane(4):V*a*2 
            for ppp = 1:length(FCCLattice)
                if FCCLattice(ppp,1)*Plane(1)+ FCCLattice(ppp,2)*Plane(2) + FCCLattice(ppp,3)*Plane(3) + (Plane(4) + m) == 0
                    p = p+1;
                    SlipPlanes(p,w:w+3) = FCCLattice(ppp,:);
                    StuffToAdd = [StuffToAdd (Plane(4) + m)];
                    StuffToAdd = unique(StuffToAdd);
                end
            end
            w = w+4;
        end
        another = 1;
        PlaneFix = [];
        Dims = size(SlipPlanes);
        for l = 1:Dims(2)/4
            if max(SlipPlanes(:,4*l)) > 0
                PlaneFix(:,another*4-3:another*4) = SlipPlanes(:,4*l-3:4*l);
                another = another + 1;
            end
        end
        SlipPlanes = PlaneFix;
        Dims = size(SlipPlanes);
        NumsPerPlane = zeros(1,3*Dims(2)/4);
        w3 = 0;
        w1 = 1;
        w2 = 1;
        w4 = 4;
        for ii= 1:Dims(2)/4 
            w3 = w3 + 4;
            CurrentPlane = [];
            for o = 1:Dims(1)
                if SlipPlanes(o,w4)== 1 | SlipPlanes(o,w4)== 2
                    CurrentPlane = [CurrentPlane; SlipPlanes(o,w4-3:w4)];
                end
            end
            w4 = w4+4;
            numOpp = 0;
            numSam1 = 0;
            numSam2 = 0;
            DimsN = size(CurrentPlane);
            for iii = 1:DimsN(1)
                for iiii = 1:DimsN(1)
                    dist = sqrt((0-0)^2+(a/2 - 0)^2+(a/2 - a)^2)+ a/10;
                    if sqrt((CurrentPlane(iii,1)-CurrentPlane(iiii,1))^2+(CurrentPlane(iii,2)-CurrentPlane(iiii,2))^2+(CurrentPlane(iii,3)-CurrentPlane(iiii,3))^2) <= dist
                        if sqrt((CurrentPlane(iii,1)-CurrentPlane(iiii,1))^2+(CurrentPlane(iii,2)-CurrentPlane(iiii,2))^2+(CurrentPlane(iii,3)-CurrentPlane(iiii,3))^2) ~= 0
                            if CurrentPlane(iii,4) == 1 & CurrentPlane(iiii,4) == 1
                                numSam1 = numSam1 + 1;
                            elseif CurrentPlane(iii,4) == 2 & CurrentPlane(iiii,4) == 2
                                numSam2 = numSam2 + 1;
                            else 
                                numOpp = numOpp + 1;
                            end
                        end
                    end
                end
            end
                NumsPerPlane(w1,w2) = NumsPerPlane(w1,w2) + numSam1;
                NumsPerPlane(w1,w2+1) = NumsPerPlane(w2+1) +  numSam2;
                NumsPerPlane(w1,w2+2) = NumsPerPlane(w2+2) + numOpp;
                w2 = w2 + 3;
        end
        R = [];
        jj = 1;
            for j = 1:Dims(2)/4
                Naa = NumsPerPlane(jj);
                Nbb = NumsPerPlane(jj+1);
                Nab = NumsPerPlane(jj+2);
                R(j) = exp(-.0001*(Naa+Nbb*(1-.001)+Nab*(1-.001/2))/1);
                jj = jj+3;
            end
            
            Sorted = sort(R);
            Maxes = Sorted(end-1:end);
            Plane1 = 1;
            Plane2 = 1;
            num = 1;
            for jjj = 1:Dims(2)/4
                if  R(jjj) == Maxes(1)
                    Plane1 = num;
                elseif R(jjj) == Maxes(2)
                    Plane2 = num;
                end
                num = num + 1;
            end
            
            NewLattice = FCCLattice;
            ww = 2;
            Random = randi([1 3]);
            for jjjj = 1:length(StuffToAdd)-2
                for kkkk = 1:Dims(1)
                    if FCCLattice(kkkk,1)*Plane(1) + FCCLattice(kkkk,2)*Plane(2) + FCCLattice(kkkk,3)*Plane(3) + StuffToAdd(ww) == 0
                        if Random == 1
                            NewLattice(kkkk,1) = NewLattice(kkkk,1) + 0;
                            NewLattice(kkkk,2) = NewLattice(kkkk,2) + a/2;
                            NewLattice(kkkk,3) = NewLattice(kkkk,3) - a/2;
                        else
                            NewLattice(kkkk,1) = NewLattice(kkkk,1) + 0;
                            NewLattice(kkkk,2) = NewLattice(kkkk,2) - a/2;
                            NewLattice(kkkk,3) = NewLattice(kkkk,3) + a/2;
                        end
                    end
                end
                ww = ww +1;
            end
            
            clf
            for i = 1:length(NewLattice)
                if FCCLattice(i,4) == 1
                    plot3(NewLattice(i,1),NewLattice(i,2),NewLattice(i,3), 'o', 'MarkerSize',20, 'MarkerEdgeColor', 'red', 'MarkerFaceColor', 'red')
                    hold on
                elseif FCCLattice(i,4) == 2 
                    plot3(NewLattice(i,1),NewLattice(i,2),NewLattice(i,3), 'o', 'MarkerSize',20, 'MarkerEdgeColor', 'blue', 'MarkerFaceColor', 'blue')
                    hold on
                else
                    plot3(NewLattice(i,1),NewLattice(i,2),NewLattice(i,3), 'o', 'MarkerSize',20, 'MarkerEdgeColor', 'yellow', 'MarkerFaceColor', 'yellow')
                    hold on
                end 
            end
            axis equal
    elseif sysnum == 6
        x = [-a 0 0];
        y = [0 a 0];
        z = [0 0 a];
        Vec1 = y-x;
        Vec2 = z-x;
        C = cross(Vec1,Vec2);
        Plane = [C(1) C(2) C(3) C(1)*-x(1)+C(2)*-x(2)+C(3)*-x(3)];
        StuffToAdd = [];
        for m = -V*a*2:Plane(4):V*a*2 
            for ppp = 1:length(FCCLattice)
                if FCCLattice(ppp,1)*Plane(1)+ FCCLattice(ppp,2)*Plane(2) + FCCLattice(ppp,3)*Plane(3) + (Plane(4) + m) == 0
                    p = p+1;
                    SlipPlanes(p,w:w+3) = FCCLattice(ppp,:);
                    StuffToAdd = [StuffToAdd (Plane(4) + m)];
                    StuffToAdd = unique(StuffToAdd);
                end
            end
            w = w+4;
        end
        another = 1;
        PlaneFix = [];
        Dims = size(SlipPlanes);
        for l = 1:Dims(2)/4
            if max(SlipPlanes(:,4*l)) > 0
                PlaneFix(:,another*4-3:another*4) = SlipPlanes(:,4*l-3:4*l);
                another = another + 1;
            end
        end
        SlipPlanes = PlaneFix;
        Dims = size(SlipPlanes);
        NumsPerPlane = zeros(1,3*Dims(2)/4);
        w3 = 0;
        w1 = 1;
        w2 = 1;
        w4 = 4;
        for ii= 1:Dims(2)/4 
            w3 = w3 + 4;
            CurrentPlane = [];
            for o = 1:Dims(1)
                if SlipPlanes(o,w4)== 1 | SlipPlanes(o,w4)== 2
                    CurrentPlane = [CurrentPlane; SlipPlanes(o,w4-3:w4)];
                end
            end
            w4 = w4+4;
            numOpp = 0;
            numSam1 = 0;
            numSam2 = 0;
            DimsN = size(CurrentPlane);
            for iii = 1:DimsN(1)
                for iiii = 1:DimsN(1)
                    dist = sqrt((0-0)^2+(a/2 - 0)^2+(a/2 - a)^2)+ a/10;
                    if sqrt((CurrentPlane(iii,1)-CurrentPlane(iiii,1))^2+(CurrentPlane(iii,2)-CurrentPlane(iiii,2))^2+(CurrentPlane(iii,3)-CurrentPlane(iiii,3))^2) <= dist
                        if sqrt((CurrentPlane(iii,1)-CurrentPlane(iiii,1))^2+(CurrentPlane(iii,2)-CurrentPlane(iiii,2))^2+(CurrentPlane(iii,3)-CurrentPlane(iiii,3))^2) ~= 0
                            if CurrentPlane(iii,4) == 1 & CurrentPlane(iiii,4) == 1
                                numSam1 = numSam1 + 1;
                            elseif CurrentPlane(iii,4) == 2 & CurrentPlane(iiii,4) == 2
                                numSam2 = numSam2 + 1;
                            else 
                                numOpp = numOpp + 1;
                            end
                        end
                    end
                end
            end
                NumsPerPlane(w1,w2) = NumsPerPlane(w1,w2) + numSam1;
                NumsPerPlane(w1,w2+1) = NumsPerPlane(w2+1) +  numSam2;
                NumsPerPlane(w1,w2+2) = NumsPerPlane(w2+2) + numOpp;
                w2 = w2 + 3;
        end
        R = [];
        jj = 1;
            for j = 1:Dims(2)/4
                Naa = NumsPerPlane(jj);
                Nbb = NumsPerPlane(jj+1);
                Nab = NumsPerPlane(jj+2);
                R(j) = exp(-.0001*(Naa+Nbb*(1-.001)+Nab*(1-.001/2))/1);
                jj = jj+3;
            end
            
            Sorted = sort(R);
            Maxes = Sorted(end-1:end);
            Plane1 = 1;
            Plane2 = 1;
            num = 1;
            for jjj = 1:Dims(2)/4
                if  R(jjj) == Maxes(1)
                    Plane1 = num;
                elseif R(jjj) == Maxes(2)
                    Plane2 = num;
                end
                num = num + 1;
            end
            
            NewLattice = FCCLattice;
            ww = 2;
            Random = randi([1 3]);
            for jjjj = 1:length(StuffToAdd)-2
                for kkkk = 1:Dims(1)
                    if FCCLattice(kkkk,1)*Plane(1) + FCCLattice(kkkk,2)*Plane(2) + FCCLattice(kkkk,3)*Plane(3) + StuffToAdd(ww) == 0
                        if Random == 1
                            NewLattice(kkkk,1) = NewLattice(kkkk,1) + a/2;
                            NewLattice(kkkk,2) = NewLattice(kkkk,2) + 0;
                            NewLattice(kkkk,3) = NewLattice(kkkk,3) + a/2;
                        else
                            NewLattice(kkkk,1) = NewLattice(kkkk,1) - a/2;
                            NewLattice(kkkk,2) = NewLattice(kkkk,2) + 0;
                            NewLattice(kkkk,3) = NewLattice(kkkk,3) - a/2;
                        end

                    end
                end
                ww = ww +1;
            end
            clf
            for i = 1:length(NewLattice)
                if FCCLattice(i,4) == 1
                    plot3(NewLattice(i,1),NewLattice(i,2),NewLattice(i,3), 'o', 'MarkerSize',20, 'MarkerEdgeColor', 'red', 'MarkerFaceColor', 'red')
                    hold on
                elseif FCCLattice(i,4) == 2 
                    plot3(NewLattice(i,1),NewLattice(i,2),NewLattice(i,3), 'o', 'MarkerSize',20, 'MarkerEdgeColor', 'blue', 'MarkerFaceColor', 'blue')
                    hold on
                else
                    plot3(NewLattice(i,1),NewLattice(i,2),NewLattice(i,3), 'o', 'MarkerSize',20, 'MarkerEdgeColor', 'yellow', 'MarkerFaceColor', 'yellow')
                    hold on
                end 
            end
            axis equal
    elseif sysnum == 7
        x = [a 0 0];
        y = [0 -a 0];
        z = [0 0 a];
        Vec1 = y-x;
        Vec2 = z-x;
        C = cross(Vec1,Vec2);
        Plane = [C(1) C(2) C(3) C(1)*-x(1)+C(2)*-x(2)+C(3)*-x(3)];
        StuffToAdd = [];
        for m = -V*a*2:Plane(4):V*a*2 
            for ppp = 1:length(FCCLattice)
                if FCCLattice(ppp,1)*Plane(1)+ FCCLattice(ppp,2)*Plane(2) + FCCLattice(ppp,3)*Plane(3) + (Plane(4) + m) == 0
                    p = p+1;
                    SlipPlanes(p,w:w+3) = FCCLattice(ppp,:);
                    StuffToAdd = [StuffToAdd (Plane(4) + m)];
                    StuffToAdd = unique(StuffToAdd);
                end
            end
            w = w+4;
        end
        another = 1;
        PlaneFix = [];
        Dims = size(SlipPlanes);
        for l = 1:Dims(2)/4
            if max(SlipPlanes(:,4*l)) > 0
                PlaneFix(:,another*4-3:another*4) = SlipPlanes(:,4*l-3:4*l);
                another = another + 1;
            end
        end
        SlipPlanes = PlaneFix;
        Dims = size(SlipPlanes);
        NumsPerPlane = zeros(1,3*Dims(2)/4);
        w3 = 0;
        w1 = 1;
        w2 = 1;
        w4 = 4;
        for ii= 1:Dims(2)/4 
            w3 = w3 + 4;
            CurrentPlane = [];
            for o = 1:Dims(1)
                if SlipPlanes(o,w4)== 1 | SlipPlanes(o,w4)== 2
                    CurrentPlane = [CurrentPlane; SlipPlanes(o,w4-3:w4)];
                end
            end
            w4 = w4+4;
            numOpp = 0;
            numSam1 = 0;
            numSam2 = 0;
            DimsN = size(CurrentPlane);
            for iii = 1:DimsN(1)
                for iiii = 1:DimsN(1)
                    dist = sqrt((0-0)^2+(a/2 - 0)^2+(a/2 - a)^2)+ a/10;
                    if sqrt((CurrentPlane(iii,1)-CurrentPlane(iiii,1))^2+(CurrentPlane(iii,2)-CurrentPlane(iiii,2))^2+(CurrentPlane(iii,3)-CurrentPlane(iiii,3))^2) <= dist
                        if sqrt((CurrentPlane(iii,1)-CurrentPlane(iiii,1))^2+(CurrentPlane(iii,2)-CurrentPlane(iiii,2))^2+(CurrentPlane(iii,3)-CurrentPlane(iiii,3))^2) ~= 0
                            if CurrentPlane(iii,4) == 1 & CurrentPlane(iiii,4) == 1
                                numSam1 = numSam1 + 1;
                            elseif CurrentPlane(iii,4) == 2 & CurrentPlane(iiii,4) == 2
                                numSam2 = numSam2 + 1;
                            else 
                                numOpp = numOpp + 1;
                            end
                        end
                    end
                end
            end
                NumsPerPlane(w1,w2) = NumsPerPlane(w1,w2) + numSam1;
                NumsPerPlane(w1,w2+1) = NumsPerPlane(w2+1) +  numSam2;
                NumsPerPlane(w1,w2+2) = NumsPerPlane(w2+2) + numOpp;
                w2 = w2 + 3;
        end
        R = [];
        jj = 1;
            for j = 1:Dims(2)/4
                Naa = NumsPerPlane(jj);
                Nbb = NumsPerPlane(jj+1);
                Nab = NumsPerPlane(jj+2);
                R(j) = exp(-.0001*(Naa+Nbb*(1-.001)+Nab*(1-.001/2))/1);
                jj = jj+3;
            end
            
            Sorted = sort(R);
            Maxes = Sorted(end-1:end);
            Plane1 = 1;
            Plane2 = 1;
            num = 1;
            for jjj = 1:Dims(2)/4
                if  R(jjj) == Maxes(1)
                    Plane1 = num;
                elseif R(jjj) == Maxes(2)
                    Plane2 = num;
                end
                num = num + 1;
            end
            
            NewLattice = FCCLattice;
            ww = 2;
            Random = randi([1 3]);
            for jjjj = 1:length(StuffToAdd)-2
                for kkkk = 1:Dims(1)
                    if FCCLattice(kkkk,1)*Plane(1) + FCCLattice(kkkk,2)*Plane(2) + FCCLattice(kkkk,3)*Plane(3) + StuffToAdd(ww) == 0
                        if Random == 1
                            NewLattice(kkkk,1) = NewLattice(kkkk,1) + a/2;
                            NewLattice(kkkk,2) = NewLattice(kkkk,2) + a/2;
                            NewLattice(kkkk,3) = NewLattice(kkkk,3) + 0;
                        else
                            NewLattice(kkkk,1) = NewLattice(kkkk,1) - a/2;
                            NewLattice(kkkk,2) = NewLattice(kkkk,2) - a/2;
                            NewLattice(kkkk,3) = NewLattice(kkkk,3) + 0;
                        end
                        
                    end
                end
                ww = ww +1;
            end
            clf
            for i = 1:length(NewLattice)
                if FCCLattice(i,4) == 1
                    plot3(NewLattice(i,1),NewLattice(i,2),NewLattice(i,3), 'o', 'MarkerSize',20, 'MarkerEdgeColor', 'red', 'MarkerFaceColor', 'red')
                    hold on
                elseif FCCLattice(i,4) == 2 
                    plot3(NewLattice(i,1),NewLattice(i,2),NewLattice(i,3), 'o', 'MarkerSize',20, 'MarkerEdgeColor', 'blue', 'MarkerFaceColor', 'blue')
                    hold on
                else
                    plot3(NewLattice(i,1),NewLattice(i,2),NewLattice(i,3), 'o', 'MarkerSize',20, 'MarkerEdgeColor', 'yellow', 'MarkerFaceColor', 'yellow')
                    hold on
                end 
            end
            axis equal
    elseif sysnum == 8
        x = [a 0 0];
        y = [0 -a 0];
        z = [0 0 a];
        Vec1 = y-x;
        Vec2 = z-x;
        C = cross(Vec1,Vec2);
        Plane = [C(1) C(2) C(3) C(1)*-x(1)+C(2)*-x(2)+C(3)*-x(3)];
        StuffToAdd = [];
        for m = -V*a*2:Plane(4):V*a*2 
            for ppp = 1:length(FCCLattice)
                if FCCLattice(ppp,1)*Plane(1)+ FCCLattice(ppp,2)*Plane(2) + FCCLattice(ppp,3)*Plane(3) + (Plane(4) + m) == 0
                    p = p+1;
                    SlipPlanes(p,w:w+3) = FCCLattice(ppp,:);
                    StuffToAdd = [StuffToAdd (Plane(4) + m)];
                    StuffToAdd = unique(StuffToAdd);
                end
            end
            w = w+4;
        end
        another = 1;
        PlaneFix = [];
        Dims = size(SlipPlanes);
        for l = 1:Dims(2)/4
            if max(SlipPlanes(:,4*l)) > 0
                PlaneFix(:,another*4-3:another*4) = SlipPlanes(:,4*l-3:4*l);
                another = another + 1;
            end
        end
        SlipPlanes = PlaneFix;
        Dims = size(SlipPlanes);
        NumsPerPlane = zeros(1,3*Dims(2)/4);
        w3 = 0;
        w1 = 1;
        w2 = 1;
        w4 = 4;
        for ii= 1:Dims(2)/4 
            w3 = w3 + 4;
            CurrentPlane = [];
            for o = 1:Dims(1)
                if SlipPlanes(o,w4)== 1 | SlipPlanes(o,w4)== 2
                    CurrentPlane = [CurrentPlane; SlipPlanes(o,w4-3:w4)];
                end
            end
            w4 = w4+4;
            numOpp = 0;
            numSam1 = 0;
            numSam2 = 0;
            DimsN = size(CurrentPlane);
            for iii = 1:DimsN(1)
                for iiii = 1:DimsN(1)
                    dist = sqrt((0-0)^2+(a/2 - 0)^2+(a/2 - a)^2)+ a/10;
                    if sqrt((CurrentPlane(iii,1)-CurrentPlane(iiii,1))^2+(CurrentPlane(iii,2)-CurrentPlane(iiii,2))^2+(CurrentPlane(iii,3)-CurrentPlane(iiii,3))^2) <= dist
                        if sqrt((CurrentPlane(iii,1)-CurrentPlane(iiii,1))^2+(CurrentPlane(iii,2)-CurrentPlane(iiii,2))^2+(CurrentPlane(iii,3)-CurrentPlane(iiii,3))^2) ~= 0
                            if CurrentPlane(iii,4) == 1 & CurrentPlane(iiii,4) == 1
                                numSam1 = numSam1 + 1;
                            elseif CurrentPlane(iii,4) == 2 & CurrentPlane(iiii,4) == 2
                                numSam2 = numSam2 + 1;
                            else 
                                numOpp = numOpp + 1;
                            end
                        end
                    end
                end
            end
                NumsPerPlane(w1,w2) = NumsPerPlane(w1,w2) + numSam1;
                NumsPerPlane(w1,w2+1) = NumsPerPlane(w2+1) +  numSam2;
                NumsPerPlane(w1,w2+2) = NumsPerPlane(w2+2) + numOpp;
                w2 = w2 + 3;
        end
        R = [];
        jj = 1;
            for j = 1:Dims(2)/4
                Naa = NumsPerPlane(jj);
                Nbb = NumsPerPlane(jj+1);
                Nab = NumsPerPlane(jj+2);
                R(j) = exp(-.0001*(Naa+Nbb*(1-.001)+Nab*(1-.001/2))/1);
                jj = jj+3;
            end
            
            Sorted = sort(R);
            Maxes = Sorted(end-1:end);
            Plane1 = 1;
            Plane2 = 1;
            num = 1;
            for jjj = 1:Dims(2)/4
                if  R(jjj) == Maxes(1)
                    Plane1 = num;
                elseif R(jjj) == Maxes(2)
                    Plane2 = num;
                end
                num = num + 1;
            end
            
            NewLattice = FCCLattice;
            ww = 2;
            Random = randi([1 3]);
            for jjjj = 1:length(StuffToAdd)-2
                for kkkk = 1:Dims(1)
                    if FCCLattice(kkkk,1)*Plane(1) + FCCLattice(kkkk,2)*Plane(2) + FCCLattice(kkkk,3)*Plane(3) + StuffToAdd(ww) == 0
                        if Random == 1
                            NewLattice(kkkk,1) = NewLattice(kkkk,1) + 0;
                            NewLattice(kkkk,2) = NewLattice(kkkk,2) + a/2;
                            NewLattice(kkkk,3) = NewLattice(kkkk,3) + a/2;
                        else
                            NewLattice(kkkk,1) = NewLattice(kkkk,1) + 0;
                            NewLattice(kkkk,2) = NewLattice(kkkk,2) - a/2;
                            NewLattice(kkkk,3) = NewLattice(kkkk,3) - a/2;
                        end
                        
                    end
                end
                ww = ww +1;
            end
            clf
            for i = 1:length(NewLattice)
                if FCCLattice(i,4) == 1
                    plot3(NewLattice(i,1),NewLattice(i,2),NewLattice(i,3), 'o', 'MarkerSize',20, 'MarkerEdgeColor', 'red', 'MarkerFaceColor', 'red')
                    hold on
                elseif FCCLattice(i,4) == 2 
                    plot3(NewLattice(i,1),NewLattice(i,2),NewLattice(i,3), 'o', 'MarkerSize',20, 'MarkerEdgeColor', 'blue', 'MarkerFaceColor', 'blue')
                    hold on
                else
                    plot3(NewLattice(i,1),NewLattice(i,2),NewLattice(i,3), 'o', 'MarkerSize',20, 'MarkerEdgeColor', 'yellow', 'MarkerFaceColor', 'yellow')
                    hold on
                end 
            end
            axis equal
    elseif sysnum == 9 
        x = [a 0 0];
        y = [0 -a 0];
        z = [0 0 a];
        Vec1 = y-x;
        Vec2 = z-x;
        C = cross(Vec1,Vec2);
        Plane = [C(1) C(2) C(3) C(1)*-x(1)+C(2)*-x(2)+C(3)*-x(3)];
        StuffToAdd = [];
        for m = -V*a*2:Plane(4):V*a*2 
            for ppp = 1:length(FCCLattice)
                if FCCLattice(ppp,1)*Plane(1)+ FCCLattice(ppp,2)*Plane(2) + FCCLattice(ppp,3)*Plane(3) + (Plane(4) + m) == 0
                    p = p+1;
                    SlipPlanes(p,w:w+3) = FCCLattice(ppp,:);
                    StuffToAdd = [StuffToAdd (Plane(4) + m)];
                    StuffToAdd = unique(StuffToAdd);
                end
            end
            w = w+4;
        end
        another = 1;
        PlaneFix = [];
        Dims = size(SlipPlanes);
        for l = 1:Dims(2)/4
            if max(SlipPlanes(:,4*l)) > 0
                PlaneFix(:,another*4-3:another*4) = SlipPlanes(:,4*l-3:4*l);
                another = another + 1;
            end
        end
        SlipPlanes = PlaneFix;
        Dims = size(SlipPlanes);
        NumsPerPlane = zeros(1,3*Dims(2)/4);
        w3 = 0;
        w1 = 1;
        w2 = 1;
        w4 = 4;
        for ii= 1:Dims(2)/4 
            w3 = w3 + 4;
            CurrentPlane = [];
            for o = 1:Dims(1)
                if SlipPlanes(o,w4)== 1 | SlipPlanes(o,w4)== 2
                    CurrentPlane = [CurrentPlane; SlipPlanes(o,w4-3:w4)];
                end
            end
            w4 = w4+4;
            numOpp = 0;
            numSam1 = 0;
            numSam2 = 0;
            DimsN = size(CurrentPlane);
            for iii = 1:DimsN(1)
                for iiii = 1:DimsN(1)
                    dist = sqrt((0-0)^2+(a/2 - 0)^2+(a/2 - a)^2)+ a/10;
                    if sqrt((CurrentPlane(iii,1)-CurrentPlane(iiii,1))^2+(CurrentPlane(iii,2)-CurrentPlane(iiii,2))^2+(CurrentPlane(iii,3)-CurrentPlane(iiii,3))^2) <= dist
                        if sqrt((CurrentPlane(iii,1)-CurrentPlane(iiii,1))^2+(CurrentPlane(iii,2)-CurrentPlane(iiii,2))^2+(CurrentPlane(iii,3)-CurrentPlane(iiii,3))^2) ~= 0
                            if CurrentPlane(iii,4) == 1 & CurrentPlane(iiii,4) == 1
                                numSam1 = numSam1 + 1;
                            elseif CurrentPlane(iii,4) == 2 & CurrentPlane(iiii,4) == 2
                                numSam2 = numSam2 + 1;
                            else 
                                numOpp = numOpp + 1;
                            end
                        end
                    end
                end
            end
                NumsPerPlane(w1,w2) = NumsPerPlane(w1,w2) + numSam1;
                NumsPerPlane(w1,w2+1) = NumsPerPlane(w2+1) +  numSam2;
                NumsPerPlane(w1,w2+2) = NumsPerPlane(w2+2) + numOpp;
                w2 = w2 + 3;
        end
        R = [];
        jj = 1;
            for j = 1:Dims(2)/4
                Naa = NumsPerPlane(jj);
                Nbb = NumsPerPlane(jj+1);
                Nab = NumsPerPlane(jj+2);
                R(j) = exp(-.0001*(Naa+Nbb*(1-.001)+Nab*(1-.001/2))/1);
                jj = jj+3;
            end
            
            Sorted = sort(R);
            Maxes = Sorted(end-1:end);
            Plane1 = 1;
            Plane2 = 1;
            num = 1;
            for jjj = 1:Dims(2)/4
                if  R(jjj) == Maxes(1)
                    Plane1 = num;
                elseif R(jjj) == Maxes(2)
                    Plane2 = num;
                end
                num = num + 1;
            end
            
            NewLattice = FCCLattice;
            ww = 2;
            Random = randi([1 3]);
            for jjjj = 1:length(StuffToAdd)-2
                for kkkk = 1:Dims(1)
                    if FCCLattice(kkkk,1)*Plane(1) + FCCLattice(kkkk,2)*Plane(2) + FCCLattice(kkkk,3)*Plane(3) + StuffToAdd(ww) == 0
                        if Random == 1
                            NewLattice(kkkk,1) = NewLattice(kkkk,1) + a/2;
                            NewLattice(kkkk,2) = NewLattice(kkkk,2) + 0;
                            NewLattice(kkkk,3) = NewLattice(kkkk,3) - a/2;
                        else
                            NewLattice(kkkk,1) = NewLattice(kkkk,1) - a/2;
                            NewLattice(kkkk,2) = NewLattice(kkkk,2) + 0;
                            NewLattice(kkkk,3) = NewLattice(kkkk,3) + a/2;
                        end
                       
                    end
                end
                ww = ww +1;
            end
            clf
            for i = 1:length(NewLattice)
                if FCCLattice(i,4) == 1
                    plot3(NewLattice(i,1),NewLattice(i,2),NewLattice(i,3), 'o', 'MarkerSize',20, 'MarkerEdgeColor', 'red', 'MarkerFaceColor', 'red')
                    hold on
                elseif FCCLattice(i,4) == 2 
                    plot3(NewLattice(i,1),NewLattice(i,2),NewLattice(i,3), 'o', 'MarkerSize',20, 'MarkerEdgeColor', 'blue', 'MarkerFaceColor', 'blue')
                    hold on
                else
                    plot3(NewLattice(i,1),NewLattice(i,2),NewLattice(i,3), 'o', 'MarkerSize',20, 'MarkerEdgeColor', 'yellow', 'MarkerFaceColor', 'yellow')
                    hold on
                end 
            end
            axis equal
    elseif sysnum == 10
        x = [a 0 0];
        y = [0 a 0];
        z = [0 0 -a];
        Vec1 = y-x;
        Vec2 = z-x;
        C = cross(Vec1,Vec2);
        Plane = [C(1) C(2) C(3) C(1)*-x(1)+C(2)*-x(2)+C(3)*-x(3)];
        StuffToAdd = [];
        for m = -V*a*2:Plane(4):V*a*2 
            for ppp = 1:length(FCCLattice)
                if FCCLattice(ppp,1)*Plane(1)+ FCCLattice(ppp,2)*Plane(2) + FCCLattice(ppp,3)*Plane(3) + (Plane(4) + m) == 0
                    p = p+1;
                    SlipPlanes(p,w:w+3) = FCCLattice(ppp,:);
                    StuffToAdd = [StuffToAdd (Plane(4) + m)];
                    StuffToAdd = unique(StuffToAdd);
                end
            end
            w = w+4;
        end
        another = 1;
        PlaneFix = [];
        Dims = size(SlipPlanes);
        for l = 1:Dims(2)/4
            if max(SlipPlanes(:,4*l)) > 0
                PlaneFix(:,another*4-3:another*4) = SlipPlanes(:,4*l-3:4*l);
                another = another + 1;
            end
        end
        SlipPlanes = PlaneFix;
        Dims = size(SlipPlanes);
        NumsPerPlane = zeros(1,3*Dims(2)/4);
        w3 = 0;
        w1 = 1;
        w2 = 1;
        w4 = 4;
        for ii= 1:Dims(2)/4 
            w3 = w3 + 4;
            CurrentPlane = [];
            for o = 1:Dims(1)
                if SlipPlanes(o,w4)== 1 | SlipPlanes(o,w4)== 2
                    CurrentPlane = [CurrentPlane; SlipPlanes(o,w4-3:w4)];
                end
            end
            w4 = w4+4;
            numOpp = 0;
            numSam1 = 0;
            numSam2 = 0;
            DimsN = size(CurrentPlane);
            for iii = 1:DimsN(1)
                for iiii = 1:DimsN(1)
                    dist = sqrt((0-0)^2+(a/2 - 0)^2+(a/2 - a)^2)+ a/10;
                    if sqrt((CurrentPlane(iii,1)-CurrentPlane(iiii,1))^2+(CurrentPlane(iii,2)-CurrentPlane(iiii,2))^2+(CurrentPlane(iii,3)-CurrentPlane(iiii,3))^2) <= dist
                        if sqrt((CurrentPlane(iii,1)-CurrentPlane(iiii,1))^2+(CurrentPlane(iii,2)-CurrentPlane(iiii,2))^2+(CurrentPlane(iii,3)-CurrentPlane(iiii,3))^2) ~= 0
                            if CurrentPlane(iii,4) == 1 & CurrentPlane(iiii,4) == 1
                                numSam1 = numSam1 + 1;
                            elseif CurrentPlane(iii,4) == 2 & CurrentPlane(iiii,4) == 2
                                numSam2 = numSam2 + 1;
                            else 
                                numOpp = numOpp + 1;
                            end
                        end
                    end
                end
            end
                NumsPerPlane(w1,w2) = NumsPerPlane(w1,w2) + numSam1;
                NumsPerPlane(w1,w2+1) = NumsPerPlane(w2+1) +  numSam2;
                NumsPerPlane(w1,w2+2) = NumsPerPlane(w2+2) + numOpp;
                w2 = w2 + 3;
        end
        R = [];
        jj = 1;
            for j = 1:Dims(2)/4
                Naa = NumsPerPlane(jj);
                Nbb = NumsPerPlane(jj+1);
                Nab = NumsPerPlane(jj+2);
                R(j) = exp(-.0001*(Naa+Nbb*(1-.001)+Nab*(1-.001/2))/1);
                jj = jj+3;
            end
            
            Sorted = sort(R);
            Maxes = Sorted(end-1:end);
            Plane1 = 1;
            Plane2 = 1;
            num = 1;
            for jjj = 1:Dims(2)/4
                if  R(jjj) == Maxes(1)
                    Plane1 = num;
                elseif R(jjj) == Maxes(2)
                    Plane2 = num;
                end
                num = num + 1;
            end
            
            NewLattice = FCCLattice;
            ww = 2;
            Random = randi([1 3]);
            for jjjj = 1:length(StuffToAdd)-2
                for kkkk = 1:Dims(1)
                    if FCCLattice(kkkk,1)*Plane(1) + FCCLattice(kkkk,2)*Plane(2) + FCCLattice(kkkk,3)*Plane(3) + StuffToAdd(ww) == 0
                        if Random == 1
                            NewLattice(kkkk,1) = NewLattice(kkkk,1) + a/2;
                            NewLattice(kkkk,2) = NewLattice(kkkk,2) - a/2;
                            NewLattice(kkkk,3) = NewLattice(kkkk,3) + 0;
                        else
                            NewLattice(kkkk,1) = NewLattice(kkkk,1) - a/2;
                            NewLattice(kkkk,2) = NewLattice(kkkk,2) + a/2;
                            NewLattice(kkkk,3) = NewLattice(kkkk,3) + 0;
                        end
                        
                    end
                end
                ww = ww +1;
            end
            clf
            for i = 1:length(NewLattice)
                if FCCLattice(i,4) == 1
                    plot3(NewLattice(i,1),NewLattice(i,2),NewLattice(i,3), 'o', 'MarkerSize',20, 'MarkerEdgeColor', 'red', 'MarkerFaceColor', 'red')
                    hold on
                elseif FCCLattice(i,4) == 2 
                    plot3(NewLattice(i,1),NewLattice(i,2),NewLattice(i,3), 'o', 'MarkerSize',20, 'MarkerEdgeColor', 'blue', 'MarkerFaceColor', 'blue')
                    hold on
                else
                    plot3(NewLattice(i,1),NewLattice(i,2),NewLattice(i,3), 'o', 'MarkerSize',20, 'MarkerEdgeColor', 'yellow', 'MarkerFaceColor', 'yellow')
                    hold on
                end 
            end
            axis equal
    elseif sysnum == 11
        x = [a 0 0];
        y = [0 a 0];
        z = [0 0 -a];
        Vec1 = y-x;
        Vec2 = z-x;
        C = cross(Vec1,Vec2);
        Plane = [C(1) C(2) C(3) C(1)*-x(1)+C(2)*-x(2)+C(3)*-x(3)];
        StuffToAdd = [];
        for m = -V*a*2:Plane(4):V*a*2 
            for ppp = 1:length(FCCLattice)
                if FCCLattice(ppp,1)*Plane(1)+ FCCLattice(ppp,2)*Plane(2) + FCCLattice(ppp,3)*Plane(3) + (Plane(4) + m) == 0
                    p = p+1;
                    SlipPlanes(p,w:w+3) = FCCLattice(ppp,:);
                    StuffToAdd = [StuffToAdd (Plane(4) + m)];
                    StuffToAdd = unique(StuffToAdd);
                end
            end
            w = w+4;
        end
        another = 1;
        PlaneFix = [];
        Dims = size(SlipPlanes);
        for l = 1:Dims(2)/4
            if max(SlipPlanes(:,4*l)) > 0
                PlaneFix(:,another*4-3:another*4) = SlipPlanes(:,4*l-3:4*l);
                another = another + 1;
            end
        end
        SlipPlanes = PlaneFix;
        Dims = size(SlipPlanes);
        NumsPerPlane = zeros(1,3*Dims(2)/4);
        w3 = 0;
        w1 = 1;
        w2 = 1;
        w4 = 4;
        for ii= 1:Dims(2)/4 
            w3 = w3 + 4;
            CurrentPlane = [];
            for o = 1:Dims(1)
                if SlipPlanes(o,w4)== 1 | SlipPlanes(o,w4)== 2
                    CurrentPlane = [CurrentPlane; SlipPlanes(o,w4-3:w4)];
                end
            end
            w4 = w4+4;
            numOpp = 0;
            numSam1 = 0;
            numSam2 = 0;
            DimsN = size(CurrentPlane);
            for iii = 1:DimsN(1)
                for iiii = 1:DimsN(1)
                    dist = sqrt((0-0)^2+(a/2 - 0)^2+(a/2 - a)^2)+ a/10;
                    if sqrt((CurrentPlane(iii,1)-CurrentPlane(iiii,1))^2+(CurrentPlane(iii,2)-CurrentPlane(iiii,2))^2+(CurrentPlane(iii,3)-CurrentPlane(iiii,3))^2) <= dist
                        if sqrt((CurrentPlane(iii,1)-CurrentPlane(iiii,1))^2+(CurrentPlane(iii,2)-CurrentPlane(iiii,2))^2+(CurrentPlane(iii,3)-CurrentPlane(iiii,3))^2) ~= 0
                            if CurrentPlane(iii,4) == 1 & CurrentPlane(iiii,4) == 1
                                numSam1 = numSam1 + 1;
                            elseif CurrentPlane(iii,4) == 2 & CurrentPlane(iiii,4) == 2
                                numSam2 = numSam2 + 1;
                            else 
                                numOpp = numOpp + 1;
                            end
                        end
                    end
                end
            end
                NumsPerPlane(w1,w2) = NumsPerPlane(w1,w2) + numSam1;
                NumsPerPlane(w1,w2+1) = NumsPerPlane(w2+1) +  numSam2;
                NumsPerPlane(w1,w2+2) = NumsPerPlane(w2+2) + numOpp;
                w2 = w2 + 3;
        end
        R = [];
        jj = 1;
            for j = 1:Dims(2)/4
                Naa = NumsPerPlane(jj);
                Nbb = NumsPerPlane(jj+1);
                Nab = NumsPerPlane(jj+2);
                R(j) = exp(-.0001*(Naa+Nbb*(1-.001)+Nab*(1-.001/2))/1);
                jj = jj+3;
            end
            
            Sorted = sort(R);
            Maxes = Sorted(end-1:end);
            Plane1 = 1;
            Plane2 = 1;
            num = 1;
            for jjj = 1:Dims(2)/4
                if  R(jjj) == Maxes(1)
                    Plane1 = num;
                elseif R(jjj) == Maxes(2)
                    Plane2 = num;
                end
                num = num + 1;
            end
            
            NewLattice = FCCLattice;
            ww = 2;
            Random = randi([1 3]);
            for jjjj = 1:length(StuffToAdd)-2
                for kkkk = 1:Dims(1)
                    if FCCLattice(kkkk,1)*Plane(1) + FCCLattice(kkkk,2)*Plane(2) + FCCLattice(kkkk,3)*Plane(3) + StuffToAdd(ww) == 0
                        if Random == 1
                            NewLattice(kkkk,1) = NewLattice(kkkk,1) + 0;
                            NewLattice(kkkk,2) = NewLattice(kkkk,2) + a/2;
                            NewLattice(kkkk,3) = NewLattice(kkkk,3) + a/2;
                        else
                            NewLattice(kkkk,1) = NewLattice(kkkk,1) + 0;
                            NewLattice(kkkk,2) = NewLattice(kkkk,2) - a/2;
                            NewLattice(kkkk,3) = NewLattice(kkkk,3) - a/2;
                        end
                       
                    end
                end
                ww = ww +1;
            end
            clf
            for i = 1:length(NewLattice)
                if FCCLattice(i,4) == 1
                    plot3(NewLattice(i,1),NewLattice(i,2),NewLattice(i,3), 'o', 'MarkerSize',20, 'MarkerEdgeColor', 'red', 'MarkerFaceColor', 'red')
                    hold on
                elseif FCCLattice(i,4) == 2 
                    plot3(NewLattice(i,1),NewLattice(i,2),NewLattice(i,3), 'o', 'MarkerSize',20, 'MarkerEdgeColor', 'blue', 'MarkerFaceColor', 'blue')
                    hold on
                else
                    plot3(NewLattice(i,1),NewLattice(i,2),NewLattice(i,3), 'o', 'MarkerSize',20, 'MarkerEdgeColor', 'yellow', 'MarkerFaceColor', 'yellow')
                    hold on
                end 
            end
            axis equal
    else
        x = [a 0 0];
        y = [0 a 0];
        z = [0 0 -a];
        Vec1 = y-x;
        Vec2 = z-x;
        C = cross(Vec1,Vec2);
        Plane = [C(1) C(2) C(3) C(1)*-x(1)+C(2)*-x(2)+C(3)*-x(3)];
        StuffToAdd = [];
        for m = -V*a*2:Plane(4):V*a*2 
            for ppp = 1:length(FCCLattice)
                if FCCLattice(ppp,1)*Plane(1)+ FCCLattice(ppp,2)*Plane(2) + FCCLattice(ppp,3)*Plane(3) + (Plane(4) + m) == 0
                    p = p+1;
                    SlipPlanes(p,w:w+3) = FCCLattice(ppp,:);
                    StuffToAdd = [StuffToAdd (Plane(4) + m)];
                    StuffToAdd = unique(StuffToAdd);
                end
            end
            w = w+4;
        end
        another = 1;
        PlaneFix = [];
        Dims = size(SlipPlanes);
        for l = 1:Dims(2)/4
            if max(SlipPlanes(:,4*l)) > 0
                PlaneFix(:,another*4-3:another*4) = SlipPlanes(:,4*l-3:4*l);
                another = another + 1;
            end
        end
        SlipPlanes = PlaneFix;
        Dims = size(SlipPlanes);
        NumsPerPlane = zeros(1,3*Dims(2)/4);
        w3 = 0;
        w1 = 1;
        w2 = 1;
        w4 = 4;
        for ii= 1:Dims(2)/4 
            w3 = w3 + 4;
            CurrentPlane = [];
            for o = 1:Dims(1)
                if SlipPlanes(o,w4)== 1 | SlipPlanes(o,w4)== 2
                    CurrentPlane = [CurrentPlane; SlipPlanes(o,w4-3:w4)];
                end
            end
            w4 = w4+4;
            numOpp = 0;
            numSam1 = 0;
            numSam2 = 0;
            DimsN = size(CurrentPlane);
            for iii = 1:DimsN(1)
                for iiii = 1:DimsN(1)
                    dist = sqrt((0-0)^2+(a/2 - 0)^2+(a/2 - a)^2)+ a/10;
                    if sqrt((CurrentPlane(iii,1)-CurrentPlane(iiii,1))^2+(CurrentPlane(iii,2)-CurrentPlane(iiii,2))^2+(CurrentPlane(iii,3)-CurrentPlane(iiii,3))^2) <= dist
                        if sqrt((CurrentPlane(iii,1)-CurrentPlane(iiii,1))^2+(CurrentPlane(iii,2)-CurrentPlane(iiii,2))^2+(CurrentPlane(iii,3)-CurrentPlane(iiii,3))^2) ~= 0
                            if CurrentPlane(iii,4) == 1 & CurrentPlane(iiii,4) == 1
                                numSam1 = numSam1 + 1;
                            elseif CurrentPlane(iii,4) == 2 & CurrentPlane(iiii,4) == 2
                                numSam2 = numSam2 + 1;
                            else 
                                numOpp = numOpp + 1;
                            end
                        end
                    end
                end
            end
                NumsPerPlane(w1,w2) = NumsPerPlane(w1,w2) + numSam1;
                NumsPerPlane(w1,w2+1) = NumsPerPlane(w2+1) +  numSam2;
                NumsPerPlane(w1,w2+2) = NumsPerPlane(w2+2) + numOpp;
                w2 = w2 + 3;
        end
        R = [];
        jj = 1;
            for j = 1:Dims(2)/4
                Naa = NumsPerPlane(jj);
                Nbb = NumsPerPlane(jj+1);
                Nab = NumsPerPlane(jj+2);
                R(j) = exp(-.0001*(Naa+Nbb*(1-.001)+Nab*(1-.001/2))/1);
                jj = jj+3;
            end
            
            Sorted = sort(R);
            Maxes = Sorted(end-1:end);
            Plane1 = 1;
            Plane2 = 1;
            num = 1;
            for jjj = 1:Dims(2)/4
                if  R(jjj) == Maxes(1)
                    Plane1 = num;
                elseif R(jjj) == Maxes(2)
                    Plane2 = num;
                end
                num = num + 1;
            end
            
            NewLattice = FCCLattice;
            ww = 2;
            Random = randi([1 3]);
            for jjjj = 1:length(StuffToAdd)-2
                for kkkk = 1:Dims(1)
                    if FCCLattice(kkkk,1)*Plane(1) + FCCLattice(kkkk,2)*Plane(2) + FCCLattice(kkkk,3)*Plane(3) + StuffToAdd(ww) == 0
                        if Random == 1
                            NewLattice(kkkk,1) = NewLattice(kkkk,1) + 0;
                            NewLattice(kkkk,2) = NewLattice(kkkk,2) + a/2;
                            NewLattice(kkkk,3) = NewLattice(kkkk,3) + a/2;
                        else
                            NewLattice(kkkk,1) = NewLattice(kkkk,1) + 0;
                            NewLattice(kkkk,2) = NewLattice(kkkk,2) - a/2;
                            NewLattice(kkkk,3) = NewLattice(kkkk,3) - a/2;
                        end
                        
                    end
                end
                ww = ww +1;
            end
            clf
            for i = 1:length(NewLattice)
                if FCCLattice(i,4) == 1
                    plot3(NewLattice(i,1),NewLattice(i,2),NewLattice(i,3), 'o', 'MarkerSize',20, 'MarkerEdgeColor', 'red', 'MarkerFaceColor', 'red')
                    hold on
                elseif FCCLattice(i,4) == 2 
                    plot3(NewLattice(i,1),NewLattice(i,2),NewLattice(i,3), 'o', 'MarkerSize',20, 'MarkerEdgeColor', 'blue', 'MarkerFaceColor', 'blue')
                    hold on
                else
                    plot3(NewLattice(i,1),NewLattice(i,2),NewLattice(i,3), 'o', 'MarkerSize',20, 'MarkerEdgeColor', 'yellow', 'MarkerFaceColor', 'yellow')
                    hold on
                end 
            end
            axis equal
    end
end