% dimensions
n = 8;

%Initial field configuration
Phi_old_r = 2*pi*rand(8,8,8);
Phi_old_I = 2*pi*rand(8,8,8);

% Generating an array for Noise
eta = zeros(n,n,n);

% beta value
b = 0.7;

% Total steps
totalsteps = 10000;

% Step size
targetedstepsize = 0.005;

% Creating a file to input data
filename = append('b = ',sprintf('%3.1f',b) ,'.txt');
fileID = fopen(filename,'w');
fprintf(fileID,'%s\t%s\t%s\t%s\t%s\t%s\t%s\n','langevin time','Phi_r(0,0,0)','Phi_i(0,0,0)','Action_r','Action_i','Kmax','Noise(0,0,0)');

for uu = 0:0.05:0.5
    u = uu^(0.5);
    time = 0;
    count = 0;
    imagphi = 0;
    imagphisq = 0;
    Sr = 0;
    Si = 0;
    KMAX = 0;
    for sweeps = 1:totalsteps
        eta = normrnd(0,sqrt(2),8,8,8);
        
        gradSR = complexgradientreal(Phi_old_r, Phi_old_I, u, b,n);
        gradSI = complexgradientimag(Phi_old_r, Phi_old_I, u, b,n);
        
        maxgrad = power((power(gradSR,2) + power(gradSI,2)),0.5);
        highestgrad = max(maxgrad,[],'all');
        
        epsilon = targetedstepsize;
        time = time + epsilon;
        
        Phi_NEW_r = complexlangevinreal(Phi_old_r, gradSR, epsilon,eta,n);
        Phi_old_r = mod(Phi_NEW_r, 2*pi);
        
        Phi_NEW_I = complexlangevinimag( Phi_old_I, gradSI, epsilon,n);
        Phi_old_I = mod(Phi_NEW_I, 2*pi); 
        
        if sweeps >= 4000 && mod(sweeps,10) == 0
            Realaction = (realpartofaction(Phi_old_r, Phi_old_I, u, b,n) / power(n,3));
            Imagaction = (imagpartofaction(Phi_old_r, Phi_old_I, u, b,n) / power(n,3));
            Sr = Sr + Realaction;
            Si = Si + Imagaction;
            
            imagphi = imagphi + sum(Phi_old_I,"all");
            imagphisq = imagphisq + sum(power(Phi_old_I,2),"all");
            
            KMAX = KMAX + highestgrad;
            
            count = count + 1;
            
            fprintf(fileID,'%d\t%d\t%d\t%d\t%d\t%d\t%d\n',time,Phi_old_r(1,1,1),Phi_old_I(1,1,1),Realaction,Imagaction,highestgrad,eta(1,1,1));  
        end
    end
    fprintf(fileID,'\n');
        
    avgactionstring = append('Average value of Action density = ',sprintf('%d',Sr/count),'\n');
    avgkmaxstring = append('Average value of Kmax = ',sprintf('%d',KMAX/count),'\n');
    avgimagphistring = append('sum of imag phi values = ',sprintf('%d',imagphi),'\n');
    avgimagphisqstring = append('sum of square of imag phi values = ',sprintf('%d',imagphisq),'\n');
    disp(avgactionstring)
    disp(avgkmaxstring)
    disp(avgimagphistring)
    disp(avgimagphisqstring)
end


%Functions are defined below
function gradSr = complexgradientreal(Phi_old_r,Phi_old_I,u,b,n)
    gradSr = zeros(n,n,n);

    for x = 1:n
        for y = 1:n
            for z = 1:n
                aax = (sin(Phi_old_r(x,y,z) - Phi_old_r(modul((x + 1) ,n),y,z))) * (cosh(Phi_old_I(x,y,z) - Phi_old_I(modul((x + 1) ,n),y,z)));
                bbx = (sin(Phi_old_r(x,y,z) - Phi_old_r(modul((x - 1) ,n),y,z))) * (cosh(Phi_old_I(x,y,z) - Phi_old_I(modul((x - 1) ,n),y,z)));

                aay = (sin(Phi_old_r(x,y,z) - Phi_old_r(x,modul((y + 1) ,n),z))) * (cosh(Phi_old_I(x,y,z) - Phi_old_I(x,modul((y + 1) ,n),z)));
                bby = (sin(Phi_old_r(x,y,z) - Phi_old_r(x,modul((y - 1) ,n),z))) * (cosh(Phi_old_I(x,y,z) - Phi_old_I(x,modul((y - 1) ,n),z)));

                aaz = (sin(Phi_old_r(x,y,z) - Phi_old_r(x,y,modul((z + 1) ,n)))) * (cosh(Phi_old_I(x,y,z) - (u) - Phi_old_I(x,y,modul((z + 1) ,n))));
                bbz = (sin(Phi_old_r(x,y,z) - Phi_old_r(x,y,modul((z - 1) ,n)))) * (cosh(Phi_old_I(x,y,z) + (u) - Phi_old_I(x,y,modul((z - 1) ,n))));
                
                gradSr(x,y,z) = (-b) * (aax + bbx + aay + bby + aaz + bbz);
            end
        end
    end
end

function gradSi = complexgradientimag(Phi_old_r,Phi_old_I,u,b,n)
    gradSi = zeros(n,n,n);

    for x = 1:n
        for y = 1:n
            for z = 1:n
                ccx = (cos(Phi_old_r(x,y,z) - Phi_old_r(modul((x + 1) ,n),y,z))) * (sinh(Phi_old_I(x,y,z) - Phi_old_I(modul((x + 1) ,n),y,z)));
                ddx = (cos(Phi_old_r(x,y,z) - Phi_old_r(modul((x - 1) ,n),y,z))) * (sinh(Phi_old_I(x,y,z) - Phi_old_I(modul((x - 1) ,n),y,z)));

                ccy = (cos(Phi_old_r(x,y,z) - Phi_old_r(x,modul((y + 1) ,n),z))) * (sinh(Phi_old_I(x,y,z) - Phi_old_I(x,modul((y + 1) ,n),z)));
                ddy = (cos(Phi_old_r(x,y,z) - Phi_old_r(x,modul((y - 1) ,n),z))) * (sinh(Phi_old_I(x,y,z) - Phi_old_I(x,modul((y - 1) ,n),z)));

                ccz = (cos(Phi_old_r(x,y,z) - Phi_old_r(x,y,modul((z + 1), n)))) * (sinh(Phi_old_I(x,y,z) - Phi_old_I(x,y,modul((z + 1), n)) - u));
                ddz = (cos(Phi_old_r(x,y,z) - Phi_old_r(x,y,modul((z - 1), n)))) * (sinh(Phi_old_I(x,y,z) - Phi_old_I(x,y,modul((z - 1), n)) + u));
                
                gradSi(x,y,z) = (-b) * (ccx + ddx + ccy + ddy + ccz + ddz);
            end
        end
    end
end

function S_r = realpartofaction(Phi_old_r, Phi_old_I, u, b,n)
    S_r = 0;
    for x = 1:n
        for y = 1:n
            for z = 1:n
                S_r = S_r +  (-b) * ((cos(Phi_old_r(x,y,z) - Phi_old_r(modul((x + 1) ,n),y,z))) * (cosh(Phi_old_I(x,y,z) - Phi_old_I(modul((x + 1) ,n),y,z))));                                  
                S_r = S_r +  (-b) * ((cos(Phi_old_r(x,y,z) - Phi_old_r(x,modul((y + 1) ,n),z))) * (cosh(Phi_old_I(x,y,z) - Phi_old_I(x,modul((y+1),n),z))));
                S_r = S_r +  (-b) * ((cos(Phi_old_r(x,y,z) - Phi_old_r(x,y,modul((z + 1) ,n)))) * (cosh(Phi_old_I(x,y,z) - Phi_old_I(x,y,modul((z + 1) ,n)) - u)));
            end
        end
    end
end


function S_I = imagpartofaction(Phi_old_r, Phi_old_I, u, b,n)
    S_I = 0;
    for x = 1:n
        for y = 1:n
            for z = 1:n
                S_I = S_I +  (b) * ((sin(Phi_old_r(x,y,z) - Phi_old_r(modul((x + 1) ,n),y,z))) * (sinh(Phi_old_I(x,y,z) - Phi_old_I(modul((x + 1) ,n),y,z))));                                  
                S_I = S_I +  (b) * ((sin(Phi_old_r(x,y,z) - Phi_old_r(x,modul((y + 1) ,n),z))) * (sinh(Phi_old_I(x,y,z) - Phi_old_I(x,modul((y+1),n),z))));
                S_I = S_I +  (b) * ((sin(Phi_old_r(x,y,z) - Phi_old_r(x,y,modul((z + 1) ,n)))) * (sinh(Phi_old_I(x,y,z) - Phi_old_I(x,y,modul((z + 1) ,n)) - u)));
            end
        end
    end
end

function Phi_new_r = complexlangevinreal(Phi_old_r,gradSr,epsilon,eta,n)
    Phi_new_r = zeros(n,n,n);
    Phi_new_r = Phi_old_r + (epsilon * gradSr) + power(epsilon,0.5) * eta;    
end

function Phi_new_I = complexlangevinimag(Phi_old_I,gradSi,epsilon,n)
    Phi_new_I = zeros(n,n,n);
    Phi_new_I = Phi_old_I + (epsilon * gradSi);    
end

function epsilon = adaptivestepsize(avgKmax, Kmax , targetedstepsize)
    epsilon = 0;
    if Kmax > avgKmax
        epsilon = targetedstepsize * (avgKmax/Kmax);
    else
        epsilon = targetedstepsize;
    end
end

function modss = modul(x,y)
    if mod(x,y) == 0
        modss = y;
    else
        modss = mod(x,y);
    end
end
