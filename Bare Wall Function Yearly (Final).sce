//Base Case: Bare Wall
//Initialize script
close;
clear;
clc;

r = 1
q = 1
J = 1
while J < 366 || J == 366 then
Date = [r q 2022]

//General input data 
/*
Layers
A -> External Environment
B -> Green Layer
C -> Aeroponic Wall 1
D -> Aeroponic Chamber
E -> Aeroponic Wall 2
F -> Air Gap
G -> Wall
H -> Internal Environment 

Temperatures 
1 -> External Environment
2  -> External wall
3 -> Internal wall
4 -> Internal Environment
*/

//Cross-sectional area
Gl = 1 //Wall length [m]
Gh = 1 //Wall height [m]
Gn = 0.22 //Wall thickness [m] 
AG = Gl*Gh //Wall face area [m^2]

//Thermal Conductivity 
kG = 1 //Wall [W/mK]

//Additional inpute
T4 = 21 + 273.15 //Internal room temperature [degC]
Gsc = 0.0820        //(MJ/(m^2.min))
deg = -33.96        //Degrees of latitude of location (-ve for Southern hemisphere)


//Calculating J
hr = Date(1)
Month = Date(2)
Year = Date(3)

Months = [31 28 31 30 31 30 31 31 30 31 30 31]
if int(Year/4) == Year/4  then
    Months(2) = 29
end

J = sum(Months(1:Month-1))+hr

//Access Data from Excel
[fd,SST,Sheetnames,Sheetpos] = xls_open('Model Inputs.xls')
[Value,TextInd] = xls_read(fd,Sheetpos(5))

X = J*24-22 //Date placeholder in excel
Y = X+23

Time = Value(X:Y,4)
T1s = Value(X:Y,5)
RHs = Value(X:Y,6)
Rss = Value(X:Y,7)
Rsos = Value(X:Y,8)
Wss = Value(X:Y,9)

//Access physical properties of air from excel
[Value,TextInd] = xls_read(fd,Sheetpos(3))
mclose();
Temps = Value(22:36,1)
ps = Value(22:36,2)
Cps = Value(22:36,3)
vs = Value(22:36,4)
ks = Value(22:36,6)
Prs = Value(22:36,8)
gBvs = Value(22:36,9)

i = 1
for i = 1:length(T1s)
    
    T1 = T1s(i) + 273.15 //Outside temperature [degC]
    
//Radiation
//Extracting Azimuth Angle
[fd,SST,Sheetnames,Sheetpos] = xls_open('Orientation.xls')
[Value,TextInd] = xls_read(fd,Sheetpos(2))
mclose()
Az = Value(J+1,2:25)'
Az = Az(:)

//For a north-facing wall - Azimuth between 90 & 270 deg
for w = 1:length(Az)
    if (Az(w) > 0 & Az(w) < 90) | (Az(w) > 270 & Az(w) < 360)
    Rss(w) = Rss(w)
    else 
    Rss(w) = 0
    end
end
    //Calculation of Solar Radiation (Rs) - [W/m^2] on an hourly basis
    //Given constants
    as = 0.25           //regression constant
    bs = 0.5
    
    //Calculation of Net Shortwave Radiation (Rns) - [W/m^2]
    albedo = 0.85          //fraction of light reflected by surface
    Rs = Rss(i)
    Rns = (1-albedo)*Rs
    
    //Clear-sky solar radiation (Rso) - [W/m^2]
    Rso = Rsos(i)
    
    //Calculation of Net Longwave Radiation (Rnl) - [W/m^2]
    //Given constants for solving 
    sbc = 4.903*10^(-9)             //Stefan-Boltzmann Constant [MJ/K^4.m^2.hr]
    TmaxK = 35 + 273.15                      //max absolute temperature during 24-hr period [K]
    TminK = 20 + 273.15                       //min absolute temperature during 24-hr period [K]
    
    //Solving for vapour pressure using dew point temp
    a = 17.625          //magnus coefficient - set
    b = 243.04          //magnus coefficient - set
    RH = RHs(i)           //Annual average for Cape Town
    
    Tdew = (b*(log(RH/100)+((a*(T1-273.15))/(b+(T1-273.15)))))/(a-(log(RH/100)+((a*(T1-273.15))/(b+(T1-273.15)))))      //Dewpoint temperature calculated using Magnus-Tetens formula [deg C]
    
    ea = 0.6108*exp((17.27*Tdew)/(Tdew+237.3))      //Actual vapour pressure [kPa] - remember that Tdew is in deg. C
    
    Rnl = sbc*T1*(0.34-0.14*sqrt(ea))*(Rs/Rso)
    
    //Calculation of Net Radiation (Rn) - [W/m^2]
    Rnet = Rns-Rnl
    Rn(i) = (Rnet/2)*AG
    
    //Getting rid of Nan (from 0/0)
    nanIndices = find(isnan(Rn));
    Rn(nanIndices) = 0
    
    function f = BW(c)
        T2 = c(1)
        T3 = c(2)
        
        //Convection 1
        Tf = 0.5*(T2+T1) //Calculate Tf
        //Interpolate to find Pr and gBv value for Tf
        p = interp1(Temps,ps,Tf,'linear')
        v = interp1(Temps,vs,Tf,'linear')
        Pr = interp1(Temps,Prs,Tf,'linear')
        gBv = interp1(Temps,gBvs,Tf,'linear')
        uA = Wss(i)
       
        //Determining air frontal velocity
        h = Gh/2         //Half of wall height
        ho = 2              //Wind measurement height - 2 m
        WSE = 0.25          //Wind Shear Exponent - 0.25 for area surrounded with buildings
        V = uA*(h/ho)^WSE   //Air frontal velocity [m/s]
        
        //Reynolds Number
        Re = Gl*V*p/v
        
        if Re < 2*10^5 then
            Nu = 0.664*(Re^0.5)*(Pr^(1/3))
        else Nu = 0.0360*(Re^0.8)*(Pr^(1/3))
        end
        
        k1 = interp1(Temps,ks,Tf,'linear')
        
        h1 = Nu*k1/Gl
        
        CV1 = h1*AG*(T1-T2)
        
        //Convection 2
        Tf1 = 0.5*(T3+T4) //Calculate Tf
        //Interpolate to find Pr and gBv value for Tf
        Pr = interp1(Temps,Prs,Tf1,'linear')
        gBv = interp1(Temps,gBvs,Tf1,'linear')
        
        //Calculate h
        Gr = gBv*(Gl^3)*abs(T3-T4) //Check that its wl
        Ra = Gr*Pr
        if Ra < 10^9 then
            Nu = 0.68 + (0.670*Ra^(1/4))/((1 + (0.492/Pr)^(9/16))^4/9)
        else Nu = (0.825 + (0.387*Ra^(1/6))/((1 + (0.492/Pr)^(9/16))^(8/27)))^2
        end
        
        k3 = interp1(Temps,ks,Tf1,'linear')
        
        h2 = Nu*k3/Gl
        
        CV2 = h2*(T3-T4)
        
        CD = (kG/Gn)*(T2-T3)
        f(1) = Rn(i) + CV1 - CD
        f(2) = CD - CV2
    endfunction
    
        T20 = [T1+20;T4+20]
    [c,v,info] = fsolve(T20,BW)
    
    T2 = c(1)
    T3 = c(2)
    
    //Resolve for heat transfer
        //Convection 1
        Tf = 0.5*(T2+T1) //Calculate Tf
        //Interpolate to find Pr and gBv value for Tf
        p = interp1(Temps,ps,Tf,'linear')
        v = interp1(Temps,vs,Tf,'linear')
        Pr = interp1(Temps,Prs,Tf,'linear')
        gBv = interp1(Temps,gBvs,Tf,'linear')
        uA = Wss(i)
       
        //Determining air frontal velocity
        h = Gh/2         //Half of wall height
        ho = 2              //Wind measurement height - 2 m
        WSE = 0.25          //Wind Shear Exponent - 0.25 for area surrounded with buildings
        V = uA*(h/ho)^WSE   //Air frontal velocity [m/s]
        
        //Reynolds Number
        Re = Gl*V*p/v
        
        if Re < 2*10^5 then
            Nu = 0.664*(Re^0.5)*(Pr^(1/3))
        else Nu = 0.0360*(Re^0.8)*(Pr^(1/3))
        end
        
        k1 = interp1(Temps,ks,Tf,'linear')
        
        h1 = Nu*k1/Gl
        
        CV1 = h1*AG*(T1-T2)
        
        //Convection 2
        Tf1 = 0.5*(T3+T4) //Calculate Tf
        //Interpolate to find Pr and gBv value for Tf
        Pr = interp1(Temps,Prs,Tf1,'linear')
        gBv = interp1(Temps,gBvs,Tf1,'linear')
        
        //Calculate h
        Gr = gBv*(Gl^3)*abs(T3-T4) //Check that its wl
        Ra = Gr*Pr
        if Ra < 10^9 then
            Nu = 0.68 + (0.670*Ra^(1/4))/((1 + (0.492/Pr)^(9/16))^4/9)
        else Nu = (0.825 + (0.387*Ra^(1/6))/((1 + (0.492/Pr)^(9/16))^(8/27)))^2
        end
        
        k3 = interp1(Temps,ks,Tf1,'linear')
        
        h2 = Nu*k3/Gl
        
        CV2 = h2*(T3-T4)
        
        CD = (kG/Gn)*(T2-T3)
    T2x(i) = T2 - 273.15
    T3x(i) = T3 - 273.15
    T1x(i) = T1 - 273.15
    T4x(i) = T4 - 273.15
    CV1x(i) = CV1
    CDx(i) = CD
    CV2x(i) = CV2
    if CV2 > 0 
        CV2C(i) = CV2
        CV2H(i) = 0
    else 
        CV2C(i) = 0
        CV2H(i) = CV2
        end
end

mprintf("J = %1.0f\n",J)
TotalH(r,q) = sum(CV2H)
TotalC(r,q) = sum(CV2C)
r = r + 1
if r > Months(q)
    r = 1
    q = q + 1
end
end 
