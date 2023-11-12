//Base Case: Aeroponic System
//Initialize script
close;
clear;
clc;

Date = [23 01 2022]

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
2 -> Surface of the green layer
3 -> External surface of Aeroponic Wall 1
4 -> Internal surface of Aeroponic Wall 1
5 -> Internal surface of Aeroponic Wall 2
6 -> External surface of Aeroponic Wall 2
7 -> External wall
8 -> Internal wall
9 -> Internal Environment
*/

//PLANT LAYER INCLUDED
//Layer dimensions
L = 1
H = 1
A = L*H

//Green Layer (B)
AB = 0.25 //Cross-sectional ace area [m^2]
Bl = 0.5 //Length [m]
Bh = 0.5 //Height[m]
Bn = 0.2 //Thickness [m]
dleaf = 0.00371

//Aeroponic Wall 1 (C)
Cl = Bl //Length [m]
Ch = Bh //Height[m]
Cn = 0.00635 //Thickness [m]
AC = Cl*Ch //Cross-sectional ace area [m^2]

//Aeroponic Chamber (D)
Dl = Bl //Length [m]
Dh = Bh //Height[m]
Dn = 0.1 //Thickness [m]
DH = (4*Dl*Dn)/(2*Dl+2*Dn)
AD = Dl*Dh //Cross-sectional ace area [m^2]

//Aeroponic Wall 2 (E)
El = Bl //Length [m]
Eh = Bh //Height[m]
En = 0.00635 //Thickness [m]
AE = El*Eh //Cross-sectional ace area [m^2]

//Air Gap (F)
Fl = Bl //Length [m]
Fh = Bh //Height[m]
Fn = 0.05 //Thickness [m]
AF = Fl*Fh //Cross-sectional ace area [m^2]

//Wall (G)
Gl = Bl //Wall length [m]
Gh = Bh //Wall width [m]
Gn = 0.22 //Wall thickness [m] 
AG = Gl*Gh //Wall face area [m^2]

//Thermal Conductivity 
kB = 0.35 //Green layer [W/mK]
kC = 0.19 //PVC [W/mK]
kE = kC //PVC [W/mK]
kG = 1 //Wall [W/mK]

//Additional inputs
Tin = 21
T9 = Tin + 273.15  //Internal room temperature [degC]
Gsc = 0.0820        //(MJ/(m^2.min))
deg = -33.96        //Degrees of latitude of location (-ve for Southern hemisphere)
rw = 3*10^(-3)      //Radius of water droplet [m]
g = 9.81            //Graviational aceleration [m/s2]

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
Temps = Value(22:36,1)
ps = Value(22:36,2)
Cps = Value(22:36,3)
vs = Value(22:36,4)
ks = Value(22:36,6)
Prs = Value(22:36,8)
gBvs = Value(22:36,9)

//Access physical properties of moist air from excel
[Value,TextInd] = xls_read(fd,Sheetpos(6))
mclose();
Tempsm = Value(3:15,2)
ps100 = Value(3:15,3)
ps0 = Value(3:15,4)
Cps100 = Value(3:15,5)
Cps0 = Value(3:15,6)
vs100 = Value(3:15,7)
vs0 = Value(3:15,8)
ks100 = Value(3:15,9)
ks0 = Value(3:15,10)

i = 1
for i = 1:length(T1s)
    
    T1 = T1s(i) + 273.15 //Outside temperature [degC]

    //Radiation for Plants
    //Calculation of Solar Radiation (Rs) - [W/m^2] on an hourly basis
    //Given constants
    as = 0.25           //regression constant
    bs = 0.5
    
    //Calculation of Net Shortwave Radiation (Rns) - [W/m^2]
    albedo = 0.225          //fraction of light reflected by surface
    Rs = Rss(i)
    Rns = (1-albedo)*Rs
    
    //Clear-sky solar radiation (Rso) - [W/m^2]
    Rso = Rsos(i)
    
    //Calculation of Net Longwave Radiation (Rnl) - [W/m^2]
    //Given constants for solving 
    sbc = 4.903*10^(-9)             //Stefan-Boltzmann Constant [MJ/K^4.m^2.hr]
    TmaxK = 35+273.15                      //max absolute temperature during 24-hr period [K]
    TminK = 20+273.15                       //min absolute temperature during 24-hr period [K]
    
    //Solving for vapour pressure using dew point temp
    a = 17.625          //magnus coefficient - set
    b = 243.04          //magnus coefficient - set
    RH = RHs(i)           //Annual average for Cape Town
    
    Tdew = (b*(log(RH/100)+((a*(T1-273.15))/(b+(T1-273.15)))))/(a-(log(RH/100)+((a*(T1-273.15))/(b+(T1-273.15)))))      //Dewpoint temperature calculated using Magnus-Tetens formula [deg C]
    
    ea = 0.6108*exp((17.27*Tdew)/(Tdew+237.3))      //Actual vapour pressure [kPa] - remember that Tdew is in deg. C
    
    Rnl = sbc*T1*(0.34-0.14*sqrt(ea))*(Rs/Rso)
    
    //Calculation of Net Radiation (Rn) - [W/m^2]
    Rnet = Rns-Rnl
    Rn(i) = (Rnet/2)*AB
    
    //Getting rid of Nan (from 0/0)
    nanIndices = find(isnan(Rn));
    Rn(nanIndices) = 0
    
    function f = BW(c)
        T2 = c(1)
        T3 = c(2)
        T4 = c(3)
        T5 = c(4)
        T6 = c(5)
        T7 = c(6)
        T8 = c(7)

        //With plant coverage 
        //Convection A to B
        Tf = 0.5*(T1+T2) //Calculate Tf
        pA = interp1(Temps,ps,Tf,'linear')
        vA = interp1(Temps,vs,Tf,'linear')
        Pr = interp1(Temps,Prs,Tf,'linear')
        Cp = interp1(Temps,Cps,Tf,'linear')
        uA = Wss(i)
        //Determining air frontal velocity
        LAD = 10  //Assumption of between 2-10 m^2/m^3
        h = Bh/2         //Half of wall height
        ho = 2              //Wind measurement height - 2 m
        WSE = 0.25          //Wind Shear Exponent - 0.25 for area surrounded with buildings
        V = uA*(h/ho)^WSE   //Air frontal velocity [m/s]
        ReA = Bl*V*pA/vA
            if dleaf <= 0.05 then 
        NuA = 1.86*(Pr^0.33)*(ReA^0.5)
            else dleaf > 0.05 
        NuA = 1.18*(Pr^0.33)*(ReA^0.5)
        end
        kA = interp1(Temps,ks,Tf,'linear')
        hA = NuA*kA/Bl
        CVA = hA*AB*(T1-T2)
        
        //Calculating the heat associated with Evapotranspiration
        //Calculating Transpiration Rate 
        slope = (4098*(0.6108*exp((17.27*(T1-273.15))/((T1-273.15)+237.3))        ))/(((T1-273.15)+237.3)^2) //Slope of saturation vapour pressure c        urve [kPa/deg C](NB! - Tinf in deg C!!!)
        Ktime = 60*60   //Unit time conversion [s/hr]
        lhvap = 2260*1000   //Latent heat of vapourisation of water [J/kg]
        //Determining the saturation vapour pressure deficit
        TmaxC = TmaxK-273.15
        TminC = TminK-273.15
        e0Tmax = 0.6108*exp((17.27*TmaxC)/(TmaxC+237.3))
        e0Tmin = 0.6108*exp((17.27*TminC)/(TminC+237.3))
        es = (e0Tmax+e0Tmin)/2
        
        //Determining the psychrometric constant
        P = 101.325     //Atmospheric pressure [kPa]
        rmw = 0.622     //ratio molecular weight of water vapour over dry air 
        pc = (Cp*P)/(rmw*lhvap) //Psychrometric constant [kPa/deg C]
        //Determining aerodynamic resistance (Allen et al.)
        //Adjusting wind speed to be at wall height
        z = 10              //height of wind speed measurement
        u10 = Wss(i)        //wind speed measurement at 10 m
        u2 = u10*(4.87/(log((67.8*z)-5.42)))    //wind speed at 2 m above ground [m/s]
        zm = 2
        ch = 0.1            //crop height [m] - Should vary this
        d = (2/3)*ch
        zom = 0.123*ch
        zh = 2
        zoh = 0.1*zom
        k = 0.41            //von Karman's constant
        ra = (log((zm-d)/zom)*log((zh-d)/zoh))/((k^2)*u2)     //aerodynamic resistance [s/m]
        //Determining canopy resistance from Yazdanseta (2017)
        rl = 100        //Allen et al assumes this for grass species - need to find data for different plant species - use past thesis
        LAI = 1.7           //Just a general guess for now - species specific
        LAIslit = 0.95*LAI   //From past thesis
        rs = rl/LAIslit     //Canopy resistance [s/m]
        G = 0               //Soil Heat Flux - assumed to be negligible
        eThr = 0.6108*exp((17.27*(T1-273.15))/((T1-273.15)+237.3))
        
        RnET = Rn(i)*3600         //Converting net radiation to J/m^2
        RnET = RnET/(10^6)    //Converting net radiation to MJ/m^2
        
        //Transpiration rate [MJ/m^2.hr]
        ET = ((slope*(RnET-G))+(Ktime*pA*(Cp/10^6)*((eThr-ea)/ra)))/(slope+pc*(1+rs/ra)) 
        ETkg = (ET*10^6)/lhvap          //[kg/m^2.hr]
        ETW = ET*((10^6)/3600)          //[W/m^2]
        
        
        //Maximum possible cooling power of vegetation
        Qmax = ETW
        
        //Calculating Number of Transfer Units (NTU) to calculate the 
        NTU = (hA*LAD*LAIslit*Bn)/(V*pA*Cp)
        
        //Calculating effectiveness (e)
        e = 1-exp(-NTU)
        
        //Calculating actual heat released by evapotranspiration [W/m^2]
        QE = e*Qmax

        //Conduction through B
        CDB = (AB*kB/Bn)*(T2-T3)
        
        //Conduction through C
        CDC = (AC*kC/Cn)*(T3-T4)
        
        //Convection C to E (Through D) (Natural)
        Tf = 0.5*(T4+T5) //Calculate Tf
        Pr = interp1(Temps,Prs,Tf,'linear')
        gBv = interp1(Temps,gBvs,Tf,'linear')
        Gr = gBv*(Dn^3)*abs(T4 - T5)
        Ra = Gr*Pr
        Ratio = Dh/Dn
        if 1 < Ratio && Ratio < 2 || Ratio == 1 || Ratio == 2 then
            Nu = 0.18*((Pr/(0.2+Pr))*Ra)^0.29
        elseif 2 < Ratio && Ratio < 10 || Ratio == 10 then
            Nu = 0.22*(((Pr/(0.2+Pr))*Ra)^0.29)*(Ratio^(-0.25))
        elseif 10 < Ratio && Ratio < 40 then
            if Ra < 10^6.5
                Nu = 0.42*(Ra^(0.25))*(Pr^(0.012))*(Ratio^(-0.3))
            else Nu = 0.046*(Ra^(1/3))
            end 
        end
        
        //Adjusting Nu for moist air 
        p100 = interp1(Tempsm,ps100,Tf,'linear')
        p0 = interp1(Tempsm,ps0,Tf,'linear')
        Cp100 = interp1(Tempsm,Cps100,Tf,'linear')
        Cp0 = interp1(Tempsm,Cps0,Tf,'linear')
        v100 = interp1(Tempsm,vs100,Tf,'linear')
        v0 = interp1(Tempsm,vs0,Tf,'linear')
        k100 = interp1(Tempsm,ks100,Tf,'linear')
        k0 = interp1(Tempsm,ks0,Tf,'linear')
        
        n = 1/4 //Laminar Flow assumed
        Num = Nu*((p100/p0)^(2*n))*((v0/v100)^n)*((Cp100/Cp0)^n)*((k0/k100)^n)
        hD = Num*k100/Dn
        CVD = hD*AD*(T4-T5)
        
        //Conduction through E
        CDE = (AE*kE/En)*(T5-T6)

        //Convection E to G (Through F) (Natural)
        Tf = 0.5*(T6+T7)
        Pr = interp1(Temps,Prs,Tf,'linear')
        gBv = interp1(Temps,gBvs,Tf,'linear')
        Gr = gBv*(Fn^3)*abs(T6 - T7)
        Ra = Gr*Pr
        Ratio = Fh/Fn
        if 1 < Ratio && Ratio < 2 || Ratio == 1 || Ratio == 2 then
            Nu = 0.18*((Pr/(0.2+Pr))*Ra)^0.29
        elseif 2 < Ratio && Ratio < 10 || Ratio == 10 then
            Nu = 0.22*(((Pr/(0.2+Pr))*Ra)^0.29)*(Ratio^(-0.25))
        elseif 10 < Ratio && Ratio < 40 then
            if Ra < 10^6.5
                Nu = 0.42*(Ra^(0.25))*(Pr^(0.012))*(Ratio^(-0.3))
            else Nu = 0.046*(Ra^(1/3))
            end 
        end
        k = interp1(Temps,ks,Tf,'linear')
        hE = Nu*k/Fn
        CVF = hE*AE*(T6-T7)
        
        //Conduction through G
        CDG = (AG*kG/Gn)*(T7-T8)
        
        //Convection G to H
        Tf = 0.5*(T8+T9) //Calculate Tf
        //Interpolate to find Pr and gBv value for Tf
        Pr = interp1(Temps,Prs,Tf,'linear')
        gBv = interp1(Temps,gBvs,Tf,'linear')
        Gr = gBv*(Gl^3)*abs(T8-T9) //Check that its wl
        Ra = Gr*Pr
        if Ra < 10^9 then
            Nu = 0.68 + (0.670*Ra^(1/4))/((1 + (0.492/Pr)^(9/16))^4/9)
        else Nu = (0.825 + (0.387*Ra^(1/6))/((1 + (0.492/Pr)^(9/16))^(8/27)))^2
        end
        kH = interp1(Temps,ks,Tf,'linear')
        hH = Nu*kH/Gl
        CVH = hH*AG*(T8-T9)
        
        f(1) = Rn(i) + CVA - QE - CDB
        f(2) = CDB - CDC
        f(3) = CDC - CVD
        f(4) = CVD - CDE
        f(5) = CDE - CVF
        f(6) = CVF - CDG
        f(7) = CDG - CVH

    endfunction
    T0 = [T1;T1;T1;T1+10;T1+10;T9;T9]//Guess initial Ts
    [c,v,info] = fsolve(T0,BW)
    
        T2 = c(1)
        T3 = c(2)
        T4 = c(3)
        T5 = c(4)
        T6 = c(5)
        T7 = c(6)
        T8 = c(7)

        T1xP(i) = T1 - 273.15
        T2xP(i) = T2 - 273.15
        T3xP(i) = T3 - 273.15
        T4xP(i) = T4 - 273.15
        T5xP(i) = T5 - 273.15
        T6xP(i) = T6 - 273.15
        T7xP(i) = T7 - 273.15
        T8xP(i) = T8 - 273.15
        T9xP(i) = T9 - 273.15
        
        //Resolve for heat transfer
        //Convection A to B
        Tf = 0.5*(T1+T2) //Calculate Tf
        pA = interp1(Temps,ps,Tf,'linear')
        vA = interp1(Temps,vs,Tf,'linear')
        Pr = interp1(Temps,Prs,Tf,'linear')
        Cp = interp1(Temps,Cps,Tf,'linear')
        uA = Wss(i)
        //Determining air frontal velocity
        LAD = 10  //Assumption of between 2-10 m^2/m^3
        h = Bh/2         //Half of wall height
        ho = 2              //Wind measurement height - 2 m
        WSE = 0.25          //Wind Shear Exponent - 0.25 for area surrounded with buildings
        V = uA*(h/ho)^WSE   //Air frontal velocity [m/s]
        ReA = Bl*V*pA/vA
            if dleaf <= 0.05 then 
        NuA = 1.86*(Pr^0.33)*(ReA^0.5)
            else dleaf > 0.05 
        NuA = 1.18*(Pr^0.33)*(ReA^0.5)
        end
        kA = interp1(Temps,ks,Tf,'linear')
        hA = NuA*kA/Bl
        CVA = hA*AB*(T1-T2)
        
        //Calculating the heat associated with Evapotranspiration
        //Calculating Transpiration Rate 
        slope = (4098*(0.6108*exp((17.27*(T1-273.15))/((T1-273.15)+237.3))        ))/(((T1-273.15)+237.3)^2) //Slope of saturation vapour pressure c        urve [kPa/deg C](NB! - Tinf in deg C!!!)
        Ktime = 60*60   //Unit time conversion [s/hr]
        lhvap = 2260*1000   //Latent heat of vapourisation of water [J/kg]
        //Determining the saturation vapour pressure deficit
        TmaxC = TmaxK-273.15
        TminC = TminK-273.15
        e0Tmax = 0.6108*exp((17.27*TmaxC)/(TmaxC+237.3))
        e0Tmin = 0.6108*exp((17.27*TminC)/(TminC+237.3))
        es = (e0Tmax+e0Tmin)/2
        
        //Determining the psychrometric constant
        P = 101.325     //Atmospheric pressure [kPa]
        rmw = 0.622     //ratio molecular weight of water vapour over dry air 
        pc = (Cp*P)/(rmw*lhvap) //Psychrometric constant [kPa/deg C]
        //Determining aerodynamic resistance (Allen et al.)
        //Adjusting wind speed to be at wall height
        z = 10              //height of wind speed measurement
        u10 = Wss(i)        //wind speed measurement at 10 m
        u2 = u10*(4.87/(log((67.8*z)-5.42)))    //wind speed at 2 m above ground [m/s]
        zm = 2
        ch = 0.1            //crop height [m] - Should vary this
        d = (2/3)*ch
        zom = 0.123*ch
        zh = 2
        zoh = 0.1*zom
        k = 0.41            //von Karman's constant
        ra = (log((zm-d)/zom)*log((zh-d)/zoh))/((k^2)*u2)     //aerodynamic resistance [s/m]
        //Determining canopy resistance from Yazdanseta (2017)
        rl = 100        //Allen et al assumes this for grass species - need to find data for different plant species - use past thesis
        LAI = 1.7           //Just a general guess for now - species specific
        LAIslit = 0.95*LAI   //From past thesis
        rs = rl/LAIslit     //Canopy resistance [s/m]
        G = 0               //Soil Heat Flux - assumed to be negligible
        eThr = 0.6108*exp((17.27*(T1-273.15))/((T1-273.15)+237.3))
        
        RnET = Rn(i)*3600         //Converting net radiation to J/m^2
        RnET = RnET/(10^6)    //Converting net radiation to MJ/m^2
        
        //Transpiration rate [MJ/m^2.hr]
        ET = ((slope*(RnET-G))+(Ktime*pA*(Cp/10^6)*((eThr-ea)/ra)))/(slope+pc*(1+rs/ra)) 
        ETkg = (ET*10^6)/lhvap          //[kg/m^2.hr]
        ETW = ET*((10^6)/3600)          //[W/m^2]
        
        
        //Maximum possible cooling power of vegetation
        Qmax = ETW
        
        //Calculating Number of Transfer Units (NTU) to calculate the 
        NTU = (hA*LAD*LAIslit*Bn)/(V*pA*Cp)
        
        //Calculating effectiveness (e)
        e = 1-exp(-NTU)
        
        //Calculating actual heat released by evapotranspiration [W/m^2]
        QE = e*Qmax

        //Conduction through B
        CDB = (AB*kB/Bn)*(T2-T3)
        
        //Conduction through C
        CDC = (AC*kC/Cn)*(T3-T4)
        
        //Convection C to E (Through D) (Natural)
        Tf = 0.5*(T4+T5) //Calculate Tf
        Pr = interp1(Temps,Prs,Tf,'linear')
        gBv = interp1(Temps,gBvs,Tf,'linear')
        Gr = gBv*(Dn^3)*abs(T4 - T5)
        Ra = Gr*Pr
        Ratio = Dh/Dn
        if 1 < Ratio && Ratio < 2 || Ratio == 1 || Ratio == 2 then
            Nu = 0.18*((Pr/(0.2+Pr))*Ra)^0.29
        elseif 2 < Ratio && Ratio < 10 || Ratio == 10 then
            Nu = 0.22*(((Pr/(0.2+Pr))*Ra)^0.29)*(Ratio^(-0.25))
        elseif 10 < Ratio && Ratio < 40 then
            if Ra < 10^6.5
                Nu = 0.42*(Ra^(0.25))*(Pr^(0.012))*(Ratio^(-0.3))
            else Nu = 0.046*(Ra^(1/3))
            end 
        end
        
        //Adjusting Nu for moist air 
        p100 = interp1(Tempsm,ps100,Tf,'linear')
        p0 = interp1(Tempsm,ps0,Tf,'linear')
        Cp100 = interp1(Tempsm,Cps100,Tf,'linear')
        Cp0 = interp1(Tempsm,Cps0,Tf,'linear')
        v100 = interp1(Tempsm,vs100,Tf,'linear')
        v0 = interp1(Tempsm,vs0,Tf,'linear')
        k100 = interp1(Tempsm,ks100,Tf,'linear')
        k0 = interp1(Tempsm,ks0,Tf,'linear')
        
        n = 1/4 //Laminar Flow assumed
        Num = Nu*((p100/p0)^(2*n))*((v0/v100)^n)*((Cp100/Cp0)^n)*((k0/k100)^n)
        hD = Num*k100/Dn
        CVD = hD*AD*(T4-T5)
        
        //Conduction through E
        CDE = (AE*kE/En)*(T5-T6)

        //Convection E to G (Through F) (Natural)
        Tf = 0.5*(T6+T7)
        Pr = interp1(Temps,Prs,Tf,'linear')
        gBv = interp1(Temps,gBvs,Tf,'linear')
        Gr = gBv*(Fn^3)*abs(T6 - T7)
        Ra = Gr*Pr
        Ratio = Fh/Fn
        if 1 < Ratio && Ratio < 2 || Ratio == 1 || Ratio == 2 then
            Nu = 0.18*((Pr/(0.2+Pr))*Ra)^0.29
        elseif 2 < Ratio && Ratio < 10 || Ratio == 10 then
            Nu = 0.22*(((Pr/(0.2+Pr))*Ra)^0.29)*(Ratio^(-0.25))
        elseif 10 < Ratio && Ratio < 40 then
            if Ra < 10^6.5
                Nu = 0.42*(Ra^(0.25))*(Pr^(0.012))*(Ratio^(-0.3))
            else Nu = 0.046*(Ra^(1/3))
            end 
        end
        k = interp1(Temps,ks,Tf,'linear')
        hE = Nu*k/Fn
        CVF = hE*AE*(T6-T7)
        
        //Conduction through G
        CDG = (AG*kG/Gn)*(T7-T8)
        
        //Convection G to H
        Tf = 0.5*(T8+T9) //Calculate Tf
        //Interpolate to find Pr and gBv value for Tf
        Pr = interp1(Temps,Prs,Tf,'linear')
        gBv = interp1(Temps,gBvs,Tf,'linear')
        Gr = gBv*(Gl^3)*abs(T8-T9) //Check that its wl
        Ra = Gr*Pr
        if Ra < 10^9 then
            Nu = 0.68 + (0.670*Ra^(1/4))/((1 + (0.492/Pr)^(9/16))^4/9)
        else Nu = (0.825 + (0.387*Ra^(1/6))/((1 + (0.492/Pr)^(9/16))^(8/27)))^2
        end
        kH = interp1(Temps,ks,Tf,'linear')
        hH = Nu*kH/Gl
        CVH = hH*AG*(T8-T9)
        
        CVAxP(i) = CVA
        CDBxP(i) = CDB
        CDCxP(i) = CDC
        CVDxP(i) = CVD
        CDExP(i) = CDE
        CVFxP(i) = CVF
        CDGxP(i) = CDG
        CVHxP(i) = CVH/AG
        QExP(i) = QE
        
    if info ~= 1
    mprintf("info = %1.0f\n",info)
    mprintf("i = %1.0f\n",i)
    mprintf("m = %1.0f\n",m)
    end
end 

//NO PLANT LAYER
//Layer dimensions
L = 1
H = 1
A = L*H

//Green Layer (B)
AB = 0.25 //Cross-sectional ace area [m^2]
Bl = 0.5 //Length [m]
Bh = 0.5 //Height[m]
Bn = 0.2 //Thickness [m]
dleaf = 0.00371

//Aeroponic Wall 1 (C)
AC = 0.25 //Cross-sectional ace area [m^2]
Cl = 0.5 //Length [m]
Ch = 0.5 //Height[m]
Cn = 0.00635 //Thickness [m]
AC = Cl*Ch 

//Aeroponic Chamber (D)
Dl = Cl //Length [m]
Dh = Ch //Height[m]
Dn = 0.1016 //Thickness [m]
DH = (4*Dl*Dn)/(2*Dl+2*Dn)
AD = Dl*Dh //Cross-sectional ace area [m^2]

//Aeroponic Wall 2 (E)
El = Cl //Length [m]
Eh = Ch //Height[m]
En = 0.00635 //Thickness [m]
AE = El*Eh //Cross-sectional ace area [m^2]

//Air Gap (F)
Fl = Cl //Length [m]
Fh = Ch //Height[m]
Fn = 0.05 //Thickness [m]
AF = Fl*Fh //Cross-sectional ace area [m^2]

//Wall (G)
Gl = Cl //Wall length [m]
Gh = Ch //Wall width [m]
Gn = 0.22 //Wall thickness [m] 
AG = Gl*Gh //Wall face area [m^2]

//Thermal Conductivity 
kB = 0.35 //Green layer [W/mK]
kC = 0.19 //PVC [W/mK]
kE = kC //PVC [W/mK]
kG = 1 //Wall [W/mK]

//Additional inputs
Tin = 21
T9 = Tin + 273.15  //Internal room temperature [degC]
Gsc = 0.0820        //(MJ/(m^2.min))
deg = -33.96        //Degrees of latitude of location (-ve for Southern hemisphere)
rw = 3*10^(-3)      //Radius of water droplet [m]
g = 9.81            //Graviational aceleration [m/s2]

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
Temps = Value(22:36,1)
ps = Value(22:36,2)
Cps = Value(22:36,3)
vs = Value(22:36,4)
ks = Value(22:36,6)
Prs = Value(22:36,8)
gBvs = Value(22:36,9)

//Access physical properties of moist air from excel
[Value,TextInd] = xls_read(fd,Sheetpos(6))
mclose();
Tempsm = Value(3:15,2)
ps100 = Value(3:15,3)
ps0 = Value(3:15,4)
Cps100 = Value(3:15,5)
Cps0 = Value(3:15,6)
vs100 = Value(3:15,7)
vs0 = Value(3:15,8)
ks100 = Value(3:15,9)
ks0 = Value(3:15,10)

i = 1
for i = 1:length(T1s)
    
    T1 = T1s(i) + 273.15 //Outside temperature [degC]

    //Radiation without Plants
    //Calculation of Solar Radiation (Rs) - [W/m^2] on an hourly basis
    //Given constants
    as = 0.25           //regression constant
    bs = 0.5
    
    //Calculation of Net Shortwave Radiation (Rns) - [W/m^2]
    albedo = 0.44          //fraction of light reflected by surface
    Rs = Rss(i)
    Rns = (1-albedo)*Rs
    
    //Clear-sky solar radiation (Rso) - [W/m^2]
    Rso = Rsos(i)
    
    //Calculation of Net Longwave Radiation (Rnl) - [W/m^2]
    //Given constants for solving 
    sbc = 4.903*10^(-9)             //Stefan-Boltzmann Constant [MJ/K^4.m^2.hr]
    TmaxK = 35+273.15                      //max absolute temperature during 24-hr period [K]
    TminK = 20+273.15                       //min absolute temperature during 24-hr period [K]
    
    //Solving for vapour pressure using dew point temp
    a = 17.625          //magnus coefficient - set
    b = 243.04          //magnus coefficient - set
    RH = RHs(i)           //Annual average for Cape Town
    
    Tdew = (b*(log(RH/100)+((a*(T1-273.15))/(b+(T1-273.15)))))/(a-(log(RH/100)+((a*(T1-273.15))/(b+(T1-273.15)))))      //Dewpoint temperature calculated using Magnus-Tetens formula [deg C]
    
    ea = 0.6108*exp((17.27*Tdew)/(Tdew+237.3))      //Actual vapour pressure [kPa] - remember that Tdew is in deg. C
    
    Rnl = sbc*T1*(0.34-0.14*sqrt(ea))*(Rs/Rso)
    
    //Calculation of Net Radiation (Rn) - [W/m^2]
    Rnet = Rns-Rnl
    Rn(i) = (Rnet/2)*AB
    
    //Getting rid of Nan (from 0/0)
    nanIndices = find(isnan(Rn));
    Rn(nanIndices) = 0
    
    function f = BW2(c)
        T3 = c(1)
        T4 = c(2)
        T5 = c(3)
        T6 = c(4)
        T7 = c(5)
        T8 = c(6)

        //With no plant coverage 
        //Convection A to C
        Tf = 0.5*(T1+T3) //Calculate Tf
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
        
        CVA = h1*AC*(T1-T3)

        //Conduction through C
        CDC = (AC*kC/Cn)*(T3-T4)
        
        //Convection C to E (Through D) (Natural)
        Tf = 0.5*(T4+T5) //Calculate Tf
        Pr = interp1(Temps,Prs,Tf,'linear')
        gBv = interp1(Temps,gBvs,Tf,'linear')
        Gr = gBv*(Dn^3)*abs(T4 - T5)
        Ra = Gr*Pr
        Ratio = Dh/Dn
        if 1 < Ratio && Ratio < 2 || Ratio == 1 || Ratio == 2 then
            Nu = 0.18*((Pr/(0.2+Pr))*Ra)^0.29
        elseif 2 < Ratio && Ratio < 10 || Ratio == 10 then
            Nu = 0.22*(((Pr/(0.2+Pr))*Ra)^0.29)*(Ratio^(-0.25))
        elseif 10 < Ratio && Ratio < 40 then
            if Ra < 10^6.5
                Nu = 0.42*(Ra^(0.25))*(Pr^(0.012))*(Ratio^(-0.3))
            else Nu = 0.046*(Ra^(1/3))
            end 
        end
        
        //Adjusting Nu for moist air 
        p100 = interp1(Tempsm,ps100,Tf,'linear')
        p0 = interp1(Tempsm,ps0,Tf,'linear')
        Cp100 = interp1(Tempsm,Cps100,Tf,'linear')
        Cp0 = interp1(Tempsm,Cps0,Tf,'linear')
        v100 = interp1(Tempsm,vs100,Tf,'linear')
        v0 = interp1(Tempsm,vs0,Tf,'linear')
        k100 = interp1(Tempsm,ks100,Tf,'linear')
        k0 = interp1(Tempsm,ks0,Tf,'linear')
        
        n = 1/4 //Laminar Flow assumed
        Num = Nu*((p100/p0)^(2*n))*((v0/v100)^n)*((Cp100/Cp0)^n)*((k0/k100)^n)
        hD = Num*k100/Dn
        CVD = hD*AD*(T4-T5)
        
        //Conduction through E
        CDE = (AE*kE/En)*(T5-T6)

        //Convection E to G (Through F) (Natural)
        Tf = 0.5*(T6+T7)
        Pr = interp1(Temps,Prs,Tf,'linear')
        gBv = interp1(Temps,gBvs,Tf,'linear')
        Gr = gBv*(Fn^3)*abs(T6 - T7)
        Ra = Gr*Pr
        Ratio = Fh/Fn
        if 1 < Ratio && Ratio < 2 || Ratio == 1 || Ratio == 2 then
            Nu = 0.18*((Pr/(0.2+Pr))*Ra)^0.29
        elseif 2 < Ratio && Ratio < 10 || Ratio == 10 then
            Nu = 0.22*(((Pr/(0.2+Pr))*Ra)^0.29)*(Ratio^(-0.25))
        elseif 10 < Ratio && Ratio < 40 then
            if Ra < 10^6.5
                Nu = 0.42*(Ra^(0.25))*(Pr^(0.012))*(Ratio^(-0.3))
            else Nu = 0.046*(Ra^(1/3))
            end 
        end
        k = interp1(Temps,ks,Tf,'linear')
        hE = Nu*k/Fn
        CVF = hE*AE*(T6-T7)
        
        //Conduction through G
        CDG = (AG*kG/Gn)*(T7-T8)
        
        //Convection G to H
        Tf = 0.5*(T8+T9) //Calculate Tf
        //Interpolate to find Pr and gBv value for Tf
        Pr = interp1(Temps,Prs,Tf,'linear')
        gBv = interp1(Temps,gBvs,Tf,'linear')
        Gr = gBv*(Gl^3)*abs(T8-T9) //Check that its wl
        Ra = Gr*Pr
        if Ra < 10^9 then
            Nu = 0.68 + (0.670*Ra^(1/4))/((1 + (0.492/Pr)^(9/16))^4/9)
        else Nu = (0.825 + (0.387*Ra^(1/6))/((1 + (0.492/Pr)^(9/16))^(8/27)))^2
        end
        kH = interp1(Temps,ks,Tf,'linear')
        hH = Nu*kH/Gl
        CVH = hH*AG*(T8-T9)
        
        f(1) = Rn(i) + CVA - CDC
        f(2) = CDC - CVD
        f(3) = CVD - CDE
        f(4) = CDE - CVF
        f(5) = CVF - CDG
        f(6) = CDG - CVH

    endfunction
    T0 = [T1+20;T1+10;T1+10;T1+5;T9+20;T9+10]//Guess initial Ts
    [c,v,info] = fsolve(T0,BW2)

        T3 = c(1)
        T4 = c(2)
        T5 = c(3)
        T6 = c(4)
        T7 = c(5)
        T8 = c(6)

        T1x(i) = T1 - 273.15
        T3x(i) = T3 - 273.15
        T4x(i) = T4 - 273.15
        T5x(i) = T5 - 273.15
        T6x(i) = T6 - 273.15
        T7x(i) = T7 - 273.15
        T8x(i) = T8 - 273.15
        T9x(i) = T9 - 273.15
        
        //Resolve for heat transfer
        //With no plant coverage 
        //Convection A to C
        Tf = 0.5*(T3+T1) //Calculate Tf
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
        
        CVA = h1*AG*(T1-T3)

        //Conduction through C
        CDC = (AC*kC/Cn)*(T3-T4)
        
        //Convection C to E (Through D) (Natural)
        Tf = 0.5*(T4+T5) //Calculate Tf
        Pr = interp1(Temps,Prs,Tf,'linear')
        gBv = interp1(Temps,gBvs,Tf,'linear')
        Gr = gBv*(Dn^3)*abs(T4 - T5)
        Ra = Gr*Pr
        Ratio = Dh/Dn
        if 1 < Ratio && Ratio < 2 || Ratio == 1 || Ratio == 2 then
            Nu = 0.18*((Pr/(0.2+Pr))*Ra)^0.29
        elseif 2 < Ratio && Ratio < 10 || Ratio == 10 then
            Nu = 0.22*(((Pr/(0.2+Pr))*Ra)^0.29)*(Ratio^(-0.25))
        elseif 10 < Ratio && Ratio < 40 then
            if Ra < 10^6.5
                Nu = 0.42*(Ra^(0.25))*(Pr^(0.012))*(Ratio^(-0.3))
            else Nu = 0.046*(Ra^(1/3))
            end 
        end
        
        //Adjusting Nu for moist air 
        p100 = interp1(Tempsm,ps100,Tf,'linear')
        p0 = interp1(Tempsm,ps0,Tf,'linear')
        Cp100 = interp1(Tempsm,Cps100,Tf,'linear')
        Cp0 = interp1(Tempsm,Cps0,Tf,'linear')
        v100 = interp1(Tempsm,vs100,Tf,'linear')
        v0 = interp1(Tempsm,vs0,Tf,'linear')
        k100 = interp1(Tempsm,ks100,Tf,'linear')
        k0 = interp1(Tempsm,ks0,Tf,'linear')
        
        n = 1/4 //Laminar Flow assumed
        Num = Nu*((p100/p0)^(2*n))*((v0/v100)^n)*((Cp100/Cp0)^n)*((k0/k100)^n)
        hD = Num*k100/Dn
        CVD = hD*AD*(T4-T5)
        
        //Conduction through E
        CDE = (AE*kE/En)*(T5-T6)

        //Convection E to G (Through F) (Natural)
        Tf = 0.5*(T6+T7)
        Pr = interp1(Temps,Prs,Tf,'linear')
        gBv = interp1(Temps,gBvs,Tf,'linear')
        Gr = gBv*(Fn^3)*abs(T6 - T7)
        Ra = Gr*Pr
        Ratio = Fh/Fn
        if 1 < Ratio && Ratio < 2 || Ratio == 1 || Ratio == 2 then
            Nu = 0.18*((Pr/(0.2+Pr))*Ra)^0.29
        elseif 2 < Ratio && Ratio < 10 || Ratio == 10 then
            Nu = 0.22*(((Pr/(0.2+Pr))*Ra)^0.29)*(Ratio^(-0.25))
        elseif 10 < Ratio && Ratio < 40 then
            if Ra < 10^6.5
                Nu = 0.42*(Ra^(0.25))*(Pr^(0.012))*(Ratio^(-0.3))
            else Nu = 0.046*(Ra^(1/3))
            end 
        end
        k = interp1(Temps,ks,Tf,'linear')
        hE = Nu*k/Fn
        CVF = hE*AE*(T6-T7)
        
        //Conduction through G
        CDG = (AG*kG/Gn)*(T7-T8)
        
        //Convection G to H
        Tf = 0.5*(T8+T9) //Calculate Tf
        //Interpolate to find Pr and gBv value for Tf
        Pr = interp1(Temps,Prs,Tf,'linear')
        gBv = interp1(Temps,gBvs,Tf,'linear')
        Gr = gBv*(Gl^3)*abs(T8-T9) //Check that its wl
        Ra = Gr*Pr
        if Ra < 10^9 then
            Nu = 0.68 + (0.670*Ra^(1/4))/((1 + (0.492/Pr)^(9/16))^4/9)
        else Nu = (0.825 + (0.387*Ra^(1/6))/((1 + (0.492/Pr)^(9/16))^(8/27)))^2
        end
        kH = interp1(Temps,ks,Tf,'linear')
        hH = Nu*kH/Gl
        CVH = hH*AG*(T8-T9)
        
        CVAx(i) = CVA
        CDCx(i) = CDC
        CVDx(i) = CVD
        CDEx(i) = CDE
        CVFx(i) = CVF
        CDGx(i) = CDG
        CVHx(i) = CVH
        
    if info ~= 1
    mprintf("W\n")
    mprintf("info = %1.0f\n",info)
    mprintf("i = %1.0f\n",i)
    mprintf("m = %1.0f\n",m)
    end
end 
CVHxW = CVHx/AG

WCs = 0:0.1:1

for m = 1:length(WCs)
    WC = WCs(m)
CV(:,m) = (1-WC)*A*CVHxW + WC*A*CVHxP
end 

//Plot
scf(1)
clf()
x2 = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23]
plot(x2,CV')
xlabel('Time [Hour]')
ylabel('Q [W]')
legend('0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1')
