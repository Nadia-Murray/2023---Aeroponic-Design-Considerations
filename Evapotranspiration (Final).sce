close;
clear;
clc;

//Calculating the heat associated with Evapotranspiration
//Retrieving plant species properties
[fd,SST,Sheetnames,Sheetpos] = xls_open('Model Inputs.xls')
[Value,TextInd] = xls_read(fd,Sheetpos(7))
PS = (2:8)
for m = 1:length(PS)
[Value,TextInd] = xls_read(fd,Sheetpos(7))
rls = Value(PS(m),4)
LAIs = Value(PS(m),5)
ds = Value (PS(m),7)

//Given constants
Gsc = 0.0820        //(MJ/(m^2.min))
deg = -33.96        //Degrees of latitude of location (-ve for Southern hemisphere)
Bl = 1              //Green layer length [m]

//Calculating J
Date = [23 01 2022]
hr = Date(1)
Month = Date(2)
Year = Date(3)

Months = [31 28 31 30 31 30 31 31 30 31 30 31]
if int(Year/4) == Year/4  then
    Months(2) = 29
end

J = sum(Months(1:Month-1))+hr

//Access data from Excel
[Value,TextInd] = xls_read(fd,Sheetpos(5))

X = J*24-22 //Date placeholder in excel
Y = X+23

Time = Value(X:Y,4)
T1s = Value(X:Y,5)
RHs = Value(X:Y,6)
Rss = Value(X:Y,7)
Rsos = Value(X:Y,8)
Wss = Value(X:Y,9)
Tdews = Value(X:Y,10)

i = 1
for i = 1:length(T1s)
//Calculation of Solar Radiation (Rs) - [W/m^2] on an hourly basis

//Calculation of Net Shortwave Radiation (Rns) - [W/m^2]
albedo = 0.225          //fraction of light reflected by surface - Allen gives a range of 0.2-0.25 for a plant layer
Rs = Rss(i)
Rns = (1-albedo)*Rs

//Calculation of clear-sky solar radiation (Rso) - [W/m^2]
Rso = Rsos(i)

//Calculation of Net Longwave Radiation (Rnl) - [W/m^2]
//Given constants for solving 
sbc = 5.670374419*10^(-8)             //Stefan-Boltzmann Constant [W/K^4.m^2]

//Solving for vapour pressure using dew point temp
a = 17.625              //magnus coefficient - set
b = 243.04              //magnus coefficient - set
T1 = T1s(i)+273.15      //Outside temperature [K]
RH = RHs(i)             //Relative humidity per hour in CPT
Tdew = (b*(log(RH/100)+((a*(T1-273.15))/(b+(T1-273.15)))))/(a-(log(RH/100)+((a*(T1-273.15))/(b+(T1-273.15)))))      //Dewpoint temperature calculated using Magnus-Tetens formula [deg C]

ea = 0.6108*exp((17.27*Tdew)/(Tdew+237.3))      //Actual vapour pressure [kPa] - remember that Tdew is in deg. C

Rnl = sbc*(T1)*(0.34-0.14*sqrt(ea))*(Rs/Rso)
disp(Rns,Rnl)
//Calculation of Net Radiation (Rn) - [W/m^2]
Rnet = Rns-Rnl
Rn(i) = Rnet/2

//Getting rid of Nan (from 0/0)
nanIndices = find(isnan(Rn));
Rn(nanIndices) = 0


//Calculating Transpiration Rate 
slope = (4098*(0.6108*exp((17.27*(T1-273.15))/((T1-273.15)+237.3))))/(((T1-273.15)+237.3)^2) //Slope of saturation vapour pressure curve [kPa/deg C](NB! - T1 in deg C!!!)

Ktime = 60*60   //Unit time conversion [s/hr]

//Properties of air - dependent on temperature
[Value,TextInd] = xls_read(fd,Sheetpos(3))

Temps = Value(22:36,1)
ps = Value(22:36,2)
Cps = Value(22:36,3)
vs = Value(22:36,4)
ks = Value(22:36,6)
Prs = Value(22:36,8)
T = T1
pA = interp1(Temps,ps,T,'linear')       //mean density of dry air [kg/m^3]
Cp = interp1(Temps,Cps,T,'linear')      //Specific heat capacity of air [J/kg.deg C]
vA = interp1(Temps,vs,T,'linear')       //viscosity of air [Pa.s]
kA = interp1(Temps,ks,T,'linear')       //Thermal conductivity [W/m.K]
Pr = interp1(Temps,Prs,T,'linear')      //Prandtl number
lhvap = 2260*1000   //Latent heat of vapourisation of water [J/kg]


//Determining the vapour pressure deficit
es =  0.6108*exp((17.27*T1)/(T1+237.3))
VPD = es-ea

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
d = (2/3)*ch        //zero plane displacement height
zom = 0.123*ch
zh = 2
zoh = 0.1*zom
k = 0.41            //von Karman's constant

ra = (log((zm-d)/zom)*log((zh-d)/zoh))/((k^2)*u2)     //aerodynamic resistance [s/m]

//Determining canopy resistance from Yazdanseta (2017)
rl = rls
LAI = LAIs
LAIslit = 0.95*LAI   //From past thesis
rs = rl/LAIslit     //Canopy resistance [s/m]
G = 0               //Soil Heat Flux - assumed to be negligible
eThr = 0.6108*exp((17.27*(T1-273.15))/((T1-273.15)+237.3))

RnET = Rn(i)*3600         //Converting net radiation to J/m^2
RnET = RnET/(10^6)    //Converting net radiation to MJ/m^2

//Transpiration rate [MJ/m^2.hr]
ET = ((slope*(RnET-G))+(Ktime*pA*(Cp/10^6)*((eThr-ea)/ra)))/(slope+pc*(1+rs/ra)) 
ETkg = (ET*10^6)/lhvap          //[kg/m^2.hr] - this is just used to get transpiration in terms of kg - not used in next line
ETW = ET*((10^6)/3600)          //[W/m^2]


//Maximum possible cooling power of vegetation
Qmax = ETW

//Determining effectiveness and NTU to calculate Q
LAD = 10  //Assumption of between 2-10 m^2/m^3
Bn = 0.2           //Thickness of green wall (in range of 20-130)
Bh = 2           //Wall height [m]
h = Bh/2         //Half of wall height
ho = 2              //Wind measurement height - 2 m
WSE = 0.25          //Wind Shear Exponent - 0.25 for area surrounded with buildings
V = u2*(h/ho)^WSE   //Air frontal velocity [m/s]

//Calculating Dimensionless Quantities for NTU
dleaf = ds
ReNTU = (pA*V*Bl)/vA

    if dleaf <= 0.05 then 
        NuA = 1.86*(Pr^0.33)*(ReNTU^0.5)
    else dleaf > 0.05 
        NuA = 1.18*(Pr^0.33)*(ReNTU^0.5)
    end

//Determining the heat transfer coefficient for vegetation
hA = NuA*kA/Bl

//Calculating Number of Transfer Units (NTU) to calculate the 
NTU = (hA*LAD*LAIslit*Bn)/(V*pA*Cp)

//Calculating effectiveness (e)
e = 1-exp(-NTU)

//Calculating actual heat released by evapotranspiration [W/m^2]
QE = e*Qmax
Q(i) = QE
end

Qx(:,m) = Q
end

scf(1)
clf()
plot((0:23),Qx)
plot(0,24)
xlabel('Time of day [hr]')
ylabel('Heat Released from Evapotranspiration [W/m^2]')
legend('Cape Rush','Sea Rose','Carpet Daisy','Cape May','Red-stem Crassula','Spekboom','Red Grass')

