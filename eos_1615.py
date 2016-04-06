#VLE Calculation PR Equation of state
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import scipy
#P=101.325*10**3 #Atmospheric pressure 
R=8.314 # J/K-gmole
T=298.16

#Acetone
Tc1=508.1
Pc1=47*10**5
w1=.304
Z1=.232
Rc1=2.316
Tb1=329.2
Tbr1=T/Tb1
Tr1=T/Tc1
Mw1=60

 #Chloroform
Tc2=536.4
Pc2=53.7*10**5
w2=.218
Z2=.293
Rc2=3.071
Tb2=334.3
Tbr2=T/Tb2
Tr2=T/Tc2
Mw2=119.38

#Antoine equation
#Acetone
A1=7.2316	
B1=1277.03	
C1=237.23


#Chloroform
A2=6.9371	
B2=1171.2	
C2=227


def composition(A1,B1,C1,A2,B2,C2,i):
    P1sat=(10**(A1-B1/(C1+T-273.16)))*101.325*10**3/760
    #print 'P1sat',P1sat
    P2sat=(10**(A2-B2/(C2+T-273.16)))*101.325*10**3/760
    #print 'P2sat',P2sat
    x1=i*0.1
    PP1=x1*P1sat
    PP2=(1-x1)*P2sat
    Pbubble=PP1+PP2
    y1=PP1/Pbubble
    #print x1,y1
    return x1,y1,Pbubble
    
Tmix=[]
grhomix=[]
lrhomix=[]
x=[]
y=[]
st=[]
def parachor(Tc1,Pc1,Tb1,Tc2,Pc2,Tb2):
    P01=39.6431*(0.22217-2.91042*10**-3*(Rc1/Tbr1**2))*Tc1**(1.083)/Pc1**(.8333)
    #print 'P01',P01
    p1=P01*(1-Tr1)**0.37*Tr1*scipy.exp(0.3066/Tr1+.86442*Tr1**9)
    #print 'P1',p1
    P02=39.6431*(0.22217-2.91042*10**-3*(Rc2/Tbr2**2))*Tc2**(1.083)/Pc2**(.8333)
    #print 'P02',P02
    p2=P02*(1-Tr2)**0.37*Tr2*scipy.exp(0.3066/Tr2+.86442*Tr2**9)
    #print 'P2',p2
    return p1,p2
    
def density(Tc,Pc,Z,Rc,P,T):
    Tr=T/Tc
    print P
    delta=1+(0.02*(1-.92*scipy.exp(-1000*scipy.absolute(1-Tr))-0.035*(Tr-1))*(Rc-1))
    #print scipy.exp(-1000*scipy.absolute(1-Tr))
    #print 'delta',delta
    delta=1+(0.02*-0.035*(Tr-1))*(Rc-1)
    a=.42748*R**2*Tc**2.5/Pc
    #print 'a', a
    b=(.08664*R*Tc/Pc)/delta
    #print 'b', b
    a=round(a,4)
    b=round (b,9)
    def eos(rho):
        f=P-(((rho*R*T)/(1-b*rho))-(a*rho**2/(T**.5*(1+b*rho))))
        return f
    # Initial guess gas
    rho=P/(R*T)
    #print rho,1/rho
    rhog=fsolve(eos,rho)
    #print rhog,1/rhog
    #Initial guess liquid 
    rl=1/b
    #print b
    #print rl,1/b
    rhol=fsolve(eos,rl)
    print rhol,1/rhol
    return rhog,rhol
#a=density(Tc1,Pc1,Z1,R1,P,T)
#print a
    
# Mixing Rules
def mixing(Tc1,Pc1,Z1,Rc1,Tc2,Pc2,Z2,Rc2):
    
    Zc12=(Z1+Z2)/2
    #print 'Zc12',Zc12
    k12PR=1-(2*(Tc1*Tc2)**.5/(Tc1+Tc2))**Zc12
    #print 'k12PR',k12PR
    n=8*((Tc1*Tc2)/(Pc1*Pc2))**0.5
    d=((Tc1/Pc1)**(.3333)+(Tc2/Pc2)**(0.3333))**3
    #print d
    #print 'd',d
    #print 'n',n
    const1= n/d
    #const1=(8*((Tc1*Tc2)/(Pc1*Pc2))**.5)/(((Tc1/Pc1)**(1/3)+(Tc2/Pc2)**(1/3))**3)
    #print const1
    k12RM=1-(1-k12PR)*const1
    #print 'k12RM',k12RM
    Tc12=(1-k12RM)*(Tc1*Tc2)**.5
    #print'Tc12', Tc12
    Pc12=8*Tc12/((Tc1/Pc1)**(0.3333)+(Tc2/Pc2)**(0.3333))**3
    #print 'Pc12',Pc12
    R12=((Rc1)**(0.3333)+Rc2**(0.3333))**3/8
    #print R12
    for i in range (11):
        #print i
        comp=composition(A1,B1,C1,A2,B2,C2,i)
        x1=comp[0]
        x2=1-x1
        y1=comp[1]
        y2=1-y1
        
        Rcm=x1**2*Rc1+x2**2*Rc2+2*x1*x2*R12
        #print Rm
        
        #Tcm
        num=(x1**2*Tc1**2)/Pc1+(x2**2*Tc2**2)/Pc2+(2*x1*x2*Tc12**2)/Pc12
        den=(x1**2*Tc1)/Pc1+(x2**2*Tc2)/Pc2+(2*x1*x2*Tc12)/Pc12
        Tcm=num/den
        #print 'Tcm',num,den, Tcm
        Tmix.append(Tcm)
        x.append(x1)
        y.append(y1)
        
        #Pcm
        Pcm=num/den**2
        #print Pcm
        #print comp[2]
        P=comp[2]     
        #print P,T
        
        ldensitymix=density(Tcm,Pcm,Zc12,Rcm,P,T)[1]
        lrhomix.append(ldensitymix)
        
        num=(y1**2*Tc1**2)/Pc1+(y2**2*Tc2**2)/Pc2+(2*y1*y2*Tc12**2)/Pc12
        den=(y1**2*Tc1)/Pc1+(y2**2*Tc2)/Pc2+(2*y1*y2*Tc12)/Pc12
        Tcm=num/den
        #print Tcm
        Pcm=num/den**2
        gdensitymix=density(Tcm,Pcm,Zc12,Rcm,P,T)[0]
        grhomix.append(gdensitymix)
        
        Trm=T/Tcm
        #print 'Trm',Trm
        P012=((Tc12**2/Pc12)/(Pc12/Tc12)**2.333)**0.25
        #print 'P012',P012
        parameter=parachor(Tc1,Pc1,Tb1,Tc2,Pc2,Tb2)
        P1=parameter[0]
        P2=parameter[1]
        #print 'P1',P1
        #print 'P2',P2
        PL=(x1**2*(Pc1/Tc1)**(2.333)*P1**4+x2**2*(Pc2/Tc2)**(2.333)*P2**4+2*x1*x2*(Pc12/Tc12)**(2.333)*P012**4)**0.25*(x1**2*(Tc1/Pc1)+x2**2*(Tc2/Pc2)+2*x1*x2*(Tc12/Pc12))**0.5833
        PV=(y1**2*(Pc1/Tc1)**2.333*P2**4+y2**2*(Pc2/Tc2)**2.333*P2**4+2*y1*y2*(Pc12/Tc12)**2.333*P012**4)**0.25*(y1**2*(Tc1/Pc1)+y2**2*(Tc2/Pc2)+2*y1*y2*(Tc12/Pc12))**0.5833
        Mwl=x1*Mw1+x2*Mw2
        Mwg=y1*Mw1+y2*Mw2
        sigma=(((1-Trm)**0.37*Trm*scipy.exp(0.30066/Trm+0.86442*Trm**9))*(PL*ldensitymix*Mwl*10**-3-PV*gdensitymix*Mwg*10**-3))**4/(10**6)
        #print PL,PV
        st.append(sigma)
        m12=1-P012/(P1*P2)**0.5
        #print 'm12', m12
    return st
s= mixing(Tc1,Pc1,Z1,Rc1,Tc2,Pc2,Z2,Rc2)
print 'Surface Tension',st

a=mixing(Tc1,Pc1,Z1,Rc1,Tc2,Pc2,Z2,Rc2)
#for i in range(11):
    #print 'Liquid Phase Density:',x[i],'-',lrhomix[i]
    #print 'Gas Phase Density:',y[i],'-',grhomix[i]
fig=plt.figure()
ax=fig.add_subplot(111)
plt.plot(x,st,'r')#Gw
plt.plot(y,st,'g')#Gw
ax.title.set_text('Mole Fraction Vs Surface Tension')
ax.xaxis.label.set_text('Mole Fraction')
ax.yaxis.label.set_text('Surface Tension')
plt.show()


    
