import numpy as np
import astropy
import astropy.coordinates
import astropy.time
import warnings

# ERFA doesn't like ancient astronomy, but I think we can forgive
# Ptolemy for the lack of arcsecond precision it warns about
warnings.filterwarnings('ignore', 'ERFA function')

def sgn(x):
    return np.sign(x)

def cosd(x):
    return np.cos(np.radians(x))

def sind(x):
    return np.sin(np.radians(x))

nsunlong  = [  0,59, 8,17,13,12,31]
mlongsun0 = [330,45, 0, 0, 0, 0, 0]
apogeesun = [ 65,30, 0, 0, 0, 0, 0]
eccsun    = [  0, 2,30, 0, 0, 0, 0]
obliquity = [ 23,51,20, 0, 0, 0, 0]
ptropsun  = [365,14,48, 0, 0, 0, 0]
rsun      = 1210

# Moon

nlongmoon   = [ 13,10,34,58,33,30,30]
nanommoon   = [ 13, 3,53,56,17,51,59]
nlatargmoon = [ 13,13,45,39,48,56,37]
nelongmoon  = [ 12,11,26,41,20,17,59]
mlongmoon0  = [ 41,22, 0, 0, 0, 0, 0]
manommoon0  = [268,49, 0, 0, 0, 0, 0]
latargmoon0 = [354,15, 0, 0, 0, 0, 0]
melongmoon0 = [ 70,37, 0, 0, 0, 0, 0]
epimoon     = [  0, 6,20, 0, 0, 0, 0]  
eccmoon     = [  0,12,29, 0, 0, 0, 0]
incmoon     = [  5, 0, 0, 0, 0, 0, 0]
psynmoon    = [ 29,31,50, 8,20, 0, 0]

# Saturn

nlongsat    = [  0, 2, 0,33,31,28,51]
nepianomsat = [  0,57, 7,43,41,43,40]
apogeesat0  = [224,10, 0, 0, 0, 0, 0]
episat      = [  0, 6,30, 0, 0, 0, 0]  
eccsat      = [  0, 3,25, 0, 0, 0, 0]
incsat0     = [  2,30, 0, 0, 0, 0, 0]
incsat1     = [  4,30, 0, 0, 0, 0, 0]
nodesat     = [ 50, 0, 0, 0, 0, 0, 0]
rsat        = 17026

# Jupiter

nlongjup    = [  0, 4,59,14,26,46,31]
nepianomjup = [  0,54, 9, 2,46,26, 0]
apogeejup0  = [152, 9, 0, 0, 0, 0, 0]
epijup      = [  0,11,30, 0, 0, 0, 0]
eccjup      = [  0, 2,45, 0, 0, 0, 0]
incjup0     = [  1,30, 0, 0, 0, 0, 0]
incjup1     = [  2,30, 0, 0, 0, 0, 0]
nodejup     = [340, 0, 0, 0, 0, 0, 0]
rjup        = 11503.5

# Mars

nlongmar    = [  0,31,26,36,53,51,33]
nepianommar = [  0,27,41,40,19,20,58]
apogeemar0  = [106,40, 0, 0, 0, 0, 0]
epimar      = [  0,39,30, 0, 0, 0, 0]
eccmar      = [  0, 6, 0, 0, 0, 0, 0]
incmar0     = [  1, 0, 0, 0, 0, 0, 0]
incmar1     = [  2,15, 0, 0, 0, 0, 0]
nodemar     = [  0, 0, 0, 0, 0, 0, 0]
rmar        = 5040


# Venus

nepianomven = [  0,36,59,25,53,11,28]
apogeeven0  = [ 46,10, 0, 0, 0, 0, 0]
epiven      = [  0,43,10, 0, 0, 0, 0]
eccven      = [  0, 1,15, 0, 0, 0, 0]
incven0     = [  0,10, 0, 0, 0, 0, 0]
incven1     = [  2,30, 0, 0, 0, 0, 0]
incven2     = [  3,30, 0, 0, 0, 0, 0]
rven        = 622.5

# Mercury

nepianommer = [  3, 6,24, 6,59,35,50]
apogeemer0  = [181,10, 0, 0, 0, 0, 0]
epimer      = [  0,22,30, 0, 0, 0, 0]
eccmer      = [  0, 3, 0, 0, 0, 0, 0]
incmer0     = [  0,45, 0, 0, 0, 0, 0]
incmer1     = [  6,15, 0, 0, 0, 0, 0]
incmer2     = [  7, 0, 0, 0, 0, 0, 0]
rmer        = 115


pl = ["Saturn","Jupiter","Mars","Sun","Venus","Mercury","Moon"]
sign = ["Ari","Tau","Gem","Cnc","Leo","Vir","Lib","Sco","Sgr","Cap","Aqr","Psc"]
diom = [0,31,62,92,122,152,181,211,241,271,302,333]
wdn = ["Sunday","Monday","Tuesday","Wednesday","Thursday","Friday","Saturday"]
sgnrl = [2,4,5,6,3,5,4,2,1,0,0,1]
orbds = [9,9,8,15,7,7,12,12]

# day1 = 1
# month1 = 1
# year1 = -746 
# hour1 = 0
# min1 = 0
# sec1 = 0

# uses mean solar time
eotcor = 0
# tropical longitudes, based on seasons
longsys = 0

def sex2dec(s):
    d = 0
    for i in range(7):
        d = s[6-i] + d/60
    return d

def degmod(x):
    return x % 360

def radmod(x):
    return x % (2*np.pi)

def asind(x):
    return np.degrees(np.arcsin(x))

def atand(x):
    return np.degrees(np.arctan(x))

def atand2(x,y):
    return np.degrees(np.arctan2(x,y))

def dec2sex(x):
    asign = " ";
    if(x < 0):
        asign = "-"
    s = abs(x)
    s0 = np.floor(s)
    f = 60*(s-s0)
    s1 = int(np.floor(f))
    s2 = int(round(60*(f-s1)))
    if(s2 == 60):
        s1 = s1+1
        s2 = 0

    if(s1 == 60):
        s0 = s0+1
        s1 = 0

    return  f"{asign}{s0};{s1:02};{s2:02}" #  asign+str(s0)+";"+s1+","+s2

def latout(epi, ecc, inc0, inc1, node, latarg, tepianom):
    lat = 0
    rho1 = epi * cosd(tepianom)
    rho2 = epi * sind(tepianom)
    rholatmax = 1 + ecc * cosd(node)
    rholatmin = 1 - ecc * cosd(node)
    rho3 = np.sqrt((rholatmax + rho1) * (rholatmax + rho1) + rho2 * rho2)
    rho4 = np.sqrt((rholatmin + rho1) * (rholatmin + rho1) + rho2 * rho2)
    latmax = (inc0 * (rho1 + rholatmax) - inc1 * rho1) / rho3
    latmin = (inc0 * (rho1 + rholatmin) - inc1 * rho1) / rho4
    carg = cosd(latarg)
    lat = carg * ((latmax + latmin) + sgn(carg) * (latmax - latmin)) / 2
    return lat



alexandria_lat = 31.2001 * astropy.units.deg
alexandria_lon = 29.9187 * astropy.units.deg
# ptolemy age 18
# midnight on the 31st of January 118CE
obs_start_date = astropy.time.Time('0118-01-31T00:00:00.0', format="isot", scale="utc")



class Body:
    oblq = sex2dec(obliquity)

    @property
    def name(self):
        return self.__class__.__name__

    def __init__(self, observing_lat=alexandria_lat, observing_lon=alexandria_lon):
        self.observing_location = astropy.coordinates.EarthLocation(lat=observing_lat, lon=observing_lon, height=0 * astropy.units.m)
        self.local_sidereal_time = observing_lon / 15

    def radec_to_altaz(self, ra, dec, jd):
        # convert ra and dec to alt az
        observing_time = astropy.time.Time(jd, format='jd')
        ra = ra * astropy.units.deg
        dec = dec * astropy.units.deg

        aa = astropy.coordinates.AltAz(location=self.observing_location, obstime=observing_time)
        coord = astropy.coordinates.SkyCoord(ra, dec)
        altaz = coord.transform_to(aa)
        return altaz.alt.value, altaz.az.value
        
    def radec_at_date(self, jd):
        xsa, ysa, zsa = self.xyz_at_date(jd)
        rasa = degmod(atand2(ysa,xsa));
        dcsa = atand(zsa/np.sqrt(xsa*xsa+ysa*ysa))
        return rasa, dcsa

    def altaz_at_date(self, jd):
        rasa, dcsa = self.radec_at_date(jd)
        return self.radec_to_altaz(rasa, dcsa, jd)
    
    def xyz_to_altaz(self, jd, x, y, z):
        rasa = degmod(atand2(y,x));
        dcsa = atand(z/np.sqrt(x*x+y*y))
        return self.radec_to_altaz(rasa, dcsa, jd)



class Sun(Body):
    nsulong = sex2dec(nsunlong)
    apogeesu = sex2dec(apogeesun)
    eccsu = sex2dec(eccsun)
    mlongsu0 = sex2dec(mlongsun0)
    manomsu0 = degmod(mlongsu0-apogeesu)
    eqsu0 = atand(eccsu*sind(manomsu0)/(1+eccsu*cosd(manomsu0)))
    tlongsu0 = degmod(mlongsu0-eqsu0)
    xsu0 = cosd(tlongsu0)
    ysu0 = cosd(Body.oblq)*sind(tlongsu0)
    rasu0 = degmod(atand2(ysu0,xsu0))

    def get_mlongsu(self, ddays):
        return degmod(self.mlongsu0+ddays*self.nsulong)

    def xyz_at_date(self, jd):
        ddays = julian_to_egyptian(jd)
        mlongsu = self.get_mlongsu(ddays)
        manomsu = degmod(mlongsu-self.apogeesu)
        eqsu = atand(self.eccsu*sind(manomsu)/(1+self.eccsu*cosd(manomsu)))
        tlongsu = degmod(mlongsu-eqsu)
        dist60su = 60*np.sqrt(1+self.eccsu*self.eccsu+2*self.eccsu*cosd(manomsu));
        distsu = rsun*dist60su/60;
        xsu = distsu * cosd(tlongsu)
        ysu = distsu * cosd(self.oblq)*sind(tlongsu)
        zsu = distsu * sind(self.oblq)*sind(tlongsu)
        return xsu, ysu, zsu



def julian_to_egyptian(jd):
    return jd - 1448637.9044675925

class Planet(Body):
    def eqplan(self, n, ecc, epi, meccanom, mepianom):
        esin = ecc * sind(meccanom)
        ecos = ecc * cosd(meccanom)
        a = ecos + np.sqrt(1 - esin * esin)
        pros = -atand(2 * esin / a)
        b = np.sqrt(a * a + 4 * esin * esin)
        fsin = epi * sind(mepianom - pros)
        fcos = epi * cosd(mepianom - pros)
        eq = atand(fsin / (b + fcos))
        eps = asind(esin)
        px = epi * cosd(mepianom) + ecc * cosd(meccanom) + cosd(eps)
        py = epi * sind(mepianom) - 2 * ecc * sind(meccanom)
        dist = np.sqrt(px * px + py * py)
        if n == 1:
            output = pros
        if n == 2:
            output = eq
        if n == 3:
            output = dist
        return output

class InnerPlanet(Planet):
    def __init__(self, sun, *args, **kwargs):
        self.sun = sun
        super().__init__(*args, **kwargs)


class Moon(Body):
    nlongmo = sex2dec(nlongmoon)
    nanommo = sex2dec(nanommoon)
    epimo = sex2dec(epimoon)
    eccmo = sex2dec(eccmoon)
    incmo = sex2dec(incmoon)
    nlatargmo = sex2dec(nlatargmoon)
    nelongmo = sex2dec(nelongmoon)
  

    def xyz_at_date(self, jd):
        ddays = julian_to_egyptian(jd)
        # switched off equation of time correction
        ddayscor = ddays

        # Calculate the geocentric position of the Moon
        mlongmo = degmod(sex2dec(mlongmoon0)+ddayscor*self.nlongmo);
        anommo = degmod(sex2dec(manommoon0)+ddayscor*self.nanommo);
        latargmo = degmod(sex2dec(latargmoon0)+ddayscor*self.nlatargmo);
        melongmo = degmod(sex2dec(melongmoon0)+ddayscor*self.nelongmo);

        esin = self.eccmo*sind(2*melongmo);
        ecos = self.eccmo*cosd(2*melongmo);
        oc = ecos+np.sqrt(1-esin*esin);
        prosmo = atand(esin/(oc+ecos));
        tanommo = degmod(anommo+prosmo);
        fsin = self.epimo*sind(tanommo);
        fcos = self.epimo*cosd(tanommo);
        eqmo = atand(fsin/(oc+fcos));
  
        tlongmo = degmod(mlongmo-eqmo);
        latargmo = degmod(latargmo-eqmo);
        latmo = self.incmo*cosd(latargmo);
        distmo = 60*np.sqrt(1+self.epimo*self.epimo+2*self.epimo*cosd(tanommo-eqmo));
        xmo = distmo * cosd(latmo)*cosd(tlongmo);
        ymo = distmo * cosd(latmo)*sind(tlongmo)*cosd(self.oblq)-sind(latmo)*sind(self.oblq);
        zmo = distmo * cosd(latmo)*sind(tlongmo)*sind(self.oblq)+sind(latmo)*cosd(self.oblq);
        return xmo, ymo, zmo


class Mercury(InnerPlanet):

    nepianomme = sex2dec(nepianommer);
    epime = sex2dec(epimer);
    eccme = sex2dec(eccmer);
    incme0 = -sex2dec(incmer0);
    incme1 = sex2dec(incmer1);
    incme2 = -sex2dec(incmer2);


    def eqme(self, n, ecc, epi, meccanom, mepianom):
        ecos = ecc * cosd(meccanom)
        esin = ecc * sind(meccanom)
        ecoscos = 2 * ecc * cosd(meccanom / 2) * cosd(3 * meccanom / 2)
        ecossin = 2 * ecc * cosd(meccanom / 2) * sind(3 * meccanom / 2)
        a = ecos + ecoscos + np.sqrt(1 - ecossin * ecossin)
        pros = -atand(esin / a)
        b = np.sqrt(a * a + esin * esin)
        fcos = epi * cosd(mepianom - pros)
        fsin = epi * sind(mepianom - pros)
        eq = atand(fsin / (b + fcos))
        gcos = ecc * (cosd(meccanom) + cosd(2 * meccanom))
        gsin = ecc * (sind(meccanom) + sind(2 * meccanom))
        pp = np.sqrt(1 - gsin * gsin) + gcos
        qq = np.sqrt(pp * pp + ecc * ecc + 2 * ecc * pp * cosd(meccanom))
        px = qq + fcos
        py = fsin
        dist = np.sqrt(px * px + py * py)
        if n == 1:
            output = pros
        if n == 2:
            output = eq
        if n == 3:
            output = dist
        return output


    def xyz_at_date(self, jd):
        ddays = julian_to_egyptian(jd)
        ddayscor = ddays
        cent = 36525
        prectab = ddayscor / cent

        apogeeme = degmod(sex2dec(apogeemer0) + prectab)
        mepime = self.sun.get_mlongsu(ddays)
        mepianomme = degmod(21 + 55 / 60 + ddayscor * self.nepianomme)
        meccanomme = degmod(mepime - apogeeme)
        
        prosme = self.eqme(1, self.eccme, self.epime, meccanomme, mepianomme)
        teccanomme = degmod(meccanomme + prosme)
        tepianomme = degmod(mepianomme - prosme)
        eqame = self.eqme(2, self.eccme, self.epime, meccanomme, mepianomme)
        tlongme = degmod(mepime + prosme + eqame)

        etame = abs(tepianomme - 180)
        pprime = abs(self.epime * cosd(etame) * sind(self.incme1))
        xprime = 0.94444 - self.epime * cosd(etame) * cosd(self.incme1)
        yprime = self.epime * sind(etame)
        oprime = np.sqrt(xprime * xprime + yprime * yprime)
        c3me = atand2(pprime, oprime)
        c6me = abs(atand2(self.epime * sind(tepianomme), 1 + self.epime * cosd(tepianomme)))
        c4me = 6.8 * c6me / 60
        xkappa0p = degmod(teccanomme + 270)
        latme1 = -sgn(cosd(tepianomme)) * c3me * cosd(xkappa0p)
        xkappa0pp = degmod(teccanomme + 180)
        w = np.where(cosd(teccanomme) > 0)
        latme2 = np.zeros_like(tepianomme)
        latme2[w] = (0.9 * sgn(sind(tepianomme)) * c4me * cosd(xkappa0pp))[w]
        w = np.where(cosd(teccanomme) < 0)
        latme2[w] = (1.1 * sgn(sind(tepianomme)) * c4me * cosd(xkappa0pp))[w]
        latme3 = self.incme0 * cosd(teccanomme) * cosd(teccanomme)
        latme = latme1 + latme2 + latme3
        dist60me = 60*self.eqme(3,self.eccme,self.epime,meccanomme,mepianomme);
        distme = rmer*dist60me/60;

        xme = distme * cosd(latme) * cosd(tlongme)
        yme = distme * cosd(latme) * sind(tlongme) * cosd(self.oblq) - sind(latme) * sind(self.oblq)
        zme = distme * cosd(latme) * sind(tlongme) * sind(self.oblq) + sind(latme) * cosd(self.oblq)
        return xme, yme, zme


class Venus(InnerPlanet):
    nepianomve = sex2dec(nepianomven);
    epive = sex2dec(epiven)
    eccve = sex2dec(eccven)
    incve0 = sex2dec(incven0)
    incve1 = -sex2dec(incven1)
    incve2 = sex2dec(incven2)

    def xyz_at_date(self, jd):
        ddays = julian_to_egyptian(jd)
        ddayscor = ddays
        cent = 36525
        prectab = ddayscor / cent

        apogeeve = degmod(sex2dec(apogeeven0) + prectab)
        mepive = self.sun.get_mlongsu(ddays)
        mepianomve = degmod(71 + 7 / 60 + ddayscor * self.nepianomve)
        meccanomve = degmod(mepive - apogeeve)

        prosve = self.eqplan(1, self.eccve, self.epive, meccanomve, mepianomve)
        teccanomve = degmod(meccanomve + prosve)
        tepianomve = degmod(mepianomve - prosve)
        eqave = self.eqplan(2, self.eccve, self.epive, meccanomve, mepianomve)
        tlongve = degmod(mepive + prosve + eqave)

        etave = abs(tepianomve - 180.)
        pprime = abs(self.epive * cosd(etave) * sind(self.incve1))
        xprime = 0.999782 - self.epive * cosd(etave) * cosd(self.incve1)
        yprime = self.epive * sind(etave)
        oprime = np.sqrt(xprime * xprime + yprime * yprime)
        c3ve = atand2(pprime, oprime)
        c6ve = abs(atand2(self.epive * sind(tepianomve), 1 + self.epive * cosd(tepianomve)))
        c4ve = 3.25 * c6ve / 60
        xkappa0p = degmod(teccanomve + 90)
        latve1 = -sgn(cosd(tepianomve)) * c3ve * cosd(xkappa0p)
        latve2 = sgn(sind(tepianomve)) * c4ve * cosd(teccanomve)
        latve3 = self.incve0 * cosd(teccanomve) * cosd(teccanomve)
        latve = latve1 + latve2 + latve3
        dist60ve = 60*self.eqplan(3,self.eccve,self.epive,meccanomve,mepianomve);
        distve = rven*dist60ve/60;

        xve = distve * cosd(latve) * cosd(tlongve)
        yve = distve * cosd(latve) * sind(tlongve) * cosd(self.sun.oblq) - sind(latve) * sind(self.sun.oblq)
        zve = distve * cosd(latve) * sind(tlongve) * sind(self.sun.oblq) + sind(latve) * cosd(self.sun.oblq)

        return xve, yve, zve


class Mars(InnerPlanet):
    nlongma = sex2dec(nlongmar)
    nepianomma = sex2dec(nepianommar)
    epima = sex2dec(epimar)
    eccma = sex2dec(eccmar)
    incma0 = sex2dec(incmar0)
    incma1 = sex2dec(incmar1)
    nodema = sex2dec(nodemar)


    def xyz_at_date(self, jd):
        ddays = julian_to_egyptian(jd)
        ddayscor = ddays
        cent = 36525
        prectab = ddayscor / cent

        apogeema = degmod(sex2dec(apogeemar0)+prectab);
        mepima = degmod(3+32/60+ddayscor*self.nlongma);
        mepianomma = degmod(327+13/60+ddayscor*self.nepianomma);  
        meccanomma = degmod(mepima-apogeema);

        prosma = self.eqplan(1,self.eccma,self.epima,meccanomma,mepianomma);
        teccanomma = degmod(meccanomma+prosma);
        tepianomma = degmod(mepianomma-prosma);
        eqama = self.eqplan(2,self.eccma,self.epima,meccanomma,mepianomma);
        tlongma = degmod(mepima+prosma+eqama);

        latargma = degmod(teccanomma+self.nodema);
        latma = latout(self.epima,self.eccma,self.incma0,self.incma1,self.nodema,latargma,tepianomma);
        dist60ma = 60*self.eqplan(3,self.eccma,self.epima,meccanomma,mepianomma);
        distma = rmar*dist60ma/60;
        xma = distma * cosd(latma)*cosd(tlongma);
        yma = distma * cosd(latma)*sind(tlongma)*cosd(self.oblq)-sind(latma)*sind(self.oblq);
        zma = distma * cosd(latma)*sind(tlongma)*sind(self.oblq)+sind(latma)*cosd(self.oblq);
        return xma, yma, zma


class Jupiter(Planet):
    nlongju = sex2dec(nlongjup);
    nepianomju = sex2dec(nepianomjup);
    epiju = sex2dec(epijup);
    eccju = sex2dec(eccjup);
    incju0 = sex2dec(incjup0);
    incju1 = sex2dec(incjup1);
    nodeju = sex2dec(nodejup);
  
    def xyz_at_date(self, jd):
        ddays = julian_to_egyptian(jd)
        ddayscor = ddays
        cent = 36525
        prectab = ddayscor / cent
        apogeeju = degmod(sex2dec(apogeejup0)+prectab);
        mepiju = degmod(184+41/60+ddayscor*self.nlongju);
        mepianomju = degmod(146+4/60+ddayscor*self.nepianomju);  
        meccanomju = degmod(mepiju-apogeeju);

        prosju = self.eqplan(1,self.eccju,self.epiju,meccanomju,mepianomju);
        teccanomju = degmod(meccanomju+prosju);
        tepianomju = degmod(mepianomju-prosju);
        eqaju = self.eqplan(2,self.eccju,self.epiju,meccanomju,mepianomju);
        tlongju = degmod(mepiju+prosju+eqaju);

        latargju = degmod(teccanomju+self.nodeju);
        latju = latout(self.epiju,self.eccju,self.incju0,self.incju1,self.nodeju,latargju,tepianomju);
        dist60ju = 60*self.eqplan(3,self.eccju,self.epiju,meccanomju,mepianomju);
        distju = rjup*dist60ju/60;


        xju = distju * cosd(latju)*cosd(tlongju);
        yju = distju * cosd(latju)*sind(tlongju)*cosd(self.oblq)-sind(latju)*sind(self.oblq);
        zju = distju * cosd(latju)*sind(tlongju)*sind(self.oblq)+sind(latju)*cosd(self.oblq);

        return xju, yju, zju

class Saturn(Planet):
    nlongsa = sex2dec(nlongsat);
    nepianomsa = sex2dec(nepianomsat);
    episa = sex2dec(episat);  
    eccsa = sex2dec(eccsat);
    incsa0 = sex2dec(incsat0);
    incsa1 = sex2dec(incsat1);
    nodesa = sex2dec(nodesat);

    def xyz_at_date(self,jd):
        ddays = julian_to_egyptian(jd)
        ddayscor = ddays
        cent = 36525
        prectab = ddayscor / cent
        apogeesa = degmod(sex2dec(apogeesat0)+prectab);
        mepisa = degmod(296+43/60+ddayscor*self.nlongsa);
        mepianomsa = degmod(34+2/60+ddayscor*self.nepianomsa);  
        meccanomsa = degmod(mepisa-apogeesa);

        prossa = self.eqplan(1,self.eccsa,self.episa,meccanomsa,mepianomsa);
        teccanomsa = degmod(meccanomsa+prossa);
        tepianomsa = degmod(mepianomsa-prossa);
        eqasa = self.eqplan(2,self.eccsa,self.episa,meccanomsa,mepianomsa);
        tlongsa = degmod(mepisa+prossa+eqasa);

        latargsa = degmod(teccanomsa+self.nodesa);
        latsa = latout(self.episa,self.eccsa,self.incsa0,self.incsa1,self.nodesa,latargsa,tepianomsa);
        dist60sa = 60*self.eqplan(3,self.eccsa,self.episa,meccanomsa,mepianomsa);
        distsa = rsat*dist60sa/60;
        xsa = distsa * cosd(latsa)*cosd(tlongsa);
        ysa = distsa * cosd(latsa)*sind(tlongsa)*cosd(self.oblq)-sind(latsa)*sind(self.oblq);
        zsa = distsa * cosd(latsa)*sind(tlongsa)*sind(self.oblq)+sind(latsa)*cosd(self.oblq);
        return xsa, ysa, zsa
    
        
