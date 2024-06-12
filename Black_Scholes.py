import math
from datetime import date
import Volatility as volobj

class Options:

    def __init__(self):
        pass

    def get_normalCDF(self, x_value):
        a1 = 0.31938153
        a2 = -0.356563782
        a3 = 1.781477937
        a4 = -1.821255978
        a5 = 1.330274429

        L = abs( x_value );
        K = 1 / (1 + 0.2316419 * L);
        w = 1 - 1 / (math.pow( 2.0 * math.pi, 0.5 )) * math.exp( -L * L / 2.0 ) * (
                    a1 * K + a2 * K * K + a3 * math.pow( K, 3.0 ) + a4 * math.pow( K, 4.0 ) + a5 * math.pow( K, 5.0 ))

        if x_value < 0:
            w = 1 - w;

        return w

    def get_call_price(self, S0, K, T, vol, q, r):
        d1 = (math.log(S0/K)+T*(r-q+vol*vol/2.0))/(vol*math.sqrt(T))
        d2 = d1 - vol * math.sqrt(T)
        call = S0*math.exp(-q*T)*self.get_normalCDF(d1)-K*math.exp(-r*T)*self.get_normalCDF(d2)
        return call

    def get_put_price(self,S0,K,T, vol, q ,r):
        return self.get_call_price(S0,K,T,vol,q,r)-S0+K*math.exp(-r*T)

    def get_down_and_in_put_price(self,S0,K,B,T,vol,q,r):
        phi = -1
        nu = 1
        b = r-q  # cost of carry
        mu = (b - vol * vol / 2.0) / (vol * vol)
        lambda1 = math.sqrt( mu * mu + 2.0 * r / (vol * vol) )
        x1 = math.log( S0 / K ) / (vol * math.sqrt( T )) + (1 + mu) * vol * math.sqrt( T )
        x2 = math.log( S0 / B ) / (vol * math.sqrt( T )) + (1 + mu) * vol * math.sqrt( T )
        y1 = math.log( B * B / (S0 * K) ) / (vol * math.sqrt( T )) + (1 + mu) * vol * math.sqrt( T )
        y2 = math.log( B / S0 ) / (vol * math.sqrt( T )) + (1 + mu) * vol * math.sqrt( T )
        z = math.log( B / K ) / (vol * math.sqrt( T )) + lambda1 * vol * math.sqrt( T )
        A = phi * S0 * math.exp( (q - r) * T ) * self.get_normalCDF( phi * x1 ) - phi * K * math.exp(
            -r * T ) * self.get_normalCDF( phi * x1 - phi * vol * math.sqrt( T ) )
        B = phi * S0 * math.exp( (q - r) * T ) * self.get_normalCDF( phi * x2 ) - phi * K * math.exp(
            -r * T ) * self.get_normalCDF( phi * x2 - phi * vol * math.sqrt( T ) )
        C = phi * S0 * math.exp( (q - r) * T ) * math.pow( B / S0, 2.0 * (mu + 1.0) ) * self.get_normalCDF(
            nu * y1 ) - phi * K * math.exp( -r * T ) * math.pow( B / S0, 2.0 * mu ) * self.get_normalCDF(
            nu * y1 - nu * vol * math.sqrt( T ) )
        D = phi * S0 * math.exp( (q - r) * T ) * math.pow( B / S0, 2.0 * (mu + 1.0) ) * self.get_normalCDF(
            nu * y2 ) - phi * K * math.exp( -r * T ) * math.pow( B / S0, 2.0 * mu ) * self.get_normalCDF(
            nu * y2 - nu * vol * math.sqrt( T ) )
        R = 0
        E = 0
        F = R * (math.pow( B / S0, mu + lambda1 )) * self.get_normalCDF( nu * z ) + math.pow( B / S0, mu - lambda1 ) * self.get_normalCDF(
            nu * z - 2 * nu * lambda1 * vol * math.sqrt( T ) )
        if K < B:
            return A + E
        elif K > B:
            return B - C + D + E

    def get_down_and_out_put_price(self,S0,K,B,T,vol,q, r):
        # reference in (Haug, 1997)
        phi = -1
        nu = 1
        b = r-q #cost of carry
        mu = (b -vol*vol/2.0)/(vol*vol)
        lambda1 = math.sqrt(mu*mu+2.0*r/(vol*vol))
        x1 = math.log(S0/K)/(vol*math.sqrt(T))+(1+mu)*vol*math.sqrt(T)
        x2 = math.log(S0/B)/(vol*math.sqrt(T))+(1+mu)*vol*math.sqrt(T)
        y1 = math.log(B*B/(S0*K))/(vol*math.sqrt(T))+ (1+mu)*vol*math.sqrt(T)
        y2 = math.log(B/S0)/(vol*math.sqrt(T))+(1+mu)*vol*math.sqrt(T)
        z = math.log(B/K)/(vol*math.sqrt(T))+lambda1*vol*math.sqrt(T)
        A = phi*S0*math.exp((q-r)*T)*self.get_normalCDF(phi*x1)-phi*K*math.exp(-r*T)*self.get_normalCDF(phi*x1-phi*vol*math.sqrt(T))
        B = phi*S0*math.exp((q-r)*T)*self.get_normalCDF(phi*x2)-phi*K*math.exp(-r*T)*self.get_normalCDF(phi*x2-phi*vol*math.sqrt(T))
        C = phi*S0*math.exp((q-r)*T)*math.pow(B/S0, 2.0*(mu+1))*self.get_normalCDF(nu*y1)-phi*K*math.exp(-r*T)*math.pow(B/S0,2.0*mu)*self.get_normalCDF(nu*y1-nu*vol*math.sqrt(T))
        D = phi*S0*math.exp((q-r)*T)*math.pow(B/S0, 2.0*(mu+1))*self.get_normalCDF(nu*y2)-phi*K*math.exp(-r*T)*math.pow(B/S0, 2.0*mu)*self.get_normalCDF(nu*y2-nu*vol*math.sqrt(T))
        R = 0
        F = R*(math.pow(B/S0,mu+lambda1))*self.get_normalCDF(nu*z)+math.pow(B/S0,mu-lambda1)*self.get_normalCDF(nu*z-2.0*nu*lambda1*vol*math.sqrt(T))
        if K < B:
            return F
        elif K > B:
            return A - B + C - D + F

    def get_binarycall_up_and_in_price(self,S0,K,B,payoff,T,vol,r,q):
        d1 = (math.log(S0/K)+(r-q+vol*vol/2.0)*T)/(vol*math.sqrt(T))
        d2 = d1-vol*math.sqrt(T)
        d3 = (math.log(S0/B)+(r-q+vol*vol/2.0)*T)/(vol*math.sqrt(T))
        d4 = d3 - vol*math.sqrt(T)
        P = 2* self.get_normalCDF(d4)
        call=payoff*math.exp(-r*T)*self.get_normalCDF(d2)*P
        return call

    def get_dailyoption_callprice(self,S0,K,vol):
        T = 1.0/365
        return self.get_call_price(S0,K,T,vol,0,0)

    def get_dailyoption_putprice(self,S0,K,vol):
        T = 1.0/365
        return self.get_put_price(S0,K,T,vol,0,0)

def main():
    print( "this is the main" )
    myobj = Options()
    myvol = volobj.parseVol( 'clean.csv' )
    myvol.parsefile()
    myvol.initialize()
    S = 1174.75
    K = 1200
    B = 1100
    sigma = myvol.interpolate( K )/100.0
    rate = 5.48/100.0
    print("vol= {0}".format(sigma))
    expiry = myvol.get_expiry()
    tokens = expiry.split( "/" )
    d0 = date(2024,6,12)
    d1 = date(int(tokens[2]),int(tokens[1]),int(tokens[0]))
    T = ((d1-d0).days)/365.0
    print("T= {0}".format(T))
    print("Put= {0}".format(myobj.get_put_price( S,K, T, sigma, 0, rate )))
    print("Down and out Put= {0}".format(myobj.get_down_and_out_put_price(S,K,B,T,sigma,0,rate)))
    print("Down and in Put= {0}".format(myobj.get_down_and_in_put_price( S, K, B, T, sigma, 0, rate )))
    print("Binary up and in call= {0}".format(myobj.get_binarycall_up_and_in_price(S,K,1350,20,T,sigma,rate,0)))
    print("DailyOptionCall= {0}".format(myobj.get_dailyoption_callprice(S,K,sigma)))
    print( "DailyOptionPut= {0}".format( myobj.get_dailyoption_putprice( S, K, sigma ) ) )


if __name__ == "__main__":
    main()
