class complex
{
public:
	complex()
	{real=imag=0.0;}
	complex(double r)
	{real=r;imag=0.0;}
	complex(double r,double i)
	{real=r;imag=i;}
	friend complex conj(complex &a);

	friend double mod(complex &m)
	{
		return sqrt(m.real*m.real+m.imag*m.imag);
	}
	friend double mod2(complex &m)
	{
		return m.real*m.real+m.imag*m.imag;
	}
    friend double Real(complex &n)
	{
		return n.real;
	}
    friend double Imag(complex &n)
	{
		return n.imag;
	}
	friend complex operator+(const complex &c1,const complex &c2);
    friend complex operator-(const complex &c1,const complex &c2);
    friend complex operator*(const complex &c1,const complex &c2);
    friend complex operator/(const complex &c1,const complex &c2);
	friend void print(const complex &c);
private:
	double real,imag;
};
 complex conj(complex &a)
 {
	return complex(a.real,-1*a.imag);
 }
 complex operator+(const complex &c1,const complex &c2)
 {return complex(c1.real+c2.real,c1.imag+c2.imag);}
 complex operator-(const complex &c1,const complex &c2)
 {return complex(c1.real-c2.real,c1.imag-c2.imag);}
 complex operator*(const complex &c1,const complex &c2)
 {return complex(c1.real*c2.real-c1.imag*c2.imag,c1.real*c2.imag+c1.imag*c2.real);}
 complex operator/(const complex &c1,const complex &c2)
 {return complex((c1.real*c2.real+c1.imag*c2.imag)/
 (c2.real*c2.real+c2.imag*c2.imag),(c1.imag*c2.real-c1.real*c2.imag)/
 (c2.real*c2.real+c2.imag*c2.imag));
 }
 void print(const complex &c)
 {
	 using namespace std;
	 if(c.imag<0)
		 cout<<c.real<<c.imag<<"i";
	 else
		 cout<<c.real<<"+"<<c.imag<<"i";
 }
