using System;
using System.IO;
using System.Numerics;

namespace Praktyka
{
    //zmienne globalne
    public class Zmienne
    {
        public Complex temp { get; set; }
        public Complex chem { get; set; }
        public double pair { get; set; }
        public Complex g0{ get; set; }
        public double band { get; set; }
        public double bandfilling { get; set; }
        public double p0 { get; set; }
        public double p1 { get; set; }
        public double p2 { get; set; }
    }

    class Program
    {

        static void Main(string[] args)
        {
            //inicjalizacja nowej klasy zmiennych
            Zmienne a = new Zmienne();

            double st = 0.002; int nstep = 20; double dd = 1.0E-06; double zdd = 1.0E-06; double dr = 0.001; // parametry
            //Complex temp; Complex chem; double pair; Complex g0; double band; double bandfilling; double p0; double p1; double p2; //zmienne
            double gap; int occupancy;
            /*  gap is the singlet Tc equation,
	            e(x,y) is the singlet order parameter function, 
	            en(x,y) is the dispersion function, 
	            g(x,y) is the spin-orbit coupling function absolute value, 
	            g0 is the spin-orbit coupling rate, 
                temp=Tc, must take temp>0,
                pair is the pair potential,
	            bandfilling is the band filling,
	            chem is the chemical potential,
	            band=t2/temp is the ratio of the next nearest and the nearest-neighbor hopping terms.*/

            double pi = Math.Acos(-1.0);
            a.p0 = 0.0;
            a.p1 = pi;
            a.p2 = 2.0 * pi;

            a.pair = 1.0; a.g0 = 0.0; a.bandfilling = 1.2; // przypisywanie wartosci startowych
            double xt = 1.0E-06; double yt = 2.01; double xc = 0.0; double yc = 10.0; // jakies zmienne xy

            //szukanie miejsc zerowych
            Complex chp = 0.0;
            a.chem = chp;
            a.temp = zbrent(Program.gap, xt, yt, zdd);
            a.chem = zbrent(Program.occupancy, xc, yc, zdd);

            while (Math.Abs(a.chem.Real - chp.Real) > dd) //petla na szukanie miejsc zerowych
            {
                chp = a.chem.Real;
                a.temp = zbrent(Program.gap, xt, yt, zdd);
                a.chem = zbrent(Program.occupancy, xc, yc, zdd);
            }

            chp = a.chem;
            Complex t0 = a.temp;
            double o0 = a.g0.Real;
            double o1 = (a.temp / t0).Real;
            double o2 = a.temp.Real;
            double o3 = a.chem.Real;

            //zapis do pliku
            string zapis = o0 + "\t" + o1 + "\t" + o2 + "\t" + o3;
            //petla na znalezienie nieuzywanej nazwy pliku
            String lokacjadanych = "";
            int numerdanych = 0;
            for (int l = 0; l < 1000; l++)
            {
                String nazwadanych = @"Dane" + l + ".txt";
                if (!File.Exists(nazwadanych))
                {
                    lokacjadanych = nazwadanych;
                    numerdanych = l;
                    break;
                }
            }
            using (StreamWriter file = new StreamWriter(lokacjadanych, true))
            {
                file.WriteLine("{0}", zapis);
                Console.WriteLine("{0}", zapis);
            }

            a.g0 = 0.23999999;
            a.chem = 0.42988643;

            //petla
            int n = 1;
            do
            {
                a.g0 = a.g0 + st;
                a.temp = zbrent(Program.gap, xt, yt, zdd);
                a.chem = zbrent(Program.occupancy, xc, yc, zdd);
                while (Math.Abs(a.chem.Real - chp.Real) > dd)
                {
                    chp = a.chem;
                    a.temp = zbrent(Program.gap, xt, yt, zdd);
                    a.chem = zbrent(Program.occupancy, xc, yc, zdd);
                }
                chp = a.chem;
                o0 = a.g0.Real;
                o1 = (a.temp / t0).Real;
                o2 = a.temp.Real;
                o3 = a.chem.Real;

                //zapis do pliku
                zapis = o0 + "\t" + o1 + "\t" + o2 + "\t" + o3;
                using (StreamWriter file = new StreamWriter(lokacjadanych, true))
                {
                    file.WriteLine("{0}", zapis);
                    Console.WriteLine("{0}", zapis);
                }

            } while (o1 > dr);
        }

        //gap
        static double gap(double z)
        {
            //nowa klasa zmiennych
            Zmienne zm = new Zmienne();
            
            double o = calka2(csinglet, zm.p0, zm.p1, z);
            o = 4.0 * zm.pair * o / Math.Abs(zm.p2*zm.p2);
            double gap = 1.0 - o;
            return gap;
        }
        
        static double csinglet(double y, double z)
        {
            //nowa klasa zmiennych
            Zmienne zm = new Zmienne();

            double csinglet = calka1(singlet, zm.p0, zm.p1, y, z);
            return csinglet;
        }

        static double singlet(double x, double y, double z)
        {
            //nowa klasa zmiennych
            Zmienne zm = new Zmienne();

            double dda = 1.0E-06; double s1; double yt; double s2;
            double e1 = e(x, y);
            double e2 = e1 * e1;
            double x1 = en(x, y) - zm.chem.Real;
            double x2 = g(x, y);
            double y1 = x1 + x2;
            double y2 = x1 - x2;
            double t2 = 2.0 * z;
            if (Math.Abs(y1) < dda)
            {
                s1 = 1.0 / (2.0 * t2);
            }
            else
            {
                yt = y1 / t2;
                s1 = th(yt) / (2.0 * y1);
            }
            if (Math.Abs(y2) < dda)
            {
                s2 = 1.0 / (2.0 * t2);
            }
            else
            {
                yt = y2 / t2;
                s2 = th(yt) / (2.0 * y2);
            }
            double singlet = e2 * (s1 + s2);
            return singlet;
        }

        static double occupancy(double z)
        {
            //nowa klasa zmiennych
            Zmienne zm = new Zmienne();

            double o = calka2(cfermiso, zm.p0, zm.p1, z);
            o = 4.0 * o / Math.Abs(zm.p2 * zm.p2);
            double occupancy = zm.bandfilling - o;
            return occupancy;
        }

        static double cfermiso(double y, double z)
        {
            //nowa klasa zmiennych
            Zmienne zm = new Zmienne();

            double cfermiso = calka1(fermiso, zm.p0, zm.p1, y, z);
            return cfermiso;
        }

        static double fermiso(double x, double y, double z)
        {
            //nowa klasa zmiennych
            Zmienne zm = new Zmienne();

            double x1 = en(x, y) - z;
            double x2 = g(x, y);
            double y1 = (x1 + x2) / zm.temp.Real;
            double y2 = (x1 - x2) / zm.temp.Real;
            double f1 = fermi(y1);
            double f2 = fermi(y2);
            double fermiso = f1 + f2;
            return fermiso;
        }

        static double fermi(double x)
        {
            double f;
            if(Math.Abs(x) < 1.0E+02)
            {
                f = 1.0 / (Math.Exp(x) + 1.0);
            }
            else
            {
                f = 0.5 * (1.0 - x / Math.Abs(x));
            }
            double fermi = f;
            return fermi;
        }

        static double th(double x)
        {
            double th1; double th; double z;
            if(Math.Abs(x) < 1.0E+02)
            {
                z = Math.Exp(x);
                th1 = (z - 1.0 / z) / (z + 1.0 / z);
            }
            else
            {
                th1 = x / Math.Abs(x);
            }
            th = th1;
            return th;
        }

        static double e(double x, double y)
        {
            double e;
            //e = Math.Cos(y) - Math.Cos(x);
            e = 1.0;
            return e;
        }

        static double en(double x, double y)
        {
            //nowa klasa zmiennych
            Zmienne zm = new Zmienne();

            double en;
            en = -2.0 * (Math.Cos(x) + Math.Cos(y));
            //en = -2.0 * (Math.Cos(x) + Math.Cos(y)) + 4.0 * zm.band * Math.Cos(x) * Math.Cos(y);
            return en;
        }

        static double g(double x, double y)
        {
            //nowa klasa zmiennych
            Zmienne zm = new Zmienne();

            double g1 = zm.g0.Real;
            double g;
            g = g1 * Math.Sqrt(Math.Sin(x) * Math.Sin(x) + Math.Sin(y) * Math.Sin(y));
            return g;
        }

        static double calka2(Func<double, double, double> func, double a, double b, double t)
        {
            //nowa klasa zmiennych
            Zmienne zm = new Zmienne();

            //wartosci hi i fi
            double[] hi1 = new double[5]; double[] fi1 = new double[5];
            for(int zmienna1 = 0; zmienna1<5; zmienna1++)
            {
                hi1[zmienna1] = hi(zmienna1);
                fi1[zmienna1] = fi(zmienna1);
            }

            double nn = 1000; double cal = 0.0; double z; int n = 1; int i = 0; double x;
            z = (b - a) / 2.0;
            do
            {
                do
                {
                    x = a + (2.0 * n - 1.0 + fi1[i]) * z / nn;
                    cal = cal + hi1[i] * func(x, t);
                    i++;
                }
                while (i < 5);
                n++;
                i = 0;
            }
            while (n < nn + 1);
            double calka2 = z * cal / nn;
            return calka2;
        }

        static double calka1(Func<double, double, double, double> func, double a, double b, double y, double t) // bierze jako zmienna funkcje
        {
            //nowa klasa zmiennych
            Zmienne zm = new Zmienne();

            //wartosci hi i fi
            double[] hi1 = new double[5]; double[] fi1 = new double[5];
            for (int zmienna1 = 0; zmienna1 < 5; zmienna1++)
            {
                hi1[zmienna1] = hi(zmienna1);
                fi1[zmienna1] = fi(zmienna1);
            }

            double nn = 1000; double cal = 0.0; double z; int n = 1; int i = 0; double x;
            z = (b - a) / 2.0;
            do
            {
                do
                {
                    x = a + (2.0 * n - 1.0 + fi1[i]) * z / nn;
                    cal = cal + hi1[i] * func(x, y, t);
                    i++;
                }
                while (i < 5);
                n++;
                i = 0;
            }
            while (n < nn + 1);
            double calka1 = z * cal / nn;
            return calka1;
        }

        static double zbrent(Func<double, double> func, double x1, double x2, double tol)
        {
            Console.WriteLine("wezwano zbrent");
            //zmienne
            int i = 1; int itmax = 101; double EPS = 3.0E-08;
            double zbrent;
            double a = x1;
            double b = x2;
            double c = 0;
            double d = 0;
            double e = 0;
            double fa = func(a);
            double fb = func(b);
            double fc;
            double toll;
            double xm;
            double s;
            double p; double q; double r;

            //granice
            if ((fa > 0 && fb > 0) || (fa < 0 && fb < 0))
            {
                Console.WriteLine("Root must be bracketed for zbrent");
                zbrent = 0.0;
                return zbrent;
            }

            c = b;
            fc = fb;
            do
            {

                if ((fb > 0 && fc > 0) || (fb < 0 && fc < 0))
                {
                    c = a;
                    fc = fa;
                    d = b - a;
                    e = d;
                }

                if(Math.Abs(fc) < Math.Abs(fb))
                {
                    a = b;
                    b = c;
                    c = a;
                    fa = fb;
                    fb = fc;
                    fc = fa;
                }

                toll = 2.0 * EPS * Math.Abs(b) + 0.5 * tol;
                xm = 0.5 * (c - b);

                if (Math.Abs(xm) <= toll || fb == 0)
                {
                    zbrent = b;
                    return zbrent;
                }

                if (Math.Abs(e) >= toll && Math.Abs(fa) > Math.Abs(fb))
                {
                    s = fb / fa;
                    if (a == c)
                    {
                        p = 2.0 * xm * s;
                        q = 1.0 - s;
                    }
                    else
                    {
                        q = fa / fc;
                        r = fb / fc;
                        p = s * (2.0 * xm * q * (q - r) - (b - a) * (r - 1.0));
                        q = (q - 1.0) * (r - 1.0) * (s - 1.0);
                    }
                    if (p > 0)
                    {
                        q = -q;
                        p = Math.Abs(p);
                    }
                    if (2.0 * p < Math.Min(3.0 * xm * q - Math.Abs(toll * q), Math.Abs(e * q)))
                    {
                        e = d;
                        d = p / q;
                    }
                    else
                    {
                        d = xm;
                        e = d;
                    }
                }
                else
                {
                    d = xm;
                    e = d;
                }
                a = b;
                fa = fb;
                if (Math.Abs(d) > toll)
                {
                    b = b + d;
                }
                else
                {
                    if (xm >= 0)
                    {
                        b = b + Math.Abs(toll);
                    }
                    else
                    {
                        b = b - Math.Abs(toll);
                    }
                }
                fb = func(b);
            } while (i < itmax);

            //maksymalna liczba iteracji
            Console.WriteLine("Zbrent exceeding maximum iterations");
            zbrent = b;
            return zbrent;
        }

        //wartosci fi
        static double fi(int a)
        {
            switch (a)
            {
                case 0:
                    return -0.90617985;
                case 1:
                    return -0.53846931;
                case 2:
                    return 0.0;
                case 3:
                    return 0.53846931;
                case 4:
                    return 0.90617985;
                default:
                    return 0;
            }
        }

        //wartosci hi
        static double hi(int a)
        {
            switch (a)
            {
                case 0:
                    return 0.23692689;
                case 1:
                    return 0.47862867;
                case 2:
                    return 0.56888889;
                case 3:
                    return 0.47862867;
                case 4:
                    return 0.23692689;
                default:
                    return 0;
            }
        }

    }
}
