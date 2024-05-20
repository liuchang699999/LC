N = 20;
h = 10;
a = 20;
M = 20;
d = 1;
p0 = newp;
Point(p0) = {a, 0, 0, d};
For i In {1:N}
 p1 = newp;
 Point(p1) = {a*Cos(i/N*Pi/2), 0, h*Sin(i/N*Pi/2),d};
 ll = newl;
 Line(ll) = {p0, p1};
 out[i] = ll;
 p0 = p1;
EndFor 

For j In {1:M}
out = Extrude {{0, 0, 1}, {0, 0, 0}, 2*Pi/M} {
   Line{out[]}; 
};
EndFor
Coherence;
//+
//+
//+
