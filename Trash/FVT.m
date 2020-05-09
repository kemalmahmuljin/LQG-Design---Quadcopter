syms z
f2 = ((5.677e-06* z^4 + 2.271e-05* z^3 + 3.406e-05 *z^2 + 2.271e-05* z + 5.677e-06)*(z-1))/(z*(z^4 - 3.975* z^3 + 5.926* z^2 - 3.926 *z + 0.9753));
fvf2 = subs(f2,1);

f6T = (3.704e-05* z^2 + 7.407e-05* z + 3.704e-05)*(z-1)/z;
f6N = z^2 - 1.975 *z + 0.9753;
f3 = f3t/f3n;
fvf3n = subs(f3t,1)
fvf3 = subs(f3,1)

%f6T = (6.25e-05*z^2 + 0.000125*z + 6.25e-05)*(z-1)/z;
%f6N = z^2 - 2*z + 1;

f6Td = diff(f6T)
f6Nd = diff(f6N)

f6Tdd = diff(f6Td);
f6Ndd = diff(f6Nd);

f6 = subs(f6Td, 1)/subs(f6Nd, 1)


f1n = z^4 - 3.975* z^3 + 5.926* z^2 - 3.926 *z + 0.9753;
fvf1n = subs(f1n,1)
f1t = 5.677e-06 *z^4 + 2.271e-05 *z^3 + 3.406e-05 *z^2 + 2.271e-05* z + 5.677e-06;
fvf1t = subs(f1t,1)