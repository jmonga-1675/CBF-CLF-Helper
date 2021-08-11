function val = psi_gen(x,y)
a = 0.9;
val = -sign(y).*abs(y).^a-sign(phi(x,y,a)).*abs(phi(x,y,a)).^(a/(2-a));

end

function val = phi(x,y,a)

val = x + 1/(2-a)*sign(y).*abs(y).^(2-a);

end