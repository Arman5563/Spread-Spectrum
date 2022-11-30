function z=sine_noise(M,n,omega0)
A = 1.2*M; % can be 1100 for constant A
z=A.*sin(omega0.*n);
end