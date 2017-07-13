%
% basic JPL shift function from the L1C ATBD
%

function delta_tb = jpl_shift_l1b(Tb_in, v_in, v_nom);

persistent d1 a b

if isempty(d1)
   d1 = load('umbc_shift_1b');
   a = d1.a; b = d1.b;
end

dv = v_nom - v_in;

% Just uncomment lines under Yibo or L Strow to pick interpolation approach

%% Yibo's approach
% Tb_spline = interp1(v_in, Tb_in, v_nom, 'spline');
% 
% Tb_resamp = Tb_in + (a .* (Tb_spline - Tb_in) ./ dv + b) .* dv;

% L Strow approach
pp = csape(v_in, Tb_in);
fprime = fnder(pp);
btderiv = fnval(fprime,v_nom);
delta_tb =  (a.*btderiv + b).*dv;
