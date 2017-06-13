xin_2017_m;

% set_param('xin_2017','LoadInitialState','off');
% sim('xin_2017');
% 
% xInitial = xFinal;

load('xInitial');
Kp2 = 0.8;

set_param('xin_2017','LoadInitialState','on');
sim('xin_2017');