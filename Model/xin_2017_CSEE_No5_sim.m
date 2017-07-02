xin_2017_CSEE_No5_m;

set_param('xin_2017_CSEE_No5','LoadInitialState','off');
sim('xin_2017_CSEE_No5');

xInitial = xFinal;

load('xInitial');
Kp2 = 0.8;

set_param('xin_2017_CSEE_No5','LoadInitialState','on');
sim('xin_2017');