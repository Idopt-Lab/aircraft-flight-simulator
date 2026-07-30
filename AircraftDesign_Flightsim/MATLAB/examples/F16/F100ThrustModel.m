function [T,info]=F100ThrustModel(thr,M,alt_m,varargin)
if ~(isnumeric(thr)&&isscalar(thr)&&isfinite(thr)),error("F100ThrustModel:InvalidThrottle","thr must be finite scalar.");end
if ~(isnumeric(M)&&isscalar(M)&&isfinite(M)),error("F100ThrustModel:InvalidMach","M must be finite scalar.");end
if ~(isnumeric(alt_m)&&isscalar(alt_m)&&isfinite(alt_m)),error("F100ThrustModel:InvalidAltitude","alt_m must be finite scalar.");end
persistent mi ai ti mm am tm ma aa ta
if isempty(mi)
mi=[0 0.2 0.4 0.6 0.8 1.0];ai=[-10000 0 10000 20000 30000 40000 50000 60000];
ti=[0.0430 0.0488 0.0528 0.0694 0.0899 0.1183 0.1467 0;0.0500 0.0501 0.0335 0.0544 0.0797 0.1049 0.1342 0;0.0040 0.0047 0.0020 0.0272 0.0595 0.0891 0.1203 0;-0.0804 -0.0804 -0.0560 -0.0237 0.0276 0.0718 0.1073 0;-0.2129 -0.2129 -0.1498 -0.1025 0.0474 0.0868 0.0900 0;-0.2839 -0.2839 -0.1104 -0.0469 -0.0270 0.0552 0.0800 0];
mm=[0 0.2 0.4 0.6 0.8 1.0 1.2 1.4];am=ai;
tm=[1.2600 1.0000 0.7400 0.5340 0.3720 0.2410 0.1490 0;1.1710 0.9340 0.6970 0.5060 0.3550 0.2310 0.1430 0;1.1500 0.9210 0.6920 0.5060 0.3570 0.2330 0.1450 0;1.1810 0.9510 0.7210 0.5320 0.3780 0.2480 0.1540 0;1.2580 1.0200 0.7820 0.5820 0.4170 0.2750 0.1700 0;1.3690 1.1200 0.8710 0.6510 0.4750 0.3150 0.1950 0;1.4850 1.2300 0.9750 0.7440 0.5450 0.3640 0.2250 0;1.5941 1.3400 1.0860 0.8450 0.6280 0.4240 0.2630 0];
ma=[0 0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6 1.8 2.0 2.2 2.4 2.6];aa=ai;
ta=[1.1816 1.0000 0.8184 0.6627 0.5280 0.3756 0.2327 0;1.1308 0.9599 0.7890 0.6406 0.5116 0.3645 0.2258 0;1.1150 0.9474 0.7798 0.6340 0.5070 0.3615 0.2240 0;1.1284 0.9589 0.7894 0.6420 0.5134 0.3661 0.2268 0;1.1707 0.9942 0.8177 0.6647 0.5309 0.3784 0.2345 0;1.2411 1.0529 0.8648 0.7017 0.5596 0.3983 0.2467 0;1.3287 1.1254 0.9221 0.7462 0.5936 0.4219 0.2614 0;1.4365 1.2149 0.9933 0.8021 0.6360 0.4509 0.2794 0;1.5711 1.3260 1.0809 0.8700 0.6874 0.4860 0.3011 0;1.7301 1.4579 1.1857 0.9512 0.7495 0.5289 0.3277 0;1.8314 1.5700 1.3086 1.0474 0.8216 0.5786 0.3585 0;1.9700 1.6900 1.4100 1.2400 0.9100 0.6359 0.3940 0;2.0700 1.8000 1.5300 1.3400 1.0000 0.7200 0.4600 0;2.2000 1.9200 1.6400 1.4400 1.1000 0.8000 0.5200 0];
end
tf=min(max(thr,0),1);hf=alt_m/0.3048;hu=min(max(hf,-10000),60000);Mn=max(M,0);Mi=min(max(Mn,mi(1)),mi(end));Mm=min(max(Mn,mm(1)),mm(end));Ma=min(max(Mn,ma(1)),ma(end));
Fi=interp2(ai,mi,ti,hu,Mi,"linear");Fm=interp2(am,mm,tm,hu,Mm,"linear");Fa=interp2(aa,ma,ta,hu,Ma,"linear");
Ti=Fi*17800*4.44822;Tm=Fm*17800*4.44822;Ta=Fa*29000*4.44822;
if tf<=0.5,wi=1-2*tf;wm=2*tf;wa=0;else,wi=0;wm=2-2*tf;wa=2*tf-1;end
T=wi*Ti+wm*Tm+wa*Ta;ec=[-M;wi*(M-mi(end));wm*(M-mm(end));wa*(M-ma(end));-10000-hf;hf-60000;-thr;thr-1];ev=max([0;ec]);
mx=ma(end);if wi>1e-12,mx=min(mx,mi(end));end;if wm>1e-12,mx=min(mx,mm(end));end
info=struct("valid",ev<=1e-12,"envelope_c",ec,"envelope_violation",ev,"active_max_mach",mx,"mach_requested",M,"altitude_ft_requested",hf, ...
"altitude_ft_used",hu,"throttle_requested",thr,"throttle_used",tf,"weights",[wi wm wa],"T_idle_N",Ti,"T_mil_N",Tm,"T_aug_N",Ta);
end
