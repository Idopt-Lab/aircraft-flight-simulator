function C=F16Lookup(x,u,g)
if ~isnumeric(x)||numel(x)<12,error("F16Lookup:InvalidState","x must contain at least 12 numeric states.");end
if ~isnumeric(u),error("F16Lookup:InvalidControl","u must be numeric.");end
persistent ar er CLT CDT CMT ClrA CnpA CYpA CYrA Mm Cma Mw CDw
if isempty(ar)
ar=[-0.175 -0.087 0 0.087 0.175 0.262 0.349 0.436 0.524 0.611 0.698 0.785]; er=[-0.436 -0.218 0 0.218 0.436];
CLT=[-0.659 -0.709 -0.754 -0.792 -0.825;-0.151 -0.196 -0.238 -0.278 -0.316;0.183 0.141 0.100 0.059 0.017;0.491 0.454 0.414 0.371 0.326;0.797 0.763 0.725 0.680 0.630;1.108 1.079 1.041 0.993 0.940;1.395 1.366 1.327 1.274 1.214;1.615 1.587 1.547 1.490 1.427;1.804 1.777 1.737 1.674 1.610;1.900 1.872 1.829 1.766 1.699;1.898 1.869 1.822 1.757 1.689;1.753 1.724 1.674 1.612 1.546];
CDT=[0.217 0.174 0.156 0.181 0.230;0.094 0.055 0.041 0.062 0.101;0.081 0.040 0.021 0.039 0.076;0.106 0.061 0.040 0.057 0.101;0.166 0.119 0.096 0.114 0.158;0.252 0.203 0.182 0.202 0.240;0.404 0.362 0.347 0.371 0.416;0.628 0.588 0.577 0.601 0.637;0.875 0.840 0.826 0.852 0.880;1.127 1.095 1.084 1.102 1.125;1.365 1.334 1.326 1.338 1.356;1.517 1.487 1.478 1.482 1.489];
CMT=[0.205 0.081 -0.046 -0.174 -0.259;0.168 0.077 -0.020 -0.145 -0.202;0.186 0.107 -0.009 -0.121 -0.184;0.196 0.110 -0.005 -0.127 -0.193;0.213 0.111 -0.006 -0.129 -0.199;0.251 0.141 0.010 -0.102 -0.150;0.245 0.127 0.006 -0.097 -0.160;0.238 0.119 -0.001 -0.113 -0.167;0.252 0.133 0.014 -0.087 -0.104;0.231 0.108 0 -0.084 -0.076;0.198 0.081 -0.013 -0.069 -0.041;0.192 0.093 0.032 -0.006 -0.005];
ClrA=[-0.126 -0.026 0.063 0.113 0.208 0.230 0.319 0.437 0.680 0.100 0.447 -0.330];
CnpA=[-0.061 -0.052 -0.052 0.012 0.013 0.024 -0.050 -0.150 -0.130 -0.158 -0.240 -0.150];
CYpA=[-0.108 -0.108 -0.188 0.110 0.258 0.226 0.344 0.362 0.611 0.529 0.298 -0.227];
CYrA=[0.882 0.852 0.876 0.958 0.962 0.974 0.819 0.483 0.590 1.210 -0.493 -1.040];
Mm=[0.6 0.8 0.9 1.0 1.1 1.2]; Cma=[0 0.0974 -0.3323 -0.9626 -0.8480 -0.7907]; Mw=[0 0.81 1.1 1.8]; CDw=[0 0 0.023 0.015];
end
x=x(:);u=u(:);v=x(4:6);w=x(10:12);V=max(norm(v),1e-6);alpha=atan2(v(3),v(1));beta=atan2(v(2),hypot(v(1),v(3)));alt=max(-x(3),0);[~,a,~]=atmosisa(alt);M=V/max(a,1e-6);
da=0;de=0;dr=0;if numel(u)>=1,da=u(1);end;if numel(u)>=2,de=u(2);end;if numel(u)>=3,dr=u(3);end
[S,b,c]=g.get_reference_geometry();if ~(isfinite(S)&&S>0&&isfinite(b)&&b>0&&isfinite(c)&&c>0),error("F16Lookup:InvalidGeometry","Reference geometry must be positive.");end
amax=g.aero_mach_max;blim=g.aero_beta_limit_rad;alim=g.aileron_limit_rad;rlim=g.rudder_limit_rad;aq=min(max(alpha,ar(1)),ar(end));eq=min(max(de,er(1)),er(end));daq=min(max(da,-alim),alim);drq=min(max(dr,-rlim),rlim);
ec=[ar(1)-alpha;alpha-ar(end);er(1)-de;de-er(end);-alim-da;da-alim;-rlim-dr;dr-rlim;-blim-beta;beta-blim;-M;M-amax];ev=max([0;ec]);valid=ev<=1e-12;
CL=interp2(er,ar,CLT,eq,aq,"linear");CD=interp2(er,ar,CDT,eq,aq,"linear")+interp1(Mw,CDw,min(max(M,Mw(1)),Mw(end)),"linear");Cm=interp2(er,ar,CMT,eq,aq,"linear");
if M<=Mm(1),dCma=0;elseif M>=Mm(end),dCma=Cma(end);else,dCma=interp1(Mm,Cma,M,"linear");end
Cm=Cm+dCma*aq;CY=-1.146*beta;Cl=-0.186*beta;Cn=0.241*beta;
if V>1
ph=w(1)*b/(2*V);qh=w(2)*c/(2*V);rh=w(3)*b/(2*V);CL=CL+28.9*qh;Cm=Cm-5.23*qh;
Cl=Cl-0.443*ph+interp1(ar,ClrA,aq,"linear")*rh;Cn=Cn+interp1(ar,CnpA,aq,"linear")*ph-0.378*rh;CY=CY+interp1(ar,CYpA,aq,"linear")*ph+interp1(ar,CYrA,aq,"linear")*rh;
end
Cl=Cl+0.05*daq;Cn=Cn-0.08*drq;
C=struct("CL",CL,"CD",CD,"CY",CY,"Cl",Cl,"Cm",Cm,"Cn",Cn,"alpha_rad",alpha,"beta_rad",beta,"mach",M,"airspeed_mps",V,"altitude_m",alt, ...
"alpha_used_rad",aq,"elevator_used_rad",eq,"valid",valid,"lookup_clamped",any(abs([aq-alpha,eq-de,daq-da,drq-dr])>1e-12)||M>Mw(end), ...
"envelope_c",ec,"envelope_violation",ev,"aero_mach_max",amax);
end
