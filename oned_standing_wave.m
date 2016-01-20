%
Nl=100;
Nt=100;
NtCycl=10;
jNM1=1;
jN=2;
jNP1=3;
uVals=zeros(3,Nl);
dt2=(1./Nt)^2;
dx2=(1./Nl)^2;
c2=1;

% Set initial values
j=0;
for i=0:Nl-1
 uVals(jNM1,i+1)=sin(2*pi/Nl*(i+0.5))*cos(2*pi/Nt*(j));
end
j=1;
for i=0:Nl-1
 uVals(jN,i+1)=sin(2*pi/Nl*(i+0.5))*cos(2*pi/Nt*(j));
end
uVals

% Set reference solution
uRef=zeros(Nt,Nl);
for j=0:Nt-1
 for i=0:Nl-1
  uRef(j+1,i+1)=sin(2*pi/Nl*(i+0.5))*cos(2*pi/Nt*(j));
 end
end

% for j=1:Nt
%  plot(uRef(j,:));
%  pnam=sprintf('foo%d.jpeg',j);
%  axis([0 Nl -1 1]);
%  print('-djpeg',pnam);
% end

% Calculate numerical solution
for j=3:Nt*NtCycl
 for i=1:Nl
  iM1=i-1;
  if iM1 == 0
   iM1=Nl;
  end
  iP1=i+1;
  if iP1>Nl
   iP1=1;
  end
  v1=2*uVals(jN,i)-uVals(jNM1,i);
  v2=dt2*c2*(uVals(jN,iP1)+uVals(jN,iM1)-2*uVals(jN,i))/dx2;
  uVals(jNP1,i)=v1+v2;
 end
 jNM1orig=jNM1;
 jNM1=jN;
 jN=jNP1;
 jNP1=jNM1orig;
end

clf
subplot(2,1,1);
plot(uVals(jN,:))
hold on
plot(uRef(Nt,:),'r')
pnam=sprintf('fooNum.jpeg',j);
axis([0 Nl -1 1]);
subplot(2,1,2);
plot(uRef(Nt,:)-uVals(jN,:))
print('-djpeg',pnam);
err=uRef(Nt,:)-uVals(jN,:);
err1=norm(err,1);
err2=norm(err,2);
erri=norm(err,Inf);
fprintf(1,'err[1,2,Inf] = %d, %d, %d\n',err1,err2,erri);
