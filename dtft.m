function [fax,fdata]=dtft(tdata,dt,sfre,tfre,finterval);

tdata=tdata(:).';
N=length(tdata);
fdata=[];
fax=(sfre:finterval:tfre);

%wint=blackman(N);
wint=ones(N,1);
win=transpose(wint);

for nn=fax
  
    fdata=[fdata sum(tdata.*win.*exp(-j*(0:N-1)*dt*2*pi*nn))];
end

