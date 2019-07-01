% genmatlabtests
%
% Generate data using Brian Hargreaves' code for comparison. His code can be
% found at http://mrsrl.stanford.edu/~brian/bloch/. His code was modified
% slightly to save data for plotting (rather than actually plotting) and to
% account for left-handed rotations (which is what BlochSim.jl uses).

gentestA5b;
gentestB2c;
gentestB2d;
gentestF1a;
gentestF1b;
gentestF1c;
gentestF2c;
gentestF3a;
gentestF3c;
gentestF3d;
gentestF3f;




% http://mrsrl.stanford.edu/~brian/bloch/a5b.m
function gentestA5b
    
    % Bloch Equation Simulation, Excercise A-5b
    % -----------------------------------------
    % 

    dT = 1;		% 1ms delta-time.
    T = 1000;	% total duration
    N = ceil(T/dT)+1; % number of time steps.
    df = 10;	% Hz off-resonance.
    T1 = 600;	% ms.
    T2 = 100;	% ms.

    % ===== Get the Propagation Matrix ======

    [A,B] = freeprecess(dT,T1,T2,df);


    % ===== Simulate the Decay ======

    M = zeros(3,N);	% Keep track of magnetization at all time points.
    M(:,1)=[1;0;0];	% Starting magnetization.

    for k=2:N
        M(:,k) = A*M(:,k-1)+B;
    end;


    % ===== Save the Results ======

    save('matlabtestdata/testA5b.mat', 'M', '-v7.3');
    
end

function gentestB2c
    
    sig = sesignal(600, 100, 50, 1000, 0);
    
    save('matlabtestdata/testB2c.mat', 'sig', '-v7.3');
    
end

function gentestB2d
    
    sig = fsesignal(600, 100, 50, 1000, 0, 8);
    
    save('matlabtestdata/testB2d.mat', 'sig', '-v7.3');
    
end

% http://mrsrl.stanford.edu/~brian/bloch/f1a.m
function gentestF1a
    
    % Bloch Equation Simulation, Excercise F-1a
    % -----------------------------------------
    % 

    rf = [pi/4 pi/4];	% Complex sequence of RF rotations.
                % This will make more sense in later exercises.

    T1 = 600;       % ms.
    T2 = 100;       % ms.

    freq = [-500:5:500];	% Hz.
    sig = 0*freq;		% Allocate space.

    % Starting magnetization.

    for f=1:length(freq)
        M = [0;0;1];
        dT = 2.3;	% ms between pulses.
        M = throt(abs(rf(1)),angle(rf(1))) * M;	% RF Rotation.

        for k = 2:length(rf)
        [A,B] = freeprecess(dT,T1,T2,freq(f));
        M = A*M+B;				% Propagate to next pulse.
            M = throt(abs(rf(k)),angle(rf(k))) * M;	% RF Rotation.
        end;
        sig(f) = M(1)+i*M(2);

    end;
    
    save('matlabtestdata/testF1a.mat', 'sig', '-v7.3');
    
end

% http://mrsrl.stanford.edu/~brian/bloch/f1b.m
function gentestF1b
    
    % Bloch Equation Simulation, Excercise F-1b
    % -----------------------------------------
    % 

    rf = [pi/4 pi/4];	% Complex sequence of RF rotations.
                % This will make more sense in later exercises.

    T1 = 600;       % ms.
    T2 = 100;       % ms.

    xpos = [-20:.1:20];	% mm.
    sig = 0*xpos;		% Allocate space.
    grad = .1;		% G/cm.

    % Starting magnetization.

    for x=1:length(xpos)
        M = [0;0;1];
        dT = 2.3;	% ms between pulses.
        M = throt(abs(rf(1)),angle(rf(1))) * M;	% RF Rotation.

        for k = 2:length(rf)
        [A,B] = freeprecess(dT,T1,T2,0);	% No off-resonance.
        M = A*M+B;				% Propagate to next pulse.
        gradrot = 4258*2*pi*(dT/1000)*grad*(xpos(x)/10);  
                            % Rotation from gradient.
        M = zrot(gradrot)*M;
            M = throt(abs(rf(k)),angle(rf(k))) * M;	% RF Rotation.
        end;
        sig(x) = M(1)+i*M(2);

    end;
    
    save('matlabtestdata/testF1b.mat', 'sig', '-v7.3');
    
end

% http://mrsrl.stanford.edu/~brian/bloch/f1c.m
function gentestF1c
    
    % Bloch Equation Simulation, Excercise F-1c
    % -----------------------------------------
    % 

    rf = [pi/4 pi/4];	% Complex sequence of RF rotations.
                % This will make more sense in later exercises.

    T1 = 600;       % ms.
    T2 = 100;       % ms.

    xpos = [-20:.1:20];	% mm.
    sig = 0*xpos;		% Allocate space.
    grad = .1;		% G/cm.

    % Starting magnetization.

    for x=1:length(xpos)
        M = [0;0;1];
        dT = 2.3;	% ms between pulses.
        M = throt(abs(rf(1)),angle(rf(1))) * M;	% RF Rotation.

        for k = 2:length(rf)
        [A,B] = freeprecess(dT,T1,T2,0);	% No off-resonance.
        M = A*M+B;				% Propagate to next pulse.

        gradrot = 4258*2*pi*(dT/1000)*grad*(xpos(x)/10);  
                            % Rotation from gradient.
        M = zrot(gradrot)*M;
            M = throt(abs(rf(k)),angle(rf(k))) * M;	% RF Rotation.
        end;

        % --- These lines added from f1b ---

        [A,B] = freeprecess(dT,T1,T2,0);	% No off-resonance.
        M = A*M+B;				% Propagate .

        gradrot = 4258*2*pi*(dT/1000)*(-0.5*grad)*(xpos(x)/10);  
                                % Rotation from gradient.
        M = zrot(gradrot)*M;
        sig(x) = M(1)+i*M(2);

    end;
    
    save('matlabtestdata/testF1c.mat', 'sig', '-v7.3');
    
end

% http://mrsrl.stanford.edu/~brian/bloch/f2c.m
function gentestF2c
    
    % Bloch Equation Simulation, Excercise F-1a
    % -----------------------------------------
    % 


    t = 0:.04:6;           % in ms.
    dT = t(2)-t(1);         % in ms.
    b1 = .05*sinc(t-3);     % in G.
    rf = 2*pi*4258*b1*(dT/1000); % Rotation in radians.

    T1 = 600;       % ms.
    T2 = 100;       % ms.

    freq = [-1000:20:1000];	% Hz.
    sig = 0*freq;		% Allocate space.

    % Starting magnetization.

    for f=1:length(freq)
        M = [0;0;1];
        M = throt(abs(rf(1)),angle(rf(1))) * M;	% RF Rotation.

        for k = 2:length(rf)
        [A,B] = freeprecess(dT,T1,T2,freq(f));
        M = A*M+B;				% Propagate to next pulse.
            M = throt(abs(rf(k)),angle(rf(k))) * M;	% RF Rotation.
        end;
        [A,B] = freeprecess(dT/2,T1,T2,freq(f));
        M = A*M+B;
        sig(f) = M(1)+i*M(2);

    end;
    
    save('matlabtestdata/testF2c.mat', 'sig', '-v7.3');
    
end

function gentestF3a
    
    % Bloch Equation Simulation, Excercise F-1a
    % -----------------------------------------
    % 
    %function [msig,m]=sliceprofile(rf,grad,t,T1,T2,pos,df)

    t = [0:.0001:.006];
    x= [-20:.1:20];
    [sig,~]=sliceprofile(.05*sinc(1000*t-3),0.1*ones(size(t)),t,600,100,x,0);
    
    save('matlabtestdata/testF3a.mat', 'sig', '-v7.3');
    
end

% http://mrsrl.stanford.edu/~brian/bloch/f3c.m
function gentestF3c
    
    % Bloch Equation Simulation, Excercise F-3c
    % -----------------------------------------
    % 
    %	Almost the same as F-3b:


    t = .0001:.00005:.006;           % in s.
    rf = .05*sinc(1000*t-3);	
    grad = .1*ones(size(t));

    % Now add refocussing gradient, and append to rf and t.

    refocratio = -.52;
    grad = [grad refocratio*grad];
    t = [t t+.006];
    rf = [rf 0*rf];	

    T1 = 600;       % ms.
    T2 = 100;       % ms.
    x = [-20:1:20];
    df = 0;
    [sig,m] = sliceprofile(rf,grad,t,T1,T2,x,df);
    mz = m(3,:);
    
    save('matlabtestdata/testF3c.mat', 'sig', 'mz', '-v7.3');
    
end

% http://mrsrl.stanford.edu/~brian/bloch/f3d.m
function gentestF3d
    
    % Bloch Equation Simulation, Excercise F-3c
    % -----------------------------------------
    % 
    %	Almost the same as F-3b:


    t = .0001:.00005:.006;           % in s.
    rf = .05*sinc(1000*t-3);	
    grad = .1*ones(size(t));

    % Now add refocussing gradient, and append to rf and t.

    refocratio = -.52;
    grad = [grad refocratio*grad];
    t = [t t+.006];
    rf = [rf 0*rf];	

    T1 = 600;       % ms.
    T2 = 100;       % ms.
    x = [-20:1:20];
    df = 100;
    [sig,m] = sliceprofile(rf,grad,t,T1,T2,x,df);
    mz = m(3,:);
    
    save('matlabtestdata/testF3d.mat', 'sig', 'mz', '-v7.3');
    
end

% http://mrsrl.stanford.edu/~brian/bloch/f3f.m
function gentestF3f
    
    % Bloch Equation Simulation, Excercise F-3c
    % -----------------------------------------
    % 
    %	Almost the same as F-3b:


    t = .0001:.00005:.006;           % in s.
    rf = .05*sinc(1000*t-3);	
    rfmod = exp(2*pi*900*i*t);
    rf = [rf.*rfmod + rf.*conj(rfmod)];
    grad = .1*ones(size(t));

    % Now add refocussing gradient, and append to rf and t.

    refocratio = -.52;
    grad = [grad refocratio*grad];
    t = [t t+.006];
    rf = [rf 0*rf];	

    T1 = 600;       % ms.
    T2 = 100;       % ms.
    x = [-50:1:50];
    df = 0;
    [sig,m] = sliceprofile(rf,grad,t,T1,T2,x,df);
    mz = m(3,:);
    
    save('matlabtestdata/testF3f.mat', 'sig', 'mz', '-v7.3');
    
end

% http://mrsrl.stanford.edu/~brian/bloch/zrot.m
function Rz=zrot(phi)

    Rz = [cos(phi) sin(phi) 0;-sin(phi) cos(phi) 0; 0 0 1];

end

% http://mrsrl.stanford.edu/~brian/bloch/xrot.m
function Rx=xrot(phi)

    Rx = [1 0 0; 0 cos(phi) sin(phi);0 -sin(phi) cos(phi)];

end

% http://mrsrl.stanford.edu/~brian/bloch/yrot.m
function Ry=yrot(phi)

    Ry = [cos(phi) 0 -sin(phi);0 1 0;sin(phi) 0 cos(phi)];

end

% http://mrsrl.stanford.edu/~brian/bloch/throt.m
function Rth=throt(phi,theta)

    Rz = zrot(-theta);
    Ry = yrot(-phi);
    Rth = inv(Rz)*Ry*Rz;

end

% http://mrsrl.stanford.edu/~brian/bloch/freeprecess.m
function [Afp,Bfp]=freeprecess(T,T1,T2,df)
%
%	Function simulates free precession and decay
%	over a time interval T, given relaxation times T1 and T2
%	and off-resonance df.  Times in ms, off-resonance in Hz.

    phi = 2*pi*df*T/1000;	% Resonant precession, radians.
    E1 = exp(-T/T1);	
    E2 = exp(-T/T2);

    Afp = [E2 0 0;0 E2 0;0 0 E1]*zrot(phi);
    Bfp = [0 0 1-E1]';
end

% http://mrsrl.stanford.edu/~brian/bloch/sesignal.m
function [Msig,Mss] = sesignal(T1,T2,TE,TR,dfreq)
    % 
    %	function [Msig,Mss] = sesignal(T1,T2,TE,TR,dfreq)
    % 
    %	Calculate the steady state signal at TE for a spin-echo
    %	sequence, given T1,T2,TR,TE in ms.  Force the
    %	transverse magnetization to zero before each excitation.
    %	dfreq is the resonant frequency in Hz.  flip is in radians.
    %

    Rflip = yrot(pi/2);	% Rotation from excitation pulse (90)
    Rrefoc = xrot(pi);	% Rotation from refocusing pulse (usually 180)

    [Atr,Btr] = freeprecess(TR-TE,T1,T2,dfreq);	% Propagation TE to TR
    [Ate2,Bte2] = freeprecess(TE/2,T1,T2,dfreq);	% Propagation 0 to TE/2
                            % (same as TE/2 to TE)

    % Neglect residual transverse magnetization prior to excitation.
    Atr = [0 0 0;0 0 0;0 0 1]*Atr;		% (Just keep Mz component)

    % Let 	M1 be the magnetization just before the 90.
    %	M2 be just before the 180.
    %	M3 be at TE.
    %	M4 = M1
    %
    % then
    %	M2 = Ate2*Rflip*M1 + Bte2
    %	M3 = Ate2*Rrefoc*M2 + Bte2
    %	M4 = Atr * M3 + Btr
    %
    % Solve for M3... (Magnetization at TE)
    %

    Mss = inv(eye(3)-Ate2*Rrefoc*Ate2*Rflip*Atr) * (Bte2+Ate2*Rrefoc*(Bte2+Ate2*Rflip*Btr));

    Msig = Mss(1)+i*Mss(2);

end

% http://mrsrl.stanford.edu/~brian/bloch/fsesignal.m
function [Msig,Mss] = fsesignal(T1,T2,TE,TR,dfreq,ETL)
    
    % 
    %	function [Msig,Mss] = sesignal(T1,T2,TE,TR,dfreq,ETL)
    % 
    %	Calculate the steady state signal at TE for a multi-echo spin-echo
    %	sequence, given T1,T2,TR,TE in ms.  Force the
    %	transverse magnetization to zero before each excitation.
    %	dfreq is the resonant frequency in Hz.  flip is in radians.
    %

    Rflip = yrot(pi/2);	% Rotation from Excitation  (usually 90)
    Rrefoc = xrot(pi);	% Rotation from Refocusing (usually 180)

    [Atr,Btr] = freeprecess(TR-ETL*TE,T1,T2,dfreq);	% Propagation last echo to TR
    [Ate2,Bte2] = freeprecess(TE/2,T1,T2,dfreq);	% Propagation over TE/2

    % Neglect residual transverse magnetization prior to excitation.
    Atr = [0 0 0;0 0 0;0 0 1]*Atr;	% Retain only Mz component.


    % Since ETL varies, let's keep a "running" A and B.  We'll
    % calculate the steady-state signal just after the tip, Rflip.

    % Initial.
    A=eye(3);
    B=[0 0 0]';


    % For each echo, we "propagate" A and B:
    for k=1:ETL
        A = Ate2*Rrefoc*Ate2*A;			% TE/2 -- Refoc -- TE/2
        B = Bte2+Ate2*Rrefoc*(Bte2+Ate2*B);
    end;


    % Propagate A and B through to just after flip, and calculate steady-state.
    A = Rflip*Atr*A;
    B = Rflip*(Btr+Atr*B);

    Mss = inv(eye(3)-A)*B;	% -- Steady state is right after 90 pulse!
    M = Mss;


    % Calculate signal on each echo.
    for k=1:ETL
        M = Ate2*Rrefoc*Ate2*M + Bte2+Ate2*Rrefoc*Bte2;
        Mss(:,k)=M;
        Msig(k)=M(1)+i*M(2);
    end;

end

% http://mrsrl.stanford.edu/~brian/bloch/sliceprofile.m
function [msig,m]=sliceprofile(rf,grad,t,T1,T2,pos,df)

    gamma = 4258;
    dT = t(2)-t(1);         % s.
    rfrot = 2*pi*gamma*rf*dT; % Rotation in radians.

    pos = pos(:).';		% Make 1xN.
    msig = 0*pos;
    m = [msig;msig;msig];

    for x=1:length(pos)

        M = [0;0;1];
        [A,B] = freeprecess(1000*dT/2,T1,T2,df);

        for k = 1:length(rf)
        M = A*M+B;
        grot = 2*pi*gamma*(pos(x)/10)*grad(k)*dT/2;
        M = zrot(grot)*M;

            M = throt(abs(rfrot(k)),angle(rfrot(k))) * M;	% RF Rotation.

        M = A*M+B;
        grot = 2*pi*gamma*(pos(x)/10)*grad(k)*dT/2;
        M = zrot(grot)*M;
        end;
        m(:,x) = M;
        msig(x) = M(1)+i*M(2);

    end;
   
end
