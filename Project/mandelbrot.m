%{
z=0
loop ...
  z=z^2+c

Treating the real and imaginary parts of C as image coordinates 
on the complex plane, pixels may then be coloured according to how
soon the sequence crosses an arbitrarily chosen threshold, with a 
special color (usually black) used for the values of C for which the 
sequence is convergent (this is, it has not crossed the threshold 
after the predetermined number of iterations, niter ). 
%}

clear
close all

n=5;
niter=10;

x0 = -2;   x1 = 0.47;
y0 = -1.12; y1 = 1.12;

[x,y] = meshgrid(linspace(x0, x1, n), linspace(y0, y1, n));


loop1=0;
loop2=0;
loop3=1;
loop4=0;

% Original Matlab loop with complex numbers
if loop1==1
    tic
    k=zeros(n,n);
    for px=1:n
        for py=1:n
            c=x(px,py)+1i*y(px,py);
            z=0;
            k(px,py)=0; % assume it is not diverging
            for ii=1:niter
                z=z^2+c;
                if abs(z)>2 
                    k(px,py)=niter-ii;
                    break;
                end
            end
        end
    end
    toc
    %sum(k(:))
end

% Without complex numbers
if loop2==1
    tic
    k=zeros(n,n);    
    for px=1:n
        for py=1:n
            cx=x(px,py);
            cy=y(px,py);
            zx=0;
            zy=0;
            k(px,py)=0; % assume it is not diverging
            for ii=1:niter
                zx2 = cx + zx*zx-zy*zy;
                zy2 = cy + zx*zy*2;
                zx=zx2;
                zy=zy2;
                if zx*zx+zy*zy > 4 
                    k(px,py)=niter-ii;
                    break;
                end
            end
        end
    end
    toc
    %sum(k(:))
end

% loop3: with vectors of coordinates instead of matrices

cxv = linspace(x0, x1, n);
cyv = linspace(y0, y1, n);
if loop3==1
    tic
    k=zeros(n,n);    
    for px=1:n
        for py=1:n
            cx=cxv(px);
            cy=cyv(py);
            zx=0;
            zy=0;
            k(py,px)=0; % assume it is not diverging
            for ii=1:niter
                zx2 = cx + zx*zx-zy*zy;
                zy2 = cy + zx*zy*2;
                zx=zx2;
                zy=zy2;
                if zx*zx+zy*zy > 4 
                    k(py,px)=niter-ii;
                    break;
                end
            end
        end
    end
    toc
    %sum(k(:))
end
k=k';
% vectorized loop: all the points of the set are evaluated together
if loop4==1
    tic
    k=zeros(n,n);
    
    c = x + 1i * y;
    z = zeros(size(c));
    k = zeros(size(c));
    z = zeros(size(c));

    for ii = 1:niter
        z   = z.^2 + c;
        % If abs(z)>2: it has diverged so we assign it a value that decreases with iteration number
        % For the next iteration, k will be != 0 so the value assigned won't
        % change
        k(abs(z) > 2 & k == 0) = niter - ii; 
    end
    toc
    %sum(k(:))
end



figure,
imagesc(k),
colormap hot
axis square