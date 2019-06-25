function [B, cv] = direct(xx, y, h, nd)
%Searching the effective dimension reduction subspace of model
%  			y = g(B^TX) + e
%Useage: directions(x, y, h, d)
%Input: 
%     x --- expanaltory variables
%     y --- response
%     h --- bandwidth
%     d --- working dimension of the space
%Output:
%     Directions B
%Reference: Y. Xia, H. Tong, W.K. Li and L.X. Zhu, 
% "adaptive estimation of the effictive dimension space",(2001)  
%-------------------------------------------------------------
% Example
%	n = 100;
%	x = randn(n,4);
%	beta1 = [1 2 3 0]';
%	beta1 = beta1/sqrt(beta1'*beta1);
%	beta2 = [-2 1 0 1]';
%	beta2 = beta2/sqrt(beta2'*beta2);
%	y = (x*beta1).^2 + x*beta2 + 0.2*randn(n,1);
%	B = rMAVE(x, y, 0.5, 2)
%  % estimation errors
%	B0 = [beta1 beta2];
%	error = B'*(eye(4)-B0*B0')*B

[n,p] = size(xx);
[ny, py] = size(y); 
mm = mean(xx);
xx = xx - repmat(mm,n,1);
ss = inv(xx'*xx/n)^0.5;
xx = xx*ss;


if (ny ~= n)
   disp('Error: matrix of x and y dont match')
   return
end
if (p <= 1)
   disp('Error: please check your matrix x')
   return
end
if (py > 1)
   disp('Error: please check your matrix y')
   return
end
if (sum(abs(mean(xx,1))) > 0) | (abs(mean(std(xx,1))-1)>0) 
%   disp('Warning: covariates must be standardized')
end;   
b = (1:p)';             
b = b/sqrt((b')*b);
c = b;
Jzf1 = ones(n, 1);
txt = Jzf1;
J1 = ones(n,1);
zfdd = zeros(1,p);
zfa0 = zeros(n, p*p);
zfa1 = zeros(n, p);
zfb0 = zeros(n, p);
zfb1 = ones(n, 1);
zfb2 = ones(n, 1);
nref = 10;
niter = 10;
zfh0 = mean(std(xx))/n^(1/(p+4));	
for zfia = 1:n;
   	zfxy = xx - repmat(xx(zfia,:), n,1);
      zfx = zfxy.*zfxy*ones(p,1);
      zfK = exp(- zfx /(2*zfh0*zfh0*p) );
      K = zfK/sum(zfK);
      Kv = repmat(K,1,p).*zfxy;
      zfa0(zfia,:) = reshape((Kv')*zfxy, 1, p*p);
      zfa1(zfia,:) = ones(1,n,1)*Kv;
      zfb0(zfia,:) = (y')*Kv;
      zfb1(zfia) = (K')*y;
   end;
   
	zfbb0 = [];
	zfcc0 = zeros(p, nd);
	for ik = 1 : min(nd+1,p);
      if ik > 1; 
         b = (eye(p) - zfbb0*(zfbb0'))*ones(p,1);
   	end;
   	for iter = 1 : niter;
        b0 = b;
        zfbb = [zfbb0 b];
        zfAs = 0.0;
        zfBs = 0.0;
     	for zfia = 1 : n;
        zfaa = reshape(zfa0(zfia,:), p, p);
        zfat12 = zfa1(zfia,:)*zfbb ;
        zfat22 = (zfbb')*zfaa*zfbb;
        zfat1 = [1
           (zfat12')];
        zfat2 = [zfa1(zfia,:)*zfbb 
           (zfbb')*zfaa*zfbb ];
        zfat = [zfat1 zfat2]  + (1.0E-5)*eye(1+ik);
        zfbt = [zfb1(zfia)   zfb0(zfia,:) * zfbb ]';
        zfbt = inv(zfat)*zfbt;
        a = zfbt(1);

        zfbe1 = 0;
        if ik >1 ;
        		zfbe1 = zfbt(2:ik);
        end;
        d = zfbt(1+ik);

		  if ik == 1 
        		zfbb0 = zeros(p,1);
        end
        zfBs = zfBs + d*( (zfb0(zfia,:)- zfa1(zfia,:)*a)' - zfaa*zfbb0*zfbe1  );
        zfAs = zfAs + d*d*zfaa;
        if ik == 1 
           zfbb0=[];
        end
     	end;

      zfAs = zfAs + 1.0e-10*eye(p);
      b = inv(zfAs)*zfBs;
      if ik > 1 
      		b = (eye(p) - zfbb0*(zfbb0'))*b;
      end;
      b = b/sqrt((b')*b);

      if  ( (b0')*(eye(p) - b*(b'))*b0 < 1.0E-5 )
                iter = niter + 1;
      end;
   	end;
		if ik == 1 
         zfbb0 = b;
      end
      if ik > 1 
         zfbb0 = [zfbb0 b];
      end
   end;
   error1 = 0.0;
   error2 = 0.0;
   for kk = 1 : nref;
        	zfbbi = zfbb0(:,1:nd);
        	xxk = xx*zfbb0;
        	[n1,nd1] = size(xxk);  
 		  	zfh = h;
        	for zfia = 1 : n;
         	zfxyk = xxk - repmat(xxk(zfia,:), n,1);
         	zfxy = xx - repmat(xx(zfia,:), n, 1);
         	zfx = zfxyk.*zfxyk*ones(nd1,1);
         	zfK = exp(- zfx /(2*zfh*zfh) );
            K = zfK/sum(zfK);
         	Kv = repmat(K,1,p).*zfxy;
         	zfa0(zfia,:) = reshape((Kv')*zfxy, 1, p*p);
         	zfa1(zfia,:) = ones(1,n)*Kv;
         	zfb0(zfia,:) = (y')*Kv;
         	zfb1(zfia) = (K')*y;
        	end;
    	  	zfbb0 = zfbb0(:,1:nd);
    	  	for ik = 1 : nd;
    	  		b = zfbb0(:,ik);
    			for iter = 1 : niter;
        			b0 = b;
        			zfAs = 0.0;
        			zfBs = 0.0;
        			for zfia = 1 : n;
        				zfaa = reshape(zfa0(zfia,:), p, p);
        				zfat12 = zfa1(zfia,:)*zfbb0 ;
        				zfat22 = (zfbb0')*zfaa*zfbb0;
        
        				zfat1 = [1 zfat12]';
        				zfat2 = [zfa1(zfia,:)*zfbb0 
           				(zfbb0')*zfaa*zfbb0] ;
        				zfat = [zfat1 zfat2]  + (1.0E-5)*eye(1+nd);
        				zfbt = [zfb1(zfia)   zfb0(zfia,:)*zfbb0]';

        				zfbt = inv(zfat)*zfbt;
        				a = zfbt(1);
        				zfbe1 = zfbt(2:(nd+1));
        				d = zfbe1(ik);
                        DD(zfia, ik) = d;
        				zfbe1(ik) = 0.0;
        				zfBs = zfBs + d*( (zfb0(zfia,:)- zfa1(zfia,:)*a)' - zfaa*zfbb0*zfbe1  );
        				zfAs = zfAs + d*d*zfaa;
        			end;
        			zfAs = zfAs + 1.0e-10*eye(p);
        			b = inv(zfAs)*zfBs;
        			zfbb0(:,ik) = zeros(p,1);
        			b = (eye(p) - zfbb0*(zfbb0'))*b;
        			b = b/sqrt((b')*b);
        			zfbb0(:,ik) = b;
        		end;
       end;
       error1 =  sum( diag( zfbb0'*(eye(p)-zfbbi*zfbbi')*zfbb0 ) );   
       if error2 + error1 < 0.0001 
          break;
       end;
       error2 = error1;
   end;
   B = zfbb0;

xB = xx*B;
ye = y;
for i = 1:n;
    xi = xB - repmat(xB(i,:), n,1);
    kx = exp(-sum(xi.^2,2)/(2*h*h));
    kx(i) = 0;
    xi1 = [ones(n,1) xi];
    kxi1 = xi1.*repmat(kx, 1, nd+1);
    beta = inv( kxi1'*xi1 )*kxi1'*y;
    ye(i) = beta(1);
end
cv = (y-ye)'*(y-ye)/n;

DD = DD-repmat(mean(DD), n, 1);
[v, d] = eig(DD'*DD);
d = diag(d);
[d, I] = sort(d);
v = v(:,I);
B = B*v;

B = ss*B;
for i = 1:size(B,2);
    B(:,i) = B(:,i)/sqrt(B(:,i)'*B(:,i));
end

