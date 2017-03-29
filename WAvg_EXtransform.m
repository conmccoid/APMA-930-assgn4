%% Function that transforms a scaler field over eta and xi to a scalar field over x and y
% The input is a matrix, whose elements are understood to be ordered
% according to their location in the eta and xi. Then it returns a (of
% equal size) whose elements are the values of the scalar field over the
% transformed variables. (Theta ranges from 0 to pi, and r ranges from (1
% to R, where R needs to be provided). 


function    FinalPsi =  EXtransform(psiEX,R) 

     % Get the origianl size of the matrix
     [NN,MM] = size(psiEX) ;
     
     % Pick off the left and right pats
     RightPsi = psiEX(:,MM/2+1:MM) ;
     LeftPsi = psiEX(:,1:MM/2) ;
     
     % Determine the size of the two matrices
     [NR,MR] = size(RightPsi) ;
     [NL,ML] = size(LeftPsi) ;
     
     % Divide the matrix in two
     
     
     % Build an equally sized matrix for psiXY   
     RightpsiXY = zeros(NR,MR) ;
     LeftpsiXY = zeros(NL,ML) ;
     
     
 %%   
     RightPsi = flipud(RightPsi) ;
     % Loop through the matrix RightpsiXY, and fill the entries with the closest
     % corresponding entry in RightPsi
      for ii = 1:NR
         
         for jj = 1:MR
             
             % Determine the corresponding x,y values
             [x,y] = index2XY(jj,ii,NR,MR,R,'r') ;
             
           
             % Determine the corresponding eta and xi values
             [eta,xi] = XY2EX(x,y)  ;
             
             
             % Determine the corresponding n and m inidices in psiEX
             [n,m] = EX2index(eta,xi,NR,MR,R,'r') ;
             
             
             % Check for permissable indices (make sure you're not inside
             % the small circle or outside of the large one)
             
              if 1 < n & n < NR & 1 < m & m < MR
                 RightpsiXY(ii,jj) = WeiAvg(eta,xi,NR,MR,R,'r',RightPsi) ;
                 %RightpsiXY(ii,jj) = RightPsi(n,m) ;
             elseif n == 1 & 0 < m & m < MR+1 | m == 1 & 0 < n & n < NR + 1
                 RightpsiXY(ii,jj) = RightPsi(n,m) ;
             elseif n == NR & 0 < m & m < MR + 1 | m == MR & 0 < n & n < NR + 1
                 RightpsiXY(ii,jj) = RightPsi(n,m) ;
             else
                 RightpsiXY(ii,jj) = NaN ;
             
                      
             end
             
             
         end
         
       end
      
%%       
       
       % Loop through the matrix LeftpsiXY and fill in with the closest
       % value of LeftPsi
          for ii = 1:NL
         
          for jj = 1:ML
             
             % Determine the corresponding x,y values
             [x,y] = index2XY(jj,ii,NL,ML,R,'l') ;
             
           
             % Determine the corresponding eta and xi values
             [eta,xi] = XY2EX(x,y)  ;
             
             
             % Determine the corresponding n and m inidices in psiEX
             [n,m] = EX2index(eta,xi,NL,ML,R,'l') ;
             
             
             % Check for permissable indices (make sure you're not inside
             % the small circle or outside of the large one)
             
             if 1 < n & n < NL & 1 < m & m < ML
                 LeftpsiXY(ii,jj) = WeiAvg(eta,xi,NL,ML,R,'l',LeftPsi) ;
                 %LeftpsiXY(ii,jj) = LeftPsi(n,m) ;
             elseif n == 1 & 0 < m & m < ML+1 | m == 1 & 0 < n & n < NL + 1
                 LeftpsiXY(ii,jj) = LeftPsi(n,m) ;
             elseif n == NL & 0 <m & m < ML + 1 | m == ML & 0 < n & n < NL + 1
                 LeftpsiXY(ii,jj) = LeftPsi(n,m) ;
             else
                 LeftpsiXY(ii,jj) = NaN ;
             
                      
             end
             
             
         end
         
     end
       
       
       
       
       

% Produce a plot (this can be commented out if the plot should be produced
% outside of the function).
xx = linspace(-R,R,MR) ; 
yy = linspace(R,0,NR) ;

FinalPsi = [RightpsiXY ; LeftpsiXY] ;
FinalPsi = flipud(FinalPsi) ;

%contourf(xx,yy,RightpsiXY,60)
    
end

%%


% Maps indices of psiXY to actual x,y values. Takes on value (-R,0) in
% bottom left corner of the matrix. 
function  [x,y] = index2XY(j,i,N,M,R,side)

     % Since the domain is basically a semicircle we will take the x y
     % domain to be size 2R in x and size R in y
     
     % Determine the grid spacing
     dx = 2*R/M ;
     dy = R/N ;

     
     % Check if we are dealing with the left or right matrix in the eta xi
     % plane.
     if side == 'r'
     % Map the indices.
     x = (j-1)*dx - R ;
     y = (N-i)*dy ;
     
     elseif side == 'l'
     
     x = (j-1)*dx - R ;
     y = (1-i)*dy ;
     
     end
     
end 

%%

% Maps x and y values to eta and xi values
function  [eta,xi] = XY2EX(x,y)    
     
     % Map the numerical values
%      eta = atan(y./x)*pi ;
%      if eta < 0
%          eta = eta + pi^2 ;
%      end
       eta = atan2(y,x)*pi ;
       xi  = 1/(2*pi)*log(x.^2 + y.^2) ;
end


%%

% Maps the values of eta, xi to their respective values in the psiEX
% matrix.
function [n,m] = EX2index(eta,xi,N,M,R,side)
               
     if side == 'r'
     % Map the values to the indices
     n = ceil(xi*N/(log(R)/pi + eps)) ;
     m = ceil(eta*M/(pi^2) + eps) ;
     
     elseif side == 'l'
     m = ceil( M*(pi^2 + eta)/pi^2 + eps) ;
     n = ceil( N*(log(R) - xi*pi)/log(R) + eps) ;
     end
     
end


%% Function to produce one of the four closest grid points to a given pair eta, xi

function [n,m] = closeEX2index(eta,xi,N,M,R,side,location)

% Check which side of the matrix we're working with
     if side == 'r'
       
     % Check which location is needd
         if location == 'nw'    

         n = ceil(xi*(N)/(log(R)/pi + eps)) ;
         m = floor(eta*(M)/(pi^2) + eps) ;

         elseif location == 'ne'

         n = ceil(xi*(N)/(log(R)/pi + eps)) ;
         m = ceil(eta*(M)/(pi^2) + eps) ;  

         elseif location == 'sw'

         n = floor(xi*(N)/(log(R)/pi + eps)) ;
         m = floor(eta*(M)/(pi^2) + eps) ;

         elseif location == 'se'

         n = floor(xi*(N)/(log(R)/pi + eps)) ;
         m = ceil(eta*(M)/(pi^2) + eps) ;  

         end
         
     
     elseif side == 'l'
         
         
          if location == 'nw'    

         m = floor( M*(pi^2 + eta)/pi^2 + eps) ;
         n = ceil( N*(log(R) - xi*pi)/log(R) + eps) ; 

         elseif location == 'ne'

         m = ceil( M*(pi^2 + eta)/pi^2 + eps) ;
         n = ceil( N*(log(R) - xi*pi)/log(R) + eps) ; 

         elseif location == 'sw'

         m = ceil( M*(pi^2 + eta)/pi^2 + eps) ;
         n = floor( N*(log(R) - xi*pi)/log(R) + eps) ;

         elseif location == 'se'

         m = floor( M*(pi^2 + eta)/pi^2 + eps) ;
         n = floor( N*(log(R) - xi*pi)/log(R) + eps) ;  

         end
         
         
     
     
     end

end




%%  Function to find the distance of 



function      Avg = WeiAvg(eta,xi,N,M,R,side,Matrix)

        % Find the distance between eta and the closest left grid point.
        dEta = pi^2/(M) ;
          dXi = log(R)/pi/(N) ;
         
          
       if side == 'r'
           
           % find the distances        
           [n,m] = closeEX2index(eta,xi,N,M,R,'r','sw') ;
          dL = abs(eta -  dEta*m) ;
          dH = abs(xi -  dXi*n) ;
  %         dL = dEta/2;
  %         dH = dXi/2 ;
         
          
          % Pick out the four values of the matrix that we're interested in
          
          % North west value
          [n,m] = closeEX2index(eta,xi,N,M,R,side,'nw');
          ValueNW = Matrix(n,m) ;
          
          % North east value
          [n,m] = closeEX2index(eta,xi,N,M,R,side,'ne');
          ValueNE = Matrix(n,m) ;
          
          % South west value
          [n,m] = closeEX2index(eta,xi,N,M,R,side,'sw');
          ValueSW = Matrix(n,m) ;
          
          % South east value
          [n,m] = closeEX2index(eta,xi,N,M,R,side,'se');
          ValueSE = Matrix(n,m) ;
          
          % Find the averages
          % Find the top average
          AvgTop = ((dEta-dL)*ValueNW + dL*ValueNE)/dEta;
          
          % Find the bottom average
          AvgBot = ((dEta-dL)*ValueSW + dL*ValueSE)/dEta;
          
          % Find the total average
          Avg = (AvgTop*(dH) + AvgBot*(dXi-dH))/dXi ;
          
       elseif side == 'l'
           
           % Find the distances                      
            [n,m] = closeEX2index(eta,xi,N,M,R,'l','se') ;
            dL = abs(eta +  dEta*(M-m-1)) ;     
            dH = abs(xi - dXi*(N-n-1/2)) ;
                 
 %         dL = dEta/2;
 %          dH = dXi/2 ;
 
          % Pick out the four values of the matrix that we're interested in
          
          % North west value
          [n,m] = closeEX2index(eta,xi,N,M,R,side,'nw');
          ValueNW = Matrix(n,m) ;
          
          % North east value
          [n,m] = closeEX2index(eta,xi,N,M,R,side,'ne');
          ValueNE = Matrix(n,m) ;
          
          % South west value
          [n,m] = closeEX2index(eta,xi,N,M,R,side,'sw');
          ValueSW = Matrix(n,m) ;
          
          % South east value
          [n,m] = closeEX2index(eta,xi,N,M,R,side,'se');
          ValueSE = Matrix(n,m) ;
          
          % Find the averages
          % Find the top average
          AvgTop = ((dEta-dL)*ValueNW + dL*ValueNE)/dEta;
          
          % Find the bottom average
          AvgBot = ((dEta-dL)*ValueSW + dL*ValueSE)/dEta;
          
          % Find the total average
          Avg = (AvgTop*(dH) + AvgBot*(dXi-dH))/dXi ;
          
          
       end
       
end
          
          
    



    


  