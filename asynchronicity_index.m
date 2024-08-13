function [ASI] = asynchronicity_index(r_m,k_m)

% r_m : long term mean monthly value of p [12 by 1]
% k_m : long-term mean montly value of PET  [12 by 1]

%%%  estimation of probability mass function  %%%

% p_m : probability of R_m 
% q_m : probability of k_m

for i= 1:12
    p_m(i,1)=r_m(i,1)/sum(r_m) ;
end

for i=1:12
    q_m(i,1)=k_m(i,1)/sum(k_m) ;
end


%%% estimation of JS    ####

n_m=0.5.*(p_m+q_m) ;

for i=1:12
    xx_p(i,1)=p_m(i,1)*log2(p_m(i,1)/n_m(i,1)) ;
    xx_q(i,1)=q_m(i,1)*log2(q_m(i,1)/n_m(i,1)) ;
end

D_p_n= sum(xx_p) ;
D_q_n=sum(xx_q) ;

JS_obs=(0.5*D_p_n+0.5*D_q_n)^0.5 ;

%%% minimum_js_estimate
for i=1:12
    for j=1:12
     yy=i+j-1 ;
     if yy> 12
         xx=yy-12 ;
     else
         xx=yy ;
     end
        q_m_w(i,j)=q_m(xx,1) ;
        clear yy xx 
    end
end

for i=1:12
    for j=1:12
    n_m_w(i,j)=0.5.*(p_m(i,1)+q_m_w(i,j)) ;
    end
end

for i=1:12
    for j=1:12 
     xx_p_w(i,j)=p_m(i,1)*log2(p_m(i,1)/n_m_w(i,j)) ;
    xx_q_w(i,j)=q_m_w(i,j)*log2(q_m_w(i,j)/n_m_w(i,j)) ;
    end
end

for i=1:12
D_p_n_w(1,i)= sum(xx_p_w(:,i)) ;
D_q_n_w(1,i)=sum(xx_q_w(:,i)) ;
end 

for i=1:12
JS_w(1,i)=(0.5*D_p_n_w(1,i)+0.5*D_q_n_w(1,i))^0.5 ;
end 

Js_min=min(JS_w) ;

ASI=(JS_obs-Js_min) ;

end