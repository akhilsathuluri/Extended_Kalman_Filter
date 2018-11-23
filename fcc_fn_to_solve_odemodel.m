function dydt = fcc_fn_to_solve_odemodel(t,y,u)  
       
        fcc_parameters
        
        C_rc = y(1);   O_d = y(2);  T_rg = y(3);  % controlled variables
        dFa  = u(1); dFsc  = u(2);                % inputs  

        F_a   = F_a0 + dFa;      F_sc = F_sc0 + dFsc; 
        T_oil = T_oil0 ;         kc1 = kc10 ;

        T_ri0 = ((c_po*F_oil + lamda*F_oil*c_pd)*T_oil  +  c_ps*F_rc*T_rg)/(c_po*F_oil  +  lamda*F_oil*c_pd  +  c_ps*F_rc);
        gama  = delH_f * F_oil / (T_ri0*(F_rc*c_ps + F_oil*c_po  + lamda*F_oil*c_pd));
        phi0  = 1-(m*C_rc);
        K0    = k0*exp(-E_f/(R*T_ri0));

        T_r   =  T_ri0*(1-((gama*y_f0*K0*phi0*(1-exp(-alpa*tc*COR*z_r)))/(alpa + (K0*phi0*(1-exp(-alpa*tc*COR*z_r))))));
        Kr    =  k0*exp(-E_f/(R*T_r));
        y_f1  =  y_f0*alpa/(alpa + (Kr*phi0*(1-exp(-alpa*tc*COR))));
        T_ri1 =  T_ri0*(1-((gama*y_f0*Kr*phi0*(1-exp(-alpa*tc*COR)))/(alpa + (Kr*phi0*(1-exp(-alpa*tc*COR))))));
        y_g1  = (1+R_r)*F_gi*((y_f1^I_gi) - y_f1)/(1-I_gi);

        sig   =  1.1 + sig2*(T_rg - 873);
        delH  =  -h1 - (h2*(T_rg-960)) + 0.6*(T_rg-960)^2 ; 
        k     =  k_com * exp(E_cb*((1/960)-(1/T_rg))/R);

        C_cat =  kc1 * sqrt(tc*exp(-E_cf/(R*T_ri1))/(C_rc^N));
        C_sc  =  C_rc + C_cat; 

        dydt(1,:) = (F_sc*(C_sc-C_rc)/W) - k*O_d*C_rc;
        dydt(2,:) = (Ra*(O_in - O_d)/Wa) - (((n+2+((n+4)*sig))*k*O_d*C_rc*W)/(4*M_c*Wa*(1+sig)));
        dydt(3,:) = (T_ri1*F_sc/W) + (T_a*F_a*c_pa/(W*c_ps)) - (T_rg*(F_sc*c_ps + F_a*c_pa)/(W*c_ps)) - (delH*k*O_d*C_rc/(c_ps*M_c));

end


