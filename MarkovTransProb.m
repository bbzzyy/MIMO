function Probability = MarkovTransProb( d ,scenario, density )
%UNTITLED2
%input
% d = distance
% scenario = 'highway','urban'
% density = 'low','medium' and 'high'
%output
% from\to LOS NLOSb NLOSv
% LOS
% NLOSb
% NLOSv

if strcmp(scenario,'highway')
    if strcmp(density,'low')
        dT = 70;
        ProbL2L = 6.7e-7*d^2-4.8e-4*d+0.99;
        ProbL2Nb = 4e-9*d^2-2.7e-6*d+0.018;
        ProbL2Nv = 1-ProbL2L-ProbL2Nb;
        ProbNb2L = exp(-(log(d)-5.2782)^2/1.8424)/(0.0289*d);
        ProbNb2Nb = 1-exp(-(log(d)-5.2782)^2/1.8424)/(0.0289*d);
        ProbNb2Nv = 1-ProbNb2L-ProbNb2Nb;
        if d>=0 && d<dT
            ProbNv2L = -9.8e-6*d^2+8.9e-4*d+0.97;
            ProbNv2Nb = 9.8e-6*d^2-8.9e-4*d+0.03;
        else
            ProbNv2L = -2e-6*d^2+1.6e-3*d+0.051;
            ProbNv2Nb = -1.4e-7*d^2+9.1e-5*d-0.0016;
        end
        ProbNv2Nv = 1-ProbNv2L-ProbNv2Nb;
    elseif strcmp(density,'medium')
        dT = 90;
        ProbL2L = 1.6e-6*d^2-1.2e-3*d+1;
        ProbL2Nb = -8.4e-8*d^2+3.5e-5*d+0.016;
        ProbL2Nv = 1-ProbL2L-ProbL2Nb;
        ProbNb2L = exp(-(log(d)-5.021)^2/1.5875)/(0.0346*d);
        ProbNb2Nb = 0.9132-exp(-(log(d)-4.7076)^2/0.7480)/(0.0484*d);
        ProbNb2Nv = 1-ProbNb2L-ProbNb2Nb;
        if d>=0 && d<dT
            ProbNv2L = -4.8e-5*d^2-5.62e-3*d+1.11;
            ProbNv2Nb = 4.4e-6*d^2-8.335e-4*d+0.042;
        else
            ProbNv2L = -2.286e-6*d^2+1.443e-3*d+0.1022;
            ProbNv2Nb = -2.7e-7*d^2+1.5e-4*d-0.0031;
        end
        ProbNv2Nv = 1-ProbNv2L-ProbNv2Nb;
    else
        dT = 90;
        ProbL2L = 2.1e-6*d^2-1.5e-3*d+1;
        ProbL2Nb = -1.1e-7*d^2+4.3e-5*d+0.015;
        ProbL2Nv = 1-ProbL2L-ProbL2Nb;
        ProbNb2L = exp(-(log(d)-4.927)^2/1.4876)/(0.0411*d);
        ProbNb2Nb = 0.9264-exp(-(log(d)-4.7012)^2/0.8186)/(0.056*d);
        ProbNb2Nv = 1-ProbNb2L-ProbNb2Nb;
        if d>=0 && d<dT
            ProbNv2L = -6.51e-5*d^2-1.04e-3*d+0.8706;
            ProbNv2Nb = 1.254e-7*d^2-3.775e-5*d+9.853e-3;
        else
            ProbNv2L = -1.412e-6*d^2+6.196e-4*d+0.2216;
            ProbNv2Nb = -1.4e-7*d^2+8.3e-5*d-0.0065;
        end
        ProbNv2Nv = 1-ProbNv2L-ProbNv2Nb;
    end
else
    if strcmp(density,'low')
        aL2L = 1.6e-6;
        bL2L = -1.2e-3;
        cL2L = 0.99;
        
        aL2Nb = -8.7e-7;
        bL2Nb = 6.7e-4;
        cL2Nb = -0.012;
        
        aNb2L = 1.6e-6;
        bNb2L = -1.1e-3;
        cNb2L = 0.2;
        
        aNb2Nb = -1.2e-6;
        bNb2Nb = 9.1e-4;
        cNb2Nb = 0.83;
        
        aNv2L = -1.4e-6;
        bNv2L = 6.7e-4;
        cNv2L = 0.079;
        
        aNv2Nb = -3e-7;
        bNv2Nb = 2.7e-4;
        cNv2Nb = -0.0059;
        
    elseif strcmp(density,'medium')
        aL2L = 1.5e-6;
        bL2L = -1.2e-3;
        cL2L = 0.93;
        
        aL2Nb = -5.9e-7;
        bL2Nb = 5.4e-4;
        cL2Nb = 0.0069;
        
        aNb2L = 1e-6;
        bNb2L = -7.1e-4;
        cNb2L = 0.12;
        
        aNb2Nb = -1.1e-6;
        bNb2Nb = 7.8e-4;
        cNb2Nb = 0.86;
        
        aNv2L = 8.1e-8;
        bNv2L = -2.1e-4;
        cNv2L = 0.14;
        
        aNv2Nb = -4.9e-7;
        bNv2Nb = 3.6e-4;
        cNv2Nb = -0.0046;
    else
        aL2L = 2.1e-7;
        bL2L = -6.5e-4;
        cL2L = 0.86;
        
        aL2Nb = -9e-8;
        bL2Nb = 3e-4;
        cL2Nb = 0.025;
        
        aNb2L = 7.7e-7;
        bNb2L = -5.3e-4;
        cNb2L = 0.083;
        
        aNb2Nb = -9e-7;
        bNb2Nb = 6.4e-4;
        cNb2Nb = 0.89;
        
        aNv2L = 6.8e-7;
        bNv2L = -5.7e-4;
        cNv2L = 0.14;
        
        aNv2Nb = -4e-7;
        bNv2Nb = 2.7e-4;
        cNv2Nb = 0.0058;
    end
    ProbL2L = min(1,max(0, aL2L*d^2+bL2L*d+cL2L ));
    ProbL2Nb = min(1,max(0, aL2Nb*d^2+bL2Nb*d+cL2Nb ));
    ProbL2Nv = 1-ProbL2L-ProbL2Nb;
    ProbNb2L = min(1,max(0,aNb2L*d^2+bNb2L*d+cNb2L));
    ProbNb2Nb = min(1,max(0,aNb2Nb*d^2+bNb2Nb*d+cNb2Nb));
    ProbNb2Nv = 1-ProbNb2L-ProbNb2Nb;
    ProbNv2L = min(1,max(0,aNv2L*d^2+bNv2L*d+cNv2L));
    ProbNv2Nb = min(1,max(0,aNv2Nb*d^2+bNv2Nb*d+cNv2Nb));
    ProbNv2Nv = 1-ProbNv2Nb-ProbNv2L;
end

Probability = [ProbL2L,ProbL2Nb,ProbL2Nv;
    ProbNb2L,ProbNb2Nb,ProbNb2Nv;
    ProbNv2L,ProbNv2Nb,ProbNv2Nv];

end