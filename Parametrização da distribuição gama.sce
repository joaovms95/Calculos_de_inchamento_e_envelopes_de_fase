//Aplicando Whitson, Anderson e Søreide (1990) ao caso F1 de Al Ajmi
clear
clc

function F = obj(param)
    try
        xexp = [0.04114;0.04181;0.03698;0.03412;0.02954;0.02535;0.02332;0.01958;0.01854;0.01661;0.01431;0.01306;0.01259;0.01126;0.00997;0.00940;0.00845;0.00762;0.00727;0.00625;0.00587;0.00583;0.00542;0.00528;0.00541;0.00458;0.00440;0.00399;0.00391;0.06659]
        Mexp7mais = 291
        Mexp = [95;107;121;136;149;163;176;191;207;221;237;249;261;275;289;303;317;331;345;359;373;387;400;415;429;443;457;471;485;841.21475]
        alpha = param(1)
        betha = param(2)
        eta = param(3)
        function f1 = p(x)
            if x==eta & alpha<1  then
                f1 = 0
            else
                f1 = (x-eta)^(alpha-1)*exp(-(x-eta)/betha)/(betha^alpha*gamma(alpha))
            end
        endfunction

        //Definindo matriz massas molares limite inferior e superior
        Mn(1,1) = eta
        Mn(length(Mexp),2) = 1000 //Mn será uma matriz em que a primeira coluna será limite inferior de massa de cada SCN, e segunda coluna será limite superior.
        for i = 1:length(Mexp)-1
            Mn(i,2) = (Mexp(i+1)+Mexp(i))/2
        end
        for i = 2:length(Mexp)
            Mn(i,1) = Mn(i-1,2)
        end
        Mnsup = Mn(:,2) //limite superior
        Mninf = Mn(:,1) //limite inferior

        //Integrando com função do Scilab
        function f2 = P(x)
            f2 = intg(eta,x,p)
        endfunction

        function f3 = xp(x)
            f3 = x*p(x)
        endfunction

        //Propriedade média entre limite inferior a e limite superior b
        function f4 = Pmed(a,b)
            f4 = 1/(P(b)-P(a))*(intg(eta,b,xp)-intg(eta,a,xp))
        endfunction

        /*
        //Massa molar média de uma fração SCN
        for i = 1:length(Mnsup)
            M(i) = Pmed(Mninf(i),Mnsup(i))
        end

        //Fração molar de cada fração SCN
        for i = 1:length(Mnsup)
            Z(i) = P(Mnsup(i))-P(Mninf(i))
        end
        //Normalizando
        for i = 1:length(Z)
            z(i) = Z(i)/sum(Z)
        end
        */
        wexp = xexp.*Mexp./(sum(xexp.*Mexp))
        wexp($) = []
        function [f5,M,Z] = g(x)
            for i = 1:length(x)
                Z(i) = P(x(i))-P(Mninf(i))
            end
            Z = Z./sum(Z)
            for i = 1:length(x)
                M(i) = Pmed(Mninf(i),x(i))
            end
            w = Z.*M/(eta+alpha*betha)
            f5 = w-wexp
        endfunction
        Mnsup($) = []
        [Mnsup_final,res,info] = fsolve(Mnsup,g)
        [residuo,M,Z] = g(Mnsup_final)
        ERRO = 0
        for i = 1:length(M)
            ERRO = ERRO + ((Mexp(i)-M(i))/Mexp(i))^2
        end
        F = ERRO
        // Exibe progresso
        disp([alpha,betha,eta,ERRO])
    catch
        //Se algum erro ocorrer:
        F = 1e10
        disp("Erro encontrado em alpha:" + string(alpha))
        disp("Erro encontrado em beta:" + string(betha))
        disp("Erro encontrado em eta:" + string(eta))
    end
endfunction
alpha0 = 1
betha0 = 200
eta0 = 91
chute = [alpha0;betha0;eta0]
[otim,fval] = fminsearch(obj,chute)
