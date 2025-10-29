// Trabalho de termodinâmica: artigo referência Al Ajmi 2011 F1
//Usando PC-SAFT
//m, sigma e epsilon_K dos puros foram tirados de gross e sadowski 2001
//Parâmetros da PC-SAFT dos pseudos baseados em Assareh

/*
CO2 componente 1
N2 componente 2
C1 componente 3
C2 componente 4
C3 componente 5
iC4 componente 6
nC4 componente 7
iC5 componente 8
nC5 componente 9
C6 componente 10
C7 componente 11
C8 componente 12
C9 componente 13
...
C36+ componente 40
*/

clear
clc

exec("C:\Users\jvms9\OneDrive\Área de Trabalho\Finalização de dissertação\Comparação Al-Ajmi CBTERMO versão corrigida\F1\Assareh\Função PC-SAFT para obter Z.sce")

//Frações molares do óleo inicial
x0 = [0.223;0.293;21.657;6.758;7.024;1.325;4.229;1.817;2.68;4.146;4.114;4.181;3.698;3.412;2.954;2.535;2.332;1.958;1.854;1.661;1.431;1.306;1.259;1.126;0.997;0.94;0.845;0.762;0.727;0.625;0.587;0.583;0.542;0.528;0.541;0.458;0.44;0.399;0.391;6.659]
x0 = x0/100

//Frações molares do gás de injeção
xinj = [0.5966;0;0.103;0.0769;0.095;0;0.0679;0;0.0606;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0]
injetado = [[25;75;125;200;225;275],[0;0;0;0;0;0]]//número de mols de gás injetado considerando 100 mols de óleo inicial
for h = 1:2
    for j = 1:6

        select j

        case 1
            //Primeiro ponto do swelling test
            P_bar = 110.5 //bar
            inj = injetado(1,h) //mols de gás injetado
            n_inj = inj*xinj //número de mols injetado por componente
            orig = 100 //número de mols original de óleo
            n_orig = orig*x0 //número de mols por componente óleo original
            n = n_orig+n_inj
            x = n./sum(n)

        case 2
            //Segundo ponto do swelling test
            P_bar = 143.7 //bar
            inj = injetado(2,h) //percentual injetado
            n_inj = inj*xinj //número de mols injetado por componente
            orig = 100 //número de mols original de óleo
            n_orig = orig*x0 //número de mols por componente óleo original
            n = n_orig+n_inj
            x = n./sum(n)

        case 3
            //Terceiro ponto do swelling test
            P_bar = 165.9 //bar
            inj = injetado(3,h) //percentual injetado
            n_inj = inj*xinj //número de mols injetado por componente
            orig = 100 //número de mols original de óleo
            n_orig = orig*x0 //número de mols por componente óleo original
            n = n_orig+n_inj
            x = n./sum(n)

        case 4
            //Quarto ponto do swelling test
            P_bar = 198.1 //bar
            inj = injetado(4,h) //percentual injetado
            n_inj = inj*xinj //número de mols injetado por componente
            orig = 100 //número de mols original de óleo
            n_orig = orig*x0 //número de mols por componente óleo original
            n = n_orig+n_inj
            x = n./sum(n)

        case 5
            //Quinto ponto do swelling test
            P_bar = 209.6 //bar
            inj = injetado(5,h) //percentual injetado
            n_inj = inj*xinj //número de mols injetado por componente
            orig = 100 //número de mols original de óleo
            n_orig = orig*x0 //número de mols por componente óleo original
            n = n_orig+n_inj
            x = n./sum(n)

        case 6
            //Sexto ponto do swelling test
            P_bar = 230.3 //bar
            inj = injetado(6,h) //percentual injetado
            n_inj = inj*xinj //número de mols injetado por componente
            orig = 100 //número de mols original de óleo
            n_orig = orig*x0 //número de mols por componente óleo original
            n = n_orig+n_inj
            x = n./sum(n)
        end

        //Temperatura do swelling test
        T_kelvin = 350.4 //K

        T = T_kelvin//K
        P = P_bar * 100000 //Pa
        R = 8.3144621 // Constante universal dos gases (J/(K*mol))

        K = 1.380649e-23 //m^2 kg s^-2 K^-1
        Na = 6.02214076e23 //mol^-1
        R = K*Na
        //planilha = readxls("C:\Users\jvms9\OneDrive\Área de Trabalho\Trabalhos usando correlação generalizada Riazi\Memon\Função gama, Hosseinifar, Assareh\5 pseudos\kij.xls")
        k = zeros(40,40)//= planilha(1) //Parâmetro de interação binária (kij)
        sigma = [2.7852;3.3130;3.7039;3.5206;3.6184;3.7574;3.7086;3.8296;3.7729;3.7937989;3.9358132;3.9880347;4.0167024;4.0316259;4.0416412;4.0509862;4.0592475;4.0680810;4.0775961;4.0863772;4.0967946;4.1049628;4.1133217;4.1233270;4.1335240;4.1439165;4.1545088;4.1651222;4.1755737;4.1863509;4.1972567;4.2078389;4.2177967;4.2290728;4.2404228;4.2512210;4.2620821;4.2721463;4.2830981;4.6271245]//diâmetro de segmentos (independente da temperatura) (Angstrons)
        epsilon_K = [169.21;90.96;150.03;191.42;208.11;216.53;222.88;230.75;231.20;243.55931;264.24299;271.45809;275.07996;276.08610;276.88324;277.82169;279.24174;280.68098;282.09987;283.16246;284.67554;285.72453;286.91240;288.22335;289.68311;291.28000;293.00486;294.53574;295.83643;297.54681;299.34599;300.87122;302.42211;304.10375;306.20418;307.97752;309.80719;311.25981;313.18851;388.95403]//profundidade do potencial/k (K)
        m = [2.0729;1.2053;1;1.6069;2.002;2.2616;2.3316;2.5620;2.6896;2.9362318;2.8888517;3.0753914;3.3559757;3.6945202;3.9871761;4.2985321;4.5756685;4.8941399;5.2300369;5.5221727;5.8448185;6.0844265;6.3178125;6.5874651;6.8493503;7.1034682;7.3498188;7.5966798;7.8458322;8.0798831;8.3072320;8.5387755;8.7464640;8.9856850;9.1880073;9.3972651;9.6011304;9.8147940;10.008590;12.904687]//número de segmentos

        d = sigma.*(1-0.12*exp(-3*epsilon_K/T)) //angstrons

        soma = 0
        for i = 1:length(x)
            soma = soma + x(i)*m(i)*d(i)^3 //angstrons^3
        end

        etachute = 0.5//Chute inicial de eta. De acordo com Gross e Sadowski 2001, eta = 0.5 é bom chute inicial para achar V de líquido. Já eta = 1e-10 é bom chute inicial para achar V de vapor.

        function z = f(eta)
            Z = SAFTZ(eta,T,x,m,sigma,epsilon_K,k)
            Vmolar = Na*%pi*soma/(6*eta*1e10^3)
            z = P-Z*R*T/Vmolar
        endfunction

        [etaf,v,info] = fsolve(etachute,f)
        ro = 6/%pi*etaf*soma^-1 //angstrons^-3
        Vmolar(j,h) = Na/(ro*1e-10^-3) //m3 mol-1
    end
end
disp(Vmolar)
[linh,col] = size(Vmolar)
for i = 1:linh
    SF(i) = (Vmolar(i,1)*(1+injetado(i,1)/orig))/0.000206995214
end
disp(SF)
SFexp = [1.0761;1.2314;1.3916;1.6296;1.7082;1.8656]
pe = abs((SF-SFexp)./SFexp)*100
disp(pe)
disp(Vmolar(:,1))
disp(Vmolar(:,2))
