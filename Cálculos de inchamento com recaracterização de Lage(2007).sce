// Artigo referência Al Ajmi 2011 F1
//Usando PC-SAFT e recaracterização de Lage (2007)
//m, sigma e epsilon_K dos puros foram tirados de gross e sadowski 2001
//Parâmetros da PC-SAFT dos pseudos baseados no Método 4

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
C7_1 componente 11
C7_2 componente 12
C7_3 componente 13
C7_4 componente 14
C7_5 componente 15
*/

clear
clc

exec("C:\Users\jvms9\OneDrive\Área de Trabalho\Finalização de dissertação\Comparação Al-Ajmi CBTERMO versão corrigida\Etapa 2 - alfas melhores\F1\Lage, 5 pseudos\Função PC-SAFT para obter Z.sce")

//Frações molares do óleo inicial
x0 = [0.223;0.293;21.657;6.758;7.024;1.325;4.229;1.817;2.68;4.146;14.653409628;16.245857733;8.967827490;3.318413271;6.659491879]
x0 = x0/100

//Frações molares do gás de injeção
xinj = [0.5966;0;0.103;0.0769;0.095;0;0.0679;0;0.0606;0;0;0;0;0;0]
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
        k = zeros(15,15) //Parâmetro de interação binária (kij)
        sigma = [2.7852;3.3130;3.7039;3.5206;3.6184;3.7574;3.7086;3.8296;3.7729;3.2272484;3.3442025;3.3885222;3.3912680;3.4258112;3.5825506]//diâmetro de segmentos (independente da temperatura) (Angstrons)
        epsilon_K = [169.21;90.96;150.03;191.42;208.11;216.53;222.88;230.75;231.20;168.95768;186.22237;189.09430;189.49993;195.66731;223.43914]//profundidade do potencial/k (K)
        m = [2.0729;1.2053;1;1.6069;2.002;2.2616;2.3316;2.5620;2.6896;4.1325876;4.6223710;7.3561058;12.001585;16.298212;26.274769]//número de segmentos

        d = sigma.*(1-0.12*exp(-3*epsilon_K/T)) //angstrons

        soma = 0
        for i = 1:length(x)
            soma = soma + x(i)*m(i)*d(i)^3 //angstrons^3
        end

        etachute = 0.5//Estimativa inicial de eta. De acordo com Gross e Sadowski 2001, eta = 0.5 é boa estimativa inicial para achar V de líquido. Já eta = 1e-10 é boa estimativa inicial para achar V de vapor.

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
    SF(i) = (Vmolar(i,1)*(1+injetado(i,1)/orig))/0.000211853089 //O denominador é o volume molar do óleo puro no primeiro estágio (0% de gás de injeção)
end
disp(SF)
SFexp = [1.074148128;1.223921551;1.375508038;1.600641842;1.674419504;1.821659419]//Fator de inchamento calculado pelo Método 4 na Etapa 1
pe = abs((SF-SFexp)./SFexp)*100
disp(pe)
disp(Vmolar(:,1))
disp(Vmolar(:,2))
