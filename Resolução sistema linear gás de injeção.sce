//O programa a seguir é usado para gerar pseudocomponentes para o gás de injeção de hidrocarbonetos de forma que a massa molar de seus pseudos corresponda à massa molar de determinados pseudos do óleo.

clear //Limpa a memória.
clc(); //Limpa a tela.

//Massas molares do c6 ao último pseudocomponente do fluido de reservatório selecionado para gerar sistema bem condicionado (C6 e pseudo 1 neste caso):
m = [84;118.003665809646]

//Fração molar e massa molar do c6 a c11 do gás de injeção:
x = [0.00322;0.00127;0.00036;0.00006;0.00002;0.00002]
M = [84;96;107;121;134;147]

//Determinando os momentos necessários da distribuição discreta:
mu = zeros(length(M),1) //Inicialmente o vetor mu é um vetor coluna de zeros com o mesmo n° de linhas que o vetor m.
for k = 1:length(M)
    for i = 1:length(M)
        mu(k) = mu(k) + x(i)*M(i)^(k-1)
    end
end

//Montando o sistema linear para calcular y
A = zeros(length(m), length(m))
b = zeros(length(m),1)
for k = 1:length(m)
    b(k) = - mu(k)
    for i = 1:length(m)
        A(k,i) = m(i)^(k-1)
    end
end
c = cond(A)
[y,kerA]=linsolve(A,b)  // resolve Ax+b = 0
disp("y:",y)
disp("A*y+b:",A*y+b)
disp("kerA",kerA)
disp("inv(A)",inv(A))
disp("cond(A)",c)
