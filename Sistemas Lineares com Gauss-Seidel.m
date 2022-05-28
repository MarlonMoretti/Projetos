% Faculdade Anhanguera - Engenharia Elétrica
% Aluno: Marlon Moretti
% Disciplina: Sistema Elétrico de Potência 2
% Professor: Fábio Costa
%
%
% Funcionalidade: Resolução de Sistemas de Equações Lineares pelo método de Gauss-Seidel
% Ferramenta: Matlab e Octave Online
%
%
% Problema: Encontrar a Matriz O na equação [B]*[O]=[P]
% Onde:
%   B é o barramento de um sistema de rede de distribuição de energia, matriz de susceptância
%   O (teta) é o vetor ângulo das tensões
%   P é o vetor do fluxo linearizado de potência líquida em pU


% Limpando a memória e o software
clear all; clc; close all;


% Inicializando e declarando as matrizes
B = [((1/0.1)+(1/0.2)) (-1/0.1) 0 0 0;
    (-1/0.1) ((1/0.1)+(1/0.3)+(1/0.15)) (-1/0.15) 0 (-1/0.3);
    0 (-1/0.15) (1/0.15) 0 0;
    0 0 0 (1/0.5) (-1/0.5);
    0 (-1/0.3) 0 (-1/0.5) ((1/0.5)+(1/0.3)+(1/0.8))];

P = [-1;
    -1.5;
    -2;
    -0.5;
    -0.3];

Oi = zeros(length(P),1); % Ângulo do início da iteração (Ok)
    
Of = zeros(length(P),1); % Ângulo de resultado da iteração (Ok+1)

% Tolerancia é a margem percentual de erro tolerado devido às iterações
tolerancia = 0.001;

% Loop para realizar iterações e contador de iterações
loop = true;
cont = 0;

% Loop de iterações
while(loop)
    
    % Variáveis para indexar as matrizes
    i = 1;
    j = 1;
    
    % Atualização do valor inicial a cada iteração
    for i = 1:length(Of)
        Oi(i) = Of(i);
        x(i,1) = Oi(i); % Matriz x de backup para não perder as informações de Oi
    end
    
    
    % Cálculo para encontrar a matriz Of
    for i = 1:length(P)
        somatorio = 0;
        for j = 1:length(P)
            if(i==j)
                Oi(j) = 0;
            end
            somatorio = B(i,j)*Oi(j) + somatorio;
        end
        Oi = x;
        Of(i) = (1/B(i,i))*(P(i) - somatorio);
    end
    
    
    % Cálculo para encontrar o erro acumulado nas iterações
    Oerro = zeros(length(Oi),1);
    for i = 1:length(Of)
        Oerro(i) = Of(i) - Oi(i);
        Oerro(i) = abs(Oerro(i));
    end
    erro = max(Oerro)/max(abs(Of));
    
    cont++;
    
    % Erro dentro da tolerância?
    if((erro == tolerancia) || (erro < tolerancia))
        loop = false;
    end
end



% Conversão Radianos para Graus
for i = 1:length(Of)
    Of_graus(i,1) = 180*Of(i)/pi;
end

% Resultados
disp(['===== RESOLUÇÃO DE SISTEMAS LINEARES =====';'=====      MÉTODO GAUSS-SEIDEL       =====';' ';'----- FLUXO DE POTÊNCIA LINEARIZADO  -----';' ';'Por: Marlon Moretti - Engenheiro Eletricista.';' ';' ';'Equação: B*O=P -> O=inv(B)*P';' ';' ';])
B
P
disp([' ';'=====> RESULTADO: ';' ';'Taxa de Erro: ',num2str(erro*100),'%'])
disp(['Quantidade de Iterações: ',num2str(cont)])
disp(['Matriz de Resultados para teta(O) em radianos:';' '])
disp(Of)
disp([' ';'Matriz de Resultados para teta(O) em graus°:';' '])
disp(Of_graus)
