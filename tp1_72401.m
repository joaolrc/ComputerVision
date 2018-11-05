function NumMec = tp1_72401
%% Analisador de Código de Barras
%   Retorna txt com NumMec,NumSec,NumImg,TotCB,TotQR,R0,R90,R180,R270,Q0,Q90,Q180,Q270,BadCB,BadQR,TotDigCB,CBL,CBR,CBG,StringCB
clear all; close all;
NumMec=72401;

listaF=dir('../svpi2017_TP1_img_*_*.png'); % ler lista de imagens

txt=fopen('tp1_72401.txt','wt'); % Abre um txt e escreve 

for v=1:max(size(listaF)) %Ciclo para correr a imagens da lista

A=im2double(imread(strcat('../',listaF(v).name)));

        NumImg=listaF(v).name(22:23); % numero da imagem
        NumSec=listaF(v).name(18:20); % numero da sequencia


X=false(size(A));  %criar matriz de zeros lógicos (podia ser zeros mas ocupava mais espaço)
minSize=1600; %limite minimo de pixels agrupados para guardar aresta
[L,N]=bwlabel(A);

for K=1:N;
    C=(L==K);
    if (nnz(C)>minSize);  %arestas grandes
        X=X | C;
    end
end


XX=imcomplement(X);    %inverter pixels
XXX=imfill(XX,'holes');    %pintar molduras de branco
SE2=[1 1 1 1 1 1 1
    1 0 0 0 0 0 1
    1 0 1 1 1 0 1
    1 0 1 1 1 0 1
    1 0 1 1 1 0 1
    1 0 0 0 0 0 1
    1 1 1 1 1 1 1];  %matriz representativa dos quadrados dos QR codes (1=preto na matriz SE2)
SE1=~SE2;



%Matrizes para leitura das Barras
%Delimitadores de inicio e fim para várias escalas de CB
startbit1= logical([0 0 1 0 1 1 0 1 1 1 0]);
stopbit1= logical([0 1 1 1 0 0 0 1 0 1 0 0]);
startbit2= logical([0 0 0 0 1 1 0 0 1 1 1 1 0 0 1 1 1 1 1 1 0 0]);
stopbit2= logical([0 0 1 1 1 1 1 1 0 0 0 0 0 0 1 1 0 0 1 1 0 0 0 0]);
startbit3=imresize(startbit1,3);
startbit3=startbit3(1,:);
stopbit3=imresize(stopbit1,3);
stopbit3=stopbit3(1,:);
%digitos codificação L
CODL={};
digito0L= [1 1 1 0 0 1 0];
digito1L= [1 1 0 0 1 1 0];
digito2L= [1 1 0 1 1 0 0];
digito3L= [1 0 0 0 0 1 0];
digito4L= [1 0 1 1 1 0 0];
digito5L= [1 0 0 1 1 1 0];
digito6L= [1 0 1 0 0 0 0];
digito7L= [1 0 0 0 1 0 0];
digito8L= [1 0 0 1 0 0 0];
digito9L=[1 1 1 0 1 0 0];
CODL{1}= digito0L;
CODL{2}= digito1L;
CODL{3}= digito2L;
CODL{4}= digito3L;
CODL{5}= digito4L;
CODL{6}= digito5L;
CODL{7}= digito6L;
CODL{8}= digito7L;
CODL{9}= digito8L;
CODL{10}= digito9L;

%digitos codificação R
CODR={};
digito0R= [0 0 0 1 1 0 1];
digito1R= [0 0 1 1 0 0 1];
digito2R= [0 0 1 0 0 1 1];
digito3R= [0 1 1 1 1 0 1];
digito4R= [0 1 0 0 0 1 1];
digito5R= [0 1 1 0 0 0 1];
digito6R= [0 1 0 1 1 1 1];
digito7R= [0 1 1 1 0 1 1];
digito8R= [0 1 1 0 1 1 1];
digito9R= [0 0 0 1 0 1 1];
CODR{1}= digito0R;
CODR{2}= digito1R;
CODR{3}= digito2R;
CODR{4}= digito3R;
CODR{5}= digito4R;
CODR{6}= digito5R;
CODR{7}= digito6R;
CODR{8}= digito7R;
CODR{9}= digito8R;
CODR{10}= digito9R;

%digitos cod
CODG={};
digito0G= [1 0 1 1 0 0 0];
digito1G= [1 0 0 1 1 0 0];
digito2G= [1 1 0 0 1 0 0];
digito3G= [1 0 1 1 1 1 0];
digito4G= [1 1 0 0 0 1 0];
digito5G= [1 0 0 0 1 1 0];
digito6G= [1 1 1 1 0 1 0];
digito7G= [1 1 0 1 1 1 0];
digito8G= [1 1 1 0 1 1 0];
digito9G= [1 1 0 1 0 0 0];
CODG{1}= digito0G;
CODG{2}= digito1G;
CODG{3}= digito2G;
CODG{4}= digito3G;
CODG{5}= digito4G;
CODG{6}= digito5G;
CODG{7}= digito6G;
CODG{8}= digito7G;
CODG{9}= digito8G;
CODG{10}= digito9G;

for s=1:10
    CODR3{s}=imresize(logical(CODR{s}),3);
    CODR3{s}=CODR3{s}(1,:);
    CODL3{s}=imresize(logical(CODL{s}),3);
    CODL3{s}=CODL3{s}(1,:);
    CODG3{s}=imresize(logical(CODG{s}),3);
    CODG3{s}=CODG3{s}(1,:);
    
    CODR2{s}=imresize(logical(CODR{s}),2);
    CODR2{s}=CODR2{s}(1,:);
    CODL2{s}=imresize(logical(CODL{s}),2);
    CODL2{s}=CODL2{s}(1,:);
    CODG2{s}=imresize(logical(CODG{s}),2);
    CODG2{s}=CODG2{s}(1,:);
end

%Inicialização das variáveis
TotQR=0;
TotCB=0;
Q0=0;
Q90=0;
Q180=0;
Q270=0;
R0=0;
R90=0;
R180=0;
R270=0;
cbvalido=0;
cbinvalido=0;
qrvalido=0;
qrinvalido=0;
digitosL=0;
digitosR=0;
digitosG=0;
nada=0;
CBL=0;
CBR=0;
CBG=0;
TotDigCB=0;
j=0;
BadQR=0;

[J,K]=bwlabel(XXX);   %K buracos, matriz J
ce=cell(N,1);        % criar célula Nfilas 1 col

for k=1:K;
    [fila,col]= find (J==k);
    ce{k}=A(min(fila):max(fila),min(col):max(col));
    ce{k}=ce{k}(2:end-1,2:end-1);   %retirar a moldura(preta) de cada codigo
    level{k} = graythresh(ce{k});  %detetar nivel de tresholding por codigo
    ce{k}=im2bw(ce{k}, level{k});  %binarizar todas as imagens

    for l=1:7
        ce{k}=rot90(ce{k},l);    %rodar a imagem e cortar quiet zone
        for i=1:7
            if sum(~ce{k}(:,1))==0
                ce{k}=ce{k}(:,2:end);
            end
        end
    end
    
    
    %Deteção de quadrados nos QR codes
    
    if sum(ce{k}(1:7,1))==0 && sum(ce{k}(1,1:7))==0   %Detetar quadrados superiores esquerdos
        quadsupesq{k} = 1;
    else
        quadsupesq{k} = 0;
    end
    
    if sum(ce{k}(end-7:end,1))==0 && sum(ce{k}(end,1:7))==0 %Detetar quadrado inferior esquerdo
        quadinfesq{k} = 1;
    else
        quadinfesq{k} = 0;
    end
    
    if sum(ce{k}(1:7,end))==0 && sum(ce{k}(1,end-7:end))==0 %Detetar quadrado superior direito
        quadsupdir{k} = 1;
    else
        quadsupdir{k}=0;
    end
    
    if sum(ce{k}(end-7:end,end))==0 && sum(ce{k}(end,end-7:end))==0 %Detetar quadrado inferior direito
        quadinfdir{k} = 1;
    else
        quadinfdir{k}=0;
    end
    
    % Distinção entre QRcodes  e cod barras , orientação de QRcodes e Validez
    if quadsupesq{k} || quadsupdir{k} || quadinfesq{k} || quadinfdir{k} == 1 % incremento de QR codes
        TotQR=TotQR+1;
        
        
        if ~quadinfdir{k} ==1 && quadinfesq{k}==1   %QR com orientação normal
            Q0=Q0+1;
            QR{TotQR}=ce{k};
            
            
        elseif ~quadsupdir{k} ==1 && quadinfesq{k}==1 %QR com orientação +90 contrarelogio
            Q90=Q90+1;
            QR{TotQR}=ce{k};
            
            
        elseif ~quadsupesq{k} ==1 && quadinfesq{k}==1 %QR com orientação +180 contrarelogio
            Q180=Q180+1;
            QR{TotQR}=ce{k};
            
            
        elseif ~quadinfesq{k} ==1 && quadsupesq{k}==1 %QR com orientação +270 contrarelogio
            Q270=Q270+1;
            QR{TotQR}=ce{k};
            
        end
        %Redução proporcional dos QRcodes para poder confirmar validez
        for kk=1:8
            if nnz(bwhitmiss(imresize(QR{TotQR},1/kk),SE1,SE2))==3
                QR{TotQR}=imresize(QR{TotQR},1/kk);
                break
            end
        end
        %confirmação de validez dos QR através da regra lado=17+4N
        [a b]=size(QR{TotQR});
        for N=1:10
            if (a-17)/4== N
                qrvalido=qrvalido+1;
                break
            end
        end
        
    else
        TotCB=TotCB+1;  %contador de codigos de barras
        CB{TotCB}=ce{k};  %indexação de cod barras em matriz
        
        %Identificação da orientação de cod barras e rotação, se necessário
        if sum(CB{TotCB}(1:20,1))==0 && sum(CB{TotCB}(1:20,end))==0  %identificaçao e indexaçao de barras em posiçao normal
            R0=R0+1;
            CBnorm{R0}=CB{TotCB};
            CB{TotCB}= CBnorm{R0}(1,:);
            
            
        elseif sum(CB{TotCB}(1,1:20))==0 && sum(CB{TotCB}(end,1:20))==0  %barras pos vertical com numeros á direita(norm+90º)
            R90=R90+1;
            CB_90{R90}=CB{TotCB};
            CB_90{R90}=rot90(CB_90{R90},3); %rodar para a posição correta
            CB{TotCB}= CB_90{R90}(1,:);
            
            
        elseif sum(CB{TotCB}(end-20:end,1))==0 && sum(CB{TotCB}(end-20:end,end))==0  %barras pos horizontal com numeros em cima(norm+180º) (invertido)
            R180=R180+1;
            CB_180{R180}=CB{TotCB};
            CB_180{R180}=rot90(CB_180{R180},2);
            CB{TotCB}= CB_180{R180}(1,:);       %ficar com a primeira linha
            
        elseif  sum(CB{TotCB}(1,end-20:end))==0 && sum(CB{TotCB}(end,end-20:end))==0  %barras pos vertical com numeros á esquerda (norm+270º)
            R270=R270+1;
            CB_270{R270}=CB{TotCB};
            CB_270{R270}=rot90(CB_270{R270},1);
            CB{TotCB}= CB_270{R270}(1,:);    %ficar com a primeira linha

        end
    end
              
end


    %verificação da validez dos CB
    for i=1:TotCB
        for x=1:2:20  %x é o numero de digitos que existem na mensagem
        t=23+7*x;  %comprimento do vetor
        tt=46+14*x;
        ttt=69+21*x;
        if max(size(CB{i}))== t && nnz(~CB{i}(1:11))==5
            if CB{i}(1:11)== startbit1
                if CB{i}(1,end-11:end)== stopbit1    % end-11 pq já se conta com o pixel end e se fosse end-12 contavam-se 13 pixels
                    cbvalido=cbvalido+1;
                    CBvalido{cbvalido}=CB{i};   %cbvalido é só para confirmar validez de start  e stop
                    mensagem=CBvalido{cbvalido}(1,12:end-12);
                    kkl=x;
                    kkr=x;
                    kkg=x;
                    for g=1:x
                        for s=1:10
                            if mensagem( 7*g-6:7*g)==CODL{s}
                                kkl=kkl-1;
                                if kkl==x/2+0.5
                                    j=j+1;
                                    digitoscent(:,j)=s-1;      %retirar digitos centrais
                                end
                                if kkl==0
                                    CBL=CBL+1;
                                    TotDigCB=TotDigCB+x;                           %contar a totalidade de digitos reconhecidos
                                end
                                break;
                                
                            elseif mensagem( 7*g-6:7*g)==CODR{s}
                                kkr=kkr-1;
                                if kkr==x/2+0.5
                                    j=j+1;
                                    digitoscent(:,j)=s-1;
                                end
                                if kkr==0
                                    CBR=CBR+1;
                                    TotDigCB=TotDigCB+x;
                                end
                                break;
                                
                            elseif mensagem( 7*g-6:7*g)==CODG{s}
                                kkg=kkg-1;
                                if kkg==x/2+0.5
                                    j=j+1;
                                    digitoscent(:,j)=s-1;
                                end
                                if kkg==0
                                    CBG=CBG+1;
                                    TotDigCB=TotDigCB+x;
                                end
                                break;
                            end
                        end
                    end
                    break
               end
            end
            
        elseif sum(size(CB{i}))-1== tt && nnz(~CB{i}(1:11))==6
            if CB{i}(1:22)== startbit2
                if CB{i}(end-23:end)== stopbit2
                    cbvalido=cbvalido+1;
                    CBvalido{cbvalido}=CB{i};
                    mensagem=CBvalido{cbvalido}(1,23:end-24);
                    kkl=x;
                    kkr=x;
                    kkg=x;
                    for g=1:x
                        for s=1:10
                            if mensagem(14*g-13:14*g)==CODL2{s}(1,:)
                                kkl=kkl-1;
                                if kkl==x/2+0.5
                                    j=j+1;
                                    digitoscent(:,j)=s-1;
                                end
                                if kkl==0
                                    CBL=CBL+1;
                                    TotDigCB=TotDigCB+x;
                                end
                                break;
 
                            elseif mensagem(14*g-13:14*g)==CODR2{s}(1,:)
                                kkr=kkr-1;
                                if kkr==x/2+0.5
                                    j=j+1;
                                    digitoscent(:,j)=s-1;
                                end
                                if kkr==0
                                    CBR=CBR+1;
                                    TotDigCB=TotDigCB+x;
                                end
                                break;

                            elseif mensagem(14*g-13:14*g)==CODG2{s}(1,:)
                                kkg=kkg-1;
                                if kkg==x/2+0.5
                                    j=j+1;
                                    digitoscent(:,j)=s-1;
                                end
                                if kkg==0
                                    CBG=CBG+1;
                                    TotDigCB=TotDigCB+x;
                                end
                                break;
                            end
                        end
                    end                  
                end
            end
            
       elseif sum(size(CB{i}))-1== ttt && nnz(~CB{i}(1:11))==8
            if CB{i}(1:33)== startbit3
                if CB{i}(end-35:end)== stopbit3
                    cbvalido=cbvalido+1;
                    CBvalido{cbvalido}=CB{i};
                    mensagem=CBvalido{cbvalido}(1,34:end-36);  %retirar apenas informação
                    kkl=x;
                    kkr=x;
                    kkg=x;
                    for g=1:x
                        for s=1:10
                           
                            if mensagem( 21*g-20:21*g)==CODL3{s}(1,:)
                                kkl=kkl-1;
                                if kkl==x/2+0.5
                                    j=j+1;
                                    digitoscent(:,j)=s-1;
                                end
                                if kkl==0
                                    CBL=CBL+1;
                                    TotDigCB=TotDigCB+x;
                                end
                                break;
                                
                            elseif mensagem(21*g-20:21*g)==CODR3{s}(1,:)
                                kkr=kkr-1;
                                if kkr==x/2+0.5
                                    j=j+1;
                                    digitoscent(:,j)=s-1;
                                end
                                if kkr==0
                                    CBR=CBR+1;
                                    TotDigCB=TotDigCB+x;
                                end
                                break;
                               
                            elseif mensagem(21*g-20:7*21)==CODG3{s}(1,:)
                                kkg=kkg-1;
                                if kkg==x/2+0.5
                                    j=j+1;
                                    digitoscent(:,j)=s-1;
                                end
                                if kkg==0
                                    CBG=CBG+1;
                                    TotDigCB=TotDigCB+x;
                                end
                                break;
                            end
                        end
                    end
                end
            end
        end
        
    end




end
BadCB=TotCB-CBL-CBR-CBG;  %codigos de barras invalidos
ordenados=sort(digitoscent,'descend'); % organiza a o dig_central por ordem
StringCB=num2str(ordenados); % imprime os numeros na string
StringCB = regexprep(StringCB, ' ', ''); % remove o espacos entre cada numero na string
fprintf(txt,'%d,%s,%s,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%s \r\n',NumMec,NumSec,NumImg,TotCB,TotQR,R0,R90,R180,R270,Q0,Q90,Q180,Q270,BadCB,BadQR,TotDigCB,CBL,CBR,CBG,StringCB);
end

fclose(txt); % fecha o txt

%NOTAS

%cbvalido é só para confirmar validez de start  e stop



%end
