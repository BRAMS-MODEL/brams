# Como usar o Modelo BRAMS para gerar arquivos CSV com as séries temporais.



## Opção de escrita de CSV

Há duas opções no RAMSIN para escrever as saídas CSV:

| Item | IPOS | Ação                                                   |
| ---- | ---- | ------------------------------------------------------ |
| 1    | 10   | Gera as saídas CSV concomitante com as saídas em GRADS |
| 2    | 11   | Gera apenas as saídas CSV                              |

## Tempo entre escritas no CSV

O tempo entre escritas define a série temporal. Para ajustá-lo basta alterar FRQANL para o valor requerido. Por exemplo, FRQANL = 1.0. Veja como fica um exemplo do trecho do RAMSIN_BASIC:

```bash
! History/analysis file output

IPOS = 11, ! 0-no files, 2-grads files, 3 - NetCDF
IOUTPUT = 0, ! 0-no files, 1-save ASCII, 2-save binary
HFILOUT = './2022071400/dataout/HIS/APG', ! History file prefix
AFILOUT = './2022071400/dataout/ANL/APG', ! Analysis file prefix
FRQHIS = 2160000., ! History file frequency
FRQANL = 1., ! Analysis file frequency

```

> Atente-se que a FRQANL deve ser múltiplo do DTLONG usado na rodada. Logo, para saídas a cada 1 seg (caso acima), o DTLONG deve ser de 1 seg. 

## Escolhas das estações

Para escolher as estações (ou apenas uma) edite o arquivo **estacoes.csv** e proceda a colocação da lista de estações com seu nome e coordenadas geográficas.  A primeira linha é um cabeçalho. As separações são com ";". Ex:

```csv
Nome do ponto;Latitude;Longitude
DSB_E;-16.79597;-49.23139
Lugar_em_GO;-17.0;-50.0
```

É necessário que as estações definidas estejam dentro do domínio do modelo.

## Escolha das variáveis a serem plotadas

As variáveis a serem escritas (tanto no CSV quanto no grads) devem estar definidas no RAMSIN_BASIC, na seção POST.  Exemplo:

```bash
$POST
NVP = 8, !#Number of POST variables - Open variables.csv to see the availables
VP = 'maguv',
'diruv',
'tempc',
'rh',
'press',
'rshort',
'rlong',
'precip',

GPREFIX = './2022071400/dataout/POST/APG5KM',

```

Estando as variáveis definidas em POST, aquelas escolhidas para serem escritas no arquivo CSV devem ser listadas em um arquivo chamado **estvars.dat**. Basta apenas colocá-las listadas no arquivo, linha a linha. Exemplo abaixo.

> Atenção: Nesse caso as variáveis devem estar escritas em **caixa alta**. 

```bash
MAGUV
DIRUV
TEMPC
PRECIP
RH
RSHORT
PRESS
```

Observe que no arquivo a lista não contém a variável RLONG que aparece no RAMSIN. O RAMSIN pode conter mais variáveis que as escolhidas em estvars.dat.

## Escolhendo os níveis a serem colocados no arquivo CSV

Os níveis verticais são escolhidos da mesma maneira que em rodadas normais do modelo. Esses níveis são definidos na seção MODEL_GRIDS do RAMSIN_BASIC e podem tanto ser níveis definidos pelo usuário como níveis calculados a partir de definições de NNZP, DZRAT e DELTAZ. Se DELTAZ for zero os níveis serão obtidos da lista determinada por ZZ. Veja exemplo abaixo:

```bash
NNZP     = 51, ! Number of z gridpoints
...   
DELTAZ   = 0.0,   ! Z grid spacing (set to 0. to use ZZ)
DZRAT    = 1.09,   ! Vertical grid stretch ratio
DZMAX = 750., ! Maximum delta Z for vertical stretch
FIXLEVELS = 0,

! Vertical levels if DELTAZ = 0
ZZ = 0.000,
50.00,
70.00,
90.00,
110.0,
130.0,
165.5,
232.05,
...

```

Observe que são 51 níveis no exemplo e são mostrados apenas alguns que foram ajustados. 

Em geral no arquivo de saída espera-se mostrar apenas alguns níveis de interesse e não todos. Para isso basta alterar o RAMSIN e fazer IPRESSLEV = 3 e escolher quantos níveis se deseja em INPLEVS. Por exemplo, para 8 níveis:

```bash
IPRESSLEV = 3,
INPLEVS = 8,

```

## Gerando saídas Grads em conjunto com arquivos CSV

Para gerar as saídas em GRADS deve-se usar IPOS=10 conforme mostrado nesse documento. As variáveis que serão plotadas são as mostradas na seção POST (todas elas). É similar a uma rodada com IPOS=2. 

Quando IPOS=10 o tempo entre saídas em grads **não é determinado** por FRQANL. O modelo busca o valor de **meteogramFreq** na seção METEOGRAM do RAMSIN_ADVANCED. Por exemplo, para gerar saídas a cada hora o RAMSIN_ADVANCED fica assim:

```bash
$METEOGRAM
applyMeteogram = .false.,
meteogramFreq = 3600.,
...
```

As saídas são escritas com os mesmos níveis escolhidos para os arquivos CSV e são colocadas nos arquivos apontados por GPREFIX.



## Como são as saídas geradas



As saídas produzidas são arquivos CSV para cada estação determinada no arquivo **estacoes.csv** contendo INPLEVS níveis verticais (ou todos se IPRESSLEV=0). Cada arquivo  trata de uma variável (sinal) específica. 

> Nota: todas as saídas são produzidas no final da rodada do modelo.

Cada arquivo é composto pela justaposição no nome da estação, nome da variável e a data de início da simulação (com a hora). Esses campos são separados por um sublinhado. Um nome exemplo é: DSB_E_TEMPC_2022071400.csv onde,

- DSB_E - nome da estação

- TEMPC - Nome do sinal (variável)

- 2022071400 - Data de início da simulação (condições iniciais) - 00 é a hora!

O arquivo contém 3 linhas de cabeçalho.

As duas primeira linhas do cabeçalho informam o nome da estação, sua posição geográfica, o símbolo da variável (sinal), a unidade física e uma descrição da variável.

A terceira linha é um cabeçalho de colunas. Contém a indicação do tempo, seguida do valor dos níveis em metros. Para o caso de variáveis de superfície devem ser desconsiderados os níveis maiores que zero.

O tempo em cada linha de dados é mostrado no formato AAAA-MM-DDTHH:MM:SS Exemplo: 2022-07-14T00:00:01

Esse é um exemplo do arquivo de saída:

```csv
# Local: DSB_E, Lat=-16.79597, Lon= -49.23139
# Sinal: PRESS ,Unidade: 'mb' ,Desc: 'pressure'
              Tempo,   0.0000,  47.7675,  66.8745,  85.9815, 105.0885, 124.1955, 158.1105, 221.6891,
2022-07-14T00:00:01, 920.0346, 912.0775, 909.0271, 907.1441, 905.1518, 903.2356, 900.6769, 895.6725,
2022-07-14T00:00:02, 920.0748, 912.1176, 909.0668, 907.1837, 905.1918, 903.2763, 900.7194, 895.7213,
2022-07-14T00:00:03, 920.0935, 912.1360, 909.0853, 907.2023, 905.2104, 903.2949, 900.7382, 895.7404,
2022-07-14T00:00:04, 920.0930, 912.1351, 909.0841, 907.2010, 905.2090, 903.2933, 900.7360, 895.7379,
2022-07-14T00:00:05, 920.0654, 912.1077, 909.0569, 907.1741, 905.1820, 903.2668, 900.7097, 895.7119,
2022-07-14T00:00:06, 920.0212, 912.0646, 909.0142, 907.1315, 905.1399, 903.2249, 900.6686, 895.6723,
2022-07-14T00:00:07, 919.9823, 912.0263, 908.9763, 907.0936, 905.1024, 903.1876, 900.6317, 895.6366,
....


```




