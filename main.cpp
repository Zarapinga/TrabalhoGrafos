#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <queue>
#include <limits>
#include <stack>
#include <set>
#include <algorithm>
#include <random>
#include <chrono>

using namespace std;

int n;
int quantidade;

vector<pair<double, double>>  extrairDados(const string& nomeArquivo) {
  ifstream arquivo(nomeArquivo);
  string linha;
  vector<pair<double, double>> coordenadas;

  if (arquivo.is_open()) {
    // Pula a primeira linha
    getline(arquivo, linha);

    // Lê a segunda linha e extrai o número
    getline(arquivo, linha);
    stringstream ss;
    for (char c : linha) {
      if (isdigit(c)) {
        ss << c;
      }
    }
    int numero = stoi(ss.str());

    // Lê a terceira linha e extrai o texto após ":"
    getline(arquivo, linha);
    size_t pos = linha.rfind(':');
    string texto_extraido;
    if (pos != string::npos) {
      texto_extraido = linha.substr(pos + 1);
      texto_extraido.erase(0, texto_extraido.find_first_not_of(" "));
      texto_extraido.erase(texto_extraido.find_last_not_of(" ") + 1);
    } else {
      cout << "Caractere ':' não encontrado na linha." << endl;
    }

    // Pula linhas até "NODE_COORD_SECTION"
    while (getline(arquivo, linha)) {
      if (linha.find("NODE_COORD_SECTION") != string::npos) {
        break;
      }
    }

    // Pula a linha "NODE_COORD_SECTION"

    // Vetores para armazenar os dados
    vector<int> inteiros;
    vector<float> floats1;
    vector<float> floats2;



    // Lê n linhas
    for (int i = 0; i < numero; ++i) {
      if (!getline(arquivo, linha)) {
        break; // Sai do loop se não houver mais linhas
      }

      stringstream ss(linha);
      int inteiro;
      float float1, float2;

      // Extrai os valores da linha
      ss >> inteiro >> float1 >> float2;

      // Armazena os valores nos vetores
      coordenadas.push_back(make_pair(float1, float2));
    }

    // Imprime os dados extraídos
    cout << "Quantidade de vertice: " << numero << endl;
    cout << "Tipo: " << texto_extraido << endl;

    arquivo.close();
    n = coordenadas.size();
    return coordenadas;
  } else {
    cerr << "Erro ao abrir o arquivo." << endl;
    return coordenadas;
  }
  return coordenadas;
}

// Função para calcular a distância euclidiana entre dois pontos
double calcular_distancia(pair<double, double>& p1, pair<double, double>& p2) {
    double dx = p1.first - p2.first;
    double dy = p1.second - p2.second;
    return sqrt(dx * dx + dy * dy);
}

// Função para encontrar a MST usando o algoritmo de Prim
vector<pair<int, int>> encontrar_mst(vector<pair<double, double>>& coordenadas) {
  int n = coordenadas.size();
  vector<bool> mstSet(n, false);
  vector<double> key(n, numeric_limits<double>::infinity());
  vector<int> parent(n, -1);

  priority_queue<pair<double, int>, vector<pair<double, int>>, greater<pair<double, int>>> pq;
  key[0] = 0;
  pq.push({0, 0});

  while (!pq.empty()) {
    int u = pq.top().second;
    pq.pop();

    mstSet[u] = true;

    for (int v = 0; v < n; ++v) {
      if (u != v && !mstSet[v]) {
        double distancia = calcular_distancia(coordenadas[u], coordenadas[v]);
        if (distancia < key[v]) {
          key[v] = distancia;
          parent[v] = u;
          pq.push({key[v], v});
        }
      }
    }
  }

  vector<pair<int, int>> mst;
  for (int i = 1; i < n; ++i) {
    mst.push_back({parent[i], i});
  }

  return mst;
}

// Função para encontrar um circuito Euleriano usando o algoritmo de Hierholzer
vector<int> encontrar_circuito_euleriano(vector<vector<int>>& grafo_euleriano, int inicio) {
  int n = grafo_euleriano.size();
  vector<int> circuito;
  stack<int> pilha;
  pilha.push(inicio);

  while (!pilha.empty()) {
    int u = pilha.top();

    if (!grafo_euleriano[u].empty()) {
      int v = grafo_euleriano[u].back();
      grafo_euleriano[u].pop_back();
      pilha.push(v);
    } else {
      circuito.push_back(u);
      pilha.pop();
    }
  }

  reverse(circuito.begin(), circuito.end()); // Inverte o circuito para a ordem correta
  return circuito;
}

// Função para remover repetições de vértices em um circuito Euleriano
vector<int> remover_repeticoes(const vector<int>& circuito_euleriano) {
  int n = circuito_euleriano.size();
  vector<int> ciclo_hamiltoniano;
  set<int> visitados;

  for (int u : circuito_euleriano) {
    if (visitados.find(u) == visitados.end()) {
      ciclo_hamiltoniano.push_back(u);
      visitados.insert(u);
    }
  }

  ciclo_hamiltoniano.push_back(ciclo_hamiltoniano[0]); // Fecha o ciclo, voltando ao ponto inicial
  return ciclo_hamiltoniano;
}

// Função para executar o algoritmo Twice-Around
vector<int> twice_around(vector<pair<double, double>>& coordenadas) {
  // 1. Encontrar a MST
  vector<pair<int, int>> mst = encontrar_mst(coordenadas);

  // 2. Construir o grafo Euleriano
  vector<vector<int>> grafo_euleriano(n);
  for (auto aresta : mst) {
    grafo_euleriano[aresta.first].push_back(aresta.second);
    grafo_euleriano[aresta.second].push_back(aresta.first);
  }

  // 3. Encontrar o circuito Euleriano
  vector<int> circuito_euleriano = encontrar_circuito_euleriano(grafo_euleriano, 0);

  // 4. Remover repetições
  vector<int> ciclo_hamiltoniano = remover_repeticoes(circuito_euleriano);

  return ciclo_hamiltoniano;
}

double calcular_distancia_maxima( vector<int>& ciclo_hamiltoniano, vector<pair<double, double>>& coordenadas) {
  double distancia_maxima = 0;
  int n = ciclo_hamiltoniano.size();
  for (int i = 0; i < n - 1; ++i) {
    int ponto1 = ciclo_hamiltoniano[i];
    int ponto2 = ciclo_hamiltoniano[i + 1];
    double distancia = calcular_distancia(coordenadas[ponto1], coordenadas[ponto2]);
    if (distancia > distancia_maxima) {
      distancia_maxima = distancia;
    }
  }
  return distancia_maxima;
}

vector<int> otimizar_2opt_randomizado( vector<int>& ciclo_inicial,  vector<pair<double, double>>& coordenadas) {
  vector<int> ciclo = ciclo_inicial;
  int n = ciclo.size();
  random_device rd;  // Obtem uma semente aleatória do hardware
  mt19937 gerador(rd()); // Cria um gerador Mersenne Twister com a semente
  uniform_int_distribution<> distribuicao(1, n - 2); // Define a distribuição uniforme entre 1 e n-2

  for (int iteracao = 0; iteracao < quantidade; ++iteracao) {
    // Gera posições aleatórias i e j
    int i = distribuicao(gerador);
    int j = distribuicao(gerador);
    if (i >= j - 1 || j >= n - 1) continue; // Ignora se i >= j-1 ou j >= n-1

    // Calcula o comprimento das arestas atuais
    double comprimento_atual = calcular_distancia(coordenadas[ciclo[i - 1]], coordenadas[ciclo[i]]) +
                               calcular_distancia(coordenadas[ciclo[j]], coordenadas[ciclo[j + 1]]);

    // Calcula o comprimento das novas arestas após a troca
    double comprimento_novo = calcular_distancia(coordenadas[ciclo[i - 1]], coordenadas[ciclo[j]]) +
                              calcular_distancia(coordenadas[ciclo[i]], coordenadas[ciclo[j + 1]]);

    // Se a troca resultar em um ciclo menor, realiza a troca
    if (comprimento_novo < comprimento_atual) {
      reverse(ciclo.begin() + i, ciclo.begin() + j + 1); // Inverte a parte do ciclo entre i e j
    }
  }

  return ciclo;
}

int main(int argc, char* argv[]) {

  if (argc < 4) {
    cerr << "Erro: Argumentos insuficientes." << endl;
    cerr << "Uso: " << argv[0] << " <nome_da_entrada> <num_iteracoes> <nome_da_saida>" << endl;
    return 1;
  }
  string nomeArquivo = argv[1];
  quantidade = stoi(argv[2]);
  string saida = argv[3];

  auto inicio = chrono::high_resolution_clock::now();
  // Testando o algoritmo
  vector<pair<double, double>> coordenadas = extrairDados(nomeArquivo);
  vector<int> rota_inicial = twice_around(coordenadas);

  cout << "Rota inicial: ";
  for (int ponto : rota_inicial) {
    cout << ponto << " ";
  }
  cout << endl << endl;

  double comprimento_inicial = calcular_distancia_maxima(rota_inicial, coordenadas);
  cout << "Distância máxima inicial: " << comprimento_inicial << endl << endl;

  vector<int> rota_otimizada = otimizar_2opt_randomizado(rota_inicial, coordenadas);

  cout << "Rota otimizada: ";
  for (int ponto : rota_otimizada) {
    cout << ponto << " ";
  }
  cout << endl << endl;;

  double comprimento_otimizado = calcular_distancia_maxima(rota_otimizada, coordenadas);
  cout << "Distância máxima otimizada: " << comprimento_otimizado << endl;

  auto fim = chrono::high_resolution_clock::now();
  auto duracao = chrono::duration_cast<chrono::milliseconds>(fim - inicio);
  // Imprime o tempo de execução
  cout << "Tempo de execução: " << duracao.count() << " milissegundos" << endl;

  ofstream arquivo_saida(saida);
  if (arquivo_saida.is_open()) {
    // Escreve a rota no arquivo
    arquivo_saida << "Rota otimizada: ";
    for (int ponto : rota_otimizada) {
      arquivo_saida << ponto << " ";
    }
    arquivo_saida << endl;

    // Escreve a distância máxima no arquivo
    arquivo_saida << "Distância máxima entre dois pontos: " << comprimento_otimizado << endl;

    arquivo_saida.close(); // Fecha o arquivo
  }

  return 0;
}