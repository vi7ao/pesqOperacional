# guerreiros da br 277
# Vitor Teles Moreira Passos
# Renan Augusto Mattos Nutse
# Andre Negreli Conrado Bini

import cplex
import sys

# para rodar digite no terminal o seguinte comando:
# python3 arenales.py "nome-do-arquivo"

def getInfoFromCplexFile():
    filename = sys.argv[-1]
    problem = cplex.Cplex(filename, "LP")
    nomeVariaveis = problem.variables.get_names()
    numeroVariaveis = problem.variables.get_num()
    nomeRestricoes = problem.linear_constraints.get_names()
    vetorB = problem.linear_constraints.get_rhs()
    sentidoRestricoes = problem.linear_constraints.get_senses()
    coeficientesRestricoes = problem.linear_constraints.get_rows()
    coeficientesVariaveis = problem.variables.get_cols()
    
    print(nomeVariaveis)
    print(numeroVariaveis)
    print(nomeRestricoes)
    print(vetorB)
    print(sentidoRestricoes)
    print(coeficientesRestricoes)
    print(coeficientesVariaveis)
    
# essas tres funções são pra calcular matriz inversa utilizando o método gauss-jordan.
def gerarMatrizIdentidade(n):
    identidade = [[0 for x in range(n)] for y in range(n)]
    for i in range (n):
        identidade[i][i] = 1
    return identidade # ok

def gerarMatrizExtendida(matrizOriginal, identidade):
    matrizExtendida = []
    for i in range(len(matrizOriginal)):
        matrizExtendida.append(matrizOriginal[i] + identidade[i])
    print(matrizExtendida)
    return matrizExtendida # ok

def inverterMatriz(matrizOriginal): #gauss jordan
    nLinhas = len(matrizOriginal)
    identidade = gerarMatrizIdentidade(nLinhas)
    matrizExtendida = gerarMatrizExtendida(matrizOriginal, identidade)
    for i in range (nLinhas):
        maiorIndice = i #posicao
        for j in range(i, nLinhas):
            if abs(matrizExtendida[j][i]) > abs(matrizExtendida[maiorIndice][i]):
                maiorIndice = j
                matrizExtendida[i], matrizExtendida[maiorIndice] = matrizExtendida[maiorIndice], matrizExtendida[i]
        pivo = matrizExtendida[i][i]
        for j in range((nLinhas * 2)):
            matrizExtendida[i][j] = matrizExtendida[i][j] / pivo
        for j in range(nLinhas):
            if j != i:
                multiplicador = matrizExtendida[j][i]
                for k in range((nLinhas * 2)):
                    matrizExtendida[j][k] = matrizExtendida[j][k] - multiplicador * matrizExtendida[i][k]
    matrizInvertida = [linha[nLinhas:] for linha in matrizExtendida]
    return matrizInvertida # ok

def main():
    getInfoFromCplexFile()
    
    
main()