# guerreiros da br 277
# Vitor Teles Moreira Passos
# Renan Augusto Mattos Nutse
# Andre Negreli Conrado Bini

import cplex
import sys #linha incluida para ler o arquivo no terminal

# para rodar digite no terminal o seguinte comando:
# python3 arenales.py "nome-do-arquivo"

def getInfoFromCplexFile():
    filename = sys.argv[-1]
    arquivo = cplex.Cplex(filename, "LP")
    precisaFaseUm = False
    variaveisArtificiais = 0
    numeroVariaveis = arquivo.variables.get_num()
    for i in range(numeroVariaveis):
        if arquivo.variables.get_lower_bounds(i) != 0:
            arquivo.linear_constraints.add(lin_expr=[cplex.SparsePair(ind=[arquivo.variables.get_names()[i]],
                val=[1])],
                senses=["G"],
                rhs=[arquivo.variables.get_lower_bounds(i)])
        if arquivo.variables.get_upper_bounds(i) != 1e+20:
            arquivo.linear_constraints.add(lin_expr=[cplex.SparsePair(ind=[arquivo.variables.get_names()[i]],
                val=[1])],
                senses=["L"],
                rhs=[arquivo.variables.get_upper_bounds(i)])
    nomeVariaveis = arquivo.variables.get_names()
    numeroRestricoes = arquivo.linear_constraints.get_num()
    coeficientes = arquivo.linear_constraints.get_rows()
    funcaoObj = arquivo.objective.get_linear()
    # verificacao se precisa fase 1
    if arquivo.objective.get_sense() == arquivo.objective.sense.maximize:
        sentidoOriginal = "maximize"
        for i in range (numeroVariaveis):
            funcaoObj[i] = funcaoObj[i] * -1
    else:
        sentidoOriginal = "minimize"
    matrizNaoBasica = [[0.0] * numeroVariaveis for i in range(numeroRestricoes)]
    for i in range(numeroRestricoes):
        for j, var in enumerate(coeficientes[i].ind):
            matrizNaoBasica[i][var] = coeficientes[i].val[j]
    constraint_senses = arquivo.linear_constraints.get_senses()
    vetorB = [0.0] * numeroRestricoes
    for i in range(numeroRestricoes):
        vetorB[i] = arquivo.linear_constraints.get_rhs(i)
        if vetorB[i] < 0:
            if constraint_senses[i] == "L":
                constraint_senses[i] = "G"
            elif constraint_senses[i] == "G":
                constraint_senses[i] = "L"
            vetorB[i] = vetorB[i] * -1
            for j in range(len(matrizNaoBasica[i])):
                matrizNaoBasica[i][j] *= -1          
        if arquivo.linear_constraints.get_senses(i) == "E" :
            precisaFaseUm = True
            variaveisArtificiais += 1
        if arquivo.linear_constraints.get_senses(i) == "G":
            precisaFaseUm = True
    cn = arquivo.objective.get_linear()
    cb = [0.0] * numeroRestricoes
    matrizBasica = []
    
    for i in range(numeroRestricoes):
        matrizBasica.append([0] * numeroRestricoes)
    for i in range(numeroRestricoes):
        if constraint_senses[i] == "L":
            matrizBasica[i][i] = 1
        elif constraint_senses[i] == "G":
            matrizBasica[i][i] = -1
        else:
            matrizBasica[i][i] = 0    
    return matrizNaoBasica, matrizBasica, vetorB, cn, cb, precisaFaseUm, nomeVariaveis, sentidoOriginal, variaveisArtificiais

def gerarMatrizIdentidade(n):
    identidade = [[0 for x in range(n)] for y in range(n)]
    for i in range (n):
        identidade[i][i] = 1
    return identidade # ok

def gerarMatrizExtendida(matrizOriginal, identidade):
    matrizExtendida = []
    for i in range(len(matrizOriginal)):
        matrizExtendida.append(matrizOriginal[i] + identidade[i])
    return matrizExtendida # ok
    
def inverterMatriz(matrizOriginal):
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

def multiplicaMatrizPorVetor(matriz, vetor):
    vetorResultante = [0.0] * len(matriz)
    for i in range(len(matriz)):
        for j in range(len(vetor)):
            vetorResultante[i] += matriz[i][j] * vetor[j]
    return vetorResultante

def calculoVetorMultiplicadorSimplex(vetorCustoB, matrizInvertida):
    vetorMultiplicador = [0.0] * len(vetorCustoB)
    for i in range(len(vetorCustoB)):
        for j in range(len(vetorCustoB)):
            vetorMultiplicador[i] += vetorCustoB[j] * matrizInvertida[j][i]
    return vetorMultiplicador

def calculoCustosRelativos(vetorMultiplicador, matrizNaoBasica, vetorCustoN):
    nVariaveis = len(vetorCustoN)
    nRestricoes = len(matrizNaoBasica)
    custosRelativos = [0.0] * nVariaveis
    for i in range(nVariaveis):
        col_i = [matrizNaoBasica[j][i] for j in range(nRestricoes)]
        custosRelativos[i] = vetorCustoN[i] - (vetorMultiplicador[i] * col_i[i])
    return custosRelativos
        
def indiceDoValorMinimo(vetor):
    indiceMinimo = vetor.index(min(vetor))
    return indiceMinimo

def checarSolucaoOtima(custosRelativos):
    for i in range(len(custosRelativos)):
        if custosRelativos[i] < 0:
            return False
    return True

def calculoSolucaoOtima(xb, custo_b):
    z = 0
    for i in range(len(xb)):
        z += xb[i] * custo_b[i]
    return z

def calculoDirecaoSimplex(matrizInvertida, matrizNaoBasica, indicePraEntrarBase):
    coluna = [row[indicePraEntrarBase] for row in matrizNaoBasica]
    direcao = multiplicaMatrizPorVetor(matrizInvertida, coluna)
    return direcao
    
def printSolution(z, xb, nomeVariaveis):
    print("z: " + str(z))
    for i in range(len(nomeVariaveis)):
        print(nomeVariaveis[i] + ": " + str(xb[i]))
    return

def calculoTheta(xb, y):
    theta = [0.0] * len(xb)
    for i in range(len(xb)):
        if y[i] > 0:
            theta[i] = xb[i] / y[i]
    return theta

def indiceVariavelSair(theta):
    indice = -1
    thetaMinimo = float("inf")
    for i in range(len(theta)):
        if theta[i] > 0 and theta[i] < thetaMinimo:
            thetaMinimo = theta[i]
            indice = i
    return indice

def checkUnbound(y):
    for i in range(len(y)):
        if y[i] > 0:
            return False
    return True

def simplexFase1(matrizNaoBasica, vetor_b, custo_n, custo_b, nomeVariaveis):
    
    return

def simplexFase2(matrizNaoBasica, matrizBasica, vetor_b, custo_n, custo_b, nomeVariaveis, sentidoOriginal):
        solucaoOtima = False
        iteracao = 1
        while(not solucaoOtima):
            matrizInvertida = inverterMatriz(matrizBasica)
            xb = multiplicaMatrizPorVetor(matrizInvertida, vetor_b)
            vetorMultiplicador = calculoVetorMultiplicadorSimplex(custo_b, matrizInvertida)
            custosRelativos = calculoCustosRelativos(vetorMultiplicador, matrizNaoBasica, custo_n)
            indicePraEntrarBase = indiceDoValorMinimo(custosRelativos) # começa do 0, diferente do material
            if checarSolucaoOtima(custosRelativos):
                solucaoOtima = True
                z = calculoSolucaoOtima(xb, custo_b)
                if sentidoOriginal == "maximize":
                    z = z * -1.0
                printSolution(z, xb, nomeVariaveis)
                exit()
            y = calculoDirecaoSimplex(matrizInvertida, matrizNaoBasica, indicePraEntrarBase)
            if checkUnbound(y):
                print("Solution unbounded")
                exit()
            theta = calculoTheta(xb, y)
            indiceVariavelPraSair = indiceVariavelSair(theta)
            if (indiceVariavelPraSair == -1):
                print("Problem unfeasible")
                exit()
            else:
                for i in range(len(matrizBasica)):
                    matrizBasica[i][indiceVariavelPraSair], matrizNaoBasica[i][indicePraEntrarBase] = matrizNaoBasica[i][indicePraEntrarBase], matrizBasica[i][indiceVariavelPraSair]
                custo_b[indiceVariavelPraSair], custo_n[indicePraEntrarBase] = custo_n[indicePraEntrarBase], custo_b[indiceVariavelPraSair]
                iteracao += 1
            
def main():
    matrizNaoBasica, matrizBasica, vetorB, cn, cb, precisaFaseUm, nomeVariaveis, sentidoOriginal, numeroArtificiais = getInfoFromCplexFile()
    print (matrizNaoBasica)
    print (matrizBasica)
    print (vetorB)
    print (cn)
    print (cb)
    print (precisaFaseUm)
    print (nomeVariaveis)
    print (sentidoOriginal)
    print (numeroArtificiais)
    if precisaFaseUm == True:
        print("precisa fase 1")
        #simplexFase1(matrizNaoBasica, matrizBasica, vetorB, cn, cb, nomeVariaveis)
    else:
        print("nao precisa fase 1")
        simplexFase2(matrizNaoBasica, matrizBasica, vetorB, cn, cb, nomeVariaveis, sentidoOriginal)


if __name__ == "__main__":
    main()
