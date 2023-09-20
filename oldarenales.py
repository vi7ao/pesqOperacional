# guerreiros da br 277
# Vitor Teles Moreira Passos
# Renan Augusto Mattos Nutse
# Andre Negreli Conrado Bini

import cplex
import sys
import datetime

# para rodar digite no terminal o seguinte comando:
# python3 oldarenales.py "nome-do-arquivo"


# === tratamento da leitura do arquivo pelo terminal e do log de saída ===

logger = True # se true salva o log
inputFileName = sys.argv[-1]
dt = datetime.datetime.now()
strdt = dt.strftime("%Y%m%d%H%M%S")
file_name = "output" + strdt + ".txt"


def writeOutputToFile(text):
    if logger:
        with open(file_name, 'a') as file:
            file.write(text)

def getInfoFromCplexFile():
    arquivo = cplex.Cplex(inputFileName, "LP")
    nomeVariaveis = arquivo.variables.get_names()
    numeroVariaveis = arquivo.variables.get_num()
    numeroRestricoes = arquivo.linear_constraints.get_num()
    coeficientes = arquivo.linear_constraints.get_rows()
    variaveisArtificiais = 0
    precisaFaseUm = False
    # verificacao se precisa fase 1
    matrizNaoBasica = [[0.0] * numeroVariaveis for i in range(numeroRestricoes)]
    for i in range(numeroRestricoes):
        for j, var in enumerate(coeficientes[i].ind):
            matrizNaoBasica[i][var] = coeficientes[i].val[j]
    sinalRestricoes = arquivo.linear_constraints.get_senses()
    vetorB = [0.0] * numeroRestricoes
    flagMudadoSinal = [False] * numeroRestricoes
    for i in range(numeroRestricoes):
        vetorB[i] = arquivo.linear_constraints.get_rhs(i)
        if vetorB[i] < 0:
            if sinalRestricoes[i] == "L":
                sinalRestricoes[i] = "G"
            elif sinalRestricoes[i] == "G":
                sinalRestricoes[i] = "L"
            vetorB[i] = vetorB[i] * -1
            for j in range(len(matrizNaoBasica[i])):
                print("antes: " + str(matrizNaoBasica[i][j]))
                matrizNaoBasica[i][j] *= -1 
                print("depois: " + str(matrizNaoBasica[i][j]))
                flagMudadoSinal[i] = True
        if arquivo.linear_constraints.get_senses(i) == "E" :
            precisaFaseUm = True
            variaveisArtificiais += 1
        if arquivo.linear_constraints.get_senses(i) == "G":
            precisaFaseUm = True
            variaveisArtificiais += 1
    for i in range(numeroRestricoes):
        for j, var in enumerate(coeficientes[i].ind):
            print("antes: " + str(matrizNaoBasica[i][var]))
            if not flagMudadoSinal[i]:
                matrizNaoBasica[i][var] = coeficientes[i].val[j]
            print("depois: " + str(matrizNaoBasica[i][var]))
    cn = arquivo.objective.get_linear()
    cb = [0.0] * numeroRestricoes
    sinalRestricoes = arquivo.linear_constraints.get_senses()
    matrizBasica = []
    #variaveis folga/excesso
    for i in range(numeroRestricoes):
        matrizBasica.append([0] * numeroRestricoes)
    for i in range(numeroRestricoes):
        if sinalRestricoes[i] == "L":
            matrizBasica[i][i] = 1
        elif sinalRestricoes[i] == "G":
            matrizBasica[i][i] = -1
        else:
            matrizBasica[i][i] = 0        
    if arquivo.objective.get_sense() == -1:
        maxOrMin = "max"
    else:
        maxOrMin = "min"
    #criacao do vetor auxiliar das variaveis
    variaveisNaoBasicas = []
    variaveisBasicas = []
    auxCounter = 1
    for i in range (nomeVariaveis.__len__()):
        variaveisNaoBasicas.append(nomeVariaveis[i])
        auxCounter += 1
    for i in range (numeroRestricoes):
        variaveisBasicas.append("x" + str(auxCounter))
        auxCounter += 1
    writeOutputToFile("Matriz não básica: \n")
    for i in range(len(matrizNaoBasica)):
        writeOutputToFile(str(matrizNaoBasica[i]) + "\n")
    writeOutputToFile("Matriz básica: \n")
    for i in range(len(matrizBasica)):
        writeOutputToFile(str(matrizBasica[i]) + "\n")
    writeOutputToFile("Vetor B: \n")
    writeOutputToFile(str(vetorB) + "\n")
    writeOutputToFile("Vetor Cn: \n")
    writeOutputToFile(str(cn) + "\n")
    writeOutputToFile("Vetor Cb: \n")
    writeOutputToFile(str(cb) + "\n")
    writeOutputToFile("Max or Min: \n")
    writeOutputToFile(str(maxOrMin) + "\n")
    writeOutputToFile("Precisa Fase 1: \n")
    writeOutputToFile(str(precisaFaseUm) + "\n")
    writeOutputToFile("Variaveis Artificiais: \n")
    writeOutputToFile(str(variaveisArtificiais) + "\n")
    writeOutputToFile("Fim leitura arquivo \n")
    return matrizNaoBasica, matrizBasica, vetorB, cn, cb, maxOrMin, nomeVariaveis, precisaFaseUm, variaveisArtificiais, variaveisBasicas, variaveisNaoBasicas

# === funções da inversão de matriz ===

def gerarMatrizIdentidade(n):
    identidade = [[0 for x in range(n)] for y in range(n)]
    for i in range (n):
        identidade[i][i] = 1
    return identidade # ok

def gerarMatrizExtendida(matrizOriginal, identidade):
    matrizExtendida = []
    for i in range(len(matrizOriginal)):
        matrizExtendida.append(matrizOriginal[i] + identidade[i])
    writeOutputToFile("Matriz Extendida: \n")
    for(i) in range(len(matrizExtendida)):
        writeOutputToFile(str(matrizExtendida[i]) + "\n")
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
    writeOutputToFile("Matriz Invertida: \n")
    for(i) in range(len(matrizInvertida)):
        writeOutputToFile(str(matrizInvertida[i]) + "\n")
    return matrizInvertida # ok

# === funções do algoritmo simplex ===
def multiplicaMatrizPorVetor(matriz, vetor):
    vetorResultante = [0.0] * len(matriz)
    for i in range(len(matriz)):
        for j in range(len(vetor)):
            vetorResultante[i] += matriz[i][j] * vetor[j]
    writeOutputToFile("Vetor Resultante: \n")
    writeOutputToFile(str(vetorResultante) + "\n")
    return vetorResultante

def calculoVetorMultiplicadorSimplex(vetorCustoB, matrizInvertida):
    vetorMultiplicador = [0.0] * len(vetorCustoB)
    for i in range(len(vetorCustoB)):
        for j in range(len(vetorCustoB)):
            vetorMultiplicador[i] += vetorCustoB[j] * matrizInvertida[j][i]
    writeOutputToFile("Vetor Multiplicador: " + str(vetorMultiplicador) + "\n")
    return vetorMultiplicador

def calculoCustosRelativos(vetorMultiplicador, matrizNaoBasica, vetorCustoN):
    nVariaveis = len(vetorCustoN)
    nRestricoes = len(matrizNaoBasica)
    custosRelativos = [0.0] * nVariaveis
    for i in range(nVariaveis):
        col_i = [matrizNaoBasica[j][i] for j in range(nRestricoes)]
        custosRelativos[i] = vetorCustoN[i] - (vetorMultiplicador[i] * col_i[i])
    writeOutputToFile("Custos Relativos: " + str(custosRelativos) + "\n")
    return custosRelativos
        
def indiceDoValorMinimo(vetor):
    indiceMinimo = vetor.index(min(vetor))
    writeOutputToFile("Indice do Valor Minimo: " + str(indiceMinimo) + "\n")
    return indiceMinimo

def checarSolucaoOtima(custosRelativos):
    for i in range(len(custosRelativos)):
        if custosRelativos[i] < 0:
            writeOutputToFile("Custos relativos menor que zero, solucao nao otima \n")
            return False
    writeOutputToFile("Custos relativos maior que zero, solucao otima \n")
    return True

def calculoSolucaoOtima(xb, custo_b):
    z = 0
    for i in range(len(xb)):
        z += xb[i] * custo_b[i]
    return z

def calculoDirecaoSimplex(matrizInvertida, matrizNaoBasica, indicePraEntrarBase):
    coluna = [row[indicePraEntrarBase] for row in matrizNaoBasica]
    direcao = multiplicaMatrizPorVetor(matrizInvertida, coluna)
    writeOutputToFile("Direcao simples: " + str(direcao) + "\n")
    return direcao
    
def printSolution(z, xb, nomeVariaveis, variaveisBasicas, variaveisBasicasOriginais, iteracao):
    print("z: " + str(z))
    if (iteracao != 1):
        for i in range(len(variaveisBasicas)):
            print(variaveisBasicas[i] + ": " + str(xb[i]))
        for i in range(len(variaveisBasicasOriginais)):
            if(variaveisBasicasOriginais[i] not in variaveisBasicas):
                print(variaveisBasicasOriginais[i] + ": 0.0")
        print("As variaveis básicas são: ")
        for i in range(len(nomeVariaveis)):
            print(nomeVariaveis[i])
    else:
        for i in range(len(nomeVariaveis)):
            print(nomeVariaveis[i] + ": 0.0")
        print("As variaveis básicas são: ")
        for i in range(len(nomeVariaveis)):
            print(nomeVariaveis[i])
    return

def calculoTheta(xb, y):
    theta = [0.0] * len(xb)
    for i in range(len(xb)):
        if y[i] > 0:
            theta[i] = xb[i] / y[i]
    writeOutputToFile("Theta: " + str(theta) + "\n")
    return theta

def indiceVariavelSair(theta):
    indice = -1
    thetaMinimo = float("inf")
    for i in range(len(theta)):
        if theta[i] > 0 and theta[i] < thetaMinimo:
            thetaMinimo = theta[i]
            indice = i
    writeOutputToFile("Indice Variavel Sair: " + str(indice) + "\n")
    return indice

def checkUnbound(y):
    for i in range(len(y)):
        if y[i] > 0:
            return False
    return True

# === funções de variáveis artificiais ===
def adicionarVariaveisArtificiais(matrizNaoBasica, custo_n, numeroArtificiais):
    for i in range(numeroArtificiais):
        matrizNaoBasica = [row + [0.0] for row in matrizNaoBasica] 
        novaVariavelArtificial = [0.0] * len(custo_n) #ou aqui ta errado
        novaVariavelArtificial[len(custo_n) - numeroArtificiais + i] = 1.0 #acho que ta errado aqui 
        custo_n.extend(novaVariavelArtificial) #aqui ta errado provavelmente
    writeOutputToFile("Após adicionar variáveis artificiais: \n")
    writeOutputToFile("Matriz Não Básica: \n")
    for(i) in range(len(matrizNaoBasica)):
        writeOutputToFile(str(matrizNaoBasica[i]) + "\n")
    writeOutputToFile("Custo N: \n")
    writeOutputToFile(str(custo_n) + "\n")
    return matrizNaoBasica, custo_n

def verificaVariaveisArtificiaisNaBase(matrizBasica, numeroArtificiais):
    variaveisArtificiaisNaBase = 0
    for i in range(len(matrizBasica)):
        for j in range(len(matrizBasica[i])):
            if j >= len(matrizBasica[i]) - numeroArtificiais and matrizBasica[i][j] == 1:
                variaveisArtificiaisNaBase += 1
    writeOutputToFile("Variaveis Artificiais na Base: " + str(variaveisArtificiaisNaBase) + "\n")
    return variaveisArtificiaisNaBase

# === funções fase 1 do simplex ===
def otimalidadeFase1(indicePraEntrarBase, numeroArtificiais, matrizBasica):
    if indicePraEntrarBase >= 0:
        variaveisArtificiaisNaBase = verificaVariaveisArtificiaisNaBase(matrizBasica, numeroArtificiais)
        if variaveisArtificiaisNaBase > 0:
            writeOutputToFile("Problem unfeasible\n")
        else:
            return True
    
def fase1(matrizNaoBasica, matrizBasica, vetor_b, custo_n, custo_b, maxOrMin, nomeVariaveis, numeroArtificiais):
    iteracao = 1
    writeOutputToFile("Inicio Fase 1 \n")
    writeOutputToFile("Iteração: " + str(iteracao) + "\n")
    solucaoOtima = False
    while(not solucaoOtima):
        matrizInvertida = inverterMatriz(matrizBasica)
        xb = multiplicaMatrizPorVetor(matrizInvertida, vetor_b)
        vetorMultiplicador = calculoVetorMultiplicadorSimplex(custo_b, matrizInvertida)
        custosRelativos = calculoCustosRelativos(vetorMultiplicador, matrizNaoBasica, custo_n)
        indicePraEntrarBase = indiceDoValorMinimo(custosRelativos) # começa do 0, diferente do material
        if otimalidadeFase1(indicePraEntrarBase, numeroArtificiais, matrizBasica):
            solucaoOtima = True
            return matrizNaoBasica, matrizBasica, vetor_b, custo_n, custo_b, maxOrMin, nomeVariaveis
        y = calculoDirecaoSimplex(matrizInvertida, matrizNaoBasica, indicePraEntrarBase)
        if checkUnbound(y):
            writeOutputToFile("Problem unfeasible\n")
            exit()
        theta = calculoTheta(xb, y)
        indiceVariavelPraSair = indiceVariavelSair(theta)
        if (indiceVariavelPraSair == -1):
            writeOutputToFile("Problem unfeasible\n")
            exit()
        else:
            for i in range(len(matrizBasica)):
                matrizBasica[i][indiceVariavelPraSair], matrizNaoBasica[i][indicePraEntrarBase] = matrizNaoBasica[i][indicePraEntrarBase], matrizBasica[i][indiceVariavelPraSair]
                custo_b[indiceVariavelPraSair], custo_n[indicePraEntrarBase] = custo_n[indicePraEntrarBase], custo_b[indiceVariavelPraSair]
                iteracao += 1
                writeOutputToFile("Iteracao: " + str(iteracao) + "\n")
        variaveisArtificiaisNaBase = verificaVariaveisArtificiaisNaBase(matrizBasica, numeroArtificiais)
        if variaveisArtificiaisNaBase == 0:
            return matrizNaoBasica, matrizBasica, vetor_b, custo_n, custo_b, maxOrMin, nomeVariaveis

# === simplex fase 2 ===
def simplex(matrizNaoBasica, matrizBasica, vetor_b, custo_n, custo_b, maxOrMin, nomeVariaveis, variaveisBasicas, variaveisNaoBasicas):
        if maxOrMin == "max":
            custo_n = [-x for x in custo_n]
        solucaoOtima = False
        iteracao = 1
        varBasOriginais = variaveisBasicas.copy()
        varNaoBasOriginais = variaveisNaoBasicas.copy()
        writeOutputToFile("Inicio Simplex \n")
        writeOutputToFile("Iteracao: " + str(iteracao) + "\n")
        while(not solucaoOtima):
            matrizInvertida = inverterMatriz(matrizBasica)
            xb = multiplicaMatrizPorVetor(matrizInvertida, vetor_b)
            vetorMultiplicador = calculoVetorMultiplicadorSimplex(custo_b, matrizInvertida)
            custosRelativos = calculoCustosRelativos(vetorMultiplicador, matrizNaoBasica, custo_n)
            indicePraEntrarBase = indiceDoValorMinimo(custosRelativos) # começa do 0, diferente do material
            if checarSolucaoOtima(custosRelativos):
                solucaoOtima = True
                z = calculoSolucaoOtima(xb, custo_b)
                if maxOrMin == "max":
                    z = -z
                printSolution(z, xb, nomeVariaveis, variaveisBasicas, varBasOriginais, iteracao)
                writeOutputToFile("z: " + str(z) + "\n")
                for i in range(len(nomeVariaveis)):
                    writeOutputToFile("ainda nao tratei aqui k \n")
                    writeOutputToFile(nomeVariaveis[i] + ": " + str(xb[i]) + "\n")
                exit()
            y = calculoDirecaoSimplex(matrizInvertida, matrizNaoBasica, indicePraEntrarBase)
            if checkUnbound(y):
                writeOutputToFile("Solution unbounded \n")
                exit()
            theta = calculoTheta(xb, y)
            indiceVariavelPraSair = indiceVariavelSair(theta)
            if (indiceVariavelPraSair == -1):
                writeOutputToFile("Problem unfeasible\n")
                exit()
            else:
                for i in range(len(matrizBasica)):
                    matrizBasica[i][indiceVariavelPraSair], matrizNaoBasica[i][indicePraEntrarBase] = matrizNaoBasica[i][indicePraEntrarBase], matrizBasica[i][indiceVariavelPraSair]
                custo_b[indiceVariavelPraSair], custo_n[indicePraEntrarBase] = custo_n[indicePraEntrarBase], custo_b[indiceVariavelPraSair]
                writeOutputToFile("Variavel " + variaveisBasicas[indiceVariavelPraSair] + " sai da base e variavel " + variaveisNaoBasicas[indicePraEntrarBase] + " entra na base\n")
                variaveisBasicas[indiceVariavelPraSair], variaveisNaoBasicas[indicePraEntrarBase] = variaveisNaoBasicas[indicePraEntrarBase], variaveisBasicas[indiceVariavelPraSair]
                iteracao += 1
                writeOutputToFile("Iteracao: " + str(iteracao) + "\n")
            
def main():
    matrizNaoBasica, matrizBasica, vetorB, cn, cb, maxOrMin, nomeVariaveis, precisaFaseUm, numeroArtificiais, variaveisBasicas, variaveisNaoBasicas = getInfoFromCplexFile()
    if precisaFaseUm:
        matrizNaoBasica, cn = adicionarVariaveisArtificiais(matrizNaoBasica, cn, numeroArtificiais)
        matrizNaoBasica, matrizBasica, vetorB, cn, cb, maxOrMin, nomeVariaveis = fase1(matrizNaoBasica, matrizBasica, vetorB, cn, cb, maxOrMin, nomeVariaveis, numeroArtificiais)
    simplex(matrizNaoBasica, matrizBasica, vetorB, cn, cb, maxOrMin, nomeVariaveis, variaveisBasicas, variaveisNaoBasicas)

if __name__ == "__main__":
    main()