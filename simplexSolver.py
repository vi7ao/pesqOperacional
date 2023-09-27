import cplex
import sys
import datetime

# para rodar digite no terminal o seguinte comando:
# python3 simplexSolver.py "./testCases/nome-do-arquivo.cplex.lp"

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

def readCplexLPFile():
    arquivo = cplex.Cplex(inputFileName, "LP")
    nomeVariaveis = arquivo.variables.get_names()
    numeroVariaveis = arquivo.variables.get_num()
    numeroRestricoes = arquivo.linear_constraints.get_num()
    coeficientes = arquivo.linear_constraints.get_rows()
    precisaFaseUm = False
    matrizNaoBasica = [[0.0] * numeroVariaveis for i in range(numeroRestricoes)]
    for i in range(numeroRestricoes):
        for j, var in enumerate(coeficientes[i].ind):
            matrizNaoBasica[i][var] = coeficientes[i].val[j] #inicializa matriz naoBasica
    sinalRestricoes = arquivo.linear_constraints.get_senses()
    vetorB = [0.0] * numeroRestricoes
    for i in range(numeroRestricoes):
        vetorB[i] = arquivo.linear_constraints.get_rhs(i) #inicializa vetor b
        if vetorB[i] < 0: #tratamento de trocar sinal de restrições caso b < 0
            if sinalRestricoes[i] == "L":
                sinalRestricoes[i] = "G"
            elif sinalRestricoes[i] == "G":
                sinalRestricoes[i] = "L"
            vetorB[i] = vetorB[i] * -1
            for j in range(len(matrizNaoBasica[i])):
                matrizNaoBasica[i][j] *= -1 
        if sinalRestricoes[i] == "E" :
            precisaFaseUm = True
        if sinalRestricoes[i] == "G":
            precisaFaseUm = True
    cn = arquivo.objective.get_linear()
    cb = [0.0] * numeroRestricoes
    matrizBasica = []
    #variaveis folga/excesso
    for i in range(numeroRestricoes):
        matrizBasica.append([0] * numeroRestricoes)
    for i in range(numeroRestricoes):
        if sinalRestricoes[i] == "L": #add variavel de folga
            matrizBasica[i][i] = 1
        elif sinalRestricoes[i] == "G": #add variavel de excesso
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
    writeOutputToFile("Fim leitura arquivo \n")
    if precisaFaseUm:
        formulacaoProblemaArtificial(matrizNaoBasica, matrizBasica, vetorB, nomeVariaveis, variaveisBasicas, variaveisNaoBasicas, numeroRestricoes, cn, cb, maxOrMin)
    else:
        simplex(matrizNaoBasica, matrizBasica, vetorB, cn, cb, maxOrMin, nomeVariaveis, variaveisBasicas, variaveisNaoBasicas, False)

def formulacaoProblemaArtificial(matrizNaoBasica, matrizBasica, vetorB, nomeVariaveis, variaveisBasicas, variaveisNaoBasicas, numeroRestricoes, cnOriginal, cbOriginal, maxOrMin):
    writeOutputToFile("Inicio Fase 1 \n")
    matrizNaoBasicaFase1 = []
    for i in range(len(matrizNaoBasica)):
        matrizNaoBasicaFase1.append(matrizNaoBasica[i] + matrizBasica[i])
    writeOutputToFile("Matriz não básica fase 1: \n")
    for i in range(len(matrizNaoBasicaFase1)):
        writeOutputToFile(str(matrizNaoBasicaFase1[i]) + "\n")
    matrizBasicaFase1 = []
    writeOutputToFile("Matriz básica fase 1: \n")
    for i in range(numeroRestricoes):
        matrizBasicaFase1.append([0] * numeroRestricoes)
    for i in range(numeroRestricoes):
        matrizBasicaFase1[i][i] = 1
    for i in range(len(matrizBasicaFase1)):
        writeOutputToFile(str(matrizBasicaFase1[i]) + "\n")
    cnAux = []
    cn = [0.0] * len(matrizNaoBasicaFase1[0])
    writeOutputToFile("Vetor Cn fase 1: \n")
    writeOutputToFile(str(cn) + "\n")
    cb = [1.0] * numeroRestricoes
    writeOutputToFile("Vetor Cb fase 1: \n")
    writeOutputToFile(str(cb) + "\n")
    variaveisNaoBasicasFase1 = variaveisNaoBasicas.copy()
    for i in range(variaveisBasicas.__len__()):
        variaveisNaoBasicasFase1.append(variaveisBasicas[i])
    variaveisBasicasFase1 = []
    for i in range(numeroRestricoes):
        variaveisBasicasFase1.append("artf" + str(i+1))
    cnAux.append(cnOriginal + cbOriginal)
    cnFim = cnAux[0]
    cbAux = variaveisBasicasFase1.copy()
    writeOutputToFile("Variaveis basicas fase 1: \n")
    writeOutputToFile(str(variaveisBasicasFase1) + "\n")
    writeOutputToFile("Variaveis nao basicas fase 1: \n")
    writeOutputToFile(str(variaveisNaoBasicasFase1) + "\n")
    writeOutputToFile("Fim formulaçao problema artificial \n")
    iteracao = 1
    solucaoOtima = False   
    writeOutputToFile("Inicio simplex fase 1 \n")
    writeOutputToFile("Iteracao: " + str(iteracao) + "\n")
    while(not solucaoOtima):
        matrizInvertida = inverterMatriz(matrizBasicaFase1)
        xb = multiplicaMatrizPorVetor(matrizInvertida, vetorB)
        vetorMultiplicador = calculoVetorMultiplicadorSimplex(cb, matrizInvertida)
        custosRelativos = calculoCustosRelativos(vetorMultiplicador, matrizNaoBasicaFase1, cn)
        indicePraEntrarBase = indiceDoValorMinimo(custosRelativos)
        if custosRelativos[indicePraEntrarBase] >= 0:
            has_artificial = any(variable.startswith('artf') for variable in variaveisBasicasFase1)
            if has_artificial:
                writeOutputToFile
                print("Problema infactivel")
                writeOutputToFile("Problema infactivel \n")
                exit()
            else:
                solucaoOtima = True
                writeOutputToFile("goes to fase 2 \n")
                break
        y = calculoDirecaoSimplex(matrizInvertida, matrizNaoBasicaFase1, indicePraEntrarBase)
        if checkUnbound(y):
            print("Problema infactivel")
            writeOutputToFile("Problema infactivel \n")
            exit()
        theta = calculoTheta(xb, y)
        indiceVariavelPraSair = indiceVariavelSair(theta)
        for i in range(len(matrizBasicaFase1)):
            matrizBasicaFase1[i][indiceVariavelPraSair], matrizNaoBasicaFase1[i][indicePraEntrarBase] = matrizNaoBasicaFase1[i][indicePraEntrarBase], matrizBasicaFase1[i][indiceVariavelPraSair]
        cb[indiceVariavelPraSair], cn[indicePraEntrarBase] = cn[indicePraEntrarBase], cb[indiceVariavelPraSair]
        cbAux[indiceVariavelPraSair], cnFim[indicePraEntrarBase] = cnFim[indicePraEntrarBase], cbAux[indiceVariavelPraSair]
        writeOutputToFile("Variavel " + variaveisBasicasFase1[indiceVariavelPraSair] + " sai da base e variavel " + variaveisNaoBasicasFase1[indicePraEntrarBase] + " entra na base\n")
        variaveisBasicasFase1[indiceVariavelPraSair], variaveisNaoBasicasFase1[indicePraEntrarBase] = variaveisNaoBasicasFase1[indicePraEntrarBase], variaveisBasicasFase1[indiceVariavelPraSair]
        has_artificial = any(variable.startswith('artf') for variable in variaveisBasicasFase1)
        if has_artificial:
            iteracao += 1
            writeOutputToFile("Iteracao: " + str(iteracao) + "\n")
        else:
            writeOutputToFile("Fim simplex fase 1 \n")
            solucaoOtima = True    
    indicesVarArtf = [i for i, variable in enumerate(variaveisNaoBasicasFase1) if variable.startswith('artf')]
    matrizNaoBasicaFinal = [
        [matrizNaoBasicaFase1[i][j] for j in range(len(variaveisNaoBasicasFase1)) if j not in indicesVarArtf]
        for i in range(len(matrizNaoBasicaFase1))
    ]
    variaveisNaoBasicasFinal = [variable for i, variable in enumerate(variaveisNaoBasicasFase1) if i not in indicesVarArtf]
    cnFim = [cnFim[i] for i in range(len(cnFim)) if i not in indicesVarArtf]
    writeOutputToFile("Indo para a fase 2 com as seguintes infos: \n")
    writeOutputToFile("Posicao variaveis nao basicas: \n")
    writeOutputToFile(str(variaveisNaoBasicasFinal) + "\n")
    writeOutputToFile("Matriz não básica: \n")
    for i in range(len(matrizNaoBasicaFinal)):
        writeOutputToFile(str(matrizNaoBasicaFinal[i]) + "\n")
    writeOutputToFile("Vetor Cn: \n")
    writeOutputToFile(str(cnFim) + "\n")
    writeOutputToFile("Posicao variaveis basicas: \n")
    writeOutputToFile(str(variaveisBasicasFase1) + "\n")
    writeOutputToFile("Matriz básica: \n")
    for i in range(len(matrizBasicaFase1)):
        writeOutputToFile(str(matrizBasicaFase1[i]) + "\n")
    writeOutputToFile("Vetor B: \n")
    writeOutputToFile(str(vetorB) + "\n")
    writeOutputToFile("Vetor Cb: \n")
    writeOutputToFile(str(cbAux) + "\n")
    
    simplex(matrizNaoBasicaFinal, matrizBasicaFase1, vetorB, cnFim, cbAux, maxOrMin, nomeVariaveis, variaveisBasicasFase1, variaveisNaoBasicasFinal, True)

      
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
        custosRelativos[i] = vetorCustoN[i] - sum(col_i[j] * vetorMultiplicador[j] for j in range(nRestricoes))
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

# === simplex fase 2 ===
def simplex(matrizNaoBasica, matrizBasica, vetor_b, custo_n, custo_b, maxOrMin, nomeVariaveis, variaveisBasicas, variaveisNaoBasicas, veioFase1):
        if maxOrMin == "max":
            custo_n = [-x for x in custo_n]
            if veioFase1:
                custo_b = [-x for x in custo_b]
        solucaoOtima = False
        iteracao = 1
        varBasOriginais = variaveisBasicas.copy()
        varNaoBasOriginais = variaveisNaoBasicas.copy()
        writeOutputToFile("Inicio Simplex Fase 2 \n")
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
                print("Solution unbounded")
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
    readCplexLPFile()

if __name__ == "__main__":
    main()